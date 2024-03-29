//
// Run multiblock Euler amr
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>


// Headers for basic SAMRAI objects

#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h" 
#include "tbox/InputDatabase.h" 
#include "tbox/InputManager.h" 
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h" 
#include "tbox/Utilities.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "VariableDatabase.h"

// Headers for major algorithm/data structure objects

#include "BergerRigoutsos.h"
#include "BlockGridGeometry.h"
#include "MultiblockGriddingAlgorithm.h"
#include "MblkHyperbolicLevelIntegrator.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "MultiblockPatchHierarchy.h"
#include "StandardTagAndInitialize.h"
#include "TimeRefinementIntegrator.h"
#include "TimeRefinementLevelStrategy.h"
#include "VisItDataWriter.h"

// Header for application-specific algorithm/data structure object

#include "MblkEuler.h"

using namespace SAMRAI;

void setupHierarchy(
   tbox::Pointer<tbox::Database> main_input_db,
   tbox::Pointer<hier::MultiblockGridGeometry<NDIM> >& geometry,
   tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> >& mblk_hierarchy);


//
// ===================================== The main code =======================
//

int main( int argc, char *argv[])
{
   //
   // initialize startup
   //
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   string input_filename;
   string restart_read_dirname;
   int restore_num = 0;

   bool is_from_restart = false;

   if ( (argc != 2) && (argc != 4) ) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
		 << "<restart dir> <restore number> [options]\n"
		 << "  options:\n"
		 << "  none at this time"
		 << endl;
      tbox::SAMRAI_MPI::abort();
      return (-1);
   } else {
      input_filename = argv[1];
      if (argc == 4) {
	 restart_read_dirname = argv[2];
	 restore_num = atoi(argv[3]);
	 
	 is_from_restart = true;
      }
   }

   //
   // fire up the log file
   //
   string log_file_name = "MblkEuler.log";
   bool log_all_nodes = false;
   if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
   } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
   }

   tbox::plog << "input_filename       = " << input_filename << endl;
   tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
   tbox::plog << "restore_num          = " << restore_num << endl;

   //
   // Create input database and parse all data in input file.
   //
   
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   tbox::plog << "---- done parsing input file" << endl << endl;

   //
   // Retrieve "GlobalInputs" section of the input database and set
   // values accordingly.
   //
   if (input_db->keyExists("GlobalInputs")) {
      tbox::Pointer<tbox::Database> global_db =
	 input_db->getDatabase("GlobalInputs");
      if (global_db->keyExists("tag_clustering_method")) {
	 string tag_clustering_method =
	    global_db->getString("tag_clustering_method");
	 mesh::BergerRigoutsos<NDIM>::setClusteringOption(tag_clustering_method);
      }
      if (global_db->keyExists("refine_schedule_generation_method")) {
	 string refine_schedule_generation_method =
	    global_db->getString("refine_schedule_generation_method");
	 xfer::RefineSchedule<NDIM>::setScheduleGenerationMethod(
								 refine_schedule_generation_method);
      }
      if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
	 bool flag = global_db->
	    getBool("call_abort_in_serial_instead_of_exit");
	 tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }
   
   //
   // Retrieve "Main" section of the input database.  First, read dump
   // information, which is used for writing plot files.  Second, if 
   // proper restart information was given on command line, and the 
   // restart interval is non-zero, create a restart database.
   //
   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

   //
   //..Initialize VisIt
   //
   int viz_dump_interval = 0;
   if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
   }
   
   string viz_dump_dirname = "";
   string visit_dump_dirname = "";
   int visit_number_procs_per_file = 1;

   if ( viz_dump_interval > 0 ) {
      if (main_db->keyExists("viz_dump_dirname")) {
	 viz_dump_dirname = main_db->getString("viz_dump_dirname");
      }
      visit_dump_dirname = viz_dump_dirname;         
      if (viz_dump_dirname.empty()) {
	 TBOX_ERROR("main(): "
		    << "\nviz_dump_dirname is null ... "
		    << "\nThis must be specified for use with VisIt"
		    << endl);
      }
      if (main_db->keyExists("visit_number_procs_per_file")) {
	 visit_number_procs_per_file = 
	    main_db->getInteger("visit_number_procs_per_file");
      }
   }

   const bool viz_dump_data = (viz_dump_interval > 0);
   
   int restart_interval = 0;
   if (main_db->keyExists("restart_interval")) {
      restart_interval = main_db->getInteger("restart_interval");
   } 
   
   string restart_write_dirname;
   if ( restart_interval > 0 ) {
      if (main_db->keyExists("restart_write_dirname")) {
	 restart_write_dirname = main_db->getString("restart_write_dirname");
      } else {
	 TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname'"
		    << " not specifed in input file");
      }
   }
   
   bool use_refined_timestepping = true;
   if ( main_db->keyExists("timestepping")) {
      string timestepping_method = main_db->getString("timestepping");
      if (timestepping_method == "SYNCHRONIZED") {
	 use_refined_timestepping = false;
      }
   }
   
   const bool write_restart = (restart_interval > 0)
      && !(restart_write_dirname.empty());
   
   //
   // Get restart manager and root restart database.  If run is from
   // restart, open the restart file.
   //
   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
   
   if (is_from_restart) {
      restart_manager->
	 openRestartFile(restart_read_dirname, restore_num, 
			 tbox::SAMRAI_MPI::getNodes() );
   }

   //
   // Setup the timer manager to trace timing statistics during execution
   // of the code.  The list of timers is given in the TimerManager
   // section of the input file.  Timing information is stored in the
   // restart file.  Timers will automatically be initialized to their
   // previous state if the run is restarted, unless they are explicitly
   // reset using the TimerManager::resetAllTimers() routine. 
   //

   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   //
   // CREATE THE MULTIBLOCK HIERARCHY
   //
   tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > mblk_patch_hierarchy;
   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > geom;
   
   setupHierarchy(input_db, 
                  geom,
                  mblk_patch_hierarchy);
   
   //
   // -------------------- the patch operations --------------
   //
   MblkEuler* euler_model = new MblkEuler( "MblkEuler",
					   input_db,
					   geom->getBlockGeometryArray() );

   //
   // -------------------- the multiphase level operations --------------
   //
   tbox::Pointer<MblkHyperbolicLevelIntegrator> mblk_hyp_level_integrator =
      new MblkHyperbolicLevelIntegrator(
					"HyperbolicLevelIntegrator",
					input_db->getDatabase("HyperbolicLevelIntegrator"),
					euler_model,
					mblk_patch_hierarchy,
					true, 
					use_refined_timestepping);
   
   //
   // -------------------- the mesh refinement operations --------------
   //
   tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = 
      new mesh::StandardTagAndInitialize<NDIM>( "StandardTagAndInitialize",
						mblk_hyp_level_integrator,
						input_db->getDatabase("StandardTagAndInitialize"));
   
   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();
   
   tbox::Pointer<mesh::LoadBalancer<NDIM> >  load_balancer = 
      new mesh::LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));
   
   tbox::Pointer< mesh::MultiblockGriddingAlgorithm<NDIM> > 
      mblk_gridding_algorithm = 
      new mesh::MultiblockGriddingAlgorithm<NDIM>( "MultiblockGriddingAlgorithm",
						   input_db->getDatabase("GriddingAlgorithm"),
						   mblk_patch_hierarchy,
						   error_detector,
						   box_generator,
						   load_balancer );
   
   tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > time_integrator =
      new algs::TimeRefinementIntegrator<NDIM>( "TimeRefinementIntegrator",
						input_db->getDatabase("TimeRefinementIntegrator"),
						mblk_patch_hierarchy,
						mblk_hyp_level_integrator,
						mblk_gridding_algorithm);

   //
   // ----------------------------- Set up Visualization writer(s).
   // 
   bool is_multiblock = true;
   tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer = 
      new appu::VisItDataWriter<NDIM>("MblkEuler VisIt Writer",
                                      visit_dump_dirname,
                                      visit_number_procs_per_file,
                                      is_multiblock);
   euler_model->
      registerVisItDataWriter(visit_data_writer);
   
   //
   // Initialize hierarchy configuration and data on all patches.
   // Then, close restart file and write initial state for visualization.
   //
   double dt_now = time_integrator->initializeHierarchy();
   
   tbox::RestartManager::getManager()->closeRestartFile();
   
   //
   // After creating all objects and initializing their state, we
   // print the input database and variable database contents
   // to the log file.
   //
   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

   tbox::plog << "\nPrinting a summary of model input... " << endl;
   euler_model->printClassData(tbox::plog);
   
   if ( viz_dump_data ) {
      visit_data_writer->writePlotData(
				       mblk_patch_hierarchy,
				       time_integrator->getIntegratorStep(),
				       time_integrator->getIntegratorTime()); 
   }

   //
   // ==============================================================
   // Time step loop.  Note that the step count and integration
   // time are maintained by TimeRefinementIntegrator<NDIM>.
   // ==============================================================
   //

   double loop_time = time_integrator->getIntegratorTime();
   double loop_time_end = time_integrator->getEndTime();

   int iteration_num = time_integrator->getIntegratorStep();

   int old_log_style = 1;

   while ( (loop_time < loop_time_end) && time_integrator->stepsRemaining() ) {
      
      iteration_num = time_integrator->getIntegratorStep() + 1;
      
      if ( old_log_style ) {
	 tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	 tbox::pout << "At begining of timestep # " <<  iteration_num - 1 << endl;
	 tbox::pout << "Simulation time is " << loop_time << endl;
	 tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      }

      //
      // advance the heirarchy a timestep
      //
      double dt_new = time_integrator->advanceHierarchy(dt_now); 
      
      loop_time += dt_now;
      dt_now = dt_new;
      
      if ( !old_log_style ) {
	 //
	 // write out the timestep header
	 //
	 char my_line[256];
	 sprintf( my_line, "%4d time: %9.5e dt: %10.6e  ",
		  iteration_num,
		  loop_time,
		  dt_new );
	 
	 tbox::pout << my_line << endl;
      }

      if ( old_log_style ) {
	 tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	 tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
	 tbox::pout << "Simulation time is " << loop_time << endl;
	 tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      }

      //
      // At specified intervals, write restart files.
      //
      if ( write_restart ) {
	 
	 if ( (iteration_num % restart_interval) == 0 ) {
	    tbox::RestartManager::getManager()->
	       writeRestartFile(restart_write_dirname,
				iteration_num);
	 }
      }
      
      //
      // At specified intervals, write out data files for plotting.
      //
      if ( viz_dump_data ) {
	 if ( (iteration_num % viz_dump_interval) == 0 ) {
	    visit_data_writer->writePlotData(mblk_patch_hierarchy,
					     iteration_num,
					     loop_time);
	 }
      }
      
   }   //-----------------------------------------------END TIME STEPPING LOOP

   //
   // Output timer results.
   //
   tbox::TimerManager::getManager()->print(tbox::pout);
   
   //
   // At conclusion of simulation, deallocate objects.
   //
   visit_data_writer.setNull();
   time_integrator.setNull();
   mblk_gridding_algorithm.setNull();
   load_balancer.setNull();
   box_generator.setNull(); 
   error_detector.setNull();
   mblk_hyp_level_integrator.setNull();

   if (euler_model) delete euler_model;

   mblk_patch_hierarchy.setNull();
   geom.setNull();

   input_db.setNull();
   main_db.setNull();

   tbox::pout << "\nPASSED:  MblkEuler" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}

// ----------------------------------------------------------------

//
// this function builds the skeleton grid geometry 
//
void setupHierarchy(
   tbox::Pointer<tbox::Database> main_input_db,
   tbox::Pointer<hier::MultiblockGridGeometry<NDIM> >& geometry,
   tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> >& mblk_hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!main_input_db.isNull());
#endif
   
   tbox::Pointer<tbox::Database> mult_db =
      main_input_db->getDatabase("Multiblock");

   int num_blocks = mult_db->getInteger("num_blocks");

   tbox::Array< tbox::Pointer< hier::GridGeometry<NDIM> > > geoms(num_blocks);

   /*
    * Read the geometry information and build array of geometries
    */

   char geom_name[32];

   for (int nb = 0; nb < num_blocks; nb++) {

      sprintf(geom_name, "BlockGeometry%d", nb);
      if (mult_db->keyExists(geom_name)) {
         geoms[nb] =  
            new geom::BlockGridGeometry<NDIM>(
               geom_name,
               mult_db->getDatabase(geom_name),
               nb);
      } else {
         TBOX_ERROR("main::setupHierarchy(): could not find entry `"
                    << geom_name << "' in input.");
      }
   }

   geometry = new hier::MultiblockGridGeometry<NDIM>(geoms);

   mblk_hierarchy =
      new hier::MultiblockPatchHierarchy<NDIM>("MultiblockPatchHierarchy",
                                               mult_db, geometry, true);
   
}

