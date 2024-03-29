//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/Euler/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Main program for SAMRAI Euler gas dynamics sample application
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
#include "tbox/Array.h"
#include "BoxList.h"
#include "tbox/Database.h"
#include "Index.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "PatchLevel.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "VariableDatabase.h"

// Headers for major algorithm/data structure objects from SAMRAI

#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "GriddingAlgorithm.h"
#include "LoadBalancer.h"
#include "HyperbolicLevelIntegrator.h"
#include "PatchHierarchy.h"
#include "StandardTagAndInitialize.h"
#include "TimeRefinementIntegrator.h"
#include "TimeRefinementLevelStrategy.h"
#include "VisItDataWriter.h"

// Header for application-specific algorithm/data structure object

#include "Euler.h"

// Classes for autotesting.

#if (TESTING == 1)
#include "AutoTester.h"
#endif

using namespace SAMRAI;

/*
 ************************************************************************
 *                                                                      *
 * This is the main program for an AMR Euler gas dynamics application   *
 * built using SAMRAI.   The application program is constructed by      *
 * composing a variety of algorithm objects found in SAMRAI plus some   *
 * others that are specific to this application.   The following brief  *
 * discussion summarizes these objects.                                 *
 *                                                                      *
 *    hier::PatchHierarchy<NDIM> - A container for the AMR patch hierarchy and      *
 *       the data on the grid.                                          *
 *                                                                      *
 *    geom::CartesianGridGeometry<NDIM> - Defines and maintains the Cartesian       *
 *       coordinate system on the grid.  The hier::PatchHierarchy<NDIM>             *
 *       maintains a reference to this object.                          *
 *                                                                      *
 * A single overarching algorithm object drives the time integration    *
 * and adaptive gridding processes:                                     *
 *                                                                      *
 *    algs::TimeRefinementIntegrator<NDIM> - Coordinates time integration and       *
 *       adaptive gridding procedures for the various levels            *
 *       in the AMR patch hierarchy.  Local time refinement is          *
 *       employed during hierarchy integration; i.e., finer             *
 *       levels are advanced using smaller time increments than         *
 *       coarser level.  Thus, this object also invokes data            *
 *       synchronization procedures which couple the solution on        *
 *       different patch hierarchy levels.                              *
 *                                                                      *
 * The time refinement integrator is not specific to the numerical      *
 * methods used and the problem being solved.   It maintains references *
 * to two other finer grain algorithmic objects, more specific to       *
 * the problem at hand, with which it is configured when they are       *
 * passed into its constructor.   They are:                             *
 *                                                                      *
 *    algs::HyperbolicLevelIntegrator<NDIM> - Defines data management procedures    *
 *       for level integration, data synchronization between levels,    *
 *       and tagging cells for refinement.  These operations are        *
 *       tailored to explicit time integration algorithms used for      *
 *       hyperbolic systems of conservation laws, such as the Euler     *
 *       equations.  This integrator manages data for numerical         *
 *       routines that treat individual patches in the AMR patch        *
 *       hierarchy.  In this particular application, it maintains a     *
 *       pointer to the Euler object that defines variables and         *
 *       provides numerical routines for the Euler model.               *
 *                                                                      *
 *       Euler - Defines variables and numerical routines for the       *
 *          discrete Euler equations on each patch in the AMR           *
 *          hierarchy.                                                  *
 *                                                                      *
 *    mesh::GriddingAlgorithm<NDIM> - Drives the AMR patch hierarchy generation     *
 *       and regridding procedures.  This object maintains              *
 *       references to three other algorithmic objects with             *
 *       which it is configured when they are passed into its           *
 *       constructor.   They are:                                       *
 *                                                                      *
 *       mesh::BergerRigoutsos<NDIM> - Clusters cells tagged for refinement on a    *
 *          patch level into a collection of logically-rectangular      *
 *          box domains.                                                *
 *                                                                      *
 *       mesh::LoadBalancer<NDIM> - Processes the boxes generated by the            *
 *          mesh::BergerRigoutsos<NDIM> algorithm into a configuration from         *
 *          which patches are contructed.  The algorithm we use in this *
 *          class assumes a spatially-uniform workload distribution;    *
 *          thus, it attempts to produce a collection of boxes          *
 *          each of which contains the same number of cells.  The       *
 *          load balancer also assigns patches to processors.           *
 *                                                                      *
 *       mesh::StandardTagAndInitialize<NDIM> - Couples the gridding algorithm      *
 *          to the HyperbolicIntegrator. Selects cells for              *
 *          refinement based on either Gradient detection, Richardson   *
 *          extrapolation, or pre-defined Refine box region.  The       *
 *          object maintains a pointer to the algs::HyperbolicLevelIntegrator<NDIM>,*
 *          which is passed into its constructor, for this purpose.     *
 *                                                                      *
 ************************************************************************
 */

/* 
 *******************************************************************
 *                                                                 *
 * For each run, the input filename and restart information        *
 * (if needed) must be given on the command line.                  *
 *                                                                 *
 *      For non-restarted case, command line is:                   *
 *                                                                 *
 *          executable <input file name>                           *
 *                                                                 *
 *      For restarted run, command line is:                        *
 *                                                                 *
 *          executable <input file name> <restart directory> \     *
 *                     <restart number>                            *
 *                                                                 *
 * Accessory routines used within the main program:                *
 *                                                                 *
 *   dumpVizData1dPencil - Writes 1d pencil of Euler solution data *
 *      to plot files so that it may be viewed in MatLab.  This    *
 *      routine assumes a single patch level in 2d and 3d.  In     *
 *      other words, it only plots data on level zero.  It can     *
 *      handle AMR in 1d.                                          *
 *                                                                 *
 *******************************************************************
 */

void dumpMatlabData1dPencil(const string& dirname, 
                            const string& filename,
                            const int ext,
                            const double plot_time,
                            const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                            const int pencil_direction,
                            const bool default_pencil,
                            const tbox::Array<int>& pencil_index,
                            Euler* euler_model);

int main( int argc, char *argv[] )
{

   /*
    * Initialize tbox::MPI and SAMRAI, enable logging, and process command line.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   int num_failures = 0;

   {

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

   tbox::plog << "input_filename = " << input_filename << endl;
   tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
   tbox::plog << "restore_num = " << restore_num << endl;

   /*
    * Create input database and parse all data in input file.
    */
   
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   /*
    * Retrieve "GlobalInputs" section of the input database and set 
    * values accordingly.
    */

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

   /*
    * Retrieve "Main" section of the input database.  First, read dump
    * information, which is used for writing plot files.  Second, if 
    * proper restart information was given on command line, and the 
    * restart interval is non-zero, create a restart database.
    */

   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

   string log_file_name = "euler.log";
   if (main_db->keyExists("log_file_name")) {
      log_file_name = main_db->getString("log_file_name");
   }
   bool log_all_nodes = false;
   if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
   }
   if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
   } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
   }

   int viz_dump_interval = 0;
   if (main_db->keyExists("viz_dump_interval")){
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
   }

   tbox::Array<string> viz_writer(1);
   viz_writer[0] = "VisIt";
   string viz_dump_filename;
   string vizamrai_dump_dirname;
   string visit_dump_dirname;
   bool uses_vizamrai = false;
   bool uses_visit    = false;
   int visit_number_procs_per_file = 1;
   if ( viz_dump_interval > 0 ) {
      if (main_db->keyExists("viz_writer")) {
         viz_writer = main_db->getStringArray("viz_writer");
      }
      if (main_db->keyExists("viz_dump_filename")) {
         viz_dump_filename = main_db->getString("viz_dump_filename");
      }
      string viz_dump_dirname;
      if (main_db->keyExists("viz_dump_dirname")) {
         viz_dump_dirname = main_db->getString("viz_dump_dirname");
      }
      for (int i = 0; i < viz_writer.getSize(); i++) {
         if (viz_writer[i] == "Vizamrai") uses_vizamrai = true;
         if (viz_writer[i] == "VisIt") uses_visit = true;
      }
      if (uses_vizamrai && !uses_visit) {
         vizamrai_dump_dirname = viz_dump_dirname;
      } else if (!uses_vizamrai && uses_visit) {
         visit_dump_dirname = viz_dump_dirname;         
      } else if (uses_vizamrai && uses_visit) {
         vizamrai_dump_dirname = viz_dump_dirname + "_Vizamrai";
         visit_dump_dirname = viz_dump_dirname + "_VisIt";
      } else {
         TBOX_ERROR("main(): "
                    << "\nUnrecognized 'viz_writer' entry..."
                    << "\nOptions are 'Vizamrai' and/or 'VisIt'" 
                    << endl);
      }
      if (uses_visit) {
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
      if (uses_vizamrai) {
         if (viz_dump_filename.empty()) {
         TBOX_ERROR("main(): "
                    << "\nviz_dump_filename is null ... "
                    << "\nThis must be specified for use with Vizamrai"
                    << endl);
         }
      }
   }

   string matlab_dump_filename;
   string matlab_dump_dirname;
   int matlab_dump_interval = 0;
   int  matlab_pencil_direction = 0;
   tbox::Array<int> matlab_pencil_index(NDIM-1);
   bool matlab_default_pencil = true; 
   for (int id = 0; id < NDIM-1; id++) {
      matlab_pencil_index[id] = 0;
   }
   
   if (main_db->keyExists("matlab_dump_interval")) {
      matlab_dump_interval = main_db->getInteger("matlab_dump_interval");
   }
   if (matlab_dump_interval > 0) {
      if (main_db->keyExists("matlab_dump_filename")) {
         matlab_dump_filename = main_db->getString("matlab_dump_filename");
      }
      if (main_db->keyExists("matlab_dump_dirname")) {
         matlab_dump_dirname = main_db->getString("matlab_dump_dirname");
      }
      if (main_db->keyExists("matlab_pencil_direction")) {
         matlab_pencil_direction = 
            main_db->getInteger("matlab_pencil_direction");
      }
      if (main_db->keyExists("matlab_pencil_index")) {
         matlab_default_pencil = false;
         matlab_pencil_index = main_db->getIntegerArray("matlab_pencil_index");
         if (matlab_pencil_index.getSize() != NDIM-1) {
            TBOX_ERROR("`matlab_pencil_index' has "
                       << matlab_pencil_index.getSize() << " values in input. "
                       << NDIM-1 << " values must be specified when default"
                       << " is overridden.");
         }
      }
   }

   int restart_interval = 0;
   if (main_db->keyExists("restart_interval")) {
      restart_interval = main_db->getInteger("restart_interval");
   }

   string restart_write_dirname;
   if ( restart_interval > 0 ) {
      if (main_db->keyExists("restart_write_dirname")) {
         restart_write_dirname = main_db->getString("restart_write_dirname");
      } else {
         TBOX_ERROR("`restart_interval' > 0, but key `restart_write_dirname'"
                    << "not found in input file.");
      }
   }

   bool use_refined_timestepping = true;
   if ( main_db->keyExists("timestepping")) {
      string timestepping_method = main_db->getString("timestepping");
      if (timestepping_method == "SYNCHRONIZED") {
         use_refined_timestepping = false;
      }
   }

#if (TESTING == 1) && !(HAVE_HDF5)
   /* 
    * If we are autotesting on a system w/o HDF5, the read from
    * restart will result in an error.  We want this to happen
    * for users, so they know there is a problem with the restart,
    * but we don't want it to happen when autotesting.
    */
   is_from_restart  = false;
   restart_interval = 0;
#endif

   const bool write_restart = (restart_interval > 0)
                               && !(restart_write_dirname.empty());

   /*
    * Get restart manager and root restart database.  If run is from
    * restart, open the restart file.
    */

   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

   if (is_from_restart) {
      restart_manager->
         openRestartFile(restart_read_dirname, restore_num,
                         tbox::SAMRAI_MPI::getNodes());
   }

   /*
    * Setup the timer manager to trace timing statistics during execution
    * of the code.  The list of timers is given in the tbox::TimerManager
    * section of the input file.  Timing information is stored in the
    * restart file.  Timers will automatically be initialized to their
    * previous state if the run is restarted, unless they are explicitly
    * reset using the tbox::TimerManager::resetAllTimers() routine. 
    */

   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   /*
    * Create major algorithm and data objects which comprise application.
    * Each object is initialized either from input data or restart
    * files, or a combination of both.  Refer to each class constructor
    * for details.  For more information on the composition of objects
    * and the roles they play in this application, see comments at top of file.
    */

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
      new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                input_db->getDatabase("CartesianGeometry"));

   tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy = 
      new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   Euler* euler_model = new Euler("Euler",
                                  input_db->getDatabase("Euler"),
                                  grid_geometry);
  
   tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> > hyp_level_integrator =
      new algs::HyperbolicLevelIntegrator<NDIM>("HyperbolicLevelIntegrator",
                                    input_db->getDatabase(
                                              "HyperbolicLevelIntegrator"),
                                    euler_model, true,
                                    use_refined_timestepping);

   tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = 
      new mesh::StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",      
          hyp_level_integrator,
          input_db->getDatabase("StandardTagAndInitialize"));

   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = 
      new mesh::LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));

   tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm = 
      new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            error_detector,
                            box_generator,
                            load_balancer);

   tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > time_integrator =
      new algs::TimeRefinementIntegrator<NDIM>("TimeRefinementIntegrator",
                                   input_db->getDatabase(
                                             "TimeRefinementIntegrator"),
                                   patch_hierarchy,
                                   hyp_level_integrator,
                                   gridding_algorithm);

   /*
    * Set up Visualization writer(s).  Note that the Euler application
    * creates some derived data quantities so we register the Euler model
    * as a derived data writer.  If no derived data is written, this step
    * is not necessary.
    */
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > vizamrai_data_writer = 
      new appu::CartesianVizamraiDataWriter<NDIM>("Euler Vizamrai Writer");

#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer = 
      new appu::VisItDataWriter<NDIM>("Euler VisIt Writer",
                          visit_dump_dirname,
                          visit_number_procs_per_file);
   if (uses_visit) {
      euler_model->registerVisItDataWriter(visit_data_writer);
   }
#endif

   if (uses_vizamrai) {
      euler_model->registerVizamraiDataWriter(vizamrai_data_writer);
   }
   
   if (viz_dump_interval > 0) {
      if (uses_vizamrai) {
         vizamrai_data_writer->setDirectoryName(vizamrai_dump_dirname);
         vizamrai_data_writer->setPlotDataToFloat();
         vizamrai_data_writer->setFinestLevelToPlot(
            gridding_algorithm->getMaxLevels()-1);
         for (int ln = 1; ln < gridding_algorithm->getMaxLevels(); ln++) {
            const hier::IntVector<NDIM>& lratio =
               gridding_algorithm->getRatioToCoarserLevel(ln);
            vizamrai_data_writer->setRatioToCoarserLevel(ln, lratio);
         }
         vizamrai_data_writer->setDerivedDataWriter(euler_model);
      }
   }

   /*
    * Initialize hierarchy configuration and data on all patches.
    * Then, close restart file and write initial state for visualization.
    */

   double dt_now = time_integrator->initializeHierarchy();

   tbox::RestartManager::getManager()->closeRestartFile();


#if (TESTING == 1)
   /*
    * Create the autotesting component which will verify correctness
    * of the problem. If no automated testing is done, the object does 
    * not get used.
    */
   AutoTester autotester("AutoTester", input_db);
#endif

   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents
    * to the log file.
    */

#if 0
   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

#endif
   tbox::plog << "\nCheck Euler data... " << endl;
   euler_model->printClassData(tbox::plog);

   /*
    * Create timers for measuring I/O.
    */
   tbox::Pointer<tbox::Timer> t_write_viz = tbox::TimerManager::getManager()->
                                getTimer("apps::main::write_viz");
   tbox::Pointer<tbox::Timer> t_write_restart = tbox::TimerManager::getManager()->
                                    getTimer("apps::main::write_restart");

   t_write_viz->start();
   if (matlab_dump_interval > 0) {
      dumpMatlabData1dPencil(matlab_dump_dirname,
                             matlab_dump_filename,
                             time_integrator->getIntegratorStep(),
                             time_integrator->getIntegratorTime(),
                             patch_hierarchy,
                             matlab_pencil_direction,
                             matlab_default_pencil,
                             matlab_pencil_index,
                             euler_model);
   } 
   if (viz_dump_interval > 0) { 
      if (uses_vizamrai) {
         vizamrai_data_writer->writePlotData(
            patch_hierarchy,
            viz_dump_filename,
            time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime()); 
      }
#ifdef HAVE_HDF5
      if (uses_visit) {
         visit_data_writer->writePlotData(
            patch_hierarchy,
            time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime()); 
      }
#endif
   }
   t_write_viz->stop();

   /*
    * Time step loop.  Note that the step count and integration
    * time are maintained by algs::TimeRefinementIntegrator<NDIM>.
    */

   double loop_time = time_integrator->getIntegratorTime();
   double loop_time_end = time_integrator->getEndTime();

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
   num_failures += autotester.evalTestData(time_integrator->getIntegratorStep(),
                                           patch_hierarchy,
                                           time_integrator,
                                           hyp_level_integrator,
                                           gridding_algorithm);
#endif

   while ( (loop_time < loop_time_end) && 
          time_integrator->stepsRemaining() ) {

      int iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::plog <<endl<<endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "Current dt is " << dt_now << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog <<endl<<endl;

      double dt_new = time_integrator->advanceHierarchy(dt_now);

      loop_time += dt_now;
      dt_now = dt_new;

      tbox::plog <<endl<<endl;
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
      tbox::plog <<endl<<endl;

      /*
       * At specified intervals, write restart files.
       */
      if ( write_restart ) {

         if ( (iteration_num % restart_interval) == 0 ) { 
            t_write_restart->start();
            tbox::RestartManager::getManager()->
               writeRestartFile(restart_write_dirname,
                                iteration_num);
            t_write_restart->stop();
         }
      }

      /*
       * At specified intervals, write out data files for plotting.
       */
      t_write_viz->start();
      if ( (viz_dump_interval > 0) 
           && (iteration_num % viz_dump_interval) == 0 ) {
         if (uses_vizamrai) {
            vizamrai_data_writer->writePlotData(patch_hierarchy,
                                                viz_dump_filename,
                                                iteration_num,
                                                loop_time);
         }
#ifdef HAVE_HDF5
         if (uses_visit) {
            visit_data_writer->writePlotData(patch_hierarchy,
                                             iteration_num,
                                             loop_time);
         }
#endif
         
      }
      if ((matlab_dump_interval > 0)
           && (iteration_num % matlab_dump_interval) == 0 ) {
          dumpMatlabData1dPencil(matlab_dump_dirname,
                                 matlab_dump_filename, 
                                 iteration_num,
                                 loop_time,
                                 patch_hierarchy,
                                 matlab_pencil_direction,
                                 matlab_default_pencil,
                                 matlab_pencil_index,
                                 euler_model);
      }
      t_write_viz->stop();

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
      num_failures += autotester.evalTestData( iteration_num,
                                               patch_hierarchy,
                                               time_integrator,
                                               hyp_level_integrator,
                                               gridding_algorithm);
#endif

      /*
       * Write byte transfer information to log file.
       */
#if 0
      char num_buf[8];
      sprintf(num_buf, "%02d", iteration_num);
      tbox::plog << "Step " << num_buf 
           << " P" << tbox::SAMRAI_MPI::getRank() 
           << ": " << tbox::SAMRAI_MPI::getIncomingBytes() 
           << " bytes in" << endl;
#endif

   }

   /*
    * Output timer results.
    */
   tbox::TimerManager::getManager()->print(tbox::plog);

   /*
    * At conclusion of simulation, deallocate objects.
    */
   patch_hierarchy.setNull();
   grid_geometry.setNull();

   box_generator.setNull();
   load_balancer.setNull();
   hyp_level_integrator.setNull();
   error_detector.setNull();
   gridding_algorithm.setNull();
   time_integrator.setNull();
   vizamrai_data_writer.setNull();
#ifdef HAVE_HDF5
   visit_data_writer.setNull();
#endif

   if (euler_model) delete euler_model;

   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  Euler" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(num_failures); 
}

void dumpMatlabData1dPencil(const string& dirname,
                            const string& filename,
                            const int ext,
                            const double plot_time,
                            const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                            const int pencil_direction,
                            const bool default_pencil,
                            const tbox::Array<int>& pencil_index,
                            Euler* euler_model)
{

   /*
    * Compute the boxes to write out data at each level of the hierarchy.
    */

   int nlevels = 1;

#if (NDIM == 1)
   nlevels = hierarchy->getNumberOfLevels();
#endif

   hier::BoxList<NDIM> domain(hierarchy->getGridGeometry()->getPhysicalDomain());
   hier::Box<NDIM> pencil_box(domain.getBoundingBox());

#if (NDIM > 1)
   int indx = 0;
   int id = 0;
   tbox::Array<int> tmp(NDIM-1);
   for (id = 0; id < NDIM-1; id++) {
      tmp[id] = pencil_index[id];
   }
   if (default_pencil) {
      hier::Index<NDIM> ifirst = domain.getBoundingBox().lower();
      indx = 0;
      for (id = 0; id < NDIM; id++) {
         if (id != pencil_direction) { 
            tmp[indx] = ifirst(id);
            indx++;
         }
      } 
   }
   indx = 0;
   for (id = 0; id < NDIM; id++) {
      if (id != pencil_direction) {
         pencil_box.lower(id) = tmp[indx];
         pencil_box.upper(id) = tmp[indx];
         indx++;
      }
   }
#endif

   tbox::Array<hier::BoxList<NDIM> > outboxes(nlevels);

   for (int l1 = 0; l1 < nlevels; l1++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(l1);
      outboxes[l1] = level->getBoxes();

      if (l1 < nlevels-1) {

	 tbox::Pointer<hier::PatchLevel<NDIM> > finer_level = 
	    hierarchy->getPatchLevel(l1+1);
         hier::IntVector<NDIM> coarsen_ratio 
	    = finer_level->getRatioToCoarserLevel();
         hier::BoxList<NDIM> takeaway = finer_level->getBoxes();
         takeaway.coarsen(coarsen_ratio);
         outboxes[l1].removeIntersections(takeaway);
      }

   }

   /*
    * Create matlab filename and open the output stream.
    */
   
   string dump_filename = filename;
   if (!dirname.empty()) {
      tbox::Utilities::recursiveMkdir(dirname);
      dump_filename = dirname;
      dump_filename += "/";
      dump_filename += filename;
   }

   const int size = dump_filename.length() + 16;
   char *buffer = new char[size];

   if (tbox::SAMRAI_MPI::getNodes() > 1) {
      sprintf(buffer, "%s.%04d.dat.%05d", dump_filename.c_str(),
              ext, tbox::SAMRAI_MPI::getRank());
   }
   else {
      sprintf(buffer, "%s_%04d.dat", dump_filename.c_str(), ext);
   }

   /*
    * Open a new output file having name name of buffer character array.
    */

   ofstream outfile(buffer, ios::out);
   outfile.setf(ios::scientific);
   outfile.precision(10);

   delete [] buffer;

   /*
    * There are 7 values dumped for every cell.  Here we dump the time.
    */
   for (int i = 0; i < 6+1; i++) {
      outfile << plot_time << "  ";
   }
   outfile << endl;

   euler_model->setDataContext(
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT") );

   for (int l5 = 0; l5 < nlevels; l5++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(l5);

      hier::Box<NDIM> level_pencil_box = pencil_box;
      if (l5 > 0) {
         level_pencil_box.refine(level->getRatio());
      }
      
      for (hier::PatchLevel<NDIM>::Iterator i(level); i; i++) {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(i());
         hier::Box<NDIM> pbox = patch->getBox();

         for (hier::BoxList<NDIM>::Iterator b(outboxes[l5]); b; b++) {
            const hier::Box<NDIM> box = b() * pbox * level_pencil_box; 
            
            euler_model->writeData1dPencil(patch,
                                           box,
                                           pencil_direction,
                                           outfile);
         }

      }

   }

   euler_model->clearDataContext();

   outfile.close();

}

