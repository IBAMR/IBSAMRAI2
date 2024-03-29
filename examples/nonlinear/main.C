//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/nonlinear/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1981 $
// Modified:    $LastChangedDate: 2008-02-13 11:09:35 -0800 (Wed, 13 Feb 2008) $
// Description: Main program for modified-Bratu problem
//

#include "SAMRAI_config.h"

#include <fstream>
#include <string>
#include <sys/stat.h>
#include <string>
using namespace std;

/*
 * Headers for basic SAMRAI objects used in this sample code.
 */
#include "tbox/IOStream.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "tbox/Database.h"
#include "tbox/FileStream.h"
#include "Geometry.h"
#include "GridGeometry.h"
#include "tbox/InputManager.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/SAMRAI_MPI.h"
#define included_String
#include "tbox/Utilities.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "VariableDatabase.h"

/*
 * Headers for major algorithm/data structure objects from SAMRAI
 */
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "GriddingAlgorithm.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "StandardTagAndInitialize.h"
#include "VisItDataWriter.h"

/*
 * Headers for major application-specific and solver testing objects
 */
#include "ImplicitIntegrator.h"
#ifdef HAVE_PETSC 
#include "SNES_SAMRAIContext.h"
#endif
#ifdef HAVE_SUNDIALS
#include "KINSOL_SAMRAIContext.h"
#endif
#include "ModifiedBratuProblem.h"
#include "NonlinearSolverStrategy.h"

using namespace SAMRAI;

/*!
 * @brief Main program for nonlinear example.
 *                                                                      
 * This is the main program for a sample application that employs       
 * implicit time integration on an AMR patch hierarchy and requires     
 * a nonlinear solver to advance the solution at each time step.        
 * The problem to be solved is:                                         
 *                                                                    
 *        du/dt = div(D(x)*grad(u)) + lambda*exp(u) + f(x,u)           
 *                                                                     
 * The application program is constructed by composing a variety of    
 * algorithm objects found in SAMRAI plus numerical operations that     
 * are specific to this application.  The following discussion          
 * summarizes these objects.                                            
 *                                                                      
 *    hier::PatchHierarchy<NDIM> - A container for the AMR patch        
 *        hierarchy and the data on the grid.                           
 *                                                                      
 *    geom::CartesianGridGeometry<NDIM> - Defines and maintains the     
 *       Cartesian coordinate system on the grid.  The                  
 *       hier::PatchHierarchy<NDIM> maintains a reference to this       
 *       object.                                                        
 *                                                                      
 * The implicit integrator algorithm uses these two components to       
 * integrate data one timestep.  Time-stepping is done in the main      
 * routine. This integrator uses the same timestep across all patches,  
 * regardless of the refinement.                                        
 *                                                                     
 *    algs::ImplicitIntegrator<NDIM> - advances solution on patches.    
 *             Controls integration of quantities over entire patch     
 *             hierarchy. Takes the user-defined class that defines     
 *             patch operations (ModifiedBratuProblem in this case)     
 *           as an argument.                                            
 *                                                                      
 * The user-supplied object defines characteristics and operations to   
 * be carried out on each patch.  We pass this object to the integrator 
 * which will orchestrate the integration on all patches.               
 *                                                                      
 *    ModifiedBratuProblem - Defines variables and numerical routines   
 *            on the hierarchy.  This includes routines needed by the   
 *          implicit integrator as well as nonlinear solvers used       
 *          by the integrator.                                          
 *                                                                      
 *    mesh::GriddingAlgorithm<NDIM> - Drives the AMR patch hierarchy    
 *       generation and regridding procedures.  This object maintains   
 *       references to three other algorithmic objects with             
 *       which it is configured when they are passed into its           
 *       constructor.   They are:                                       
 *                                                                      
 *       mesh::BergerRigoutsos<NDIM> - Clusters cells tagged for        
 *          refinement on a patch level into a collection of            
 *          logically-rectangular box domains.                          
 *                                                                      
 *       mesh::LoadBalancer<NDIM> - Processes the boxes generated by    
 *          the mesh::BergerRigoutsos<NDIM> algorithm into a            
 *          configuration from                                          
 *          which patches are contructed.  The routines in this         
 *          class assume a spatially-uniform workload distribution;     
 *          thus, they attempt to produce a collection of boxes         
 *          each of which contains the same number of cells.  The       
 *          load balancer also assigns patches to processors.           
 *                                                                      
 *       mesh::StandardTagAndInitialize<NDIM> -                         
 *          Couples the gridding algorithm to                           
 *          the ModifiedBratuProblem class which selects cells for      
 *          refinement based on criteria defined specifically for       
 *          the problem.  This object maintains a pointer to the        
 *          ModifiedBratuProblem object, which is passed into its       
 *          constructor, for this purpose.                              
 *                                                                      
 */

int main( int argc, char *argv[] )
{

   /*
    * Initialize MPI, SAMRAI.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
#ifdef HAVE_PETSC
   PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
#endif
   tbox::SAMRAIManager::startup();


   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {


#if (NDIM > 1)

#if !defined(HAVE_PETSC) || !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)
      tbox::pout << "This example requires the packages PETSC, SUNDIALS, "
		 << "\nand HYPRE to work properly.  SAMRAI was not configured"
		 << "\nwith one or more of these packages." 
		 << endl;
#else


      /*
       * Process command line arguments.  For each run, the input 
       * filename must be specified.  Usage is:
       * 
       *     executable <input file name>
       */
      string input_filename;

      if (argc != 2)  {
	 TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n" 
		    << "  options:\n"
		    << "  none at this time" << endl);
      } else {
	 input_filename = argv[1];
      }


      /*
       * Create input database and parse all data in input file.
       */

      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section of the input database.  Read
       * dump information, which is used for writing vizamrai plot files.
       */

      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

      string base_name = "default";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      string log_file_name = base_name + ".log";
      tbox::PIO::logOnlyNodeZero(log_file_name);
      tbox::plog << "input_filename = " << input_filename << endl;

      int viz_dump_interval = 0;
      if (main_db->keyExists("viz_dump_interval")) {
	 viz_dump_interval = main_db->getInteger("viz_dump_interval");
      }

      tbox::Array<string> viz_writer(1);
      viz_writer[0] = "Vizamrai";
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

      int regrid_interval = 0;
      if (main_db->keyExists("regrid_interval")) {
	 regrid_interval = main_db->getInteger("regrid_interval");
      }

      string nonlinear_solver_package = "KINSOL";
      if (main_db->keyExists("nonlinear_solver_package")) {
	 nonlinear_solver_package =
	    main_db->getString("nonlinear_solver_package");
      }

      /*
       * Setup the timer manager to trace timing statistics during execution
       * of the code.  The list of timers is given in the tbox::TimerManager
       * section of the input file.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Create major algorithm and data objects which comprise application.
       * Each object is initialized from input data.
       * Refer to each class constructor
       * for details.  For more information on the composition of objects
       * and the roles they play in this application, see comments at top
       * of file.
       */

      tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer;
      if ( uses_visit ) {
	 visit_data_writer =
	    new appu::VisItDataWriter<NDIM>("Bratu VisIt Writer",
					    visit_dump_dirname,
					    visit_number_procs_per_file);
      }

      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer;
      if ( uses_vizamrai ) {
	 viz_data_writer =
	    new appu::CartesianVizamraiDataWriter<NDIM>("Bratu Viz Writer");
      }

      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
	 new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
					       input_db->getDatabase("CartesianGeometry"));

      tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
	 new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

      ModifiedBratuProblem *bratu_model = new ModifiedBratuProblem(
	 "ModifiedBratuProblem",
	 input_db->getDatabase("ModifiedBratuProblem"),
	 grid_geometry, viz_data_writer, visit_data_writer);

      solv::NonlinearSolverStrategy<NDIM>* nonlinear_solver = NULL;

      if (nonlinear_solver_package == "PETSc-SNES") {

#ifdef HAVE_PETSC
	 nonlinear_solver =
	    new solv::SNES_SAMRAIContext<NDIM>("SNESSolver",
					       input_db->getDatabase("SNESSolver"),
					       bratu_model);
#else
	 TBOX_ERROR("Cannot use PETSc-SNES option because SAMRAI was\n"
		    << "not configured to use it.");
#endif

      } else if (nonlinear_solver_package == "KINSOL") {

#ifdef HAVE_SUNDIALS
	 nonlinear_solver =
	    new solv::KINSOL_SAMRAIContext<NDIM>("KINSOLSolver",
						 input_db->getDatabase("KINSOLSolver"),
						 bratu_model);
#else
	 TBOX_ERROR("Cannot use KINSOL option because SAMRAI was\n"
		    << "not configured to use it.");
#endif

      } else {

	 TBOX_ERROR("Input key `nonlinear_solver_package' == "
		    << nonlinear_solver_package
		    << " found in input file is not recognized.");
      }

      algs::ImplicitIntegrator<NDIM>* imp_integrator = new algs::ImplicitIntegrator<NDIM>(
	 "ImplicitIntegrator",
	 input_db->getDatabase("ImplicitIntegrator"),
	 bratu_model,
	 nonlinear_solver,
	 patch_hierarchy);

      tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
	 new mesh::StandardTagAndInitialize<NDIM>("CellTaggingMethod",
						  bratu_model,
						  input_db->
						  getDatabase("StandardTagAndInitialize"));

      tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

      tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
	 new mesh::LoadBalancer<NDIM>("LoadBalancer",
				      input_db->getDatabase("LoadBalancer"));

      tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
	 new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
					   input_db->getDatabase("GriddingAlgorithm"),
					   error_detector,
					   box_generator,
					   load_balancer);

      /*
       * Tag buffers are passed to the gridding algorithm for buffering
       * tagged cells before new levels are created.
       */
      tbox::Array<int> tag_buffer(gridding_algorithm->getMaxLevels());
      for (int ln = 0; ln < gridding_algorithm->getMaxLevels(); ln++) {
	 tag_buffer[ln] = 0;
      }
      if (main_db->keyExists("tag_buffer")) {
	 tbox::Array<int> input_tags = main_db->getIntegerArray("tag_buffer");
	 if (input_tags.getSize() > 0) {
	    for (int ln0 = 0; ln0 < gridding_algorithm->getMaxLevels(); ln0++) {
	       if (input_tags.getSize() > ln0) {
		  tag_buffer[ln0] = ((input_tags[ln0] > 0) ?
				     input_tags[ln0] : 0);
	       } else {
		  tag_buffer[ln0] = ((input_tags[input_tags.getSize()-1] > 0)
				     ? input_tags[input_tags.getSize()-1] : 0);
	       }
	    }
	 }
      }

      /*
       * After creating all objects and initializing their state, we
       * print the input database and variable database contents to
       * the log file.
       */

      tbox::plog << "\nCheck input data and variables before simulation:" << endl;
      tbox::plog << "Input database..." << endl;
      input_db->printClassData(tbox::plog);
      tbox::plog << "\nVariable database..." << endl;
      hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

      if ( uses_vizamrai ) {
	 viz_data_writer->
	    setFinestLevelToPlot( gridding_algorithm->getMaxLevels()-1);
	 for ( int ln=1; ln<gridding_algorithm->getMaxLevels(); ln++ ) {
	    const hier::IntVector<NDIM> &lratio =
	       gridding_algorithm->getRatioToCoarserLevel(ln);
	    viz_data_writer->setRatioToCoarserLevel(ln, lratio);
	 }
      }

      /*
       * Initialize AMR hierarchy configuration and data on all patches
       * at initial time.
       * Also, write initial state for
       * vizualization.  Note that we also start a timer for the
       * simulation part of the main program.
       */

      tbox::Pointer<tbox::Timer> main_timer = tbox::TimerManager::getManager()->
	 getTimer("apps::main::main");
      tbox::Pointer<tbox::Timer> solve_timer = tbox::TimerManager::getManager()->
	 getTimer("apps::main::solve");

      main_timer->start();

      double sim_time = imp_integrator->getInitialTime();

      if (tbox::RestartManager::getManager()->isFromRestart()) {

	 patch_hierarchy->getFromRestart(gridding_algorithm->getMaxLevels());

	 gridding_algorithm->getTagAndInitializeStrategy()->
	    resetHierarchyConfiguration(patch_hierarchy,
					0,
					patch_hierarchy->getFinestLevelNumber());

      } else {

	 gridding_algorithm->makeCoarsestLevel(patch_hierarchy, sim_time);

	 bool done = false;
	 bool initial_time = true;
	 for (int lnum = 0;
	      gridding_algorithm->levelCanBeRefined(lnum) && !done; lnum++) {
	    gridding_algorithm->makeFinerLevel(patch_hierarchy,
					       sim_time,
					       initial_time,
					       tag_buffer[lnum]);
	    done = !(patch_hierarchy->finerLevelExists(lnum));
	 }

      }

      /*
       * Now that hierarchy is constructed, initialize vector weights.
       */

      bratu_model->setVectorWeights(patch_hierarchy);

      tbox::RestartManager::getManager()->closeRestartFile();

      if ( uses_vizamrai ) {
	 viz_data_writer->writePlotData(patch_hierarchy,
					viz_dump_filename,
					imp_integrator->getIntegratorStep(),
					sim_time);
      }
      if ( uses_visit ) {
	 visit_data_writer->writePlotData(patch_hierarchy,
					  imp_integrator->getIntegratorStep(),
					  sim_time);
      }

      /*
       * Initialize implicit integrator.  Then, loop over a sequence of time
       * steps until the final simulation time is reached, or the maximum
       * number of steps is exceeded.
       */

      imp_integrator->initialize();

      double dt       = imp_integrator->getCurrentDt();
      double end_time = imp_integrator->getFinalTime();

      bool first_step = true;

      while ( (sim_time < end_time) &&
	      imp_integrator->stepsRemaining() ) {

	 int iteration_num = imp_integrator->getIntegratorStep() + 1;
	 tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	 tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
	 tbox::pout << "Simulation time is " << sim_time << endl;
	 tbox::pout << "Time increment is " << dt << endl;
	 tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

	 solve_timer->start();
	 int solver_retcode = imp_integrator->advanceSolution(dt, first_step);
	 solve_timer->stop();

	 solv::SNES_SAMRAIContext<NDIM>* snes_solver =
	    dynamic_cast<solv::SNES_SAMRAIContext<NDIM>*>(nonlinear_solver);
	 if (snes_solver) {
	    int nonlinear_itns = snes_solver->getNonlinearIterationCount();
	    int linear_itns = snes_solver->getTotalLinearIterationCount();
	    tbox::pout << " Nonlinear iterations:  " << nonlinear_itns
		       << " Linear iterations:     " << linear_itns << endl;
	 }

	 bool good_solution = imp_integrator->checkNewSolution(solver_retcode);

	 if ( good_solution ) {
	    sim_time = imp_integrator->updateSolution();

	    tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	    tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
	    tbox::pout << "Simulation time is " << sim_time << endl;
	    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

	    /*
	     * If desired, write Vizamrai plot file.
	     */

	    if ( (iteration_num % viz_dump_interval) == 0 ) {
	       if ( uses_vizamrai ) {
		  viz_data_writer->writePlotData(patch_hierarchy,
						 viz_dump_filename,
						 iteration_num,
						 sim_time);
	       }
	       if ( uses_visit ) {
		  visit_data_writer->writePlotData(patch_hierarchy,
						   iteration_num,
						   sim_time);
	       }
	    }

	    /*
	     * If desired, regrid patch hierarchy and reset vector weights.
	     */

	    if ( (regrid_interval > 0)
		 && ((iteration_num % regrid_interval) == 0) ) {
	       gridding_algorithm->regridAllFinerLevels(patch_hierarchy, 0,
							sim_time,
							tag_buffer);
#if defined(HAVE_PETSC) || defined(HAVE_SUNDIALS)
	       bratu_model->setVectorWeights(patch_hierarchy);
#endif
	       first_step = true;
	    } else {
	       first_step = false;
	    }
	 } else {

	    tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	    tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
	    tbox::pout << "Failed to advance solution to " << sim_time << endl;
	    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
	    break;

	 }

	 dt = imp_integrator->getNextDt(good_solution, solver_retcode);

      }

      /*
       * At conclusion of simulation, stop timer and deallocate objects.
       */
      main_timer->stop();

      tbox::TimerManager::getManager()->print(tbox::pout);

      input_db.setNull();
      patch_hierarchy.setNull();
      grid_geometry.setNull();

      box_generator.setNull();
      error_detector.setNull();
      load_balancer.setNull();
      gridding_algorithm.setNull();

      if (imp_integrator) 
	 delete imp_integrator;
      if (nonlinear_solver) 
	 delete nonlinear_solver;
      if (bratu_model) 
	 delete bratu_model;

      if ( uses_visit ) {
	 visit_data_writer.setNull();
      }

      if ( uses_vizamrai ) {
	 viz_data_writer.setNull();
      }

#endif
#endif

#ifdef TESTING
      /*
       * Currently not checking on solution accuracy, so if a run
       * makes it here, it passes.  We need to find a way to check
       * the solution accuracy.
       */
      tbox::pout << "\nPASSED:  nonlinear" << endl;
#endif

   }

   tbox::SAMRAIManager::shutdown();

#ifdef HAVE_PETSC
   PetscFinalize();
#endif
   tbox::SAMRAI_MPI::finalize();

   return(0);
}
