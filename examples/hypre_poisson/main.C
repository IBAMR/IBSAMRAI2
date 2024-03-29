/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/hypre_poisson/main.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Main program for Hypre Poisson example
 */

#include "SAMRAI_config.h"

#include <string>
using namespace std;

#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "GriddingAlgorithm.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAIManager.h"
#include "StandardTagAndInitialize.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "VisItDataWriter.h"

#include "HyprePoisson.h"

using namespace SAMRAI;

/*
************************************************************************
*                                                                      *
* This is the driver program to demonstrate                            *
* how to use the Hypre Poisson solver.                                 *
*                                                                      *
* We set up the simple problem                                         *
*          u + div(grad(u)) = sin(x)*sin(y)                            *
* in the domain [0:1]x[0:1], with u=0 on the                           *
* boundary.                                                            *
*                                                                      *
* HyprePoisson is the primary object used to                           *
* set up and solve the system.  It maintains                           *
* the data for the computed solution u, the                            *
* exact solution, and the right hand side.                             *
*                                                                      *
* The hierarchy created to solve this problem                          *
* has only one level.  (The Hypre Poisson solver                       *
* is a single-level solver.)                                           *
*                                                                      *
*************************************************************************
*/

int main( int argc, char *argv[] )
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      bool converged = true;


#if !defined(HAVE_HYPRE)
      tbox::pout << "This example requires the package HYPRE"
		 << "\nto work properly.  SAMRAI was not configured"
		 << "\nwith this package." 
		 << endl;
#else

      /*
       * Process command line arguments.  For each run, the input 
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
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
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");
      string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);


      /*
       * Start logging.
       */
      const string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
      if (log_all_nodes) {
	 tbox::PIO::logAllNodes(log_file_name);
      } else {
	 tbox::PIO::logOnlyNodeZero(log_file_name);
      }


      /*
       * Create major algorithm and data objects which comprise application.
       * Each object will be initialized either from input data or restart
       * files, or a combination of both.  Refer to each class constructor
       * for details.  For more information on the composition of objects
       * for this application, see comments at top of file.
       */

      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
	 new geom::CartesianGridGeometry<NDIM>(base_name+"CartesianGeometry",
					       input_db->getDatabase("CartesianGeometry"));
      tbox::plog << "Cartesian Geometry:" << endl;
      grid_geometry->printClassData(tbox::plog);

      tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
	 new hier::PatchHierarchy<NDIM>(base_name+"::PatchHierarchy",
					grid_geometry);

      /*
       * The HyprePoisson object is the main user object specific to the
       * problem being solved.  It provides the implementations for setting
       * up the grid and plotting data.  It also wraps up the solve
       * process that includes making the initial guess, specifying the
       * boundary conditions and call the solver.
       */
      HyprePoisson hypre_poisson(base_name+"::HyprePoisson",
				 input_db->isDatabase("HyprePoisson") ?
				 input_db->getDatabase("HyprePoisson") :
				 tbox::Pointer<tbox::Database>(NULL) );

      /*
       * Create the tag-and-initializer, box-generator and load-balancer
       * object references required by the gridding_algorithm object.
       */
      tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > tag_and_initializer =
	 new mesh::StandardTagAndInitialize<NDIM>(
	    "CellTaggingMethod" ,
	    tbox::Pointer<mesh::StandardTagAndInitStrategy<NDIM> >(&hypre_poisson,false) ,
	    input_db->getDatabase("StandardTagAndInitialize")
	    );
      tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
	 new mesh::BergerRigoutsos<NDIM>();
      tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
	 new mesh::LoadBalancer<NDIM>("load balancer",
				      input_db->getDatabase("LoadBalancer"));

      /*
       * Create the gridding algorithm used to generate the SAMR grid
       * and create the grid.
       */
      tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm;
      gridding_algorithm =
	 new mesh::GriddingAlgorithm<NDIM>("Gridding Algorithm",
					   input_db->getDatabase("GriddingAlgorithm"),
					   tag_and_initializer,
					   box_generator,
					   load_balancer);
      tbox::plog << "Gridding algorithm:" << endl;
      gridding_algorithm->printClassData(tbox::plog);

      /*
       * Make the coarsest patch level where we will be solving.
       */
      gridding_algorithm->makeCoarsestLevel(patch_hierarchy,0.0);



      /*
       * Set up the plotter for the hierarchy just created.
       * The FACPoisson object handles the data and has the
       * function setupExternalPlotter to register its data
       * with the plotter.
       */
      tbox::Array<string> vis_writer(1);
      vis_writer[0] = "Vizamrai";
      if (main_db->keyExists("vis_writer")) {
	 vis_writer = main_db->getStringArray("vis_writer");
      }
      bool use_vizam = false;
      bool use_visit = false;
      for (int i = 0; i < vis_writer.getSize(); i++) {
	 if (vis_writer[i] == "Vizamrai") use_vizam = true;
	 if (vis_writer[i] == "VisIt") use_visit = true;
      }
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > vizam_writer;
      string vis_filename
	 = main_db->getStringWithDefault("vis_filename", base_name);
      if ( use_vizam ) {
	 vizam_writer = new appu::CartesianVizamraiDataWriter<NDIM>("Vizamrai Writer");
	 hypre_poisson.setupExternalPlotter(*vizam_writer);
      }
#ifdef HAVE_HDF5
      tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer;
      if ( use_visit ) {
	 visit_writer = new appu::VisItDataWriter<NDIM>("Visit Writer", vis_filename+".visit");
	 hypre_poisson.setupExternalPlotter(*visit_writer);
      }
#endif



#if 0
      /*
       * Set up the plotter for the hierarchy just created.
       * The HyprePoisson object handles the data and has the
       * function setupExternalPlotter to register its data
       * with the plotter.
       */
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer;
      string viz_filename
	 = main_db->getStringWithDefault("viz_filename", base_name);
      viz_data_writer = new appu::CartesianVizamraiDataWriter<NDIM>("Viz Writer");
      hypre_poisson.setupExternalPlotter(*viz_data_writer);
      viz_data_writer->printClassData(tbox::plog);
#endif

      /*
       * After creating all objects and initializing their state,
       * we print the input database and variable database contents
       * to the log file.
       */
      tbox::plog << "\nCheck input data and variables before simulation:" << endl;
      tbox::plog << "Input database..." << endl;
      input_db->printClassData(tbox::plog);

      /*
       * Solve.
       */
      converged = hypre_poisson.solvePoisson();

      /*
       * Plot.
       */
      if ( use_vizam ) {
	 vizam_writer->writePlotData( patch_hierarchy ,
				      vis_filename + "_vizam" );
      }
#ifdef HAVE_HDF5
      if ( use_visit ) {
	 visit_writer->writePlotData( patch_hierarchy, 0 );
      }
#endif


      /*
       * Deallocate objects when done.
       */

      tbox::TimerManager::getManager()->print(tbox::plog);

#endif

#ifdef TESTING
      if ( converged ) {
	 tbox::pout << "\nPASSED:  hypre" << endl;
      }
      else {
	 TBOX_WARNING("Hypre test did not converge.");
      }
#endif
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(0); 
}
