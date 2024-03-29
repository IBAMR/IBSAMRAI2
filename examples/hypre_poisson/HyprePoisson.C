/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/hypre_poisson/HyprePoisson.C $
 * Package:     SAMRAI application
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2043 $
 * Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
 * Description: Numerical routines for example Hypre Poisson solver
 */

#include "HyprePoisson.h"

#if defined(HAVE_HYPRE)

#include <iostream>

#include "IntVector.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "HierarchyCellDataOpsReal.h"
#include "SideData.h"
#include "PoissonSpecifications.h"
#include "tbox/Utilities.h"
#include "Variable.h"
#include "VariableDatabase.h"



extern "C" {
#if NDIM == 2
   void setexactandrhs_( const int &ifirst0,
                         const int &ilast0,
                         const int &ifirst1,
                         const int &ilast1,
                         double *exact,
                         double *rhs,
                         const double *dx,
                         const double *xlower );
#elif NDIM == 3
   void setexactandrhs_( const int &ifirst0,
                         const int &ilast0,
                         const int &ifirst1,
                         const int &ilast1,
                         const int &ifirst2,
                         const int &ilast2,
                         double *exact,
                         double *rhs,
                         const double *dx,
                         const double *xlower );
#endif
}

namespace SAMRAI {

/*
*************************************************************************
* Constructor creates a unique context for the object and register      *
* all its internal variables with the variable database.                *
*************************************************************************
*/
HyprePoisson::HyprePoisson(
  const string &object_name,
  tbox::Pointer<tbox::Database> database )
: d_object_name(object_name),
  d_hierarchy(NULL),
  d_poisson_hypre( object_name+"::poisson_hypre",
                   ( !database.isNull() &&
                     database->isDatabase("CellPoissonHypreSolver") ) ?
                   database->getDatabase("CellPoissonHypreSolver"):
                   tbox::Pointer<tbox::Database>(NULL) ),
  d_bc_coefs( object_name+"::bc_coefs",
              ( !database.isNull() &&
                database->isDatabase("bc_coefs") ) ?
              database->getDatabase("bc_coefs"):
              tbox::Pointer<tbox::Database>(NULL) ),
  d_context()
{

  hier::VariableDatabase<NDIM> *vdb = hier::VariableDatabase<NDIM>::getDatabase();

  /*
   * Get a unique context for this object.
   */
  d_context = vdb->getContext( d_object_name + ":Context" );

  /*
   * Register variables with hier::VariableDatabase<NDIM>
   * and get the descriptor indices for those variables.
   */
  tbox::Pointer<pdat::CellVariable<NDIM,double> > comp_soln =
     new pdat::CellVariable<NDIM,double>(object_name+":computed solution",1);
  d_comp_soln_id
     = vdb->registerVariableAndContext(
     comp_soln ,
     d_context ,
     hier::IntVector<NDIM>(1) /* ghost cell width is 1 for stencil widths */ );
  tbox::Pointer<pdat::CellVariable<NDIM,double> > exact_solution =
     new pdat::CellVariable<NDIM,double>(object_name+":exact solution");
  d_exact_id
     = vdb->registerVariableAndContext(
     exact_solution ,
     d_context,
     hier::IntVector<NDIM>(1) /* ghost cell width is 1 in case needed */ );
  tbox::Pointer<pdat::CellVariable<NDIM,double> > rhs_variable =
     new pdat::CellVariable<NDIM,double>(object_name+":linear system right hand side");
  d_rhs_id
     = vdb->registerVariableAndContext(
     rhs_variable ,
     d_context,
     hier::IntVector<NDIM>(0) /* ghost cell width is 0 */ );

  return;
}



/*
*************************************************************************
* Destructor does nothing interesting                                   *
*************************************************************************
*/
HyprePoisson::~HyprePoisson() 
{
} 





/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void HyprePoisson::initializeLevelData (
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy , 
   const int level_number ,
   const double init_data_time ,
   const bool can_be_refined ,
   const bool initial_time ,
   const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level ,
   const bool allocate_data )
{

   tbox::Pointer< hier::PatchHierarchy<NDIM> > patch_hierarchy = hierarchy;
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom =
      patch_hierarchy->getGridGeometry();

   tbox::Pointer<hier::PatchLevel<NDIM> > level
      = hierarchy->getPatchLevel(level_number);

   /*
    * If required, allocate all patch data on the level.
    */
   if ( allocate_data ) {
      level->allocatePatchData( d_comp_soln_id );
      level->allocatePatchData( d_rhs_id );
      level->allocatePatchData( d_exact_id );
   }

   /*
    * Initialize data in all patches in the level.
    */
   hier::PatchLevel<NDIM>::Iterator pi(*level);
   for ( pi.initialize(*level); pi; pi++ ) {

      const int pn = *pi;
      tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(pn);
      if( patch.isNull() ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot find patch.  Null patch pointer.");
      }
      hier::Box<NDIM> pbox = patch->getBox();
      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
         = patch->getPatchGeometry();

      tbox::Pointer<pdat::CellData<NDIM,double> > exact_data
         = patch->getPatchData(d_exact_id);
      tbox::Pointer<pdat::CellData<NDIM,double> > rhs_data
         = patch->getPatchData(d_rhs_id);

      /*
	Set source function and exact solution.
       */
#if NDIM == 2
      setexactandrhs_( pbox.lower()[0],
                       pbox.upper()[0],
                       pbox.lower()[1],
                       pbox.upper()[1],
                       exact_data->getPointer(),
                       rhs_data->getPointer(),
                       grid_geom->getDx(),
                       patch_geom->getXLower() );
#elif NDIM == 3
      setexactandrhs_( pbox.lower()[0],
                       pbox.upper()[0],
                       pbox.lower()[1],
                       pbox.upper()[1],
                       pbox.lower()[2],
                       pbox.upper()[2],
                       exact_data->getPointer(),
                       rhs_data->getPointer(),
                       grid_geom->getDx(),
                       patch_geom->getXLower() );
#endif

   }	// End patch loop.
}




/*
*************************************************************************
* Reset the hierarchy-dependent internal information.                   *
*************************************************************************
*/
void HyprePoisson::resetHierarchyConfiguration (
   tbox::Pointer<hier::BasePatchHierarchy<NDIM> > new_hierarchy ,
   int coarsest_level ,
   int finest_level )
{
   d_hierarchy = new_hierarchy;
   return;
}





/*
*************************************************************************
* Solve the Poisson problem.                                            *
*************************************************************************
*/
bool HyprePoisson::solvePoisson()
{

  if ( d_hierarchy.isNull() ) {
    TBOX_ERROR("Cannot solve using an uninitialized object.\n");
  }

  const int level_number=0;

  /*
   * Fill in the initial guess and Dirichlet boundary condition data.
   * For this example, we want u=0 on all boundaries.
   * The easiest way to do this is to just write 0 everywhere,
   * simultaneous setting the boundary values and initial guess.
   */
  tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
  hier::PatchLevel<NDIM>::Iterator ip(level);
  for ( ; ip; ip++ ) {
     tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(*ip);
     tbox::Pointer<pdat::CellData<NDIM,double> > data = patch->getPatchData(d_comp_soln_id);
     data->fill(0.0);
  }
  // d_poisson_hypre.setBoundaries( "Dirichlet" );
  d_poisson_hypre.setPhysicalBcCoefObject(&d_bc_coefs);

  /*
   * Set up HYPRE solver object.
   * The problem specification is set using the
   * CellPoissonSpecifications object then passed to the solver
   * for setting the coefficients.
   */
  d_poisson_hypre.initializeSolverState( d_hierarchy,
					 level_number);
  solv::PoissonSpecifications sps( "Hypre Poisson solver" );
  sps.setCZero();
  sps.setDConstant(1.0);
  d_poisson_hypre.setMatrixCoefficients( sps );

  /*
   * Solve the system.
   */
  tbox::pout << "solving..." << endl;
  int solver_ret;
  solver_ret = d_poisson_hypre.solveSystem( d_comp_soln_id ,
                                            d_rhs_id );
  /*
   * Present data on the solve.
   */
  tbox::pout << "\t" << (solver_ret?"":"NOT ") << "converged " << "\n"
       << "	iterations: " << d_poisson_hypre.getNumberOfIterations() << "\n"
       << "	residual: " << d_poisson_hypre.getRelativeResidualNorm() << "\n"
       << flush;


  /*
   * Deallocate state.
   */
  d_poisson_hypre.deallocateSolverState();

  /*
   * Return whether solver converged.
   */
  return solver_ret ? true : false;
}






/*
*************************************************************************
* Set up external plotter to plot internal data from this class.        *
* Tell the plotter about the refinement ratios.  Register variables     *
* appropriate for plotting.                                             *
*************************************************************************
*/
int HyprePoisson::setupExternalPlotter(
  appu::CartesianVizamraiDataWriter<NDIM> &viz_writer
) const {

   /*
    * This must be done once.
    */
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                 << " HyprePoisson::registerVariablesWithPlotter\n"
                 << "The hierarchy must be built before calling\n"
                 << "this function.\n");
   }
   int ln;
   for ( ln=1; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {

     tbox::Pointer<hier::PatchLevel<NDIM> > level = 
	d_hierarchy->getPatchLevel(ln);
      const hier::IntVector<NDIM> &lratio = level->getRatioToCoarserLevel();
      viz_writer.setRatioToCoarserLevel(ln, lratio);
   }
   /*
    * Register variables with plotter.
    */
   viz_writer.registerPlotScalar("Computed solution", d_comp_soln_id);
   viz_writer.registerDerivedPlotScalar("Error" ,
                                        (appu::VisDerivedDataStrategy<NDIM> *)this );
   viz_writer.registerPlotScalar("Exact solution", d_exact_id);
   viz_writer.registerPlotScalar("Poisson source", d_rhs_id);
   viz_writer.registerDerivedPlotScalar("Patch level number" ,
                                        (appu::VisDerivedDataStrategy<NDIM> *)this );

   return 0;
}




#ifdef HAVE_HDF5
/*
*************************************************************************
* Set up VisIt to plot internal data from this class.                   *
* Tell the plotter about the refinement ratios.  Register variables     *
* appropriate for plotting.                                             *
*************************************************************************
*/
int HyprePoisson::setupExternalPlotter(
  appu::VisItDataWriter<NDIM> &viz_writer
) const {

   /*
    * This must be done once.
    */
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                 << " HyprePoisson::registerVariablesWithPlotter\n"
                 << "The hierarchy must be built before calling\n"
                 << "this function.\n");
   }
   /*
    * Register variables with plotter.
    */
   viz_writer.registerPlotQuantity("Computed solution", 
                                   "SCALAR", 
                                   d_comp_soln_id);
   viz_writer.registerDerivedPlotQuantity("Error" ,
                                          "SCALAR", 
                                          (appu::VisDerivedDataStrategy<NDIM> *)this );
   viz_writer.registerPlotQuantity("Exact solution", 
                                   "SCALAR", 
                                   d_exact_id);
   viz_writer.registerPlotQuantity("Poisson source", 
                                   "SCALAR",
                                   d_rhs_id);
   viz_writer.registerDerivedPlotQuantity("Patch level number" ,
                                          "SCALAR",
                                          (appu::VisDerivedDataStrategy<NDIM> *)this );

   return 0;
}
#endif



/*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
bool HyprePoisson::packDerivedDataIntoDoubleBuffer(
   double *buffer ,
   const hier::Patch<NDIM> &patch ,
   const hier::Box<NDIM> &region ,
   const string &variable_name ,
   int depth_id ) const
{
   pdat::CellData<NDIM,double>::Iterator icell( patch.getBox() );

   if ( variable_name == "Error" ) {
      tbox::Pointer<pdat::CellData<NDIM,double> > current_solution_
        = patch.getPatchData(d_comp_soln_id);
      tbox::Pointer<pdat::CellData<NDIM,double> > exact_solution_
        = patch.getPatchData(d_exact_id);
      pdat::CellData<NDIM,double> &current_solution = *current_solution_;
      pdat::CellData<NDIM,double> &exact_solution = *exact_solution_;
      for ( ; icell; icell++ ) {
         double diff = ( current_solution(*icell) - exact_solution(*icell) );
	 *buffer = diff;
	 buffer += 1;
      }
   }
   else if ( variable_name == "Patch level number" ) {
      double pln = patch.getPatchLevelNumber();
      for ( ; icell; icell++ ) {
	 *buffer = pln;
	 buffer += 1;
      }
   }
   else {
      // Did not register this name.
      TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
                 << "HyprePoissonX::writeDerivedDataToStream");

   }
   return true;
}

}

#endif
