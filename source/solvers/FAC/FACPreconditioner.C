//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/FAC/FACPreconditioner.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: FAC algorithm for solving linear equations on a hierarchy
//


#ifndef included_solv_FACPreconditioner_C
#define included_solv_FACPreconditioner_C

#include "HierarchyDataOpsManager.h"
#include "FACPreconditioner.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


#ifdef DEBUG_NO_INLINE
#include "FACPreconditioner.I"
#endif

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Constructor checks arguments and sets uninitialized solver state.     *
*                                                                       *
*************************************************************************
*/

template<int DIM>  FACPreconditioner<DIM>::FACPreconditioner (
   const std::string& name,
   FACOperatorStrategy<DIM> &user_ops,
   tbox::Pointer<tbox::Database> database
) :
   d_object_name(name) ,
   d_fac_operator(user_ops) ,
   d_coarsest_ln(0) ,
   d_finest_ln(0) ,
   d_residual_vector() ,
   d_tmp_residual() ,
   d_error_vector() ,
   d_tmp_error() ,
   d_max_iterations(0) ,
   d_residual_tolerance(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_relative_residual_tolerance(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_presmoothing_sweeps(1) ,
   d_postsmoothing_sweeps(1) ,
   d_algorithm_choice("default") ,
   d_number_iterations(0) ,
   d_residual_norm(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_rhs_norm(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_convergence_factor() ,
   d_avg_convergence_factor(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_net_convergence_factor(tbox::MathUtilities<double>::getSignalingNaN()) ,
   d_do_log(false) ,
   d_controlled_level_ops()
{

   t_solve_system = tbox::TimerManager::getManager()->
      getTimer("solv::FACPreconditioner::solveSystem()_fac_cycling");

   /*
    * Initialize object with data read from input database.
    */
   if ( database ) {
      getFromInput(database);
   }

   return;
}

/*
*************************************************************************
*                                                                       *
* Destructor for FACPreconditioner<DIM>.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>  FACPreconditioner<DIM>::~FACPreconditioner()
{
   deallocateSolverState();
   return;
}


/*
********************************************************************
* Set state from database                                          *
********************************************************************
*/

template<int DIM> void FACPreconditioner<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> database )
{
   if ( database ) {
      if ( database->isBool("enable_logging") ) {
         bool logging = database->getBool("enable_logging");
         enableLogging(logging);
      }
      if ( database->isInteger("max_cycles") ) {
         int max_cycles = database->getInteger("max_cycles");
         setMaxCycles(max_cycles);
      }
      if ( database->isDouble("residual_tol") ) {
         double residual_tol = database->getDouble("residual_tol");
         setResidualTolerance(residual_tol, d_relative_residual_tolerance);
      }
      if ( database->isDouble("relative_residual_tol") ) {
         double relative_residual_tol = database->getDouble("relative_residual_tol");
         setResidualTolerance(d_residual_tolerance,relative_residual_tol);
      }
      if ( database->isInteger("num_pre_sweeps") ) {
         int num_pre_sweeps = database->getInteger("num_pre_sweeps");
         setPresmoothingSweeps(num_pre_sweeps);
      }
      if ( database->isInteger("num_post_sweeps") ) {
         int num_post_sweeps = database->getInteger("num_post_sweeps");
         setPostsmoothingSweeps(num_post_sweeps);
      }
   }
   return;
}



/*
*************************************************************************
*                                                                       *
* Functions for setting up and deallocating solver state.               *
*                                                                       *
*************************************************************************
*/
template<int DIM> void FACPreconditioner<DIM>::deallocateSolverState( )
{
   /*
    * Delete hierarchy-dependent state data.
    */

   if ( ! d_patch_hierarchy.isNull() ) {

      d_coarsest_ln = d_finest_ln = -1;
      d_patch_hierarchy.setNull();

      if ( d_error_vector ) {
         d_error_vector->resetLevels(d_error_vector->getCoarsestLevelNumber(), std::min(d_error_vector->getFinestLevelNumber(),d_error_vector->getPatchHierarchy()->getFinestLevelNumber()));
         d_error_vector->freeVectorComponents();
         d_error_vector.setNull();
      }
      if ( d_tmp_error ) {
         d_tmp_error->resetLevels(d_tmp_error->getCoarsestLevelNumber(), std::min(d_tmp_error->getFinestLevelNumber(),d_tmp_error->getPatchHierarchy()->getFinestLevelNumber()));
         d_tmp_error->freeVectorComponents();
         d_tmp_error.setNull();
      }
      if ( d_residual_vector ) {
         d_residual_vector->resetLevels(d_residual_vector->getCoarsestLevelNumber(), std::min(d_residual_vector->getFinestLevelNumber(),d_residual_vector->getPatchHierarchy()->getFinestLevelNumber()));
         d_residual_vector->freeVectorComponents();
         d_residual_vector.setNull();
      }
      if ( d_tmp_residual ) {
         d_tmp_residual->resetLevels(d_tmp_residual->getCoarsestLevelNumber(), std::min(d_tmp_residual->getFinestLevelNumber(),d_tmp_residual->getPatchHierarchy()->getFinestLevelNumber()));
         d_tmp_residual->freeVectorComponents();
         d_tmp_residual.setNull();
      }

      d_controlled_level_ops.setNull();
      d_fac_operator.deallocateOperatorState();
   }

   return;
}



template<int DIM> void FACPreconditioner<DIM>::initializeSolverState(
  const SAMRAIVectorReal<DIM,double> &solution ,
  const SAMRAIVectorReal<DIM,double> &rhs )
{
   /*
    * First get rid of current data.
    */
   deallocateSolverState();
   /*
    * Set hierarchy and levels to solve.
    */
   d_patch_hierarchy = solution.getPatchHierarchy();
   d_coarsest_ln = solution.getCoarsestLevelNumber();
   d_finest_ln = solution.getFinestLevelNumber();
   /*
    * Set the solution-vector-dependent scratch space.
    */
   d_error_vector = solution.cloneVector(d_object_name+"::error");
   d_error_vector->allocateVectorData();
   if ( d_algorithm_choice == "mccormick-s4.3" ) {
      d_tmp_error= solution.cloneVector(d_object_name+"::temporary_error");
      d_tmp_error->allocateVectorData();
   }
   d_residual_vector = rhs.cloneVector(d_object_name+"::residual");
   d_residual_vector->allocateVectorData();
   d_tmp_residual = rhs.cloneVector( d_object_name+"::FAC coarser residual" );
   d_tmp_residual->allocateVectorData();
   /*
    * Set the controlled level operators, which depend on the number
    * of components in the solution vector.
    */
   math::HierarchyDataOpsManager<DIM> *ops_manager
      = math::HierarchyDataOpsManager<DIM>::getManager();
   int num_components = solution.getNumberOfComponents();
   d_controlled_level_ops.resizeArray(num_components);
   for ( int i=0; i<num_components; ++i ) {
      tbox::Pointer< hier::Variable<DIM> > variable = solution.getComponentVariable(i);
      d_controlled_level_ops[i]
         = ops_manager->getOperationsDouble( variable ,
                                             d_patch_hierarchy ,
                                             true );
      /*
       * Note: the variable used above is only for the purpose of determining
       * the variable alignment on the grid.  It is not specific to any
       * instance.
       */
   }
   /*
    * Error checking.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   if(d_patch_hierarchy != rhs.getPatchHierarchy()) {
      TBOX_ERROR(d_object_name << ": vectors must have the same hierarchy.\n");
   }
   if(d_coarsest_ln < 0) {
      TBOX_ERROR(d_object_name << ": coarsest level must not be negative.\n");
   }
   if(d_coarsest_ln > d_finest_ln) {
      TBOX_ERROR(d_object_name << ": coarsest level must be <= finest"
		 << "level.\n");
   }
#endif
   for (int ln = d_coarsest_ln; ln <= d_finest_ln; ln++) {
      if (d_patch_hierarchy->getPatchLevel(ln).isNull()) {
         TBOX_ERROR("FACPreconditioner<DIM>::FACPreconditioner<DIM> error ..."
                    << "\n   object name = " << d_object_name
                    << "\n   hierarchy level " << ln 
                    << " does not exist" << std::endl);
      }
   }
   d_fac_operator.initializeOperatorState( solution, rhs );
}



template<int DIM> bool FACPreconditioner<DIM>::checkVectorStateCompatibility(
  const SAMRAIVectorReal<DIM,double> &solution ,
  const SAMRAIVectorReal<DIM,double> &rhs ) const
{
   /*
    * It is an error when the state is not initialized.
    */
   if ( ! d_patch_hierarchy ) {
      TBOX_ERROR(d_object_name << ": cannot check vector-state\n"
		 << "compatibility when the state is uninitialized.\n");
   }
   bool rvalue = true;
   const SAMRAIVectorReal<DIM,double> &error = *d_error_vector;
   if (  solution.getPatchHierarchy() != d_patch_hierarchy
      || rhs.getPatchHierarchy() != d_patch_hierarchy ) {
      rvalue = false;
   }
   if ( rvalue == true ) {
      if (  solution.getCoarsestLevelNumber() != d_coarsest_ln
         || rhs.getCoarsestLevelNumber() != d_coarsest_ln
         || solution.getFinestLevelNumber() != d_finest_ln
         || rhs.getFinestLevelNumber() != d_finest_ln ) {
         rvalue = false;
      }
   }
   if ( rvalue == true ) {
      const int ncomp = error.getNumberOfComponents();
      if (  solution.getNumberOfComponents() != ncomp
         || rhs.getNumberOfComponents() != ncomp ) {
         rvalue = false;
      }
   }
   return rvalue;
}



/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether itertion converged.        *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool FACPreconditioner<DIM>::solveSystem(
   SAMRAIVectorReal<DIM,double> &u,
   SAMRAIVectorReal<DIM,double> &f)
{

   d_residual_norm = tbox::MathUtilities<double>::getSignalingNaN();
   d_avg_convergence_factor = tbox::MathUtilities<double>::getSignalingNaN();

   /*
    * Set the solution-vector-dependent data if not preset.
    */
   bool clear_hierarchy_configuration_when_done = false;
   if ( d_patch_hierarchy.isNull() ) {
      clear_hierarchy_configuration_when_done = true;
      initializeSolverState( u ,
                             f );
   }
   else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !checkVectorStateCompatibility(u,f) ) {
	 TBOX_ERROR(d_object_name << ": Incompatible vectors for\n"
		    << "current state in solveSystem.\n");
      }
#endif
   }

   t_solve_system->start();

   d_error_vector->setToScalar(0.0,false);
   if ( d_tmp_error ) {
     d_tmp_error->setToScalar(0.0,false);
   }
   d_residual_vector->setToScalar(0.0);

   const double initial_residual_norm = d_residual_norm
     = computeFullCompositeResidual( *d_residual_vector ,
                                     u ,
                                     f );
   /*
    * Above step has the side effect of filling the residual
    * vector d_residual_vector.
    */

   double effective_residual_tolerance = d_residual_tolerance;
   if ( d_relative_residual_tolerance >= 0 ) {
      double tmp = d_fac_operator.computeResidualNorm( f ,
                                                       d_finest_ln ,
                                                       d_coarsest_ln );
      tmp *= d_relative_residual_tolerance;
      if ( effective_residual_tolerance < tmp ) effective_residual_tolerance = tmp;
   }


   if ( d_convergence_factor.getSize() < d_max_iterations )
      d_convergence_factor.resizeArray(d_max_iterations);
   d_number_iterations = 0; 
   /*
    * Use a do loop instead of a while loop until convergence.
    * It is important to go through the loop at least once
    * because the residual norm can be less than the tolerance
    * when the solution is 0 and the rhs is less than the tolerance.
    */
   do {

      /*
       * In zeroing the error vector, also zero out the ghost values.
       * This gives the FAC operator an oportunity to bypass the
       * ghost filling if it decides that the ghost values do not
       * change.
       */
      d_error_vector->setToScalar(0.0,false);
      /*
       * Both the recursive and non-recursive fac cycling functions
       * were coded to find the problem due presmoothing.  Both give
       * the same results, but the problem due to presmoothing still
       * exists.  BTNG.
       */
      if ( d_algorithm_choice == "default" ) {
        facCycle_Recursive( *d_error_vector ,
                            *d_residual_vector ,
                            u ,
                            d_finest_ln ,
                            d_coarsest_ln ,
                            d_finest_ln );
      }
      else if ( d_algorithm_choice == "mccormick-s4.3" ) {
        facCycle_McCormick( *d_error_vector ,
                            *d_residual_vector ,
                            u ,
                            d_finest_ln ,
                            d_coarsest_ln ,
                            d_finest_ln );
      }
      else if ( d_algorithm_choice == "pernice" ) {
        facCycle( *d_error_vector ,
                  *d_residual_vector ,
                  u ,
                  d_finest_ln ,
                  d_coarsest_ln );
      }

      int i, num_components = d_error_vector->getNumberOfComponents();
      /*
       * u += e
       */
      for ( i=0; i<num_components; ++i ) {
        int soln_id = u.getComponentDescriptorIndex(i);
        int err_id = d_error_vector->getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels(d_coarsest_ln,d_finest_ln);
        d_controlled_level_ops[i]->add (soln_id,
                                        soln_id,
                                        err_id );
      }


      /*
       * Synchronize solution across levels by coarsening the
       * more accurate fine-level solutions.
       */
      for ( int ln=d_finest_ln-1; ln>=d_coarsest_ln; --ln ) {
         d_fac_operator.restrictSolution( u ,
                                          u ,
                                          ln );
      }

      /*
       * Compute convergence factor and new residual norm.
       * 1. Temporarily save pre-cycle residual norm in
       *    convergence factor stack.
       * 2. Compute post-cycle residual norm.
       * 3. Set convergence factor to ratio of post-cycle to pre-cycle
       *    residual norms.
       */
      d_convergence_factor[d_number_iterations] = d_residual_norm;
      d_residual_norm = computeFullCompositeResidual( *d_residual_vector ,
                                                      u ,
                                                      f );

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
      if ( d_convergence_factor[d_number_iterations] != 0 ) {
         d_convergence_factor[d_number_iterations]
            = d_residual_norm/d_convergence_factor[d_number_iterations];
      }
      else {
         d_convergence_factor[d_number_iterations] = 0;
      }

      /*
       * Increment the iteration counter.
       * The rest of this block expects it to have the incremented value.
       * In particular, d_fac_operator.postprocessOneCycle does.
       */
      ++d_number_iterations;

      /*
       * Compute the convergence factors because they may be accessed
       * from the operator's postprocessOneCycle function.
       */
      d_net_convergence_factor = d_residual_norm/(initial_residual_norm+1e-20);
      d_avg_convergence_factor = pow( d_net_convergence_factor ,
                                      1.0/d_number_iterations );


      d_fac_operator.postprocessOneCycle( d_number_iterations-1 ,
                                          u ,
                                          *d_residual_vector );

   } while ( (d_residual_norm > effective_residual_tolerance) 
             && (d_number_iterations < d_max_iterations) );

   t_solve_system->stop();

   if ( clear_hierarchy_configuration_when_done ) {
      deallocateSolverState();
   }

   return(d_residual_norm < effective_residual_tolerance);

}


/*
*************************************************************************
*                                                                       *
* Recursive version of FAC cycling.                                     *
*                                                                       *
*************************************************************************
*/
template<int DIM> void FACPreconditioner<DIM>::facCycle_Recursive(
   SAMRAIVectorReal<DIM,double> &e,
   SAMRAIVectorReal<DIM,double> &r,
   SAMRAIVectorReal<DIM,double> &u,
   int lmax,
   int lmin,
   int ln )
{

   /*
    * The steps 1-4 in this function are similar, but not exactly
    * equivalent to McCormick's steps in his description of the
    * V-cycle form of FAC.
    *
    * 1. If on the coarsest level, solve it.
    * 2. Presmoothing.
    * 3. Recurse to next lower level.
    * 4. Postsmoothing.
    */

   /* Step 1. */
   if ( ln == lmin ) {
     d_fac_operator.solveCoarsestLevel( e ,
                                        r ,
                                        ln );
   }
   else {

     /* Step 2. */
     d_fac_operator.smoothError( e ,
                                 r ,
                                 ln ,
                                 d_presmoothing_sweeps );
     /* Step 3. */
     d_fac_operator.computeCompositeResidualOnLevel( *d_tmp_residual ,
                                                     e ,
                                                     r ,
                                                     ln ,
                                                     true );
     d_fac_operator.restrictResidual( *d_tmp_residual ,
                                      r ,
                                      ln-1 );
     facCycle_Recursive( e ,
                         r ,
                         u ,
                         lmax ,
                         lmin ,
                         ln-1 );
     d_fac_operator.prolongErrorAndCorrect( e ,
                                            e ,
                                            ln );
     /* Step 4. */
     d_fac_operator.smoothError( e ,
                                 r ,
                                 ln,
                                 d_postsmoothing_sweeps );
   }


   return;
}

/*
*************************************************************************
*                                                                       *
* Implementation of Steve McCormick's V-cycle form of FAC               *
* as described in Section 4.3 of "Multilevel Adaptive Methods for       *
* Partial Differential Equations".                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void FACPreconditioner<DIM>::facCycle_McCormick(
   SAMRAIVectorReal<DIM,double> &e,
   SAMRAIVectorReal<DIM,double> &r,
   SAMRAIVectorReal<DIM,double> &u,
   int lmax,
   int lmin,
   int ln )
{

   /*
    * The steps 1-4 in this function correspond to McCormick's steps
    * those in his description.
    *
    * 1. If on the coarsest level, solve it.
    * 2. Presmoothing.
    * 3. Recurse to next lower level.
    * 4. Postsmoothing.
    */

   /*
    * Step 1.
    */
   if ( ln == lmin ) {
     /*
      * Solve coarsest level.
      */
     d_fac_operator.solveCoarsestLevel( e ,
                                        r ,
                                        ln );
   }

   else {

      int i, num_components = e.getNumberOfComponents();

      /*
       * Step 2a.
       */
      d_fac_operator.computeCompositeResidualOnLevel( *d_tmp_residual ,
                                                      e ,
                                                      r ,
                                                      ln ,
                                                      true );
      for ( i=0; i<num_components; ++i ) {
        int tmp_id = d_tmp_error->getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels( ln ,
                                                ln );
        d_controlled_level_ops[i]->setToScalar( tmp_id ,
                                                0.0 );
      }
      /*
       * Step 2b.
       */
      d_fac_operator.smoothError( *d_tmp_error ,
                                  *d_tmp_residual ,
                                  ln ,
                                  d_presmoothing_sweeps );
      /*
       * Step 2c.
       */
      for ( i=0; i<num_components; ++i ) {
        int tmp_id = d_tmp_error->getComponentDescriptorIndex(i);
        int id = e.getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels( ln ,
                                                ln );
        d_controlled_level_ops[i]->add ( id ,
                                         id ,
                                         tmp_id );
      }
      /*
       * Step 3a.
       */
      d_fac_operator.computeCompositeResidualOnLevel( *d_tmp_residual ,
                                                      e ,
                                                      r ,
                                                      ln ,
                                                      true );
      d_fac_operator.restrictResidual( *d_tmp_residual,
                                       *d_tmp_residual ,
                                       ln-1 );
      /*
       * Step 3b.
       */
      for ( i=0; i<num_components; ++i ) {
        int tmp_id = e.getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels( ln-1 ,
                                                ln-1 );
        d_controlled_level_ops[i]->setToScalar( tmp_id ,
                                                0.0 );
      }
      /*
       * Step 3c.
       */
      facCycle_McCormick( e,
                          r,
                          u,
                          lmax,
                          lmin,
                          ln-1 );
      /*
       * Step 3d.
       */
      d_fac_operator.prolongErrorAndCorrect( e ,
                                             e ,
                                             ln );
      /*
       * Step 4a.
       */
      d_fac_operator.computeCompositeResidualOnLevel( *d_tmp_residual ,
                                                      e ,
                                                      r ,
                                                      ln ,
                                                      true );
      for ( i=0; i<num_components; ++i ) {
        int tmp_id = d_tmp_error->getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels( ln ,
                                                ln );
        d_controlled_level_ops[i]->setToScalar( tmp_id ,
                                                0.0 );
      }
      /*
       * Step 4b.
       */
      d_fac_operator.smoothError( *d_tmp_error ,
                                  *d_tmp_residual ,
                                  ln ,
                                  d_postsmoothing_sweeps );
      /*
       * Step 4c.
       */
      for ( i=0; i<num_components; ++i ) {
        int tmp_id = d_tmp_error->getComponentDescriptorIndex(i);
        int id = e.getComponentDescriptorIndex(i);
        d_controlled_level_ops[i]->resetLevels( ln ,
                                                ln );
        d_controlled_level_ops[i]->add ( id ,
                                         id ,
                                         tmp_id );
      }

   }

   return;
}



/*
*************************************************************************
*                                                                       *
* Reimplementation of the FAC cycle Michael Pernice coded up            *
* in the first version of FACPreconditioner<DIM> (the version          *
* that did not separate out the operators.                              *
*                                                                       *
*************************************************************************
*/
template<int DIM> void FACPreconditioner<DIM>::facCycle(
   SAMRAIVectorReal<DIM,double> &e,
   SAMRAIVectorReal<DIM,double> &r,
   SAMRAIVectorReal<DIM,double> &u,
   int lmax,
   int lmin )
{

   NULL_USE(u);

   int ln;

   /*
    * V-cycle descent.
    */
   for ( ln=lmax; ln>lmin; --ln ) {
      /*
       * Presmoothing.
       */
      d_fac_operator.smoothError( e ,
                                  r ,
                                  ln ,
                                  d_presmoothing_sweeps );
      /*
       * Compute residual to see how much correction is still needed:
       * d_tmp_residual <- r - A e
       */
      d_fac_operator.computeCompositeResidualOnLevel( *d_tmp_residual ,
                                                      e ,
                                                      r ,
                                                      ln ,
                                                      true );
      /*
       * Change the residual on the part of the next coarser level
       * below this level, so that the solve on the next coarser
       * level really solves for the correction for this level
       * where they overlap.
       */
      d_fac_operator.restrictResidual( *d_tmp_residual ,
                                       r ,
                                       ln-1 );
   }	// End V-cycle descent.


   /*
    * V-cycle ascent.
    */
   for ( ln=lmin; ln<=lmax; ++ln ) {
      /*
       * Update the error by prolonging and postsmoothing,
       * except on coarse level where an exact solve is performed.
       */
      if ( ln == lmin ) {
         /*
          * Solve coarsest level.
          */
         d_fac_operator.solveCoarsestLevel( e ,
                                            r ,
                                            ln );
      }
      else {
         /*
          * Apply the coarse level correction to this level.
          */
         d_fac_operator.prolongErrorAndCorrect( e ,
                                                e ,
                                                ln );
      

         /*
          * Postsmoothing on the error,
          * with inhomogeneous boundary condition.
          */
         /*
          * Postsmoothing.
          */
         d_fac_operator.smoothError( e ,
                                     r ,
                                     ln ,
                                     d_postsmoothing_sweeps );

      }


   }	// End V-cycle ascent.


   return;
}


/*
*************************************************************************
*                                                                       *
* Private member function to compute the composite residual on all      *
* levels and return residual norm.                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> double FACPreconditioner<DIM>::computeFullCompositeResidual(
   SAMRAIVectorReal<DIM,double> &r,
   SAMRAIVectorReal<DIM,double> &u,
   SAMRAIVectorReal<DIM,double> &f )
{

   d_fac_operator.computeCompositeResidualOnLevel
      ( r ,
        u ,
        f ,
        d_finest_ln ,
        false );

   for (int ln = d_finest_ln-1; ln >= d_coarsest_ln; --ln ) {

      // Bring down more accurate solution from finer level.
      d_fac_operator.restrictSolution(u, u, ln);

      d_fac_operator.computeCompositeResidualOnLevel
         ( r ,
           u ,
           f ,
           ln ,
           false );

      // Bring down more accurate residual from finer level.
      d_fac_operator.restrictResidual(r ,
                                      r ,
                                      ln );

   }

   double residual_norm =
      d_fac_operator.computeResidualNorm( r ,
                                          d_finest_ln ,
                                          d_coarsest_ln );

   return(residual_norm);

}






/*
*************************************************************************
*                                                                       *
* Print internal data, mostly for debugging.                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> void FACPreconditioner<DIM>::printClassData(std::ostream& os) const {
   os << "printing FACPreconditioner<DIM> data...\n"
      << "FACPreconditioner<DIM>: this = " << (FACPreconditioner<DIM>*)this << "\n"
      << "d_object_name = " << d_object_name << "\n"
      << "d_coarsest_ln = " << d_coarsest_ln << "\n"
      << "d_finest_ln = " << d_finest_ln << "\n"
      << "d_max_iterations = " << d_max_iterations << "\n"
      << "d_residual_tolerance = " << d_residual_tolerance << "\n"
      << "d_relative_residual_tolerance = " << d_relative_residual_tolerance << "\n"
      << "d_presmoothing_sweeps = " << d_presmoothing_sweeps << "\n"
      << "d_postsmoothing_sweeps = " << d_postsmoothing_sweeps << "\n"
      << "d_number_iterations = " << d_number_iterations << "\n"
      << "d_residual_norm = " << d_residual_norm << "\n"
      << "d_rhs_norm = " << d_rhs_norm << "\n"
      << std::endl;

}



}
}
#endif
