
//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/algorithm/method_of_lines/MethodOfLinesIntegrator.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Basic method-of-lines time integration algorithm
//

#ifndef included_algs_MethodOfLinesIntegrator_C
#define included_algs_MethodOfLinesIntegrator_C

#include "MethodOfLinesIntegrator.h"

#include <stdlib.h>
#include <fstream>

#include "Patch.h"
#include "PatchData.h"
#include "VariableDatabase.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION (2)

namespace SAMRAI {
    namespace algs {

/*
*************************************************************************
*									*
* The constructor and destructor for MethodOfLinesIntegrator<DIM>.     *
*									*
*************************************************************************
*/

template<int DIM> MethodOfLinesIntegrator<DIM>::MethodOfLinesIntegrator(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   MethodOfLinesPatchStrategy<DIM>* patch_strategy,
   bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(patch_strategy != ((MethodOfLinesPatchStrategy<DIM>*)NULL));
#endif

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }

   d_patch_strategy = patch_strategy;

   /*
    * Communication algorithms.
    */
   d_bdry_fill_advance     = new xfer::RefineAlgorithm<DIM>();
   d_fill_after_regrid     = new xfer::RefineAlgorithm<DIM>();
   d_fill_before_tagging   = new xfer::RefineAlgorithm<DIM>();
   d_coarsen_algorithm     = new xfer::CoarsenAlgorithm<DIM>();

   /*
    * hier::Variable contexts used in algorithm.
    */
   d_current = hier::VariableDatabase<DIM>::getDatabase()->getContext("CURRENT");
   d_scratch = hier::VariableDatabase<DIM>::getDatabase()->getContext("SCRATCH");
   d_patch_strategy->setInteriorContext(d_current);
   d_patch_strategy->setInteriorWithGhostsContext(d_scratch);

   /*
    * Set default to third-order SSP Runge-Kutta method.
    */
   d_order = 3;
   d_alpha_1.resizeArray(d_order);
   d_alpha_1[0] = 1.0;
   d_alpha_1[1] = 0.75;
   d_alpha_1[2] = 1.0/3.0;
   d_alpha_2.resizeArray(d_order);
   d_alpha_2[0] = 0.0;
   d_alpha_2[1] = 0.25;
   d_alpha_2[2] = 2.0/3.0;
   d_beta.resizeArray(d_order);
   d_beta[0] = 1.0;
   d_beta[1] = 0.25;
   d_beta[2] = 2.0/3.0;

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if ( is_from_restart ) {
	getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

}

/*
*************************************************************************
*									*
* Destructor tells tbox::RestartManager to remove this object from the   *
* list of restart items.                                                *
*									*
*************************************************************************
*/

template<int DIM> MethodOfLinesIntegrator<DIM>::~MethodOfLinesIntegrator()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

/*
*************************************************************************
*
   *
* Initialize integrator by:
   *
*
   *
*   (1) Setting the number of time data levels based on needs of
   *
*       the gridding algorithm
   *
*   (2) Invoking variable registration in patch strategy.
   *
*
   *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::initializeIntegrator(
   tbox::Pointer< mesh::GriddingAlgorithm<DIM> > gridding_alg)
{
   NULL_USE(gridding_alg);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!gridding_alg.isNull());
#endif

   /*
    * We may eventually need support for three (or more) time 
    * levels, information which may be accessed from the 
    * gridding algorithm.  Since we don't yet support this, 
    * this method simply registers variables with the integrator.
    */

   /*
    * Call variable registration in patch strategy.
    */
   d_patch_strategy->registerModelVariables(this);
}
/*
*************************************************************************
*									*
* Calculate the stable time increment by taking the minimum over        *
* all patches on all levels in the hierarchy.                           *
*									*
*************************************************************************
*/

template<int DIM> double MethodOfLinesIntegrator<DIM>::getTimestep(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const double time) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
#endif

   double dt = tbox::MathUtilities<double>::getMax();
   const int nlevels = hierarchy->getNumberOfLevels();

   for (int l = 0; l < nlevels; l++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->
                                             getPatchLevel(l);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!(level.isNull()));
#endif
      for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {

         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
 
         const double dt_patch = d_patch_strategy->
                                 computeStableDtOnPatch(*patch,
                                                        time);
         if (dt_patch < dt) dt = dt_patch;
      }
   }

   return(tbox::SAMRAI_MPI::minReduction(dt));

}

/*
*************************************************************************
*                                                                       *
* Advance the solution through the specified time increment using the   *
* general RK algorithm.  Each of the following steps is performed over  *
* all hierarchy levels.                                                 *
*                                                                       *
* (1) Copy solution values from current context to scratch context.     *
*                                                                       *
* (2) RK multistep loop for d(U)/dt = F(U):                             *
*                                                                       *
*    do i = 1, order                                                    *
*       U_i = U_n + alpha_i * dt/(order) * F(U_i)                       *
*    end do                                                             *
*                                                                       *
* (3) Copy last update of scratch solution to current context.          *
*                                                                       *
* Note that each update is performed by the concrete patch strategy     *
* in which the numerical routines are defined.                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::advanceHierarchy(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const double time,
   const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
#endif

   /*
    * Stamp data on all levels to current simulation time.
    */
   const int nlevels = hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < nlevels; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!(level.isNull()));
#endif

      level->setTime(time, d_current_data);
      level->setTime(time, d_scratch_data);
      level->setTime(time, d_rhs_data);

      /*
       * Allocate memory for U_scratch and rhs data
       */
      level->allocatePatchData(d_scratch_data, time);
      level->allocatePatchData(d_rhs_data, time);

      copyCurrentToScratch(level);
   }

   /*
    * Loop through Runge-Kutta steps
    */
   for (int rkstep = 0; rkstep < d_order; rkstep++) {

      /*
       * Loop through levels in the patch hierarchy and advance data on
       * each level by a single RK step.
       */
      for (int ln = 0; ln < nlevels; ln++) {

        /*
         * Fill ghost cells of all patches in level
         */
          d_bdry_sched_advance[ln]->fillData(time);

         /*
          * Loop through patches in current level and "singleStep" on each
          * patch.
          */
          tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
#ifdef DEBUG_CHECK_ASSERTIONS
          TBOX_ASSERT(!(level.isNull()));
#endif

          for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {

             tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
             d_patch_strategy->singleStep( *patch, 
                                           dt, 
                                           d_alpha_1[rkstep],
                                           d_alpha_2[rkstep],
                                           d_beta[rkstep] );

          } // patch loop

          if (ln > 0) {
             d_coarsen_schedule[ln]->coarsenData();
          }
          
      }  // levels loop
   
   }  // rksteps loop

   for (int ln = 0; ln < nlevels; ln++) {
      copyScratchToCurrent(hierarchy->getPatchLevel(ln));

      /* 
       * update timestamp to time after advance
       */
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
      level->setTime(time + dt, d_current_data);
   }

   /*
    * dallocate U_scratch and rhs data
    */
   for (int ln = 0; ln < nlevels; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
      level->deallocatePatchData(d_scratch_data);
      level->deallocatePatchData(d_rhs_data);
   }

}

/*
*************************************************************************
*                                                                       *
* Register the variables with the method of lines solution algorithm    *
* according to specified algorithm role (i.e., MOL_VAR_TYPE).           *
*                                                                       *
*       du/dt = F(u)                                                    *
*           u  - defined as "solution"                                  *
*         F(u) - defined as "right-hand-side"                           *
*                                                                       *
* Assignment of descriptor indices to variable lists, component         *
* selectors, and communication  algorithms takes place here.            *
*                                                                       *
* The different cases are:                                              *
*                                                                       *
* SOLN:                                                                 *
*            One time level of data is maintained for the current       *
*            solution and a "scratch" copy is used during the update    *
*            process.                                                   *
*                                                                       *
*            Two factories are needed: SCRATCH, CURRENT.                *
*                                                                       *
*            SCRATCH index is added to d_scratch_data.                  *
*            CURRENT index is added to d_current_data.                  *
*                                                                       *
* RHS:                                                                  *
*            Only one time level of data is stored and no scratch space *
*            is used.  Data may be set and manipulated at will in user  *
*            routines.  Data (including ghost values) is never touched  *
*            outside of user routines.                                  *
*                                                                       *
*            One factory needed: CURRENT.                               *
*                                                                       *
*            CURRENT index is added to d_current_data.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::registerVariable(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::IntVector<DIM>& ghosts,
   const MOL_VAR_TYPE m_v_type,
   const tbox::Pointer< xfer::Geometry<DIM> >& transfer_geom,
   const std::string &coarsen_name,
   const std::string &refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(variable.isNull()));
   TBOX_ASSERT(!(transfer_geom.isNull()));
#endif
   hier::VariableDatabase<DIM>* variable_db = hier::VariableDatabase<DIM>::getDatabase();

   switch(m_v_type) {

      case SOLN: {
      /*
       * Associate the current and scratch contexts with the variable in the
       * database and get the patch data identifiers.  The flag arrays will
       * be used to manage the allocation and deallocation of current and
       * scratch data.
       */
         d_soln_variables.appendItem(variable);

         const hier::IntVector<DIM> no_ghosts(0);

         const int current = variable_db->registerVariableAndContext(variable, 
                                                                     d_current,
                                                                     no_ghosts);

         const int scratch = variable_db->registerVariableAndContext(variable,
                                                                     d_scratch,
                                                                     ghosts);

         d_current_data.setFlag(current);

         d_scratch_data.setFlag(scratch);

         /*
          * Register variable and context needed for restart.
          */
         hier::VariableDatabase<DIM>::getDatabase()->
            registerPatchDataForRestart(current);

   /*
    * Ask the geometry for the appropriate refinement operator and register
    * that operator and the variables with the communication algorithms.
    * Two different communication algorithms are required by the RK method.
    * The Fillghosts algorithm is called during the normal Runge-Kutta time
    * stepping and fills the ghost cells of the scratch variables.  The regrid
    * algorithm is called after regrid and fills the current data on the new
    * level.
    */

         tbox::Pointer< xfer::RefineOperator<DIM> > refine_operator =
            transfer_geom->lookupRefineOperator(variable, refine_name);

         //  Fill ghosts for a variable using always the "scratch" context
         d_bdry_fill_advance->registerRefine(
            scratch,    // destination
            scratch,    // source
            scratch,    // temporary work space
            refine_operator);

         //  After regrid, use "current" context to communicate information
         //  to updated patches.  Use "scratch" as the temporary storage.
         d_fill_after_regrid->registerRefine(
            current,    // destination
            current,    // source
            scratch,    // temporary work space
            refine_operator);

         //  Before tagging error cells, copy data in current context to
         //  scratch context.  Note that this operation is not a simple
         //  copy - it also requires filling of ghost cells.  This is why
         //  it is designated as a refine operation.
         d_fill_before_tagging->registerRefine(
            scratch,    // destination
            current,    // source
            scratch,    // temporary work space
            refine_operator);

         tbox::Pointer< xfer::CoarsenOperator<DIM> > coarsen_operator =
            transfer_geom->lookupCoarsenOperator(variable, coarsen_name);

         //  Coarsen solution between levels during RK process so that
         //  coarser levels see the fine solution during integration.
         d_coarsen_algorithm->registerCoarsen(scratch,
                                              scratch,
                                              coarsen_operator);

         break;
      }

      case RHS: {
      /*
       * Associate the current context with the RHS variable in the
       * database.  The d_rhs_data component selector will be used to allocate and
       * de-allocate rhs data.
       * NOTE:  The d_rhs_data component selector was added 3/23/00 to facilitate 
       * allocation and de-allocation of rhs data for restarts.
       */
         d_rhs_variables.appendItem(variable);

         const int current = variable_db->registerVariableAndContext(variable,
                                                                     d_current,
                                                                     ghosts);

         d_rhs_data.setFlag(current);

         break;
      }

      default: {

         TBOX_ERROR(d_object_name << ":  "
                    << "unknown MOL_VAR_TYPE = " << m_v_type);

      }

   }

}

/*
*************************************************************************
*                                                                       *
* Allocate data for new level in hierarchy and initialize that data.    *
* If the new level replaces a pre-existing level in the hierarchy,      *
* data is copied from that level to the new level on their intersection.*
* Other data on the new level is set by interpolating from coarser      *
* levels in the hierarchy.  Then, user-defined initialization routines  *
* are called.                                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::initializeLevelData(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer< hier::BasePatchLevel<DIM> > old_level, 
   const bool allocate_data)
{

   NULL_USE(can_be_refined);
   NULL_USE(allocate_data);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number).isNull()));
   TBOX_ASSERT(level_number >= 0);
   if ( !(old_level.isNull()) ) {
     TBOX_ASSERT( level_number == old_level->getLevelNumber() );
   }
#endif

   tbox::Pointer< hier::PatchLevel<DIM> > level =
     hierarchy->getPatchLevel(level_number);

   /*
    * Allocate storage needed to initialize level and fill data from
    * coarser levels in AMR hierarchy.  
    */
   level->allocatePatchData(d_current_data, time);
   level->allocatePatchData(d_scratch_data, time);
 
   if ((level_number > 0) || !old_level.isNull()) {
     d_fill_after_regrid->createSchedule(level, 
                                         old_level,
                                         level_number-1,
                                         hierarchy,
                                         d_patch_strategy)->fillData(time);
   }

   level->deallocatePatchData(d_scratch_data);

   /*
    * Initialize current data for new level.
    */
   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      d_patch_strategy->initializeDataOnPatch(*patch, 
                                              time,
                                              initial_time);
   } 
}

/*
*************************************************************************
*									*
* Re-generate communication schedule after changes to the specified     *
* range of levels in the hierarchy.                                     *  
*									*
*************************************************************************
*/
template<int DIM> void MethodOfLinesIntegrator<DIM>::resetHierarchyConfiguration(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(finest_level);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (coarsest_level >= 0)
           && (coarsest_level <= finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()) );
#endif

   int finest_hiera_level = hierarchy->getFinestLevelNumber();

   //  If we have added or removed a level, resize the schedule arrays
   d_bdry_sched_advance.resizeArray(finest_hiera_level+1);
   d_coarsen_schedule.resizeArray(finest_hiera_level+1);

   //  Build coarsen and refine communication schedules.
   for (int ln = coarsest_level; ln <= finest_hiera_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!(level.isNull()));
#endif

      d_bdry_sched_advance[ln] =
         d_bdry_fill_advance->createSchedule(level,
                                             ln-1,
                                             hierarchy,
                                             d_patch_strategy);

      // coarsen schedule only for levels > 0
      if (ln > 0) {
         tbox::Pointer< hier::PatchLevel<DIM> > coarser_level = 
                                        hierarchy->getPatchLevel(ln-1);
         d_coarsen_schedule[ln] =
            d_coarsen_algorithm->createSchedule(coarser_level, level, NULL);
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Fill ghost cells for patches on level and call application-specific   *
* cell tagging routines.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::applyGradientDetector(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy, 
   const int level_number,
   const double time, 
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number).isNull()));
#endif

   tbox::Pointer< hier::PatchLevel<DIM> > level =
                                  hierarchy->getPatchLevel(level_number);

   level->allocatePatchData(d_scratch_data, time);

   /*
    * Transfer information from the "current" context to the "scratch"
    * context, on the current level. Note that ghosts will be filled
    * in this process.  We create and apply the schedule at the same
    * time because this routine is only called during
    * a regrid step, and the changing grid system means the schedule will
    * change since the last time it was called.
    */
   d_fill_before_tagging->createSchedule(level,
                                         level,
                                         level_number-1,
                                         hierarchy,
                                         d_patch_strategy)->fillData(time);

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      d_patch_strategy->
         tagGradientDetectorCells(*patch, 
                                  time,
                                  initial_time,
                                  tag_index,
                                  uses_richardson_extrapolation_too);
   }

   level->deallocatePatchData(d_scratch_data);

}

/*
*************************************************************************
*                                                                       *
* Writes the class version number, order, and				*
* alpha array to the database.						*
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION",
                   ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION);

   db->putInteger("d_order", d_order);
   db->putDoubleArray("d_alpha_1", d_alpha_1);
   db->putDoubleArray("d_alpha_2", d_alpha_2);
   db->putDoubleArray("d_beta", d_beta);
}

/*
*************************************************************************
*                                                                       *
* Reads in paramemters from the database overriding any values		*
* read in from the restart database. Also checks to make sure that	*
* number of alpha values specified equals order of Runga-Kutta scheme.	*
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> input_db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!input_db.isNull());
#endif

   if (is_from_restart) {

      if (input_db->keyExists("order")) {
         d_order = input_db->getInteger("order");
         if (d_order < 0) {
            TBOX_ERROR(d_object_name << ":  "
                      << "Negative `order' value specified in input.");
         }

      }

   } else {

      d_order = input_db->getIntegerWithDefault("order", d_order);


      if (input_db->keyExists("alpha_1")) {
         d_alpha_1 = input_db->getDoubleArray("alpha_1");
      } else { 
         TBOX_WARNING(d_object_name << ":  "
                      << "Key data `alpha_1' not found in input.  "
                      << "Using default values.  See class header.");
      }

      if (input_db->keyExists("alpha_2")) {
         d_alpha_2 = input_db->getDoubleArray("alpha_2");
      } else { 
         TBOX_WARNING(d_object_name << ":  "
                      << "Key data `alpha_2' not found in input.  "
                      << "Using default values.  See class header.");
      }

      if (input_db->keyExists("beta")) {
         d_beta = input_db->getDoubleArray("beta");
      } else { 
         TBOX_WARNING(d_object_name << ":  "
                      << "Key data `beta' not found in input.  "
                      << "Using default values.  See class header.");
      }

   }

   if ( d_alpha_1.getSize() != d_alpha_2.getSize() ||
       d_alpha_2.getSize() != d_beta.getSize() ) {
       TBOX_ERROR(d_object_name << ":  "
                 << "The number of alpha_1, alpha_2, and beta values "
		 << "specified in input is not consistent");
   }

   if ( d_alpha_1.getSize() != d_order ) {
      TBOX_WARNING(d_object_name << ":  "
                 << "The number of alpha values specified in input "
                 << "does not equal the Runga-Kutta order");
      d_order = d_alpha_1.getSize();
   }
   
}

/*
*************************************************************************
*                                                                       *
* Checks that class and restart file version numbers are equal.  If so, *
* reads in d_order and d_alpha from the database.  Also, does a 	*
* consistency check to make sure that the number of alpha values 	*
* specified equals the order of the Runga-Kutta scheme.			*
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> restart_db;
   if ( root_db->isDatabase(d_object_name) ) {
      restart_db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   int ver = restart_db->getInteger("ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION");
   if (ver != ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
                 << "Restart file version different than class version.");
   }

   d_order = restart_db->getInteger("d_order");
   d_alpha_1 = restart_db->getDoubleArray("d_alpha_1");
   d_alpha_2 = restart_db->getDoubleArray("d_alpha_2");
   d_beta = restart_db->getDoubleArray("d_beta");

   if ( d_alpha_1.getSize() != d_order ) {
      TBOX_ERROR(d_object_name << ":  "
                 << "The number of alpha values read from restart "
                 << "does not equal the Runga-Kutta order");
   }

}

/*
*************************************************************************
*                                                                       *
* Copy all solution data from current context to scratch context.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::copyCurrentToScratch(
   const tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(level.isNull()));
#endif

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      typename tbox::List< tbox::Pointer< hier::Variable<DIM> > >::Iterator
           soln_var = d_soln_variables.listStart();
      while (soln_var) {

         tbox::Pointer< hier::PatchData<DIM> > src_data = 
            patch->getPatchData(soln_var(), d_current);

         tbox::Pointer< hier::PatchData<DIM> > dst_data = 
            patch->getPatchData(soln_var(), d_scratch);

         dst_data->copy(*src_data);
         soln_var++;

       }

   }

}

/*
*************************************************************************
*                                                                       *
* Copy all solution data from scratch context to current context.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::copyScratchToCurrent(
   const tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(level.isNull()));
#endif

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      typename tbox::List< tbox::Pointer< hier::Variable<DIM> > >::Iterator
           soln_var = d_soln_variables.listStart();
      while (soln_var) {

          tbox::Pointer< hier::PatchData<DIM> > src_data = 
             patch->getPatchData(soln_var(), d_scratch);

          tbox::Pointer< hier::PatchData<DIM> > dst_data = 
             patch->getPatchData(soln_var(), d_current);

          dst_data->copy(*src_data);
          soln_var++;

       }

   }

}

/*
*************************************************************************
*                                                                       *
* Print all class data members for MethodOfLinesIntegrator<DIM> object.*
*                                                                       *
*************************************************************************
*/

template<int DIM> void MethodOfLinesIntegrator<DIM>::printClassData(std::ostream& os) const
{
   os << "\nMethodOfLinesIntegrator<DIM>::printClassData..." << std::endl;
   os << "\nMethodOfLinesIntegrator<DIM>: this = "
      << (MethodOfLinesIntegrator<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_order = " << d_order << std::endl;

   for (int j = 0; j < d_order; j++) {
     os << "d_alpha_1[" << j << "] = " << d_alpha_1[j] << std::endl;
     os << "d_alpha_2[" << j << "] = " << d_alpha_2[j] << std::endl;
     os << "d_beta[" << j << "] = " << d_beta[j] << std::endl;
   }

   os << "d_patch_strategy = " 
      << (MethodOfLinesPatchStrategy<DIM>*)d_patch_strategy << std::endl;
}

}
}
#endif
