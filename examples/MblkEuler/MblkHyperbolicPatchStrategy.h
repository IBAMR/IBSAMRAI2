//
// File:        MblkHyperbolicPatchStrategy.h
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Interface to patch routines for hyperbolic integration scheme.  
//

#ifndef included_MblkHyperbolicPatchStrategyXD
#define included_MblkHyperbolicPatchStrategyXD

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"
#include "RefinePatchStrategy.h"
#include "MultiblockCoarsenPatchStrategy.h"
#include "MultiblockGriddingAlgorithm.h" 
#include "MultiblockRefineSchedule.h" 
#include "MultiblockRefinePatchStrategy.h" 


/**
 * Class MblkHyperbolicPatchStrategy<DIM> is an abstract base class defining the 
 * interface between an MblkHyperbolicLevelIntegrator<DIM> object and operations 
 * applied to a single patch in a structured AMR hierarchy.   The operations
 * include patch initialization, dt calculation, flux computation, 
 * conservative differencing, and error estimation.  This class is derived 
 * from the xfer::RefinePatchStrategy<DIM> and xfer::CoarsenPatchStrategy<DIM> 
 * abstract base classes.   These base classes provide the interface for 
 * user-defined interlevel data refining and coarsening operations and the 
 * specification of physical boundary conditions.  The functions 
 * setPhysicalBoundaryConditions(), and pre/postprocessRefine() are
 * overloaded from the class xfer::RefinePatchStrategy<DIM>.  The operations
 * pre/postprocessCoarsen() are overloaded from xfer::CoarsenPatchStrategy<DIM>.
 * The pre/postprocessCoarsen/Refine() operations are given empty 
 * implementations here so that the user does not need to proovide them if
 * the operations are not needed.
 *
 * It is important to recognize that for the concrete patch strategy subclass 
 * and the MblkHyperbolicLevelIntegrator<DIM> to work together, the concrete 
 * strategy must know which patch data to operate on.  The patch data storage
 * is manipulated by the level integrator.  The set/clearDataContext() methods
 * allow the integrator to inform the patch strategy of the correct data 
 * context.  The concrete patch strategy subclass can access the appropriate
 * context via the getDataContext() method.
 *
 * @see algs::MblkHyperbolicLevelIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

using namespace SAMRAI;

class MblkHyperbolicLevelIntegrator;

class MblkHyperbolicPatchStrategy :
   public xfer::MultiblockRefinePatchStrategy<NDIM>,
   public xfer::MultiblockCoarsenPatchStrategy<NDIM> 
{
public:
   /**
    * Blank constructor for MblkHyperbolicPatchStrategy<DIM>.
    */
   MblkHyperbolicPatchStrategy();

   /**
    * Virtual destructor for MblkHyperbolicPatchStrategy<DIM>.
    */
   virtual ~MblkHyperbolicPatchStrategy();


   /**
    * Register specific variables needed in the numerical routines with the 
    * hyperbolic level integrator using the registerVariable() function in that
    * class.  The integrator manipulates storage for the data and this
    * registration defines the way in which data for each quantity will 
    * be manipulated on the patches.  Typically, the derived data quantities
    * for plotting are registered with a visualization data writer in this
    * routine as well, since the hyperbolic level integrator provides the
    * variable context for plotting (i.e., which data is available when a
    * plot file is generated).  The integrator pointer cannot be null in 
    * most cases.
    *
    * The gridding algorithm pointer is provided so that patch data objects 
    * may be registered with the load balancer object (owned by the gridding 
    * algorithm) for non-uniform load balancing, if needed.
    */
   virtual void registerModelVariables(
      MblkHyperbolicLevelIntegrator* integrator) = 0;

   /**
    * Set up parameters in the load balancer object (owned by the gridding
    * algorithm) if needed.  This function is called immediately after the
    * registerModelVariables() function is called.  The hyperbolic level 
    * integrator pointer is provided so that the integrator can be used
    * to manage data for the load balancer if needed (e.g., when using
    * non-uniform load balancing).
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual void setupLoadBalancer(
      MblkHyperbolicLevelIntegrator* integrator,
      mesh::BaseGriddingAlgorithm<NDIM>* gridding_algorithm)
   {
      NULL_USE(integrator);
      NULL_USE(gridding_algorithm);
   } 


   /**
    * Set the initial data on a patch interior only.  Note that no ghost cells 
    * need to be set in this routine regardless of whether the patch data
    * corresponding to the context requires ghost cells.  The data_time 
    * is the simulation time when the routine is called.  The boolean
    * initial_time is true if the routine is called at the initial time
    * when the hierarchy is initially constructed, otherwise it is false.
    */
   virtual void initializeDataOnPatch(
      hier::Patch<NDIM>& patch, 
      const double data_time,
      const bool initial_time) = 0;

   /**
    * Compute the stable time increment for a patch on the level with the
    * given number.  The boolean flag initial_time is true if the routine
    * is called at the initial simulation time; otherwise it is false.
    * The double argument dt_time is the simulation time.
    */
   virtual double computeStableDtOnPatch(
      hier::Patch<NDIM>& patch, 
      const bool initial_time,
      const double dt_time) = 0;

   /**
    * Compute TIME INTEGRALS of fluxes to be used in conservative difference
    * for patch integration.  That is, it is assumed that this numerical 
    * routine will compute the fluxes corresponding to the cell faces 
    * multiplied by the time increment.  Typically, the numerical flux is 
    * the normal flux at the cell face.  The flux integrals will be used in 
    * the conservative difference that updates the conserved quantities.
    *
    * Note that the numerical routines in this method generally require 
    * ghost cells.  Ghost cells data is filled before this routine is called.
    */
   virtual void computeFluxesOnPatch(hier::Patch<NDIM>& patch,
                                     const double time,
                                     const double dt) = 0;

   /**
    * Update patch data with a conservative difference (approximating
    * the divergence theorem) using the flux integrals computed in 
    * computeFluxesOnPatch() routine.  The boolean flag is true when this
    * routine is called during a flux synchronization step.   Otherwise,
    * it is false.  Note that the computeFluxesOnPatch() routine computes 
    * TIME INTEGRALs of the numerical fluxes (e.g., they have been multiplied
    * by dt).  So the conservative difference routine should be consistent
    * with this.
    */
   virtual void 
   conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                 const double time, 
                                 const double dt,
                                 bool at_syncronization) = 0;

   /**
    * This is an optional routine for user to process any application-specific 
    * patch strategy data BEFORE patches are advanced on the given level.
    * This routine is called after patch boundary data is filled 
    * (i.e., ghosts) and before computeFluxesOnPatch().  The arguments are:
    * level -- level that will be advanced, current_time -- current 
    * integration time, dt -- current time increment, first_step -- boolean 
    * flag that is true if advance is first in time step sequence on level 
    * (i.e., previous advance step was on another level, false otherwise, 
    * last_step -- boolean flag that is true if advance is last in time step 
    * sequence on level (i.e., synchronization with coarser level will occur 
    * immediately after this advance), regrid_advance -- boolean flag that 
    * is true if the advance is during a regridding phase (i.e., the advance 
    * is not used to integrate data on the hierarchy) in which case the 
    * results of the advance will be discarded.
    *
    * Note that when this routine is called, the scratch data is filled on
    * all patches (i.e., ghost cells) and that data is the same as the 
    * current level data on all patch interiors.  That is, both scratch and
    * current data correspond to current_time.
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual void 
   preprocessAdvanceLevelState(const tbox::Pointer< hier::PatchLevel<NDIM> >& level,
                               double current_time,
                               double dt,
                               bool first_step,
                               bool last_step,
                               bool regrid_advance)
   {
      NULL_USE(level);
      NULL_USE(current_time); 
      NULL_USE(dt); 
      NULL_USE(first_step);
      NULL_USE(last_step); 
      NULL_USE(regrid_advance); 
   }  

   /**
    * This is an optional routine for user to process any application-specific 
    * patch strategy data AFTER patches are advanced on the given level.
    * This routine is called after conservativeDifferenceOnPatch() is called 
    * and before computeStableDtOnPatch().  The arguments are:
    * level -- level that will be advanced, current_time -- current 
    * integration time, dt -- current time increment, first_step -- boolean 
    * flag that is true if advance is first in time step sequence on level 
    * (i.e., previous advance step was on another level, false otherwise, 
    * last_step -- boolean flag that is true if advance is last in time step 
    * sequence on level (i.e., synchronization with coarser
    * level will occur immediately after this advance), regrid_advance --
    * boolean flag that is true if the advance is during a regridding phase 
    * (i.e., the advance is not used to integrate data on the hierarchy) in 
    * which case the results of the advance will be discarded.
    * 
    * Note that when this routine is called, the scratch data is filled on
    * all patches (i.e., ghost cells) and that data is the same as the
    * new level data on all patch interiors.  That is, both scratch and
    * new data correspond to current_time + dt on patch interiors.  
    * The current data and ghost values correspond to the current_time.
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual void 
   postprocessAdvanceLevelState(const tbox::Pointer< hier::PatchLevel<NDIM> >& level,
                                double current_time,
                                double dt,
                                bool first_step,
                                bool last_step,
                                bool regrid_advance) 
   {
      NULL_USE(level); 
      NULL_USE(current_time); 
      NULL_USE(dt); 
      NULL_USE(first_step); 
      NULL_USE(last_step); 
      NULL_USE(regrid_advance); 
   }  

   /**
    * Tag cells on the given patch that require refinement based on
    * application-specific numerical quantities.  The tag index argument
    * indicates the index of the tag data on the patch data array.  The
    * boolean argument initial_error is true if tagging is being done at the
    * initial simulation time; otherwise, it is false.  The other boolean
    * flag uses_richardson_extrapolation_too is true when Richardson
    * extrapolation is used in addition to the gradient detector.  This flag
    * helps users manage multiple regridding criteria.
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual void tagGradientDetectorCells(
      hier::Patch<NDIM>& patch, 
      const double regrid_time, 
      const bool initial_error, 
      const int tag_index, 
      const bool uses_richardson_extrapolation_too);

   /**
    * Tag cells based from differences computed in the Richardson 
    * extrapolation.  The Richardson
    * extrapolation algorithm creates a coarsened version of some hierarchy
    * patch level and advances data in time on both the coarsened patch
    * level and the hierarchy level. This routine takes the data resulting
    * from the advance on both the coarse and fine levels, compares them, and
    * tags cells according to the difference.
    * \verbatim
    *              (2)
    *      n+1 ^------->x finish      (1) advanced_coarse
    *          |        ^             (2) coarsened_fine
    *  time  n -        |
    *          ^        |(1)
    *          |        |
    *          <--------o start
    *        fine     coarse
    * \endverbatim
    *
    * The patch supplied to this routine is on the coarsened level.  However,
    * the error_level_number corresponds to the actual hierarchy level 
    * from which it was coarsened.   Data resides on this patch in two 
    * contexts - ``advanced_coarse'' and ``coarsened_fine''.  Advanced 
    * coarse is data advanced on the coarsened version of the level, while
    * coarsened fine is the data advanced on the fine level and then
    * coarsened to the coarse level.  The regrid time and the time increment
    * are given for the actual hierarchy level.  The error coarsen ratio 
    * argument is the ratio between the index spaces on the hierarchy level
    * and the coarsened hierarchy level.  The boolean flag ``initial_error''
    * is true when the error estimation is performed at the initial simulation
    * time; i.e., when the hierarchy levels are being constructed for the first
    * time.  The tag index argument is the index of the tag data on the patch 
    * data array.  The other boolean flag uses_gradient_detector_too is 
    * true when a gradient detector scheme is used in addition to Richardson
    * extrapolation.  This flag helps users manage multiple regridding 
    * criteria. 
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual void tagRichardsonExtrapolationCells(
      hier::Patch<NDIM>& patch,
      const int error_level_number,
      const tbox::Pointer<hier::VariableContext> coarsened_fine,
      const tbox::Pointer<hier::VariableContext> advanced_coarse,
      const double regrid_time,
      const double deltat, 
      const int error_coarsen_ratio,
      const bool initial_error, 
      const int tag_index,
      const bool uses_gradient_detector_too);

   /**
    * Set user-defined boundary conditions at the physical domain boundary.
    */
   virtual void setPhysicalBoundaryConditions(
      hier::Patch<NDIM>& patch, 
      const double fill_time, 
      const hier::IntVector<NDIM>& ghost_width_to_fill) = 0;


   /**
    * Fill the singularity conditions for the multi-block case
    */
   virtual void fillSingularityBoundaryConditions(
      hier::Patch<NDIM>& patch,
      tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
      const double fill_time,
      const hier::Box<NDIM>& fill_box,
      const hier::BoundaryBox<NDIM>& bbox) = 0;
   
   /**
    * Return maximum stencil width needed for user-defined 
    * data interpolation operations.  Default is to return 
    * zero, assuming no user-defined operations provided.
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual hier::IntVector<NDIM> getRefineOpStencilWidth() const
   {
      return(hier::IntVector<NDIM>(0));
   } 

   /**
    * Pre- and post-processing routines for implementing user-defined
    * spatial interpolation routines applied to variables.  The 
    * interpolation routines are used in the hyperbolic AMR algorithm
    * for filling patch ghost cells before advancing data on a level
    * and after regridding a level to fill portions of the new level
    * from some coarser level.  These routines are called automatically
    * from within patch boundary filling schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch routines.  However, the routines only need to perform some 
    * operations when "USER_DEFINED_REFINE" is given as the interpolation 
    * method for some variable when the patch routines register variables
    * with the hyperbolic level integration algorithm, typically.  If the 
    * user does not provide operations that refine such variables in either 
    * of these routines, then they will not be refined.
    *
    * The order in which these operations are used in each patch 
    * boundary filling schedule is:
    * 
    * - \b (1) {Call user's preprocessRefine() routine.}
    * - \b (2) {Refine all variables with standard interpolation operators.}
    * - \b (3) {Call user's postprocessRefine() routine.}
    * 
    * Note that these functions are not pure virtual. They are given 
    * dummy implementations here so that users may ignore them when 
    * inheriting from this class.
    */
   virtual void preprocessRefine(hier::Patch<NDIM>& fine,
                                 const hier::Patch<NDIM>& coarse,
                                 const hier::Box<NDIM>& fine_box,
                                 const hier::IntVector<NDIM>& ratio) 
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   ///
   virtual void postprocessRefine(hier::Patch<NDIM>& fine,
                                  const hier::Patch<NDIM>& coarse,
                                  const hier::Box<NDIM>& fine_box,
                                  const hier::IntVector<NDIM>& ratio) 
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /**
    * Return maximum stencil width needed for user-defined
    * data coarsen operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    *
    * Note that this function is not pure virtual. It is given a 
    * dummy implementation here so that users may ignore it when 
    * inheriting from this class.
    */
   virtual hier::IntVector<NDIM> getCoarsenOpStencilWidth() const
   {
      return(hier::IntVector<NDIM>(0));
   }

   /**
    * Pre- and post-processing routines for implementing user-defined
    * spatial coarsening routines applied to variables.  The coarsening 
    * routines are used in the hyperbolic AMR algorithm synchronizing 
    * coarse and fine levels when they have been integrated to the same
    * point.  These routines are called automatically from within the 
    * data synchronization coarsen schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch routines.  However, the routines only need to perform some
    * operations when "USER_DEFINED_COARSEN" is given as the coarsening
    * method for some variable when the patch routines register variables
    * with the hyperbolic level integration algorithm, typically.  If the
    * user does not provide operations that coarsen such variables in either
    * of these routines, then they will not be coarsened.
    *
    * The order in which these operations are used in each coarsening
    * schedule is:
    * 
    * - \b (1) {Call user's preprocessCoarsen() routine.}
    * - \b (2) {Coarsen all variables with standard coarsening operators.}
    * - \b (3) {Call user's postprocessCoarsen() routine.}
    * 
    * Note that these functions are not pure virtual. They are given 
    * dummy implementations here so that users may ignore them when 
    * inheriting from this class.
    */
   virtual void preprocessCoarsen(hier::Patch<NDIM>& coarse,
                                  const hier::Patch<NDIM>& fine,
                                  const hier::Box<NDIM>& coarse_box,
                                  const hier::IntVector<NDIM>& ratio) 
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   ///
   virtual void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                                   const hier::Patch<NDIM>& fine,
                                   const hier::Box<NDIM>& coarse_box,
                                   const hier::IntVector<NDIM>& ratio)
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   /**
    * Return pointer to patch data context.
    */
   tbox::Pointer<hier::VariableContext> getDataContext() const
   {
      return(d_data_context);
   }

   /**
    * The hyperbolic integrator controls the context for the data to be used 
    * in the numerical routines implemented in the concrete patch strategy.
    * The setDataContext() allows the integrator to set the context for 
    * data on a patch on which to operate.
    */
   void setDataContext(tbox::Pointer<hier::VariableContext> context)
   {
      d_data_context = context;
   }

   /**
    * The clearDataContext() routine resets the data context to be null.
    */
   void clearDataContext()
   {
      d_data_context.setNull();
   }

private:

   tbox::Pointer<hier::VariableContext> d_data_context; 

};


#endif
