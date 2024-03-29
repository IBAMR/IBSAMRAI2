//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mesh/gridding/StandardTagAndInitStrategy.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Strategy interface for Richardson Extrapolation algorithm.
//
 
#ifndef included_mesh_StandardTagAndInitStrategy
#define included_mesh_StandardTagAndInitStrategy

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace mesh {

/**
 * Class StandardTagAndInitStrategy<DIM> is an abstract base class 
 * that defines a Strategy pattern interface for concrete cell tagging 
 * and level initialization operations that are needed by the 
 * StandardTagAndInitialize<DIM> class.  This base class insulates 
 * that algorithm class from routines for initializing a new level in 
 * the hierarchy and for tagging cells to be refined.  Generally, these 
 * operations are specific to the problem being solved and the solution 
 * methods being employed.  An object of this type is passed into the 
 * StandardTagAndInitialize<DIM> constructor.
 * 
 * This base class has two pure virtual functions:
 * initializeLevelData(), and resetHierarchyConfiguration() that must
 * be implemented for any regridding method.  The first function sets 
 * data on a new level after regridding.  The second function is called
 * at the end of the regridding process and can be used to set communication
 * schedules, for example, which depend on the configuration of the AMR
 * patch hierarchy.  Other routines are virtual here and given default 
 * implementations as they are specific to only one type of error estimation 
 * method.  Gradient detector functionality requires an implementation of 
 * the applyGradientDetector() routine.  The Richardson extrapolation method 
 * requires implementations of the methods: applyRichardsonExtrapolation(), 
 * coarsenDataForRichardsonExtrapolation(), getLevelDt(), advanceLevel(), 
 * resetTimeDependentData(), and resetDataToPreadvanceState().
 *
 * @see mesh::StandardTagAndInitialize.
 */

template<int DIM> class StandardTagAndInitStrategy
:
public virtual tbox::DescribedClass
{
public:
   /**
    * Default constructor for 
    * StandardTagAndInitStrategy<DIM>. 
    */     
   StandardTagAndInitStrategy();

   /**
    * Empty destructor for 
    * StandardTagAndInitStrategy<DIM>.
    */
   virtual ~StandardTagAndInitStrategy();

   /**
    * Determine time increment to advance data on level. The 
    * recompute_dt option specifies whether to compute 
    * the timestep using the current level data or to return the value
    * stored by the time integrator. The default true setting means
    * the timestep will be computed if no value is supplied.  
    *
    * This routine is only when Richardson extrapolation is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */
   virtual double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
              const double dt_time,
              const bool initial_time);

   /**
    * Advance data on all patches on specified patch level from current time
    * (current_time) to new time (new_time).   This routine is called only
    * during time-dependent regridding procedures, such as Richardson
    * extrapolation.  It is virtual with an empty implementation here (rather
    * than pure virtual) so that users are not required to provide an
    * implementation when the function is not needed.  The boolean arguments
    * are used to determine the state of the algorithm and the data when the
    * advance routine is called.  Note that this advance function is also
    * used during normal time integration steps.
    *
    * When this function is called, the level data required to begin the
    * advance must be allocated and be defined appropriately.  Typically,
    * this is equivalent to what is needed to initialize a new level after
    * regridding.  Upon exiting this routine, both current and new data may
    * exist on the level.  This data is needed until level synchronization
    * occurs, in general. Current and new data may be reset by calling
    * the member function resetTimeDependentData().
    *
    * This routine is called from two different points within the Richardson
    * exptrapolation process: to advance a temporary level that is coarser
    * than the hierarchy level on which error estimation is performed, and
    * to advance the hierarchy level itself.  In the first case, the values of
    * the boolean flags are:
    * 


    *    - \b  first_step 
    *        = true.
    *    - \b  last_step 
    *        = true.
    *    - \b  regrid_advance 
    *        = true.
    * 


    * In the second case, the values of the boolean flags are:
    * 


    *    - \b  first_step 
    *      (when regridding during time integration sequence)
    *        = true when the level is not coarsest level to synchronize
    *          immediately before the regridding process; else, false.
    *      (when generating initial hierarchy construction)
    *        = true, even though there may be multiple advance steps.
    *    - \b  last_step 
    *        = true when the advance is the last in the Richardson
    *          extrapolation step sequence; else false.
    *    - \b  regrid_advance 
    *        = true.
    * 


    */
   virtual double
   advanceLevel(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
                const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                const double current_time,
                const double new_time,
                const bool first_step,
                const bool last_step,
                const bool regrid_advance = false);

   /**
    * Reset time-dependent data storage for the specified patch level.
    *
    * This routine only applies when Richardson extrapolation is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */
   virtual void
   resetTimeDependentData(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
                          const double new_time,
                          const bool can_be_refined);

   /**
    * Reset data on the patch level by destroying all patch data other
    * than that which is needed to initialize the solution on that level.
    * In other words, this is the data needed to begin a time integration
    * step on the level.
    *
    * This routine is only when Richardson extrapolation is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */
   virtual void
   resetDataToPreadvanceState(const tbox::Pointer< hier::BasePatchLevel<DIM> > level);

   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old 
    * level.  The specific operations that occur when initializing level 
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
   virtual void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                       const int level_number,
                       const double init_data_time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<DIM> > old_level = 
                             tbox::Pointer< hier::BasePatchLevel<DIM> >(NULL),
		       const bool allocate_data = true) = 0;

   /**
    * After hierarchy levels have changed and data has been initialized on 
    * the new levels, this routine can be used to reset any information 
    * needed by the solution method that is particular to the hierarchy 
    * configuration.  For example, the solution procedure may cache 
    * communication schedules to amortize the cost of data movement on the 
    * AMR patch hierarchy.  This function will be called by the gridding 
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must 
    * make the solution consistent across multiple levels after the hierarchy 
    * is changed, this process may be invoked by this routine.  Of course the 
    * details of these processes are determined by the particular solution 
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void resetHierarchyConfiguration(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int coarsest_level,
      const int finest_level) = 0; 

   /**
    * Set integer tags to "one" in cells where refinement of the given
    * level should occur according to some user-supplied gradient criteria.
    * The double time argument is the regrid time.  The integer "tag_index"
    * argument is the patch descriptor index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument 
    * initial_time indicates whether the level is being subject to refinement 
    * at the initial simulation time.  If it is false, then the error 
    * estimation process is being invoked at some later time after the AMR 
    * hierarchy was initially constructed.  Typically, this information is 
    * passed to the user's patch tagging routines since the error 
    * estimator or gradient detector may be different in each case.
    * 
    * The boolean uses_richardson_extrapolation_too is true when Richardson
    * extrapolation error estimation is used in addition to the gradient
    * detector, and false otherwise.  This argument helps the user to
    * manage multiple regridding criteria.
    *
    * This routine is only when gradient detector is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */
   virtual void 
   applyGradientDetector(const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                         const int level_number,
                         const double error_data_time,
                         const int tag_index,
                         const bool initial_time,
                         const bool uses_richardson_extrapolation_too);

   /**
    * Set integer tags to "one" in cells where refinement of the given
    * level should occur according to some user-supplied Richardson
    * extrapolation criteria.  The "error_data_time" argument is the
    * regrid time.  The "deltat" argument is the time increment to advance
    * the solution on the level to be refined.  Note that that level is
    * finer than the level in the argument list, in general.  The
    * ratio between the argument level and the actual hierarchy level
    * is given by the integer "coarsen ratio".
    *
    * The integer "tag_index" argument is the patch descriptor index of
    * the cell-centered integer tag array on each patch in the hierarchy.
    *
    * The boolean argument initial_time indicates whether the level is being
    * subject to refinement at the initial simulation time.  If it is false,
    * then the error estimation process is being invoked at some later time
    * after the AMR hierarchy was initially constructed.  Typically, this
    * information is passed to the user's patch tagging routines since the
    * application of the Richardson extrapolation process may be different
    * in each case.
    *
    * The boolean uses_gradient_detector_too is true when a gradient
    * detector procedure is used in addition to Richardson extrapolation,
    * and false otherwise.  This argument helps the user to manage multiple
    * regridding criteria.
    *
    * This routine is only when Richardson extrapolation is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */
   virtual void
   applyRichardsonExtrapolation(const tbox::Pointer< hier::PatchLevel<DIM> > level,
                                const double error_data_time,
                                const int tag_index,
                                const double deltat,
                                const int error_coarsen_ratio,
                                const bool initial_time,
                                const bool uses_gradient_detector_too);

   /**
    * Coarsen solution data from level to coarse_level for Richardson
    * extrapolation.  Note that this routine will be called twice during
    * the Richardson extrapolation error estimation process, once to set
    * data on the coarser level and once to coarsen data from after 
    * advancing the fine level.  The init_coarse_level boolean argument 
    * indicates whether data is set on the coarse level by coarsening the 
    * "old" time level solution or by coarsening the "new" solution on the 
    * fine level (i.e., after it has been advanced).
    *
    * This routine is only when Richardson extrapolation is being used.
    * It is virtual with an empty implementation here (rather than pure
    * virtual) so that users are not required to provide an implementation
    * when the function is not needed.
    */ 
   virtual void
   coarsenDataForRichardsonExtrapolation(
			   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
			   const int level_number,
                           const tbox::Pointer< hier::PatchLevel<DIM> > coarser_level,
                           const double coarsen_data_time,
			   const bool before_advance);

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "StandardTagAndInitStrategy.C"
#endif
