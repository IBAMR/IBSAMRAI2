/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/poisson/CartesianRobinBcHelper.h $
 * Package:     SAMRAI application utilities
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2132 $
 * Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
 * Description: Robin boundary condition support on cartesian grids.
 */

#ifndef included_solv_CartesianRobinBcHelper
#define included_solv_CartesianRobinBcHelper


#include "SAMRAI_config.h"


/*
 * SAMRAI classes
 */

#include "BoundaryBox.h"

#include "PatchLevel.h"

#include "Patch.h"

#include "ArrayData.h"

#include "CellData.h"

#include "NodeData.h"

#include "RobinBcCoefStrategy.h"

#include "tbox/Pointer.h"

#ifndef included_String
#include <string>
#define included_String
#endif

#ifndef xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif


namespace SAMRAI {
    namespace solv {


/*!
 * @brief Helper utility for setting Robin boundary conditions.
 *
 * This class is intended as a helper for performing the tedious
 * task of setting boundary values for scalar quantities for the
 * general case of boundary conditions known as the Robin boundary
 * condition.
 *
 * It uses the Robin boundary condition coefficients specified
 * by a RobinBcCoefStrategy<DIM> object to determine the boundary
 * value to set.  The exact value set depends on the allignment of
 * the data and is derived from various discrete approximations
 * of the Robin formula.  This class currently supports cell-centered
 * alignment and will support node-centered alignment in the future.
 *
 * See RobinBcCoefStrategy<DIM> for the description of the Robin
 * boundary condition.
 *
 * This class inherits and implements virtual functions from
 * xfer::RefinePatchStrategy<DIM> so it may be used to help create
 * communication schedules if desired.
 */
template<int DIM> class CartesianRobinBcHelper :
  public xfer::RefinePatchStrategy<DIM>
{

public:

   /*!
    * @brief Constructor using.
    *
    * Requires a concrete implementation of RobinBcCoefStrategy<DIM>.
    *
    * @param object_name Name of the object, for general referencing.
    * @param coef_strategy Coefficients strategy being helped.
    */
   CartesianRobinBcHelper(
      std::string object_name=std::string(),
      RobinBcCoefStrategy<DIM> *coef_strategy=NULL );



   /*!
    * @brief Destructor.
    */
   virtual ~CartesianRobinBcHelper();



   //@{ @name xfer::RefinePatchStrategy<DIM> virtuals

   virtual void setPhysicalBoundaryConditions (
      hier::Patch<DIM> &patch ,
      const double fill_time ,
      const hier::IntVector<DIM> &ghost_width_to_fill );
   hier::IntVector<DIM> getRefineOpStencilWidth () const;
   virtual void preprocessRefineBoxes (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::BoxList<DIM> &fine_boxes ,
      const hier::IntVector<DIM> &ratio );
   virtual void preprocessRefine (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::Box<DIM> &fine_box ,
      const hier::IntVector<DIM> &ratio );
   virtual void postprocessRefineBoxes (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::BoxList<DIM> &fine_boxes ,
      const hier::IntVector<DIM> &ratio );
   virtual void postprocessRefine (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::Box<DIM> &fine_box ,
      const hier::IntVector<DIM> &ratio );

   //@}



   //@{

   /*!
    * @name Functions to set boundary condition values
    */

   /*!
    * @brief Set the physical boundary condition by setting the
    * value of the first ghost cells.
    *
    * This function has an interface similar to the virtual function
    * xfer::RefinePatchStrategy<DIM>::setPhysicalBoundaryConditions(),
    * and it may be used to help implement that function,
    * but it does not serve the same purpose.  The primary
    * differences are:
    * -# It specializes to cell-centered variables.
    * -# Only one ghost cell width is filled.  Setting a Robin
    *    boundary condition for cell-centered quantities requires
    *    only one ghost cell to be set.
    *    (More ghost cells can be filled by continuing the linear
    *    distribution of data beyond the first cell, but that is
    *    not implemented at this time.)
    * -# User must specify the index of the data whose ghost
    *    cells need to be filled.  This index is used to determine
    *    the variable for which to set the boundary coefficients
    *    and to get the data to be set.
    *
    * This function calls RobinBcStrategy::setBcCoefs() to
    * get the coefficients, then it sets the values in the first
    * ghost cell on the boundary.
    *
    * To determine the value for the ghost cell,
    * a @em linear approximation in the direction normal to the
    * boundary is assumed.  We write the following discrete
    * approximations:
    * @f[ u_b = \frac{ u_i + u_o }{2} @f]
    * @f[ [u_n]_b = \frac{ u_o - u_i }{h} @f]
    * where the subscript b stands for the the boundary,
    * i stands for the first cell inside the boundary and
    * o stands for the first cell outside the boundary
    * and h is the grid spacing normal to the boundary.
    * Applying this to the Robin formula gives
    * @f[ u_o = \frac{ h\gamma + u_i( \beta - \frac{h}{2} \alpha ) }
    * { \beta + \frac{h}{2} \alpha } @f] or equivalently
    * @f[ u_o = \frac{ hg + u_i (1-a(1+\frac{h}{2})) }{ 1-a(1-\frac{h}{2}) } @f]
    *
    * After setting the edge (face in 3D) boundary conditions,
    * linear approximations are used to set the boundary conditions
    * of higher boundary types (nodes in 2D, edges and nodes in 3D).
    *
    * In some cases, the calling function wants to set the
    * boundary condition homogeneously, with g=0.
    * This is useful in problems where the the solution of the
    * homogeneous problem is required in solving the inhomogeneous
    * problem.  This function respects such requests specified
    * through the argument @c homogeneous_bc.
    *
    * @internal To be more general to other data types,
    * there could be versions for other data types also,
    * such as ...InNodes, ...InFaces, etc.  However, I'm not
    * sure exactly how to implement the Robin boundary condition
    * on the faces and nodes when m != 1.  Should the boundary
    * value be set or should the first ghost value be set?
    *
    * @internal I have not addressed possibility of differences
    * in chosing the discrete formulation with which to compute
    * the boundary value.  The above formulation is obviously
    * one specific approximation, but there could be others.
    * If anoter approximation is required, there should be
    * another class like this or this class can offer a choice
    * to be set by the user.  I favor another class.
    *
    * @internal Since the data alignment can be found through
    * the target_data_id, these types of functions may be changed
    * to just plain setBoundaryValues or setBoundaryValuesInBoundaryBoxes
    * since it does assume boundary boxes.  This may have to be
    * expanded to later include coarse-fine boundary boxes
    * for more generality.
    *
    * @param patch hier::Patch on which to set boundary condition
    * @param fill_time Solution time corresponding to filling
    * @param ghost_width_to_fill Max ghost width requiring fill
    * @param target_data_id hier::Patch data index of data to be set.
    *        This data must be a cell-centered double.
    * @param homogeneous_bc Set a homogeneous boundary condition.
    *    This means g=0 for the boundary.
    */
   void setBoundaryValuesInCells (
      hier::Patch<DIM> &patch ,
      const double fill_time ,
      const hier::IntVector<DIM> &ghost_width_to_fill ,
      int target_data_id ,
      bool homogeneous_bc=false ) const;


   /*!
    * @brief Set ghost cells for an entire level.
    *
    * Loop through all patches on the given level and call
    * setBoundaryValuesInCells(hier::Patch<DIM> &patch,
    *                          const double fill_time ,
    *                          const hier::IntVector<DIM> &ghost_width_to_fill ,
    *                          int target_data_id ,
    *                          bool homogeneous_bc=false )
    * for each.
    *
    * @param level PatchLevel on which to set boundary condition
    * @param fill_time Solution time corresponding to filling
    * @param ghost_width_to_fill Max ghost width requiring fill
    * @param target_data_id hier::Patch data index of data to be set.
    *        This data must be a cell-centered double.
    * @param homogeneous_bc Set a homogeneous boundary condition.
    *    This means g=0 for the boundary.
    */
   void setBoundaryValuesInCells (
      hier::PatchLevel<DIM> &level ,
      const double fill_time ,
      const hier::IntVector<DIM> &ghost_width_to_fill ,
      int target_data_id ,
      bool homogeneous_bc=false ) const;


   /*!
    * @brief Set the physical boundary condition by setting the
    * value of the boundary nodes.
    *
    * This function is not yet implemented!
    *
    * There are some decisions that must be made before
    * the implementation can be written.
    * -# Do we set the values on the boundary or one cell
    *    away from the boundary?
    * -# What is the discrete formulation we should use
    *    to compute the value to be set?
    *
    * This function has an interface similar to the virtual function
    * xfer::RefinePatchStrategy<DIM>::setPhysicalBoundaryConditions(),
    * and it may be used to help implement that function,
    * but it does not serve the same purpose.  The primary
    * differences are:
    * -# It specializes to node-centered variables.
    * -# User must specify the index of the data whose ghost
    *    cells need to be filled.  This index is used to determine
    *    the variable for which to set the boundary coefficients
    *    and to get the data to be set.
    *
    * This function calls RobinBcStrategy::setBcCoefs() to get the
    * coefficients, then it sets the values at the boundary nodes.
    *
    * In some cases, the calling function wants to set the
    * boundary condition homogeneously, with g=0.
    * This is useful in problems where the the solution of the
    * homogeneous problem is required to solving the inhomogeneous
    * problem.  This function respects such requests specified
    * through the argument @c homogeneous_bc.
    *
    * @param patch hier::Patch on which to set boundary condition
    * @param fill_time Solution time corresponding to filling
    * @param target_data_id hier::Patch data index of data to be set.
    * @param homogeneous_bc Set a homogeneous boundary condition.
    *    This means g=0 for the boundary.
    */
   void setBoundaryValuesAtNodes (
      hier::Patch<DIM> &patch ,
      const double fill_time ,
      int target_data_id ,
      bool homogeneous_bc=false ) const;

   //@}



   //@{
   /*!
    * @name Ways to provide the Robin bc coefficients
    */

   /*!
    * @brief Provide an implementation of the RobinBcCoefStrategy<DIM>
    * for determining the boundary coefficients.
    *
    * Provide the implementation that can be used to set the
    * Robin bc coefficients.
    *
    * @param coef_strategy tbox::Pointer to a concrete inmplementation of
    *        the coefficient strategy.
    */
   void setCoefImplementation(
      const RobinBcCoefStrategy<DIM> *coef_strategy );


   /*!
    * @brief Set the data id that should be filled when setting
    * physical boundary conditions.
    *
    * When setPhysicalBoundaryConditions is called, the data
    * specified will be set.  This information is required because
    * the it is not passed in through the argument list of
    * setPhysicalBounaryConditions.
    */
   void setTargetDataId( int target_data_id );



   /*!
    * @brief Set whether boundary filling should assume homogeneous
    * conditions.
    *
    * In certain circumstances, only the value of a is needed, while
    * the value of g is temporarily not required and taken to be zero.
    * (An example is in setting the boundary condition for error
    * value in an iterative method.)  In such cases, use this function
    * to set a flag that will cause a null pointer to be given to
    * setBcCoefs() to indicate that fact.
    */
   void setHomogeneousBc( bool homogeneous_bc );


   //@}



private:

   /*!
    * @brief Trim a boundary box so that it does not stick out
    * past a given box.
    *
    * Certain boundary-related operations occur on patch data that
    * do not or cannot extend past the edgr or corner of a patch.
    * This function is used to trim down the parts of the boundary box
    * that extend past those points so that a suitable index range
    * is achieved.
    *
    * The boundary box trimmed must be of type 1 or 2.
    *
    * @param boundary_box Boundary box to be trimmed.
    * @param limit_box hier::Box to not stick past
    *
    * @return New trimmed boundary box.
    */
   hier::BoundaryBox<DIM> trimBoundaryBox(
      const hier::BoundaryBox<DIM> &boundary_box,
      const hier::Box<DIM> &limit_box ) const;

   /*!
    * @brief Return box describing the index space of boundary nodes
    * defined by a boundary box.
    *
    * Define a box describing the indices of the nodes corresponding
    * to the input boundary box.  These nodes lie on the boundary
    * itself.
    *
    * The input boundary_box must be of type 1
    * (see hier::BoundaryBox::getBoundaryType()).
    *
    * @param boundary_box input boundary box
    * @return a box to define the node indices corresponding to
    *   boundary_box
    */
   hier::Box<DIM> makeNodeBoundaryBox(
      const hier::BoundaryBox<DIM> &boundary_box ) const;

   /*!
    * @brief Return box describing the index space of faces
    * defined by a boundary box.
    *
    * Define a box describing the indices of the codimension 1
    * surface corresponding to the input boundary box.
    *
    * The input boundary_box must be of type 1
    * (see hier::BoundaryBox::getBoundaryType()).
    *
    * This is a utility function for working with the
    * indices coresponding to a boundary box but coincide
    * with the patch boundary.
    *
    * @param boundary_box input boundary box
    * @return a box to define the face indices corresponding to
    *    boundary_box
    */
   hier::Box<DIM> makeFaceBoundaryBox(
      const hier::BoundaryBox<DIM> &boundary_box ) const;

   std::string d_object_name;

   /*!
    * @brief Coefficient strategy giving a way to get to
    * Robin bc coefficients.
    */
   const RobinBcCoefStrategy<DIM> *d_coef_strategy;

   /*!
    * @brief hier::Index of target patch data when filling ghosts.
    */
   int d_target_data_id;

   /*!
    * @brief Whether to assumg g=0 when filling ghosts.
    */
   bool d_homogeneous_bc;



   /*!
    * @brief Timers for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_set_boundary_values_in_cells;
   tbox::Pointer<tbox::Timer> t_use_set_bc_coefs;
};

}
}

#endif	// included_solv_CartesianRobinBcHelper

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CartesianRobinBcHelper.C"
#endif
