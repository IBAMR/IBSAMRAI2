/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/poisson/CellPoissonHypreSolver.h $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 3276 $
 * Modified:    $LastChangedDate: 2009-06-17 10:19:07 -0700 (Wed, 17 Jun 2009) $
 * Description:	Hypre solver interface for diffusion-like elliptic problems.
 */

#ifndef included_solv_CellPoissonHypreSolver
#define included_solv_CellPoissonHypreSolver

#include "SAMRAI_config.h"

#include "tbox/SAMRAI_MPI.h"

#ifdef HAVE_HYPRE

#ifndef included_HYPRE_struct_ls
/*
 * This might break things if F77_FUNC_ is different for hypre vs
 * SAMRAI autoconf detection.  But then C/C++ macros are totally
 * broken due to namespace collision as this example highlights so
 * resorting to hacks are necessary.
 */
#ifdef F77_FUNC_
#undef F77_FUNC_
#endif
#include "HYPRE_struct_ls.h"
#define included_HYPRE_struct_ls
#endif

#include "BoxList.h"
#include "CoarseFineBoundary.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "VariableContext.h"
#include "CellData.h"
#include "SideData.h"
#include "OutersideVariable.h"
#include "GhostCellRobinBcCoefs.h"
#include "RobinBcCoefStrategy.h"
#include "PoissonSpecifications.h"
#include "SimpleCellRobinBcCoefs.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
    namespace solv {

/*!
 * @brief Use the HYPRE preconditioner library to solve (the cell-centered)
 * Poisson's equation on a single level in a hierarchy.
 *
 * Class CellPoissonHypreSolver<DIM> uses the HYPRE preconditioner library
 * to solve linear equations of the form
 * @f$ \nabla ( D \nabla u ) + C u = f @f$, where
 * C is a cell-centered array, D is a face-centered array,
 * and u and f are cell-centered arrays
 * (see PoissonSpecifications).
 * The discretization is the standard second order
 * finite difference stencil.
 *
 * Robin boundary conditions are used through the
 * interface class RobinBcCoefStrategy<DIM>.
 * Periodic boundary conditions are not supported yet.
 *
 * The user must perform the following steps to use
 * CellPoissonHypreSolver<DIM>:
 * - Create a CellPoissonHypreSolver<DIM> object.
 * - Initialize CellPoissonHypreSolver<DIM> object with a patch hierarchy,
 *   using the function initializeSolverState().
 * - Use the functions setPhysicalBcCoefObject() 
 *   to provide implementations of RobinBcCoefStrategy<DIM>.
 *   (For most problems you can probably find a suitable
 *   implementation to use without implementing the
 *   strategy yourself.  See for example
 *   SimpleCellRobinBcCoefs<DIM> and GhostCellRobinBcCoefs<DIM>.)
 * - Set the matrix coefficients in the linear system,
 *   using the function setMatrixCoefficients().
 * - Specify the stopping criteria using setStoppingCriteria().
 * - Solve the linear system, passing in u and f as the patch
 *   indices of the solution and the right hand side, respectively.
 *
 * Sample parameters for initialization from database (and their
 * default values):
 * @verbatim
 *     print_solver_info = FALSE      // Whether to print some data for debugging
 *     max_iterations = 10            // Max iterations used by Hypre
 *     relative_residual_tol = 1.0e-8 // Residual tolerance used by Hypre
 *     num_pre_relax_steps = 1        // # of presmoothing steps used by Hypre
 *     num_post_relax_steps = 1       // # of postsmoothing steps used by Hypre
 *     use_smg = FALSE                // Whether to use hypre's smg solver
 *                                    // (alternative is the pfmg solver)
 * @endverbatim
 */

template<int DIM> class CellPoissonHypreSolver
{
public:
   /*!
    * @brief Constructor.
    *
    * @param object_name Name of object.
    * @param database tbox::Database for input.
    */
   CellPoissonHypreSolver(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> database=NULL );

   /*!
    * The Poisson destructor releases all internally managed data.
    */
   ~CellPoissonHypreSolver();

   /*!
    * @brief Initialize to a given hierarchy.
    *
    * Initializer Poisson solver for a patch level in a hierarchy.
    *
    * @param hierarchy Hierarchy
    * @param ln Level number
    */
   void initializeSolverState(
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy ,
      int ln = 0 );

   /*!
    * @brief Reset to an uninitialized state.
    */
   void deallocateSolverState();

   /*!
    * @brief Set the matrix coefficients
    *
    * For information describing the Poisson equation parameters,
    * see the light-weight PoissonSpecifications class where
    * you set the values of C and D.
    *
    * This method must be called before solveSystem().
    */
   void setMatrixCoefficients( const PoissonSpecifications &spec );

   /*!
    * @brief Set default depth of the solution data involved in the solve.
    *
    * If the solution data has multiple depths,
    * the solver uses just one depth at a time.
    * The default depth is the first depth.
    * Use this function to change it.
    * This is not used to set the depth of the data (which is not
    * controled by this class) but the depth used in the solve.
    *
    * Changing the depth after setting up the matrix is permissible,
    * as the solution data does not affect the matrix.
    */
   void setSolnIdDepth( const int depth );

   /*!
    * @brief Set default depth of the rhs data involved in the solve.
    *
    * If the rhs data has multiple depths,
    * the solver uses just one depth at a time.
    * The default depth is the first depth.
    * Use this function to change it.
    * This is not used to set the depth of the data (which is not
    * controled by this class) but the depth used in the solve.
    *
    * Changing the depth after setting up the matrix is permissible,
    * as the rhs data does not affect the matrix.
    */
   void setRhsIdDepth( const int depth );

   /*!
    * @brief Set the stopping criteria (max iterations and residual
    * tolerance) for the linear solver.
    *
    * @param max_iterations gives the maximum number of iterations
    * @param relative_residual_tol the maximum error tolerance
    */
   void setStoppingCriteria(const int max_iterations = 10,
                            const double relative_residual_tol = 1.0e-6);

   /*!
    * @brief Solve the linear system Au=f.
    *
    * The solution u and the right hand side f are
    * specified via patch indices on the patch hierarchy.
    *
    * Member functions getNumberOfIterations() return the iterations
    * from the solver.
    * Note that the matrix coefficients and boundary condition object
    * must have been set up before this function is called.
    * As long as the matrix coefficients do not change, 
    * this routine may be called repeatedly to solve any number of linear 
    * systems (with the right-hand side varying).
    * If the boundary conditions or matrix coefficients are changed
    * then function setMatrixCoefficients() must be called again.
    *
    * When computing the matrix coefficients in setMatrixCoefficients(),
    * the inhomogeneous portion of the boundary condition (constant
    * terms, independent of u and thus having no effect on the matrix)
    * are saved and added to the source term, f,
    * before performing the matrix solve.  In some situations, it may be
    * useful to not add the inhomogeneous portion to f.  The flag argument
    * @c homoegneous_bc is used for this.  (This is a sort of optimization,
    * to avoid having to re-call setMatrixCoefficients() to change the
    * inhomogeneous portion.)
    *
    * @param u Descriptor of cell-centered unknown variable.
    * @param f Descriptor of cell-centered source variable.
    * @param homogeneous_bc Whether homogeneous boundary conditions
    *        are assumed.
    *
    * @return whether solver converged to specified level
    */
   int solveSystem( const int u ,
                    const int f ,
                    bool homogeneous_bc=false );


   /*!
    * @brief Return the number of iterations taken by the solver to converge.
    *
    * @return number of iterations taken by the solver to converge
    */
   int getNumberOfIterations() const;

   /*!
    * @brief Set the number of pre-relax steps used by the Hypre solve.
    */
   void setNumPreRelaxSteps(const int steps);

   /*!
    * @brief Set the number of post-relax steps used by the Hypre solve.
    */
   void setNumPostRelaxSteps(const int steps);

   /*!
    * @brief Return the final residual norm returned by the Hypre solve.
    * @return final residual norm returned by the Hypre solve.
    */
   double getRelativeResidualNorm() const;

   /*!
    * @brief Set whether to use Hypre's PFMG algorithm instead of the
    * SMG algorithm.
    *
    * The flag is used to select which of HYPRE's linear solver algorithms
    * to use if true, the semicoarsening multigrid algorithm is used, and if
    * false, the "PF" multigrid algorithm is used.
    * By default, the SMG algorithm is used.
    *
    * Changing the algorithm must be done before setting up the matrix
    * coefficients.
    */
   void setUseSMG( bool use_smg );



   /*!
    * @brief Specify boundary condition directly, without using
    * a RobinBcCoefStrategy<DIM> object.
    *
    * Use @em either setBoundaries() @em or setPhysicalBcCoefObject(),
    * but not both.
    *
    * A SimpleCelBcCoef object is used to interpret and implement
    * the specified boundary conditions.
    * See SimpleCellRobinBcCoefs<DIM>::setBoundaries()
    * for an explanation of the arguments.
    */
   void setBoundaries(const std::string& boundary_type,
                      const int fluxes = -1,
                      const int flags = -1,
                      int* bdry_types = NULL);

   /*!
    * @brief Specify boundary condition through the use of a
    * Robin boundary condition object.
    *
    * Use @em either setBoundaries() @em or setPhysicalBcCoefObject(),
    * but not both.
    *
    * The Robin boundary condition object is used when setting
    * the matrix coefficient and when solving the system.
    * If your boundary conditions are fixed values at ghost
    * cell centers, use the GhostCellRobinBcCoefs<DIM>
    * implementation of the RobinBcCoefStrategy<DIM> strategy.
    *
    * @param physical_bc_coef_strategy tbox::Pointer a concrete
    *        implementation of the Robin bc strategy.
    * @param variable hier::Variable pointer to be passed
    *        to RobinBcCoefStrategy<DIM>::setBcCoefs(),
    *        but otherwise unused by this class.
    */
   void setPhysicalBcCoefObject(
      const RobinBcCoefStrategy<DIM> *physical_bc_coef_strategy,
      const tbox::Pointer< hier::Variable<DIM> > variable=NULL );


   /*!
    * @brief Set the flag for printing solver information.
    *
    * This optional function is used primarily for debugging.
    *
    * If set true, it will print the HYPRE matrix information
    * to the following files:
    *
    * - mat_bA.out - before setting matrix coefficients in matrix assemble
    * - mat_aA.out - after setting matrix coefficients in matrix assemble
    * - sol0.out   - u before solve (i.e. for system Au = b)
    * - sol.out    - u after solve
    * - mat0.out   - A before solve
    * - mat.out    - A after solve
    * - rhs.out    - b before and after solve
    *
    * If this method is not called, or the flag is set false, no printing 
    * will occur.
    */
   void setPrintSolverInfo(const bool print);

private:

   /*!
    * @brief Set state using database
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * @param database Input database.  If a NULL pointer is given,
    * nothing is done.
    */
   void getFromInput( tbox::Pointer<tbox::Database> database );

   void setupHypreSolver();
   void destroyHypreSolver();
   void allocateHypreData();
   void deallocateHypreData();

   void copyToHypre(HYPRE_StructVector vector,
                    pdat::CellData<DIM,double> &src,
                    int depth,
                    const hier::Box<DIM> &box);
   void copyFromHypre(pdat::CellData<DIM,double> &dst,
                      int depth,
                      HYPRE_StructVector vector,
                      const hier::Box<DIM> box);


   /*!
    * @brief Add g*A*k0(a) from boundaries to rhs.
    *
    * Move the constant portion of the boundary condition
    * contribution to the right hand side and add it to the existing rhs.
    * This operation is done for physical as well as cf boundaries,
    * so it is placed in a function.
    *
    * The boundary boxes given must be to either the physical
    * boundary or coarse-fine boundary for the patch.  The
    * bc coefficient implementation should correspond to the
    * boundary being worked on.
    */
   void add_gAk0_toRhs( const hier::Patch<DIM> &patch,
                        const tbox::Array< hier::BoundaryBox<DIM> > &bdry_boxes,
                        const RobinBcCoefStrategy<DIM> *robin_bc_coef,
                        pdat::CellData<DIM,double> &rhs );


   //@{

   /*!
    * @name Dimension-independent functions to organize Fortran interface.
    */

   //! @brief Compute diagonal entries of the matrix when C is variable.
   void computeDiagonalEntries(
      pdat::CellData<DIM,double> &diagonal,
      const pdat::CellData<DIM,double> &C_data,
      const pdat::SideData<DIM,double> &variable_off_diagonal,
      const hier::Box<DIM> &patch_box );
   //! @brief Compute diagonal entries of the matrix when C is constant.
   void computeDiagonalEntries(
      pdat::CellData<DIM,double> &diagonal,
      const double C,
      const pdat::SideData<DIM,double> &variable_off_diagonal,
      const hier::Box<DIM> &patch_box );
   //! @brief Compute diagonal entries of the matrix when C is zero.
   void computeDiagonalEntries(
      pdat::CellData<DIM,double> &diagonal,
      const pdat::SideData<DIM,double> &variable_off_diagonal,
      const hier::Box<DIM> &patch_box );
   /*!
    * @brief Adjust boundary entries for variable off-diagonals.
    *
    * At the same time, save information that are needed to adjust
    * the rhs.
    */
   void adjustBoundaryEntries(
      pdat::CellData<DIM,double> &diagonal,
      const pdat::SideData<DIM,double> &variable_off_diagonal,
      const hier::Box<DIM> &patch_box,
      const pdat::ArrayData<DIM,double> &acoef_data,
      const pdat::ArrayData<DIM,double> &bcoef_data,
      const hier::Box<DIM> bccoef_box,
      pdat::ArrayData<DIM,double> &Ak0_data,
      const hier::BoundaryBox<DIM> &trimmed_boundary_box,
      const double h[DIM] );

   //@}

   //! @brief Free static variables at shutdown time.
   static void freeVariables();



   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   /*!
    * @brief Associated hierarchy. 
    */
   tbox::Pointer< hier::PatchHierarchy<DIM> > d_hierarchy;

   /*!
    * @brief Associated level number.
    *
    * Currently, this must be level number 0.
    */
   int d_ln;

   /*!
    * @brief Scratch context for this object.
    */
   tbox::Pointer<hier::VariableContext> d_context;

   //@{ @name Boundary condition handling

   /*!
    * @brief The coarse-fine boundary description for level d_ln.
    *
    * The coarse-fine boundary is computed when the operator
    * state is initialized.  It is used to allow solves on
    * levels that are not the coarsest in the hierarchy.
    */
   hier::CoarseFineBoundary<DIM> d_cf_boundary;

   /*!
    * @brief Robin boundary coefficient object for physical
    * boundaries.
    *
    * If d_physical_bc_coef_strategy is set, use it, otherwise,
    * use d_physical_bc_simple_case.
    */
   const RobinBcCoefStrategy<DIM> *d_physical_bc_coef_strategy;
   tbox::Pointer< hier::Variable<DIM> > d_physical_bc_variable;

   /*!
    * @brief Implementation of Robin boundary conefficients
    * for the case of simple boundary conditions.
    */
   SimpleCellRobinBcCoefs<DIM> d_physical_bc_simple_case;

   /*!
    * @brief Robin boundary coefficient object for coarse-fine
    * boundaries.
    *
    * This is a GhostCellRobinBcCoefs<DIM> object because we
    * expect the users to have the correct ghost cell values
    * in the coarse-fine boundaries before solving.
    */
   GhostCellRobinBcCoefs<DIM> d_cf_bc_coef;
   tbox::Pointer< hier::Variable<DIM> > d_coarsefine_bc_variable;

   //@}

   /*!
    * @brief hier::Patch index of A*k0(a) quantity
    *
    * A*k0(a) is the quantity that is saved for
    * later adding to the rhs.
    *
    * The Robin bc is expressed by the coefficients a and g
    * on the boundary (see RobinBcCoefStrategy).
    * This class uses a central difference approximation of
    * the Robin bc, which results in the value at a ghost cell,
    * uo, being writen as uo = g*k0(a) + k1(a)*ui, where ui is
    * the first interior cell value, k0 and k1 depend on a as
    * indicated.
    *
    * In setting up the Au=f system, the contribution of k1(a)*ui
    * is incorporated into the product Au.  The contribution of
    * A*g*k0(a) should be moved to the right hand side and saved for
    * later adding to f.  However, the value of g is not provided
    * until solve time.  Therefore, we save just A*k0(a) at the
    * patch data index d_Ak0_id.
    */
   int d_Ak0_id;

   static tbox::Pointer<pdat::OutersideVariable<DIM,double> > s_Ak0_var;

   /*!
    * @brief Depth of the solution variable.
    */
   int d_soln_depth;

   /*!
    * @brief Depth of the rhs variable.
    */
   int d_rhs_depth;

   int d_max_iterations;
   double d_relative_residual_tol;


   int d_number_iterations; // iterations in solver
   int d_num_pre_relax_steps;  // pre-relax steps in solver
   int d_num_post_relax_steps; // post-relax steps in solver
   double d_relative_residual_norm;  // norm from solver

   /*@
    * @brief Flag to use SMG or PFMG (default)
    */
   bool d_use_smg;

   //@{
   //! @name Hypre object
   //! @brief HYPRE grid
   HYPRE_StructGrid    d_grid;
   //! @brief HYPRE stencil
   HYPRE_StructStencil d_stencil;
   //! @brief HYPRE structured matrix
   HYPRE_StructMatrix  d_matrix;
   //! @brief Hypre RHS vector for linear solves
   HYPRE_StructVector  d_linear_rhs;
   //! @brief Hypre solution vector
   HYPRE_StructVector  d_linear_sol;
   //! @brief Hypre SMG solver data
   HYPRE_StructSolver  d_mg_data;
   //@}


   //@{

   //! @name Variables for debugging and analysis.

   /*!
    * @brief Flag to print solver info
    *
    * See setPrintSolverInfo().
    */
   bool d_print_solver_info;

   //@}

   /*!
    * @brief Timers for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_solve_system;
   tbox::Pointer<tbox::Timer> t_set_matrix_coefficients;
};

}
} // namespace SAMRAI

#ifndef DEBUG_NO_INLINE
#include "CellPoissonHypreSolver.I"
#endif

#endif

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellPoissonHypreSolver.C"
#endif
