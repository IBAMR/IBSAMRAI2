//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/packages/petsc/SNES_SAMRAIContext.h $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 3281 $
// Modified:    $LastChangedDate: 2009-06-17 19:10:03 -0700 (Wed, 17 Jun 2009) $
// Description: Wrapper for SNES solver for use in a SAMRAI-based application.
//

#ifndef included_solv_SNES_SAMRAIContext
#define included_solv_SNES_SAMRAIContext

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
************************************************************************
*/

#include "tbox/SAMRAI_MPI.h"

#ifdef HAVE_PETSC

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

extern "C" {
#ifdef PETSC2028
#include "snes.h"
#else
#include "petscsnes.h"
#endif
}

#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#endif

#include "SAMRAI_config.h"

#ifdef HAVE_PETSC

#include "NonlinearSolverStrategy.h"
#include "SNESAbstractFunctions.h"
#include "tbox/Database.h"
#include "tbox/Serializable.h"


namespace SAMRAI {
    namespace solv {

/*!
 * Class SNES_SAMRAIContext<DIM> provides an interface to the SNES
 * nonlinear solver capabilities in PETSc to facilitate their use with
 * SAMRAI.  While PETSc is implemented in an object-based manner, this
 * class makes it easier to use PETSc's routines with SAMRAI data structures.
 * In particular, it hides from the user some of the messy details required
 * to link the C++ class library with PETSc which is written in C and to
 * override PETSc objects, like vectors, with those found in SAMRAI.
 *
 * This class declares five private static member functions to link
 * user-defined routines for nonlinear residual calculation, Jacobian
 * evaluation, preconditioner setup and solve, and Jacobian-vector product 
 * operations.  The implementation of these functions is defined by the user 
 * in a subclass of the abstract base class SNESAbstractFunctions.  
 * The vector objects used within the solver are provided by the 
 * PETSc_SAMRAIVectorReal<DIM> wrapper class.
 *
 * If no parameters are read from input, PETSc defaults are used.  See the
 * PETSc documentation (http://www-unix.mcs.anl.gov/petsc/).
 * for more information on default parameters and SNES functionality. 
 *
 * Optional input keys and types that can be read by this class are:
 * 
 *    - absolute_tolerance
 *        double value for absolute nonlinear convergence tolerance
 *
 *    - relative_tolerance
 *        double value for relative nonlinear convergence tolerance
 *
 *    - step_tolerance
 *        double value for minimum tolerance on change in solution norm 
 *        between nonlinear iterates
 *
 *    - maximum_nonlinear_iterations
 *        integer value for maximum number of nonlinear iterations
 *
 *    - maximum_function_evals
 *        integer value for maximum number of nonlinear function evaluations 
 *
 *    - forcing_term_strategy
 *        integer value for forcing term choice for linear solvers 
 *        within the inexact Newton method.  Choices are "CONSTANT" (default),
 *        "EWCHOICE1", and "EWCHOICE2".
 *
 *    - constant_forcing_term
 *        double value for constant relative convergence tolerance in
 *        Krylov solver (default case)
 *
 *    - initial_forcing_term
 *        double value for initial relative convergence tolerance in
 *        Krylov solver (used in Eisenstat-Walker case). Value must
 *        satisfy @f$0 \le \eta_0 < 1@f$.
 *
 *    - maximum_forcing_term
 *        double value for maximum relative convergence tolerance in
 *        Krylov solver (used in Eisenstat-Walker case). Value must
 *        satisfy @f$0 \le \eta_{max} < 1@f$.
 *
 *    - EW_choice2_alpha
 *        double value for power used in Eisenstat-Walker choice 2 
 *        relative convergence tolerance computation. Value must
 *        satisfy @f$1 < \alpha \le 2@f$.
 *
 *    - EW_safeguard_exponent
 *        double value for power for safeguard used in Eisenstat-Walker 
 *        choice 2
 *
 *    - EW_choice2_gamma
 *        double value for multiplicative factor used in Eisenstat-Walker
 *        choice 2 relative convergence tolerance computation.  Value
 *        must satisfy @f$0 \le \gamma \le 1@f$.
 *
 *    - EW_safeguard_disable_threshold
 *        double value for threshold for imposing safeguard in 
 *        Eisenstat-Walker choice 2.  Value must 
 *        satisfy @f$0 < \eta_{threshold} < 1@f$. 
 *
 *    - linear_solver_type
 *        string value for type of linear solver.  See the
 *        <A HREF="http://www-unix.mcs.anl.gov/petsc/">PETSc documentation</A>
 *        for the valid types.
 *
 *    - uses_preconditioner
 *        boolean value for whether or not a preconditioner is used in the
 *        solution of the Newton equations.
 *
 *    - linear_solver_absolute_tolerance
 *        double value for absolute convergence tolerance in linear solver
 *     
 *    - linear_solver_divergence_tolerance
 *        double value for amount linear solver residual can increase 
 *        before solver concludes method is diverging
 *     
 *    - maximum_linear_iterations
 *        integer value for maximum number of linear solver iterations
 *
 *    - maximum_gmres_krylov_dimension
 *        integer value for maximum dimension of Krylov subspace before
 *        restarting.  Valid only if GMRES is used as the linear solver.
 *
 *    - gmres_orthogonalization_algorithm
 *        string value for algorithm used to incrementally construct the
 *        orthonormal basis of the Krylov subspace used by GMRES.  Valid 
 *        only if GMRES is used as the linear solver.  Valid values are:
 *
 *             modifiedgramschmidt
 *             gmres_cgs_refine_ifneeded
 *             gmres_cgs_refine_always
 * 
 *        See the
 *        <A HREF="http://www-unix.mcs.anl.gov/petsc/">PETSc documentation</A>
 *        for more information.
 *
 *    - uses_explicit_jacobian
 *        boolean value for whether or not the user provides code to 
 *        explicitly calculate Jacobian-vector products.
 *
 *    - differencing_parameter_strategy
 *        string value indicating strategy used for computing the differencing
 *        parameter when Jacobian-vector products are approximated via finite
 *        differences.  See the
 *        <A HREF="http://www-unix.mcs.anl.gov/petsc/">PETSc documentation</A>
 *        for valid values.
 *
 *    - function_evaluation_error
 *        double value that is the square root of the estimated relative error 
 *        in function evaluation.
 *
 * Note that all input values may override values read in from restart.  If
 * no new input value is given, the restart value is used.
 *
 * A sample input file entry might look like:
 *
 * @code
      absolute_tolerance             = 10.e-10
      relative_tolerance             = 10.e-6
      step_tolerance                 = 10.e-8
      maximum_nonlinear_iterations   = 200
      forcing_term_strategy          = "EWCHOICE1"
   @endcode
 *
 * Note that input values can also be set using accessor functions.
 * Values that are set via this mechanism will be cached both in the
 * solver context as well as in the corresponding PETSc object.  Thus
 * values changed on-the-fly will be written to restart.  These input
 * values can also be changed by directly accessing the corresponding
 * PETSc object and using native PETSc function calls; however such
 * settings/changes will NOT be cached in the solver context, and so
 * will not be written to restart.
 *
 * @see solv::SNESAbstractFunctions
 * @see solv::NonlinearSolverStrategy
 */

template<int DIM> class SNES_SAMRAIContext 
: 
public NonlinearSolverStrategy<DIM>,
public tbox::Serializable
{
public:
   /*!
    * Constructor for SNES_SAMRAIContext<DIM> allocates the SNES 
    * object and initializes rudimentary state associated with 
    * user-supplied solver components.  Then, it reads solver parameter
    * from input and restart which may override default values.
    *
    * When assertion checking is active, an unrecoverable exception
    * will result if the name string is empty or the pointer to the
    * user-defined SNES functions object is null. 
    */
   SNES_SAMRAIContext(const std::string& object_name,
                            tbox::Pointer<tbox::Database> input_db,
                            SNESAbstractFunctions* my_functions);

   /*!
    * Destructor for solve_SNES_SAMRAIContext destroys the SNES 
    * and the PETSc solution vector wrapper.
    */
   ~SNES_SAMRAIContext<DIM>();

   /*!
    * Return the PETSc nonlinear solver object.
    */
   SNES getSNESSolver() const;

   /*!
    * Return pointer to object providing user-defined functions for SNES.
    */
   SNESAbstractFunctions* getSNESFunctions() const;

   /*!
    * Return the PETSc linear solver object.
    */
//   SLES getSLESSolver() const;

   /*!
    * Return the PETSc Krylov solver object.
    */
   KSP getKrylovSolver() const;

   /*!
    * Return the PETSc Mat object for the Jacobian.
    */
   Mat getJacobianMatrix() const;

   /*!
    *  Get absolute tolerance for nonlinear solver.
    */
   double getAbsoluteTolerance() const;

   /*!
    *  Set absolute tolerance for nonlinear solver.
    */
   void setAbsoluteTolerance(double abs_tol);

   /*!
    *  Get relative tolerance for nonlinear solver.
    */
   double getRelativeTolerance() const;

   /*!
    *  Set relative tolerance for nonlinear solver.
    */
   void setRelativeTolerance(double rel_tol);

   /*!
    *  Get step tolerance for nonlinear solver.
    */
   double getStepTolerance() const;

   /*!
    *  Set step tolerance for nonlinear solver.
    */
   void setStepTolerance(double step_tol);

   /*!
    *  Get maximum iterations for nonlinear solver.
    */
   int getMaxNonlinearIterations() const;

   /*!
    *  Set maximum iterations for nonlinear solver.
    */
   void setMaxNonlinearIterations(int max_nli);

   /*!
    *  Get maximum function evaluations by nonlinear solver.
    */
   int getMaxFunctionEvaluations() const;

   /*!
    *  Set maximum function evaluations in nonlinear solver.
    */
   void setMaxFunctionEvaluations(int max_feval);

   /*!
    *  Get strategy for forcing term.
    */
   std::string getForcingTermStrategy() const;

   /*!
    *  Set strategy for forcing term.
    */
   void setForcingTermStrategy(std::string& strategy);

   /*!
    *  Get value of constant forcing term.
    */
   double getConstantForcingTerm() const;

   /*!
    *  Set value of constant forcing term.
    */
   void setConstantForcingTerm(double fixed_eta);

   /*!
    *  Get value of initial forcing term.
    */
   double getInitialForcingTerm() const;

   /*!
    *  Set value of initial forcing term.
    */
   void setInitialForcingTerm(double initial_eta);

   /*!
    *  Get value of maximum forcing term.
    */
   double getMaximumForcingTerm() const;

   /*!
    *  Set value of maximum forcing term.
    */
   void setMaximumForcingTerm(double max_eta);

   /*!
    *  Get value of exponent in Eisenstat-Walker Choice 2 forcing term.
    */
   double getEWChoice2Exponent() const;

   /*!
    *  Set value of exponent in Eisenstat-Walker Choice 2 forcing term.
    */
   void setEWChoice2Exponent(double alpha);

   /*!
    *  Get value of exponent in Eisenstat-Walker Choice 2 safeguard.
    */
   double getEWChoice2SafeguardExponent() const;

   /*!
    *  Set value of exponent in Eisenstat-Walker Choice 2 safeguard.
    */
   void setEWChoice2SafeguardExponent(double beta);

   /*!
    * Get value of factor used to scale Eisenstat-Walker Choice 2
    * forcing term.
    */
   double getEWChoice2ScaleFactor() const;

   /*!
    * Set value of factor used to scale Eisenstat-Walker Choice 2
    * forcing term.
    */
   void setEWChoice2ScaleFactor(double gamma);

   /*!
    *  Get value of threshold to disable safeguard in Eisenstat-Walker 
    *  forcing terms.
    */
   double getEWSafeguardThreshold() const;

   /*!
    *  Set value of threshold to disable safeguard in Eisenstat-Walker 
    *  forcing terms.
    */
   void setEWSafeguardThreshold(double threshold);

   /*!
    * Get type of linear solver.
    */
   std::string getLinearSolverType() const;

   /*!
    * Set type of linear solver.
    */
   void setLinearSolverType(std::string& type);

   /*!
    * Get whether a preconditioner is used.
    */
   bool getUsesPreconditioner() const;

   /*!
    * Set whether to use a preconditioner.
    */
   void setUsesPreconditioner(bool uses_preconditioner);

   /*!
    * Get absolute tolerance for linear solver.
    */
   double getLinearSolverAbsoluteTolerance() const;

   /*!
    * Set absolute tolerance for linear solver.
    */
   void setLinearSolverAbsoluteTolerance(double abs_tol);

   /*!
    * Get divergence tolerance for linear solver.
    */
   double getLinearSolverDivergenceTolerance() const;

   /*!
    * Set divergence tolerance for linear solver.
    */
   void setLinearSolverDivergenceTolerance(double div_tol);

   /*!
    * Get maximum linear iterations for linear solver.
    */
   int getMaximumLinearIterations() const;

   /*!
    * Set maximum linear iterations for linear solver.
    */
   void setMaximumLinearIterations(int max_li);

   /*!
    * Get maximum Krylov dimension in GMRES linear solver.
    */
   int getMaximumGMRESKrylovDimension() const;

   /*!
    * Set maximum Krylov dimension in GMRES linear solver.
    */
   void setMaximumGMRESKrylovDimension(int d);

   /*!
    * Get orthogonalization method used in GMRES linear solver.
    */
   std::string getGMRESOrthogonalizationMethod() const;

   /*!
    * Set orthogonalization method used in GMRES linear solver.
    */
   void setGMRESOrthogonalizationMethod(std::string& method);

   /*!
    * Get whether a method for explicit Jacobian-vector products is provided.
    */
   bool getUsesExplicitJacobian() const;

   /*!
    * Set whether a method for explicit Jacobian-vector products is provided.
    */
   void setUsesExplicitJacobian(bool use_jac);

   /*!
    * Get method for computing differencing parameter.
    */
   std::string getDifferencingParameterMethod() const;

   /*!
    * Set method for computing differencing parameter.
    */
   void setDifferencingParameterMethod(std::string& method);

   /*!
    * Get estimate of error in function evaluation.
    */
   double getFunctionEvaluationError() const;

   /*!
    * Set estimate of error in function evaluation.
    */
   void setFunctionEvaluationError(double evaluation_error);

   /*!
    * Initialize the state of the SNES solver based on vector argument 
    * representing the solution of the nonlinear system.  In general, this 
    * routine must be called before the solve() routine is invoked.
    */
   void initialize(tbox::Pointer< SAMRAIVectorReal<DIM,double> > solution);

   /*!
    *  Reset the state of the nonlinear solver after regridding.
    */
   void resetSolver(const int coarsest_level, const int finest_level);

   /*!
    * Solve the nonlinear problem.  In general, the initialize() routine
    * must be called before this solve function to set up the solver.
    * Returns 1 if successful, 0 otherwise.  
    */
   int solve();

   /*!
    * Obtain number of nonlinear iterations.  
    */
   int getNonlinearIterationCount() const;

   /*!
    * Obtain total number of linear iterations accumulated over all
    * nonlinear iterations.  
    */
   int getTotalLinearIterationCount() const;

   /*!
    * Report reason for termination.
    */
   void reportCompletionCode(std::ostream & os = tbox::plog) const;

   /*!
    * Write solver parameters to restart database matching object name.
    *
    * When assertion checking is active, an unrecoverable exception
    * will result if database pointer is null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   /*!
    * Print out all members of integrator instance to given output stream.
    */
   virtual void printClassData(std::ostream& os) const;

private:
   /*
    * Static member functions for linkage with PETSc routines.
    * See header file for SNESAbstractFunctions for more information.
    */
   static int SNESFuncEval(SNES snes,       // SNES context
                           Vec x,           // input vector     
                           Vec f,           // residual vector
                           void* ctx);      // user-defined context

   static int SNESJacobianSet(SNES snes,              // SNES context
                              Vec x,                  // input vector
                              Mat* A,                 // Jacobian matrix
                              Mat* B,                 // precond matrix
                              MatStructure* mstruct,  // precond matrix structure
                              void* ctx);             // user-defined context

   static int SNESJacobianTimesVector(Mat M,     //  Jacobian matrix
                                      Vec xin,   //  input vector
                                      Vec xout); //  output vector

   static int SNESsetupPreconditioner(void* ctx);  // input vector

   static int SNESapplyPreconditioner(void* ctx,  // user-defined context
                                      Vec xin,    // input vector
                                      Vec xout);  // output vector

   /*
    * Create and cache needed Petsc objects.
    */
   void createPetscObjects();

   /*
    * Initialize cached Petsc objects.
    */
   void initializePetscObjects();

   /*
    * Destroy cached Petsc objects.
    */
   void destroyPetscObjects();

   /*
    * Read input values from given database.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db); 

   /*!
    * Read solver parameters from restart database matching object name.
    */
   void getFromRestart();

   /*!
    * Internal state parameters:
    */
   std::string d_object_name;
   bool d_context_needs_initialization;

   /*
    * PETSc solver and preconditioner objects:
    *
    * d_SNES_solver ..... PETSc nonlinear solver object.
    *
    * d_linear_solver ... PETSc linear solver object, cached here so that
    *                     users may manipulate it through the
    *                     interface of this class.
    *
    * d_krylov_solver ... PETSc Krylov solver context.
    *
    * d_jacobian ........ PETSc matrix object, cached here so that users
    *                     may specify Jacobian operations through the
    *                     interface of this class without having to know 
    *                     about PETSc matrix shells or matrix-free matrices.
    *
    * d_preconditioner... PETSc preconditioner object, cached here so that
    *                     users may specify operations through the
    *                     interface of this class without having to know
    *                     about PETSc preconditioner shells.
    */

   SNES d_SNES_solver;
   KSP  d_krylov_solver;
   Mat  d_jacobian;
   PC   d_preconditioner;

   /*
    * Solution and residual vectors for nonlinear system.
    */

   Vec d_solution_vector;
   Vec d_residual_vector;

   /*
    * tbox::Pointer to object which provides user-supplied functions to SNES.
    */
   SNESAbstractFunctions* d_SNES_functions;

   /*
    * Boolean flags used during SNES initialization to provide correct
    * static function linkage with PETSc.
    */
   bool d_uses_preconditioner;
   bool d_uses_explicit_jacobian;

   /*
    * SNES state data maintained here for input/restart capabilities.
    */

   // Nonlinear solver parameters:

   int    d_maximum_nonlinear_iterations;
   int    d_maximum_function_evals;

   double d_absolute_tolerance;
   double d_relative_tolerance;
   double d_step_tolerance;

   std::string d_forcing_term_strategy;  // string is for input
   int    d_forcing_term_flag;      // int is for passing choice to PETSc

   double d_constant_forcing_term;
   double d_initial_forcing_term;
   double d_maximum_forcing_term;
   double d_EW_choice2_alpha;
   double d_EW_choice2_gamma;
   double d_EW_safeguard_exponent;
   double d_EW_safeguard_disable_threshold;

   SNESConvergedReason d_SNES_completion_code;

   // Linear solver parameters:

   std::string d_linear_solver_type;
   double d_linear_solver_absolute_tolerance;
   double d_linear_solver_divergence_tolerance;
   int    d_maximum_linear_iterations;

   int    d_maximum_gmres_krylov_dimension;
   std::string d_gmres_orthogonalization_algorithm;
 
   // "Matrix-free" parameters:

   std::string d_differencing_parameter_strategy; 
   double d_function_evaluation_error;

   // Output parameters:

   int d_nonlinear_iterations;
};

}
}
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SNES_SAMRAIContext.C"
#endif
