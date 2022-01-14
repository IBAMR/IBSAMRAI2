//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/packages/sundials/cvode/CVODESolver.h $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Wrapper class for CVODE solver function calls and data
//

#ifndef included_solv_CVODESolver
#define included_solv_CVODESolver

#include "SAMRAI_config.h"

#ifdef HAVE_SUNDIALS

#ifndef included_String
#include <string>
#define included_String
#endif

#include "tbox/IOStream.h"
#include "CVODEAbstractFunctions.h"
#include "SundialsAbstractVector.h"

// CVODE includes
#ifndef included_cvode_h
#define included_cvode_h
extern "C" {
#include "cvode/cvode.h"
}
#endif


namespace SAMRAI {
   namespace solv {

#ifndef LACKS_SSTREAM
#define CVODE_SAMRAI_ERROR(ierr) do {						\
      if (ierr != CV_SUCCESS) {                                   				\
         std::ostringstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#else
#define CVODE_SAMRAI_ERROR(ierr) do {						\
      if (ierr != CV_SUCCESS) {                                   				\
         std::ostrstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
      } 									\
} while (0)
#endif



/*!
 * @brief Class CVODESolver serves as a C++ wrapper for the CVODE 
 * ordinary differential equation solver package.  
 *
 * It is intended to be 
 * sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares one private static member function to link the
 * user-defined routine for right-hand side function evaluation and
 * two private statice member functions to link the user-defined 
 * preconditioner setup and solve routines.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base 
 * class CVODEAbstractFunctions.  The vector objects used within the 
 * solver are given in a subclass of the abstract class 
 * SundialsAbstractVector. The SundialsAbstractVector 
 * class defines the vector kernel operations required by the CVODE
 * package so that they may be easily supplied by a user who opts not 
 * to use the vector kernel supplied by the CVODE package.  (It should be
 * noted that the vector kernel used by CVODE is the same as the one
 * used by the other packages in the Sundials of solvers).
 * 
 * Note that this class provides no input or restart capabilities and 
 * relies on CVODE for output reporting.  
 *
 * CVODESolver Usage:
 * 
 *
 *    -  In order to use the CVODESolver, the user must provide a
 *           concrete subclass of CVODEAbstractFunctions abstract
 *           base class which defines the evaluateRHSFunction(),
 *           CVSpgmrPrecondSet(), and CVSpgmrPrecondSolve() methods.
 *
 *    -  Solving a system of ODEs using this CVODE C++ interface 
 *           requires four main stages.  First, a CVODESolver
 *           object is created with a user-specified name and
 *           CVODEAbstractFunctions object.  Second, the
 *           user must specify the integration parameters that s/he
 *           wishes to use.  Next, the user must call the CVODESolver
 *           method initialize(solution_vector) with the 
 *           SundialsAbstractVector that s/he wants to put the solution 
 *           in.  Finally, the solve() method is invoked to solve the 
 *           system of ODEs to the specified value of the independent 
 *           variable.
 *           
 *    -  The following is a list of integration parameters that
 *           must be specified by the user before calling the solve()
 *           method:
 *        
 * 
 *            - Either relative or absolute tolerance must
 *                  be set - setRelativeTolerance(relative_tolerance),
 *                  setAbsoluteTolerance(absolute_tolerance)
 * 
 *            - Initial value of independent variable -
 *                  setInitialValueOfIndependentVariable(init_time)
 *            - Final value of independent variable -
 *                  setFinalValueOfIndependentVariable(final_time
 *                      cvode_needs_initialization)
 *            - Initial condition vector - 
 *                  setInitialConditionVector(ic_vector)
 *
 *        
 *           
 *    -  The following is a list of default values for integration
 *           parameters:
 * 
 *        
 * 
 *           - @b Linear Multistep Method               
 *                BDF
 * 
 *           - @b Iteration Type                        
 *                FUNCTIONAL
 * 
 *           - @b Tolerance Type                        
 *                SS (scalar relative and scalar absolute tolerances)
 * 
 *           - @b Relative Tolerance                    
 *                0.0
 * 
 *           - @b Scalar Absolute Tolerance             
 *                0.0
 * 
 *           - @b Vector Absolute Tolerance             
 *                NULL
 * 
 *           - @b Stepping Method                       
 *                NORMAL
 * 
 *           - @b Maximum Order for Multistep Method    
 *                12 for ADAMS, 5 for BDF
 * 
 *           - @b Maximum Number of Internal Steps      
 *                500
 * 
 *           - @b Maximum Number of NIL Step Warnings   
 *                10
 * 
 *           - @b Initial Step Size                     
 *                determined by CVODE
 * 
 *           - @b Maximum Absolute Value of Step Size   
 *                infinity
 * 
 *           - @b Minimum Absolute Value of Step Size   
 *                0.0
 * 
 *           - @b CVSpgmr Preconditioning Type          
 *                NONE
 * 
 *           - @b CVSpgmr Gram Schmidt Algorithm        
 *                MODIFIED_GS
 * 
 *           - @b CVSpgmr Maximum Krylov Dimension      
 *                MIN(num_equations, CVSPGMR_MAXL=5)
 * 
 *           - @b CVSpgmr Tolerance Scale Factor        
 *                CVSPGMR_DELT = 0.05.
 * 
 *        
 *
 * 
 *
 * CVODE was developed in the Center for Applied Scientific Computing (CASC) 
 * at Lawrence Livermore National Laboratory (LLNL).  Many of the comments
 * in this class were taken verbatim from CVODE header files.  For more 
 * information about CVODE and a complete description of the operations 
 * and data structures used by this class, see S.D. Cohen and A.C. Hindmarsh, 
 * "CVODE User Guide", UCRL-MA-118618, Lawrence Livermore National 
 * Laboratory, 1994. 
 *
 * @see solv::CVODEAbstractFunctions
 * @see solv::SundialsAbstractVector 
 */

class CVODESolver
{
public:
   /**
    * Constructor for CVODESolver sets default CVODE parameters 
    * and initializes the solver package with user-supplied functions
    * CVODESolver parameters may be changed later using member 
    * functions described below.  
    *
    * Notes:
    * 


    *
    *    -
    *        The solution vector is not passed into the constructor.
    *        Before the solver can be used, the initialize() function must 
    *        be called.
    *
    * 


    *
    * Assertion checks: 
    * 


    *
    *    -
    *        my_functions must not be null
    * 
    *    -
    *        object_name must not be empty.
    * 
    * 


    * 
    */
   CVODESolver(const std::string& object_name,
                    CVODEAbstractFunctions* my_functions,
                    const bool uses_preconditioner);

   /**
    * Virtual destructor for CVODESolver closes the 
    * CVODE log file and frees the memory allocated for the 
    * CVODE memory record.
    */
   virtual ~CVODESolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * CVODE.  This routine must be called before the solver can be used.
    * 
    * Assertion checks:
    * 


    *
    *    -
    *        the solution vector must not be null
    * 
    *    -
    *        the solution vector must not have already been set
    * 
    * 


    */
   void initialize(SundialsAbstractVector* solution); 

   /**
    * Integrate ODE system specified t_f.  The integer return value is  
    * a termination code defined by CVODE.  The following is a table
    * of termination codes and a brief description of their meanings.
    * 
    * CVODE Termination Codes:
    * 


    *
    *    - @b SUCCESS (=0)            
    *        CVode succeeded.
    *
    *    - @b CVODE_NO_MEM (=-1)    
    *        The cvode_mem argument was NULL.
    *
    *    - @b ILL_INPUT (=-2)        
    *        One of the inputs to CVode is illegal. This    
    *        includes the situation when a component of the 
    *        error weight vectors becomes < 0 during
    *        internal time-stepping. The ILL_INPUT flag 
    *        will also be returned if the linear solver 
    *        routine CV--- (called by the user after
    *        calling CVodeMalloc) failed to set one of the
    *        linear solver-related fields in cvode_mem or 
    *        if the linear solver's init routine failed. In
    *        any case, the user should see the printed   
    *        error message for more details.
    *
    *    - @b TOO_MUCH_WORK (=-3)   
    *        The solver took maxstep internal steps but 
    *        could not reach t_f. The default value for  
    *        mxstep is MXSTEP_DEFAULT = 500.
    *
    *    - @b TOO_MUCH_ACC (=-4)    
    *        The solver could not satisfy the accuracy 
    *        demanded by the user for some internal step. 
    *
    *    - @b ERR_FAILURE (=-5)      
    *        Error test failures occurred too many times 
    *        (= MXNEF = 7) during one internal time step or 
    *        occurred with |h| = hmin.
    *
    *    - @b CONV_FAILURE (=-6)     
    *        Convergence test failures occurred too many  
    *        times (= MXNCF = 10) during one internal time
    *        step or occurred with |h| = hmin.
    *
    *    - @b SETUP_FAILURE (=-7)    
    *        The linear solver's setup routine failed in an
    *                 unrecoverable manner.
    *
    *    - @b SOLVE_FAILURE (=-8)    
    *        The linear solver's solve routine failed in an 
    *                 unrecoverable manner.
    * 
    * 


    * 
    * See cvode.h header file for more information about return values.  
    *
    * If CVODE or CVSpgmr requires re-initialization, it is 
    * automatically done before the solve.  This may be required if any 
    * of the CVODE or CVSpgmr data parameters have changed since the 
    * last call to the solver.  
    *
    * Assertion checks:
    * 


    *
    *     - 
    *        The user specified final value for the independent variable t
    *        must be greater than the specified initial value.
    *
    * 


    */
   int solve();

   /**
    * Accessor function for setting CVODE output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    *
    * If the file name string is empty the default file name "cvode.log" 
    * is used.
    */
   void setLogFileData(const std::string& log_fname = std::string());

   /**
    * Set CVODESolver to use my_functions as the concrete subclass 
    * of the CVODEAbstractFunctions class that defines the 
    * right-hand side evaluation and preconditioner functions.  The
    * uses_preconditioner argument indicates whether or not the
    * the user has defined preconditioner routines in their concrete
    * subclass of the CVODEAbstractFunctions class.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        my_function must not be a null pointer
    *
    * 


    */
   void setCVODEFunctions(CVODEAbstractFunctions* my_functions,
                          const bool uses_preconditioner);

   /**
    * Return pointer to object that provides user-defined functions for 
    * CVODE and CVSpgmr.
    */
   CVODEAbstractFunctions* getCVODEFunctions() const;

   // Methods for setting CVODE parameters.

   /**
    * Set linear multistep method.  The user can specify either
    * ADAMS or BDF (backward differentiation formula) methods 
    * The BDF method is recommended  for stiff problems, and 
    * the ADAMS method is recommended for nonstiff problems.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        linear_multistep_method must be one of ADAMS or BDF.
    *
    * 


    *
    * Note: the enumeration constants ADAMS and BDF are defined in cvode.h.
    */
   void setLinearMultistepMethod(int linear_multistep_method);

   /**
    * Set iteration type.  The user can specify either FUNCTIONAL
    * iteration, which does not require linear algebra, or a 
    * NEWTON iteration, which requires the solution of linear 
    * systems. In the NEWTON case, the user must also specify a 
    * CVODE linear solver. NEWTON is recommended in case of 
    * stiff problems.  
    *
    * Assertion checks:
    * 


    *
    *    -
    *        iteration_type must be one of FUNCTIONAL or NEWTON
    *
    * 


    *
    * Note: the enumeration constants FUNCTIONAL and NEWTON are defined 
    * in cvode.h.
    */
   void setIterationType(int iteration_type);

   /**
    * Set tolerance type.  This parameter specifies the relative 
    * and absolute tolerance types to be used. The SS tolerance type 
    * means a scalar relative and absolute tolerance, while the SV 
    * tolerance type means a scalar relative tolerance and a 
    * vector absolute tolerance (a potentially different 
    * absolute tolerance for each vector component).    
    *
    * Assertion checks:
    * 


    *
    *    -
    *        tolerance_type must be one of SS or SV
    *
    * 


    *
    * Note: the enumeration constants SS and SV are defined in cvode.h.
    */
   void setToleranceType(int tolerance_type);

   /**
    * Set the relative tolerance level.  
    *
    * Assertion checks:
    * 


    *
    *    -
    *        relative_tolerance must be greater than or equal to 0.0
    *
    * 


    *
    * Note that pure absolute tolerance can be used by
    * setting the relative tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setRelativeTolerance(double relative_tolerance);

   /**
    * Set the scalar absolute tolerance level.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        absolute_tolerance must be greater than or equal to 0.0
    *
    * 


    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setAbsoluteTolerance(double absolute_tolerance);

   /**
    * Set the vector absolute tolerance level.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        absolute_tolerance must not be a null pointer
    *
    *    -
    *        each component of absolute_tolerance must be 
    *        greater than or equal to 0.0
    *
    * 


    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setAbsoluteTolerance(SundialsAbstractVector* absolute_tolerance);

   /**
    * Set stepping method to use for integration.  There are 
    * stepping methods: NORMAL and ONE_STEP.  The NORMAL
    * method has the solver take internal steps until 
    * it has reached or just passed the user specified t_f
    * parameter. The solver then interpolates in order to 
    * return an approximate value of y(t_f). The ONE_STEP 
    * option tells the solver to just take one internal step 
    * and return the solution at the point reached by that 
    * step.                       
    *
    * Assertion checks:
    * 


    *
    *    -
    *        stepping_method must be one of NORMAL or ONE_STEP
    *
    * 


    *
    * Note: the enumeration constants NORMAL and ONE_STEP are 
    * defined in cvode.h.
    */
   void setSteppingMethod(int stepping_method);

   /**
    * Set initial value for independent variable.
    */
   void setInitialValueOfIndependentVariable(double t_0);

   /**
    * Set final value for independent variable (i.e. the value of
    * independent variable to integrate the system to).  The boolean
    * argument specifies whether CVODE should be re-initialized (i.e.
    * on first step) or if we are taking subsequent steps in a 
    * sequence, in which case it is not initialized.
    */
   void setFinalValueOfIndependentVariable(double t_f,
      bool cvode_needs_initialization);

   /**
    * Set initial condition vector.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        ic_vector must not be null
    *
    * 


    */
   void setInitialConditionVector(SundialsAbstractVector* ic_vector);

   /**
    * Set maximum order for the linear multistep method.
    * By default, this is set to 12 for ADAMS methods and 5 for BDF 
    * methods.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_order must be greater than or equal to 0
    *
    * 


    */
   void setMaximumLinearMultistepMethodOrder(int max_order);

   /**
    * Set maximum number of internal steps to be taken by 
    * the solver in its attempt to reach t_f.
    * By default, this is set to 500.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_num_internal_steps must be greater than or equal to 0
    *
    * 


    */
   void setMaximumNumberOfInternalSteps(int max_num_internal_steps);

   /**
    * Set maximum number of warning messages issued by the solver
    * that (t + h == t) on the next internal step.  By default, 
    * this is set to 10.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_num_warnings must be greater than or equal to 0
    *
    * 


    */
   void setMaximumNumberOfNilStepWarnings(int max_num_warnings);

   /**
    * Set initial step size.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        init_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setInitialStepSize(double init_step_size);

   /**
    * Set maximum absolute value of step size allowed.
    * By default, there is no upper bound on the absolute value
    * of step size.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setMaximumAbsoluteStepSize(double max_step_size);

   /**
    * Set minimum absolute value of step size allowed.
    * By default, this is set to 0.0.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        min_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setMinimumAbsoluteStepSize(double min_step_size);

   // Methods for setting CVSpgmr parameters.

   /**
    * Set the preconditioning type to be used by CVSpgmr.
    * This must be one of the four enumeration constants
    * NONE, LEFT, RIGHT, or BOTH defined in iterativ.h.
    * These correspond to no preconditioning, left preconditioning only,
    * right preconditioning only, and both left and right
    * preconditioning, respectively.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    precondition_type must be one of NONE, LEFT, RIGHT, or BOTH.
    *
    * 


    */
   void setPreconditioningType(int precondition_type);

   /**
    * Set the Gram-Schmidt orthogonalization type to be used by CVSpgmr.
    * This must be one of the two enumeration constants MODIFIED_GS
    * or CLASSICAL_GS defined in iterativ.h. These correspond to
    * using modified Gram-Schmidt and classical Gram-Schmidt, respectively.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    gs_type must be one of CLASSICAL_GS or MODIFIED_GS.
    *
    * 


    */
   void setGramSchmidtType(int gs_type);

   /**
    * Set the maximum Krylov dimension to be used by CVSpgmr.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value MIN(num_equations, CVSPGMR_MAXL=5).
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    max_krylov_dim must be nonnegative
    *
    * 


    */
   void setMaxKrylovDimension(int max_krylov_dim);

   /**
    * Set the factor by which the tolerance on the nonlinear
    * iteration is multiplied to get a tolerance on the linear iteration.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value CVSPGMR_DELT = 0.05.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    tol_scale_factor must be nonnegative
    *
    * 


    */
   void setCVSpgmrToleranceScaleFactor(double tol_scale_factor);

   /**
    * Get solution vector.
    */
   SundialsAbstractVector* getSolutionVector() const;

   /**
    * Get k-th derivative vector at the specified value of the
    * independent variable, t.  The integer return value is 
    * return code the CVODE CVodeDky() function.  The following is a table
    * of termination codes and a brief description of their meanings.
    * 
    * CVodeDky Return Codes:
    * 


    *
    *    - @b OKAY (=0)           
    *        CVodeDky succeeded.
    *
    *    - @b BAD_K (=-1)        
    *        
    *    - @b BAD_T (=-2)        
    *        
    *    - @b BAD_DKY (=-3)      
    *        
    *    - @b DKY_NO_MEM (=-4)  
    *        
    * 


    * 
    * Important Notes:
    * 


    *
    *    -
    *       t must lie in the interval [t_cur - h, t_cur] 
    *       where t_cur is the current internal time reached
    *       and h is the last internal step size successfully 
    *       used by the solver.
    *
    *    -
    *       k may take on value 0, 1, . . . q where q is the order 
    *       of the current linear multistep method being used. 
    * 
    *    -
    *       the dky vector must be allocated by the user.       
    *
    *    -
    *       it is only leagal to call this method after a 
    *       successful return from the solve() method.
    * 
    * 


    *
    */
   int getDkyVector(double t, int k, SundialsAbstractVector* dky) const;

   /**
    * Get actual value of the independent variable that CVODE integrated
    * to (i.e. the value of t that actually corresponds to the solution
    * vector y).
    */
   double getActualFinalValueOfIndependentVariable() const;

   /**
    * Print CVODE and CVSpgmr statistics.
    */
   void printStatistics(std::ostream& os) const;

   /**
    * Print CVODE statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    * 


    * 
    *    - @b lenrw            
    *       size (in double words) of memory used for doubles 
    * 
    *    - @b leniw            
    *       size (in integer words) of memory used for integers
    * 
    *    - @b nst              
    *       cumulative number of internal steps taken by solver
    * 
    *    - @b nfe              
    *       number of right-hand side function evaluations
    * 
    *    - @b nni              
    *       number of NEWTON iterations performed 
    * 
    *    - @b nsetups          
    *       number of calls made to linear solver's setup routine
    * 
    *    - @b netf             
    *       number of local error test failures
    * 
    *    - @b ncfn             
    *       number of nonlinear convergence failures
    * 
    *    - @b qu               
    *       order used during the last internal step
    * 
    *    - @b qcur             
    *       order to be used on the next internal step
    * 
    *    - @b hu               
    *       step size for the last internal step
    * 
    *    - @b hcur             
    *       step size to be attempted on the next internal step
    * 
    *    - @b tcur             
    *       current internal value of t reached by the solver
    * 
    *    - @b tolsf            
    *       suggested tolerance scaling factor
    *
    * 


    */
   void printCVODEStatistics(std::ostream& os) const;

   // CVODE optional return values.

   /**
    * Return the cumulative number of internal steps taken by
    * the solver.  
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfInternalStepsTaken() const;

   /**
    * Return the number of calls to the right-hand side function.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfRHSFunctionCalls() const;

   /**
    * Return the number of calls made to linear solver setup
    * routines.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfLinearSolverSetupCalls() const;

   /**
    * Return the number of NEWTON iterations performed.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfNewtonIterations() const;

   /**
    * Return the number of nonlinear convergence failures that have 
    * occurred.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfNonlinearConvergenceFailures() const;

   /**
    * Return the number of local error test failures.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfLocalErrorTestFailures() const;

   /**
    * Return the order of the linear multistep method used during
    * the last internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getOrderUsedDuringLastInternalStep() const;

   /**
    * Return the order of the linear multistep method to be used during
    * the next internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getOrderToBeUsedDuringNextInternalStep() const;

   /**
    * Return the size (in LLNL_REAL words) of memory used
    * for LLNL_REALS.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getCVODEMemoryUsageForDoubles() const;

   /**
    * Return the size (in integer words) of memory used
    * for integers.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getCVODEMemoryUsageForIntegers() const;

   /**
    * Return the step size for the last internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getStepSizeForLastInternalStep() const;

   /**
    * Return the step size to be used in the next internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getStepSizeForNextInternalStep() const;

   /**
    * Return the current internal value of the independent
    * variable reached by the solver.
    *
    * Note: if the solver was not set to collect statistics, 
    * the minimum double value (as defined in float.h) is 
    * returned.
    */
   double getCurrentInternalValueOfIndependentVariable() const;

   /**
    * Return the suggested tolerance scaling factor.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getCVODESuggestedToleranceScalingFactor() const;

   // CVSpgmr optional return values.

   /**
    * Print CVSpgmr statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    * 


    *
    *    - @b spgmr_lrw        
    *      size (in double words) of memory used for doubles
    *
    *    - @b spgmr_liw        
    *       size (in integer words) of memory used for integers
    *
    *    - @b nli              
    *       number of linear iterations
    *
    *    - @b ncfl             
    *       number of linear convergence failures
    *
    *    - @b npe              
    *       number of preconditioner evaluations
    *
    *    - @b nps              
    *       number of calls to CVSpgmrPrecondSolve()
    *
    * 


    */
   void printCVSpgmrStatistics(std::ostream& os) const;

   /**
    * Return the number of preconditioner evaluations.
    */
   int getNumberOfPreconditionerEvaluations() const;

   /**
    * Return the number of linear iterations.
    */
   int getNumberOfLinearIterations() const;

   /**
    * Return the number of CVSpgmrPrecondSolve() calls.
    */
   int getNumberOfPrecondSolveCalls() const;

   /**
    * Return the number of linear convergence failures.
    */
   int getNumberOfLinearConvergenceFailures() const;

   /**
    * Return the size (in double words) of memory used for doubles.
    */
   int getCVSpgmrMemoryUsageForDoubles() const;

   /**
    * Return the size (in integer words) of memory used for integers.
    */
   int getCVSpgmrMemoryUsageForIntegers() const;

   /**
    * Print out all data members for this object.
    */
   virtual void printClassData(std::ostream& os) const;

private:
   /*
    * Static member function for linkage with CVODE routines.
    */ 
   static int CVODERHSFuncEval( realtype t,
                                N_Vector y,
                                N_Vector y_dot,
                                void* my_solver);

   /*
    * Static member functions for linkage with CVSpgmr routines.
    */ 
   static int CVSpgmrPrecondSet(realtype t,
				N_Vector y,
				N_Vector fy,
				int jok,
				booleantype *jcurPtr,
				realtype gamma,
				void *my_solver,
				N_Vector vtemp1,
				N_Vector vtemp2,
				N_Vector vtemp3);

   static int CVSpgmrPrecondSolve(realtype t,
                                  N_Vector y,
                                  N_Vector fy,
                                  N_Vector r,
				  N_Vector z,
                                  realtype gamma,
                                  realtype delta,
                                  int lr,
                                  void *my_solver,
                                  N_Vector vtemp);

   /*
    * Open CVODE log file, allocate main memory for CVODE and initialize
    * CVODE memory record.  CVODE is initialized based on current state
    * of solver parameter data members.  If any solver parameters have 
    * changed since last initialization, this function will be automatically
    * invoked at next call to the solve() method.  Also, if NEWTON iteration 
    * is specified, this method also initializes the CVSpgmr linear solver. 
    *
    * Assertion checks:
    * 


    *
    *    -
    *       the solution vector must have already been set.
    *
    * 


    *
    */
   void initializeCVODE();

   std::string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the CVODESolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * CVODE may need to be re-initialized (e.g., if the linear solver
    * changes, CVODE must change its memory record).  In this case,
    * the initializeCVODE() member function is invoked in the next 
    * call to solve().
    */

   /*
    * Solution vector.
    */
   SundialsAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to CVODE
    * and CVSpgmr.
    */
   CVODEAbstractFunctions* d_cvode_functions;

   /*
    * CVODE memory record.
    */
   void *d_cvode_mem;                    // CVODE memory structure

   /*
    * CVODE log file information.
    */
   FILE*     d_cvode_log_file;               // CVODE message log file
   std::string    d_cvode_log_file_name;     // CVODE log file name

   /*
    * ODE parameters.
    */
   double d_t_0;        // initial value for independent variable
   double d_user_t_f;   // user-specified final value for independent variable
   double d_actual_t_f; // actual final value of indep. variable after a step
   SundialsAbstractVector* d_ic_vector;

   /*
    * ODE integration parameters.
    */
   int d_linear_multistep_method;
   int d_iteration_type;
   int d_tolerance_type;
   double d_relative_tolerance;
   bool d_use_scalar_absolute_tolerance;
   double d_absolute_tolerance_scalar;
   SundialsAbstractVector* d_absolute_tolerance_vector;
   int d_stepping_method;

   /*
    * Optional CVODE parameters.
    */
   int d_max_order;
   int d_max_num_internal_steps;
   int d_max_num_warnings;
   double d_init_step_size;
   double d_max_step_size;
   double d_min_step_size;
   /*
    * CVSpgmr parameters
    */
   int d_precondition_type;
   int d_gram_schmidt_type;
   int d_max_krylov_dim;
   double d_tol_scale_factor;

   /*
    * Boolean flag indicating whether CVODE needs initialization 
    * when solver is called.
    */
   bool d_CVODE_needs_initialization;

   /*
    * Boolean flag indicating whether user-supplied preconditioner
    * routines are provided in the concrete subclass of 
    * CVODEAbstractFunctions.
    */
   bool d_uses_preconditioner;
};


}
}

#endif
#endif
