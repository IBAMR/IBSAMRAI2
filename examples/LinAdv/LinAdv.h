//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/LinAdv/LinAdv.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Numerical routines for single patch in linear advection ex.
//
 
#ifndef included_LinAdvXD
#define included_LinAdvXD

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "BoundaryUtilityStrategy.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "CellVariable.h"
#include "tbox/Database.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HyperbolicLevelIntegrator.h"
#include "HyperbolicPatchStrategy.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include <string>
using namespace std;
#define included_String
#include "VariableContext.h"
#include "VisItDataWriter.h"

/**
 * The LinAdv class provides routines for a sample application code that
 * solves linear advection.  This code illustrates the manner in which 
 * a code employing the standard Berger/Oliger AMR algorithm for 
 * explicit hydrodynamics can be used in the SAMRAI framework.
 * This class is derived from the algs::HyperbolicPatchStrategy<NDIM> abstract base
 * class which defines the bulk of the interface between the hyperbolic
 * intergration algorithm provided by SAMRAI and the numerical routines
 * specific to linear advection.  In particular, this class provides routines 
 * which maybe applied to any patch in an AMR patch hierarchy.
 *
 * The linear advection problem is simply du/dt + div(a*u) = 0, where
 * "u" is a scalar-valued function and "a" is a constant vector.  The
 * primary numerical quantities are "uval" and "flux", which represent
 * "u" and "a*u", respectively.  All other variables are temporary
 * quantities used in the numerical routines.  The numerical routines 
 * use explicit timestepping and a second-order unsplit Godunov method.
 */

using namespace SAMRAI;

class LinAdv : 
  public tbox::Serializable,
  public algs::HyperbolicPatchStrategy<NDIM>,
  public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor for LinAdv sets default parameters for the linear
    * advection model.  Specifically, it creates variables that represent
    * the state of the solution.  The constructor also registers this
    * object for restart with the restart manager using the object name.
    *
    * After default values are set, this routine calls getFromRestart()
    * if execution from a restart file is specified.  Finally,
    * getFromInput() is called to read values from the given input
    * database (potentially overriding those found in the restart file).
    */
   LinAdv(const string& object_name,
          tbox::Pointer<tbox::Database> input_db,
          tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom);
 
    /**
     * The destructor for LinAdv does nothing.
     */
   ~LinAdv();

   ///
   ///  The following routines:
   ///
   ///      registerModelVariables(),
   ///      initializeDataOnPatch(),
   ///      computeStableDtOnPatch(),
   ///      computeFluxesOnPatch(),
   ///      conservativeDifferenceOnPatch(),
   ///      tagGradientDetectorCells(),
   ///      tagRichardsonExtrapolationCells()
   ///
   ///  are concrete implementations of functions declared in the
   ///  algs::HyperbolicPatchStrategy<NDIM> abstract base class.
   ///

   /**
    * Register LinAdv model variables with algs::HyperbolicLevelIntegrator<NDIM>
    * according to variable registration function provided by the integrator.
    * In other words, variables are registered according to their role
    * in the integration process (e.g., time-dependent, flux, etc.).
    * This routine also registers variables for plotting with the
    * Vis writer (Vizamrai or VisIt).  
    */
   void registerModelVariables(algs::HyperbolicLevelIntegrator<NDIM>* integrator);

   /**
    * Set up parameters in the load balancer object (owned by the gridding
    * algorithm) if needed.  The Euler model allows non-uniform load balancing
    * to be used based on the input file parameter called
    * "use_nonuniform_workload".  The default case is to use uniform
    * load balancing (i.e., use_nonuniform_workload == false).  For
    * illustrative and testing purposes, when non-uniform load balancing is
    * turned on, a weight of one will be applied to every grid cell.  This
    * should produce an identical patch configuration to the uniform load
    * balance case.
    */
   void setupLoadBalancer(algs::HyperbolicLevelIntegrator<NDIM>* integrator,
                          mesh::GriddingAlgorithm<NDIM>* gridding_algorithm);

   /**
    * Set the data on the patch interior to some initial values,
    * depending on the input parameters and numerical routines.
    * If the "initial_time" flag is false, indicating that the
    * routine is called after a regridding step, the routine does nothing.
    */
   void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                              const double data_time,
                              const bool initial_time);
 
   /**
    * Compute the stable time increment for patch using a CFL
    * condition and return the computed dt.
    */
   double computeStableDtOnPatch(hier::Patch<NDIM>& patch,
                                 const bool initial_time,
                                 const double dt_time);
 
   /**
    * Compute time integral of fluxes to be used in conservative difference
    * for patch integration.  When NDIM == 3, this function calls either
    * compute3DFluxesWithCornerTransport1(), or
    * compute3DFluxesWithCornerTransport2() depending on which
    * transverse flux correction option that is specified in input.
    * The conservative difference used to update the integrated quantities
    * is implemented in the conservativeDifferenceOnPatch() routine.
    */
   void computeFluxesOnPatch(hier::Patch<NDIM>& patch,
                             const double time,
                             const double dt);
 
   /**
    * Update linear advection solution variables by performing a conservative
    * difference with the fluxes calculated in computeFluxesOnPatch().
    */
   void conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const double dt,
                                      bool at_syncronization);

   /**
    * Tag cells for refinement using gradient detector.
    */
   void tagGradientDetectorCells(
      hier::Patch<NDIM>& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_indexindx,
      const bool uses_richardson_extrapolation_too);

   /**
    * Tag cells for refinement using Richardson extrapolation.
    */
   void tagRichardsonExtrapolationCells(
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

   ///
   ///  The following routines:
   ///
   ///      setPhysicalBoundaryConditions()
   ///
   ///  are concrete implementations of functions declared in the
   ///  RefinePatchStrategy abstract base class.
   ///

   /**
    * Set the data in ghost cells corresponding to physical boundary
    * conditions.  Specific boundary conditions are determined by
    * information specified in input file and numerical routines.
    */ 
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double fill_time,
                                      const hier::IntVector<NDIM>&
                                      ghost_width_to_fill);

   /**
    * Write state of LinAdv object to the given database for restart.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary 
    * condition applies.
    */
   void readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                       string& db_name,
                                       int bdry_location_index);

   /**
    * Register a Vizamrai data writer so this class will write
    * plot files that may be postprocessed with the Vizamrai 
    * visualization tool.
    */
   void registerVizamraiDataWriter( 
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_writer);

#ifdef HAVE_HDF5
   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt 
    * visualization tool.
    */
   void registerVisItDataWriter( 
      tbox::Pointer<appu::VisItDataWriter<NDIM> > viz_writer);
#endif

   /**
    * Reset physical boundary values in special cases, such as when
    * using symmetric (i.e., reflective) boundary conditions.
    */
   void boundaryReset(hier::Patch<NDIM>& patch,
                      pdat::FaceData<NDIM,double>& traced_left,
                      pdat::FaceData<NDIM,double>& traced_right) const;

   /**
    * Print all data members for LinAdv class.
    */
   void printClassData(ostream& os) const;

private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    *
    * An assertion results if the database pointer is null.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);

   void getFromRestart();

   void readStateDataEntry(tbox::Pointer<tbox::Database> db,
                           const string& db_name,
                           int array_indx,
                           tbox::Array<double>& uval); 

   /*
    * Private member function to check correctness of boundary data.
    */
   void checkBoundaryData(int btype,
                          const hier::Patch<NDIM>& patch,
                          const hier::IntVector<NDIM>& ghost_width_to_fill,
                          const tbox::Array<int>& scalar_bconds) const;


   /*
    * Three-dimensional flux computation routines corresponding to
    * either of the two transverse flux correction options.  These
    * routines are called from the computeFluxesOnPatch() function.
    */
   void compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                            const double dt);
   void compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                            const double dt);


   /*
    * The object name is used for error/warning reporting and also as a 
    * string label for restart database entries. 
    */
   string d_object_name;

   /*
    * We cache pointers to the grid geometry and Vizamrai data writer
    * object to set up initial data, set physical boundary conditions,
    * and register plot variables.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > d_vizamrai_writer;
#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<NDIM> > d_visit_writer;
#endif

   /*
    * Data items used for nonuniform load balance, if used.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_workload_variable;
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   /**
    * tbox::Pointer to state variable vector - [u]
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_uval;

   /**
    * tbox::Pointer to flux variable vector  - [F]
    */
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_flux;

   /** 
    * linear advection velocity vector 
    */
   double d_advection_velocity[NDIM];

   /*
    *  Parameters for numerical method:
    *
    *    d_godunov_order ....... order of Godunov slopes (1, 2, or 4)
    *
    *    d_corner_transport .... type of finite difference approximation
    *                            for 3d transverse flux correction
    *
    *    d_nghosts ............. number of ghost cells for cell-centered
    *                            and face/side-centered variables
    *
    *    d_fluxghosts .......... number of ghost cells for fluxes
    *
    */
   int d_godunov_order;
   string d_corner_transport;
   hier::IntVector<NDIM> d_nghosts;
   hier::IntVector<NDIM> d_fluxghosts;

   /*
    * Indicator for problem type and initial conditions
    */
   string d_data_problem;
   int d_data_problem_int;

   /*
    * Input for SPHERE problem
    */ 
   double d_radius;
   double d_center[NDIM];
   double d_uval_inside; 
   double d_uval_outside; 

   /*
    * Input for FRONT problem
    */
   int d_number_of_intervals;
   tbox::Array<double> d_front_position;
   tbox::Array<double> d_interval_uval; 

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   tbox::Array<int> d_scalar_bdry_edge_conds;
   tbox::Array<int> d_scalar_bdry_node_conds;
#if (NDIM == 3)
   tbox::Array<int> d_scalar_bdry_face_conds;
#endif

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
#if (NDIM ==2)
   tbox::Array<int> d_node_bdry_edge;
#endif
#if (NDIM == 3)
   tbox::Array<int> d_edge_bdry_face;
   tbox::Array<int> d_node_bdry_face;
#endif


   /*
    * Arrays of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
#if (NDIM ==2)
   tbox::Array<double> d_bdry_edge_uval;
#endif
#if (NDIM ==3)
   tbox::Array<double> d_bdry_face_uval;
#endif

   /*
    * Input for Sine problem initialization
    */
   double d_amplitude;
   double d_frequency[NDIM];

   /*
    * Refinement criteria parameters for gradient detector and
    * Richardson extrapolation. 
    */
   tbox::Array<string> d_refinement_criteria;
   tbox::Array<double> d_dev_tol;
   tbox::Array<double> d_dev;
   tbox::Array<double> d_dev_time_max;
   tbox::Array<double> d_dev_time_min;
   tbox::Array<double> d_grad_tol;
   tbox::Array<double> d_grad_time_max;
   tbox::Array<double> d_grad_time_min;
   tbox::Array<double> d_shock_onset;
   tbox::Array<double> d_shock_tol;
   tbox::Array<double> d_shock_time_max;
   tbox::Array<double> d_shock_time_min;
   tbox::Array<double> d_rich_tol;
   tbox::Array<double> d_rich_time_max;
   tbox::Array<double> d_rich_time_min;

};

#endif
