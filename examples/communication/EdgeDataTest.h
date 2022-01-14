//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/communication/EdgeDataTest.h $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: AMR communication tests for edge-centered patch data
//

#ifndef included_pdat_EdgeDataTest
#define included_pdat_EdgeDataTest

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "EdgeData.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDataTestStrategy.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#include "Variable.h"

namespace SAMRAI {

class CommTester;

/**
 * Class EdgeDataTest provides routines to test communication operations
 * for edge-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types for test:
 *
 *   NONE...
 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.  Additionally, there are two
 * optional input parameters for each edge variable.  These are:
 *
 * 


 *
 *    - \b  test_direction   edge directions to test 
 *                             (default = -1 ie, all directions)
 *    - \b  use_fine_value_at_interface   which values to use at coarse-
 *                                          fine interface (default = TRUE)
 *
 * 


 * 
 */

class EdgeDataTest : public PatchDataTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   EdgeDataTest(const string& object_name,
                tbox::Pointer<tbox::Database> main_input_db,
                bool do_refine,
                bool do_coarsen,
                const string& refine_option);

   /**
    * Virtual destructor for EdgeDataTest.
    */
   ~EdgeDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   virtual void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                              const double time,
                                              const hier::IntVector<NDIM>& gcw) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   void registerVariables(CommTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
				      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   bool verifyResults(hier::Patch<NDIM>& patch,
                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                      int level_number);

private:
   /*
    * Function for reading test data from input file.
    */
   void readTestInput(tbox::Pointer<tbox::Database> db);

   /*
    * Set constant function data for testing interpolation.
    * Note the data is piecewise constant, where constant 
    * depends on axis and NDIM.  That is, data is set to 
    * data_i = ndimfact*NDIM + axfact*i, where i is
    * axis index.
    */
   void setConstantData(tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
                        const hier::Box<NDIM>& box,
                        double ndimfact,
                        double axfact) const;

   /*
    * Set values for edges along physical boundaries; needed for zero
    * ghost cell case.
    */
   void setConstantBoundaryData(tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
                                const hier::BoundaryBox<NDIM>& bbox,
                                double ndimfact,
                                double axfact) const;

   /*
    * Set conservative linear function data for testing coarsening
    */
   void setConservativeData(tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
                            const hier::Box<NDIM>& box,
                            hier::Patch<NDIM>& patch, 
                            const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
                            int level_number) const;

   void setLinearData(tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;

   void checkPatchInteriorData(
      const tbox::Pointer< pdat::EdgeData<NDIM,double> >& data,
      const hier::Box<NDIM>& interior,
      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const;

   /*
    * Object string identifier for error reporting
    */
   string d_object_name;

   /*
    * Data members specific to this edge data test.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_cart_grid_geometry;

   tbox::Pointer<hier::PatchHierarchy<NDIM> > d_hierarchy;

   /*
    * Data members specific to this edge data test.
    */
   tbox::Array<int> d_test_direction;
   tbox::Array<bool> d_use_fine_value_at_interface;

   double d_Acoef;
   double d_Bcoef;
   double d_Ccoef;
   double d_Dcoef;

   bool d_do_refine;
   bool d_do_coarsen;
   string d_refine_option;
   int d_finest_level_number;

   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables;

};

}
#endif
