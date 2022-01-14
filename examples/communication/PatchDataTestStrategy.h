//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/examples/communication/PatchDataTestStrategy.h $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Base class for patch data test operations.
//

#ifndef included_hier_PatchDataTestStrategy
#define included_hier_PatchDataTestStrategy

#include "SAMRAI_config.h"


#include "tbox/Array.h"
#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "VariableContext.h"

using namespace std;

namespace SAMRAI {

class CommTester;

/**
 * Class PatchDataTestStrategy defines an interface for testing specific 
 * patch data transfer operations on individual patches when using 
 * the CommTester class.  This base class provides two member functions
 * for reading test input information.  These are:
 *
 * readVariableInput(): This function reads in a collection of databases,
 * each of which contains parameters for a single variable.  The names of
 * the sub-databases are arbitrary, but must be distinct.  Each variable
 * sub-database has the following input keys:
 *
 * 

 
 * 
 *    - \b  name         variable name string (required)
 *    - \b  depth        optional variable depth (default = 1)
 *    - \b  src_ghosts   optional comm source ghost width (default = 0,0,0)
 *    - \b  dst_ghosts   optional comm dest ghost width (default = 0,0,0)
 *    - \b  coarsen_operator   opt. coarsen op name (default = "NO_COARSEN")
 *    - \b  refine_operator    opt. refine op name (default = "NO_REFINE")
 *
 * 


 *
 * readRefinementInput(): This function reads in a collection of box
 * arrays, each of which describes the region to refine on each level.
 * The key names of the box arrays are arbitrary, but must be distinct. 
 * For example,
 * 
 * 


 *
 *    - \b  level0_boxes   boxes to refine on level zero.
 *    - \b  level1_boxes   boxes to refine on level one.
 *    - \b  ...            etc...
 *
 * 


 * 
 *
 * The pure virtual functions in this class that must be provided by
 * concrete test subclass:
 * \begin{enumerate}
 *    - [registerVariables(...)] register variables with CommTester.
 *    - [initializeDataOnPatch(...)] set patch data on new patch.
 *    - [verifyResults(...)] check results of communication operations.
 * \end{enumerate} 
 *
 * The following virtual functions are given default non-operations in this
 * class so that concrete test subclass can either implement them to test
 * specific functionality or simply ignore.  They are pure virtual in the
 * coarsen and refine patch strategy classes:
 * \begin{enumerate}
 *    - [setPhysicalBoundaryConditions(...)]
 *    - [preprocessRefine(...)]
 *    - [postprocessRefine(...)]
 *    - [preprocessCoarsen(...)]
 *    - [postprocessCoarsen(...)]
 * \end{enumerate}
 */

class PatchDataTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   PatchDataTestStrategy();

   /**
    * Virtual destructor for PatchDataTestStrategy.
    */
   virtual ~PatchDataTestStrategy();

   /**
    * Grid geometry access operations.
    */
   void setGridGeometry(tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom) 
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!grid_geom.isNull());
#endif
      d_grid_geometry = grid_geom;
   }

   ///
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > getGridGeometry() const
   {
      return(d_grid_geometry);
   }

   /**
    * Utility functions for managing patch data context.
    */
   void setDataContext(tbox::Pointer<hier::VariableContext> context)
   {
#ifdef DEBUG_CHECK_ASSERTIONS    
      TBOX_ASSERT(!context.isNull());
#endif 
      d_data_context = context;
   }

   ///
   tbox::Pointer<hier::VariableContext> getDataContext() const
   {
      return(d_data_context);
   }

   void clearDataContext()
   {
      d_data_context.setNull();
   }

   /**
    * Read variable parameters from input database.
    */
   void readVariableInput(tbox::Pointer<tbox::Database> db);

   /**
    * Read arrays of refinement boxes from input database.
    */
   void readRefinementInput(tbox::Pointer<tbox::Database> db); 

   /**
    * Virtual functions in interface to user-supplied boundary conditions,
    * coarsen and refine operations.
    */
   virtual void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                              const double time,
                                              const hier::IntVector<NDIM>& gcw) const;

   ///
   virtual void preprocessRefine(hier::Patch<NDIM>& fine,
                                 const hier::Patch<NDIM>& coarse,
                                 const hier::Box<NDIM>& fine_box,
                                 const hier::IntVector<NDIM>& ratio) const;

   ///
   virtual void postprocessRefine(hier::Patch<NDIM>& fine,
                                  const hier::Patch<NDIM>& coarse,
                                  const hier::Box<NDIM>& fine_box,
                                  const hier::IntVector<NDIM>& ratio) const;

   ///
   virtual void preprocessCoarsen(hier::Patch<NDIM>& coarse,
                                  const hier::Patch<NDIM>& fine,
                                  const hier::Box<NDIM>& coarse_box,
                                  const hier::IntVector<NDIM>& ratio) const;

   ///
   virtual void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                                   const hier::Patch<NDIM>& fine,
                                   const hier::Box<NDIM>& coarse_box,
                                   const hier::IntVector<NDIM>& ratio) const; 

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   virtual void registerVariables(CommTester* commtest) = 0;

   /**
    * Virtual function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
				      char src_or_dst) = 0;
   /**
    * Virtual function for checking results of communication operations.
    */
   virtual bool verifyResults(hier::Patch<NDIM>& patch,
                              const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                              int level_number) = 0;

protected:
   /*
    * Arrays of information read from input file describing test variables
    */
   tbox::Array<string>    d_variable_src_name;
   tbox::Array<string>    d_variable_dst_name;
   tbox::Array<int>       d_variable_depth;
   tbox::Array<hier::IntVector<NDIM> > d_variable_src_ghosts;
   tbox::Array<hier::IntVector<NDIM> > d_variable_dst_ghosts;
   tbox::Array<string>    d_variable_coarsen_op;
   tbox::Array<string>    d_variable_refine_op;

   /*
    * Arrays of information read from input file describing test variables
    */
   tbox::Array<hier::BoxArray<NDIM> > d_refine_level_boxes;

private:
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

   tbox::Pointer<hier::VariableContext> d_data_context;

};

}
#endif
