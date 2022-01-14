//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/transfer/datamovers/locally_active/LocallyActiveDataFillBox.h $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Routines for "smart" boxlist ops in locally-active data comm ops
//

#ifndef included_hier_LocallyActiveDataFillBox
#define included_hier_LocallyActiveDataFillBox

#include "SAMRAI_config.h"
#include "Box.h"
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "CoarsenClasses.h"
#include "RefineClasses.h"
#ifndef included_iostream
#include <iostream>
#define included_iostream
#endif

namespace SAMRAI {
   namespace xfer {

/*!
 * Class LocallyActiveDataFillBox is a utility class that is used primarily
 * by the LocallyActiveDataFillBoxSet class for boxlist operations in communication 
 * schedules that operate on "locally-active" data; i.e., where each data item may 
 * live on a different set of patches.  Specifically, this class contains a box and 
 * an associated list of either CoarseClass or RefineClass items, but not both.  
 * Each constructor accepts such a list. Once an object is constructed, it can only 
 * be used to support either refine operations or coarsen operations.
 *
 * @see hier::Box
 * @see xfer::CoarsenClasses
 * @see xfer::RefineClasses
 */

template<int DIM>
class LocallyActiveDataFillBox
{
public:
   /*!
    * @brief Construct a locally-active fill box with the given box
    * and list of refine items.  Note that this object maintains pointers
    * to the given box and refine items only.  So, an object of this class
    * becomes unusable if those items are destroyed before this object.
    *
    * @param box       Input box.
    * @param var_data  Input list of refine class data pointers.
    */
   LocallyActiveDataFillBox(
      const hier::Box<DIM>& box,
      const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data);

   /*!
    * @brief Construct a locally-active fill box with the given box
    * and list of coarsen items.  Note that this object maintains pointers
    * to the given box and coarsen items only.  So, an object of this class
    * becomes unusable if those items are destroyed before this object.
    *
    * @param box       Input box.
    * @param var_data  Input list of coarsen class data pointers.
    */
   LocallyActiveDataFillBox(
      const hier::Box<DIM>& box,
      const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data);

   /*!
    * Copy constructor to create a new locally-active fill box and copy 
    * the information from the argument locally-active fill box.
    * 
    * @param fill_box  Constant reference to fill box to be copied.
    */
   LocallyActiveDataFillBox(const LocallyActiveDataFillBox<DIM>& fill_box);

   /*!
    * The destructor releases all box and locally-active data storage.
    */
   virtual ~LocallyActiveDataFillBox<DIM>();

   /*!
    * Return constant reference to box maintained by this object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * Return constant reference to list of refine items maintained by this object.
    * 
    * Note that if this object was created using coarsen item data, an unrecoverable
    * assertion will result.
    */
   const tbox::List<const typename RefineClasses<DIM>::Data*>& 
      getActiveRefineVarData() const;

   /*!
    * Return constant reference to list of coarsen items maintained by this object.
    * 
    * Note that if this object was created using refine item data, an unrecoverable
    * assertion will result.
    */
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& 
      getActiveCoarsenVarData() const;

   /*!
    * Print all class member data for this locally-active fill box object
    * to specified output stream.
    */
   virtual void printClassData(std::ostream& os) const;

   /*!
    * Clear all class member data for this locally-active fill box object.
    */
   void clearLocallyActiveFillBoxData(); 

   /*!
    * @brief Check given box and list of refine items for equality with those maintained
    * by this class object.  Any encountered inequality will be sent to the given output
    * stream.
    * 
    * @return   Boolean true if equal, false otherwise.
    * @param    box        Input box for comparison.
    * @param    var_data   Input list of refine items for comparison.
    * @param    os         Input ostream for reporting mismatches.
    * 
    * Note that if this object was created using coarsen item data, an unrecoverable
    * assertion will result.
    */
   bool checkData(const hier::Box<DIM>& box,
                  const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data,
                  std::ostream& os) const;

   /*!
    * @brief Check given box and list of coarsen items for equality with those maintained
    * by this class object.  Any encountered inequality will be sent to the given output
    * stream.
    * 
    * @return   Boolean true if equal, false otherwise.
    * @param    box        Input box for comparison.
    * @param    var_data   Input list of coarsen items for comparison.
    * @param    os         Input ostream for reporting mismatches.
    *
    * Note that if this object was created using refine item data, an unrecoverable
    * assertion will result.
    */
   bool checkData(const hier::Box<DIM>& box,
                  const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data,
                  std::ostream& os) const;

private:
   const hier::Box<DIM>* d_box;

   tbox::List<const typename xfer::RefineClasses<DIM>::Data*> d_refine_var_data; 
   tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> d_coarsen_var_data; 

   bool d_refine_data;
};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataFillBox.C"
#endif
