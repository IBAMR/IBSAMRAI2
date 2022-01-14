//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/database/DatabaseBox.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	A box structure representing a portion of the AMR index space
//

#ifndef included_tbox_DatabaseBox
#define included_tbox_DatabaseBox

#include "SAMRAI_config.h"

#ifndef DatabaseBox_MAX_DIM
#define DatabaseBox_MAX_DIM 3
#else
Macro overloaded: DatabaseBox_MAX_DIM
#endif

namespace SAMRAI {
   namespace tbox {


/*!
 * @brief POD data for class DatabaseBox
 * 
 * The data in DatabaseBox need to reside in a POD class so that
 * HDF5's HOFFSET macro works.  (According to ANSI C++ standard,
 * it does not have to work with non-POD data.)
 */
struct DatabaseBox_POD {
   int d_dimension;
   int d_lo[DatabaseBox_MAX_DIM];
   int d_hi[DatabaseBox_MAX_DIM];
};


/**
 * Class DatabaseBox represents a one, two, or three dimensional box in the
 * AMR index space.  It is defined by lower and upper bounds given by integer
 * arrays.
 *
 * This box is an auxilliary data structure used by the database routines to
 * manipulate boxes.  This box type removes cyclic dependencies among the
 * database routines (which need a box) and the box (which needs the database
 * routines).  The box classes in the hierarchy package convert this box
 * structure into the standard SAMRAI box class used by the AMR algorithms.
 *
 * @internal This class should have @em NO data except for d_data.
 * See d_data for details.
 */

class DatabaseBox
{
public:
   /**
    * The default constructor creates a zero dimension empty box.
    */
   DatabaseBox();

   /**
    * Create a box of the specified dimension describing the index space
    * between lower and upper.  The dimension argument must be zero, one,
    * two, or three.
    */
   DatabaseBox(const int dimension, const int *lower, const int *upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   DatabaseBox(const DatabaseBox& box);

   /**
    * The assignment operator copies the index space of the argument box.
    */
   DatabaseBox& operator=(const DatabaseBox& box);

   /**
    * The destructor does nothing interesting.
    */
   ~DatabaseBox();

   /**
    * Return whether the box is empty.  A box is empty if it has dimension
    * zero or if any of the upper dimensions is less than its corresponding
    * lower dimension.
    */
   bool empty() const;

   /**
    * Return the currently defined dimension of the box.
    */
   int getDimension() const;

   /**
    * Set the current dimension of the box.  This must be in the range
    * from zero (an empty box) to three.
    */
   void setDimension(const int dimension);

   /**
    * Return the specified component (non-const) of the lower index of the box.
    */
   int& lower(const int i);

   /** 
    * Return the specified component (non-const) of the upper index of the box.
    */
   int& upper(const int i);

   /** 
    * Return the specified component (const) of the lower index of the box.
    */
   int lower(const int i) const;

   /** 
    * Return the specified component (const) of the upper index of the box.
    */
   int upper(const int i) const;

   /**
    * Check whether two boxes represent the same portion of index space.
    */
   int operator==(const DatabaseBox& box) const;

   /**
    * @brief All data members in a POD type.
    *
    * Due to the need to compute offsets for data members and that
    * offsets cannot be computed for non-POD data, we place all
    * data members in a POD struct and own an object of that
    * struct.
    *
    * Data members are public so that the HDFDatabase need not 
    * mirror this structure in defining a compound type for HDF.
    */
   DatabaseBox_POD d_data;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/DatabaseBox.I"
#endif
#endif
