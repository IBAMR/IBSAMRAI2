//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mesh/load_balance/SpatialKey.C $
// Package:	SAMRAI mesh generation
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Spatial Key used for generating space-filling curves.
//

#ifndef included_mesh_SpatialKey_C
#define included_mesh_SpatialKey_C

#include "SpatialKey.h"
#include <stdio.h>

#include <iomanip>

#define BITS_PER_BYTE 8
#define BITS_PER_HEX_CHAR 4

#ifdef DEBUG_NO_INLINE
#include "SpatialKey.I"
#endif

/* 
 * Force copy to include dir
 * template<int DIM>
 *
 */

namespace SAMRAI {
    namespace mesh {

/*
****************************************************************************
*                                                                          *
* Default Constructor.  Creates a spatial key with value 0 for all         *
* entries in d_key.                                                        *
*                                                                          *
****************************************************************************
*/
SpatialKey::SpatialKey()
{
   d_bits_per_int = BITS_PER_BYTE*sizeof(unsigned int);
   setToZero();
}

/*
****************************************************************************
*                                                                          *
* Creates a spatial key from the i,j,k coordinates by invoking             *
* setKey().                                                                *
*                                                                          *
****************************************************************************
*/
SpatialKey::SpatialKey(const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int k,
                                 const unsigned int level_num)
{
   d_bits_per_int = BITS_PER_BYTE*sizeof(unsigned int);
   setKey(i,j,k,level_num);
}

/*
****************************************************************************
*                                                                          *
* Creates a SpatialKey by copying the value from a pre-existing       *
* SpatialKey.                                                         *
*                                                                          *
****************************************************************************
*/
SpatialKey::SpatialKey(const SpatialKey& spatial_key)
{
   d_bits_per_int = BITS_PER_BYTE*sizeof(unsigned int);
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      d_key[i] = spatial_key.d_key[i];
   }
}

/*
****************************************************************************
*                                                                          *
* The destructor for a spatial key does nothing interesting.               *
*                                                                          *
****************************************************************************
*/
SpatialKey::~SpatialKey()
{
}

/*
****************************************************************************
*                                                                          *
* Less than operator for spatial keys.  Returns true if the first          *
* integer in the d_key arrays that differs is such that                    *
* this.d_key[i] < spatial_key.d_key[i].                                    *
*                                                                          *
****************************************************************************
*/
bool SpatialKey::operator<(const SpatialKey& spatial_key) const
{
   int i = NUM_COORDS_MIXED_FOR_SPATIAL_KEY - 1;

   while (i >= 0) {
      if (d_key[i] < spatial_key.d_key[i]) {
        return true;
      }
      else if (d_key[i] > spatial_key.d_key[i]) {
        return false;
      }
      i--;
   }

   // the two spatial keys are equal, so return false
   return false;
}

/*
****************************************************************************
*                                                                          *
* Less than or equal operator for spatial keys.  Returns true if           *
* this key is less or equal to the argument key.                           *
*                                                                          *
****************************************************************************
*/
bool SpatialKey::operator<=(const SpatialKey& spatial_key) const
{
   return ( ((*this) < spatial_key) || ((*this) == spatial_key) );
}

/*
****************************************************************************
*                                                                          *
* Greater than operator for spatial keys. Returns true if this key is      *
* is greater than the argument key.                                        *
*                                                                          *
****************************************************************************
*/
bool SpatialKey::operator>(const SpatialKey& spatial_key) const
{
   return ( !((*this) < spatial_key) && ((*this) != spatial_key) );
}

/*
****************************************************************************
*                                                                          *
* Greater than or equal operator for spatial keys.  Returns true if this   *
* key is greater than or equal to the argument key.                        *
*                                                                          *
****************************************************************************
*/
bool SpatialKey::operator>=(const SpatialKey& spatial_key) const
{
   return ( ((*this) > spatial_key) || ((*this) == spatial_key) );
}

/*
****************************************************************************
*                                                                          *
* Write a spatial key to an output stream.  The spatial key is             *
* output in hex to avoid the binary to decimal conversion of the           *
* extended integer key.                                                    *
*                                                                          *
* Uses snprintf() to create string because the behavior of C++ stream      *
* manipulators not standardized yet.                                       *
*                                                                          *
****************************************************************************
*/
std::ostream& operator << (std::ostream& s, const SpatialKey& spatial_key)
{
   const std::size_t size = spatial_key.d_bits_per_int / BITS_PER_HEX_CHAR *
             NUM_COORDS_MIXED_FOR_SPATIAL_KEY + 1;
   char *buf = new char[size];

   for (int i = NUM_COORDS_MIXED_FOR_SPATIAL_KEY - 1; i >= 0; i--) {
      snprintf(&(buf[spatial_key.d_bits_per_int / BITS_PER_HEX_CHAR *
                     ((NUM_COORDS_MIXED_FOR_SPATIAL_KEY - 1) - i)]), size,
               "%08x",spatial_key.d_key[i]);
   }

   s << buf;
   delete [] buf;

   return(s);
}

/*
****************************************************************************
*                                                                          *
* Blends one coordinate into the spatial key.  coord is the                *
* value of the coordinate and coord_offset is the offset for the           *
* coordinate being  blended in.                                            *
*                                                                          *
****************************************************************************
*/
void SpatialKey::blendOneCoord(const unsigned int coord, 
                                    const int coord_offset)
{
   unsigned int shifted_coord = coord;

   int bit_in_int;
   for (bit_in_int = 0; bit_in_int < d_bits_per_int; bit_in_int++) {
      if ( shifted_coord & ((unsigned int) 1) ){
         unsigned int bit_index;
         int int_index;
         int bit_offset;

         bit_index = NUM_COORDS_MIXED_FOR_SPATIAL_KEY * bit_in_int + 
                     coord_offset;
         int_index = bit_index / d_bits_per_int;
         bit_offset = bit_index & (d_bits_per_int - 1);
         d_key[int_index] |= (((unsigned int) 1) << bit_offset);

      }
      shifted_coord = shifted_coord >> 1;
   }
}

/*
****************************************************************************
*                                                                          *
* setKey() takes the index space coordinates and the level number          *
* and sets the value of the spatial key.  If the coordinates have          *
* binary representation given by                                           *
* (i32)(i31)...(i1)(i0), etc., the resulting spatial key                   *
* has the following form:                                                  *
* (i32)(j32)(k32)(ln32)(i31)(j31)(k31)(ln31)...(i0)(j0)(k0)(ln0).          *
* This result is stored as an array of four unsigned integers.             *
*                                                                          *
****************************************************************************
*/
void SpatialKey::setKey(const unsigned int i,
                             const unsigned int j,
                             const unsigned int k,
                             const unsigned int level_num)
{
   setToZero();

   /* blend in x coordinate */
   blendOneCoord(i,3);

   /* blend in y coordinate */
   blendOneCoord(j,2);

   /* blend in z coordinate */
   blendOneCoord(k,1);

   /* blend in level number */
   blendOneCoord(level_num,0);
}

}
}

#endif
