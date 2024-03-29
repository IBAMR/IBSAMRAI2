//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/mesh/load_balance/SpatialKey.I $
// Package:	SAMRAI mesh generation
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Spatial Key used for generating space-filling curves.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace mesh {

/*
**************************************************************************** 
*                                                                          *
* Copies the spatial key.                                                  *
*                                                                          *
****************************************************************************
*/
inline
SpatialKey& SpatialKey::operator=(const SpatialKey& spatial_key)
{
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      d_key[i] = spatial_key.d_key[i];
   }
   return (*this);
}

/*
**************************************************************************** 
*                                                                          *
* Compare two spatial keys for equality.  All values in d_key must be      *
* equal for this to return true.                                           *
*                                                                          *
****************************************************************************
*/
inline
bool SpatialKey::operator==(const SpatialKey& spatial_key) const
{
   bool are_equal = true;
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      if ( d_key[i] != spatial_key.d_key[i] ) are_equal = false;
   }
   return (are_equal);
}

/*
**************************************************************************** 
*                                                                          *
* Compare two spatial keys for inequality.  Returns true if the keys are   * 
* not equal.                                                               *
*                                                                          *
****************************************************************************
*/
inline
bool SpatialKey::operator!=(const SpatialKey& spatial_key) const
{
   return ( !((*this) == spatial_key) );
}

/*
**************************************************************************** 
*                                                                          *
* Zeroes key by setting all values in d_key to 0.                          *
*                                                                          *
**************************************************************************** 
*/
inline
void SpatialKey::setToZero()
{
   for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; i++) {
      d_key[i] = 0;
   }
}

}
}
#ifdef DEBUG_NO_INLINE
#undef inline
#endif
