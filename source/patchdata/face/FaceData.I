//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/face/FaceData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated face centered patch data type
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int FaceData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE* FaceData<DIM,TYPE>::getPointer(int face_normal, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[face_normal].getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* FaceData<DIM,TYPE>::getPointer(int face_normal, int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[face_normal].getPointer(depth));
}

template<int DIM, class TYPE> inline
TYPE& FaceData<DIM,TYPE>::operator()(const FaceIndex<DIM>& i, int depth)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
const TYPE&
FaceData<DIM,TYPE>::operator()(const FaceIndex<DIM>& i, int depth) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>&
FaceData<DIM,TYPE>::getArrayData(int face_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   return(d_data[face_normal]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>&
FaceData<DIM,TYPE>::getArrayData(int face_normal) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
#endif
   return(d_data[face_normal]);
}

template<int DIM, class TYPE> inline 
void FaceData<DIM,TYPE>::copyOnBox(
   const FaceData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{
   for (int axis = 0; axis < DIM; axis++) {
      const hier::Box<DIM> face_box = FaceGeometry<DIM>::toFaceBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), face_box);
   }

}
#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
