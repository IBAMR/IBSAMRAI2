//
// File:	tbox__List-int.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Description:	Template file automatically generated by make-template.pl
//

#include "tbox/List.h"
#include "tbox/List.C"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< int > *tbox::ListNode< int >::s_free_list=0;
bool tbox::ListNode< int >::s_registered_callback=false;
#endif
template class tbox::List< int >;
template class tbox::ListIterator< int >;
template class tbox::ListNode< int >;
}
}
