//
// File:	tbox__List-tbox__Pointer-tbox__Database.C
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Description:	Template file automatically generated by make-template.pl
//

#include "tbox/List.h"
#include "tbox/List.C"
#include "tbox/Pointer.h"
#include "tbox/Database.h"

namespace SAMRAI {
   namespace tbox {
#ifdef LACKS_STATIC_DATA_INSTANTIATION
tbox::ListNode< tbox::Pointer< tbox::Database > > *tbox::ListNode< tbox::Pointer< tbox::Database > >::s_free_list=0;
bool tbox::ListNode< tbox::Pointer< tbox::Database > >::s_registered_callback=false;
#endif
template class tbox::List< tbox::Pointer< tbox::Database > >;
template class tbox::ListIterator< tbox::Pointer< tbox::Database > >;
template class tbox::ListNode< tbox::Pointer< tbox::Database > >;
}
}
