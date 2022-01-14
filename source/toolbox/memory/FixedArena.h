//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/memory/FixedArena.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Fixed-size arena for efficient memory management
//

#ifndef included_tbox_FixedArena
#define included_tbox_FixedArena

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Arena.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
   namespace tbox {


/**
 * Class FixedArena implements a very simple fixed-size memory 
 * management scheme.  One large chunk of a known size is allocated and then
 * further allocation requests take memory from that pool.  Memory is never
 * returned to the local pool via free, although memory is freed when the
 * object exits.
 *
 * Note that the size of the memory arena must be calculated including
 * alignment space; all memory requests are aligned on boundaries as
 * defined by Arena::align().
 *
 * @see tbox::Arena
 */

class FixedArena : public Arena
{
public:
   /**
    * Allocate a fixed-size memory arena.  Multiple allocation requests
    * may be made to this arena via alloc() and the specially defined
    * new operators (see Arena).  The size of the fixed arena must
    * include padding space, since all memory will be allocated aligned
    * as defined in Arena.
    */
   FixedArena(const size_t bytes);

   /**
    * Release the memory associated with the fixed arena.  Objects that
    * allocate from a fixed arena must make sure that the arena lives longer
    * than the need for the buffer space.  To ensure that the fixed arena
    * lives long enough, objects might wish to cache a smart pointer to the
    * arena object so that the actual arena object is not deleted until all
    * objects using the arena are deleted.
    */
   virtual ~FixedArena();

   /**
    * Allocate an aligned portion of memory from the fixed arena.  Allocations
    * that exceed the arena space cause the program to abort with a memory
    * overrun exception.
    */
   virtual void *alloc(const size_t bytes);

   /**
    * Return allocated memory to the fixed arena.  For FixedArena,
    * the memory is not actually returned to the free memory pool until
    * the object is destroyed.  Thus, free() is a null operation.
    */
   virtual void free(void *p);

   /**
    * Allocate a fixed arena of the specified size.  Subclasses should
    * override allocateArena() to return the proper subclass type.
    */
   virtual Pointer<Arena> allocateArena(const size_t bytes);

   /**
    * Print out internal class data for debugging.
    */
   virtual void printClassData(std::ostream& os) const;

private:
   FixedArena(const FixedArena&);	// not implemented
   void operator=(const FixedArena&);	// not implemented

   Pointer<Arena> d_arena;
   void *d_pool;
   size_t d_size;
   size_t d_used;
};


}
}

#endif
