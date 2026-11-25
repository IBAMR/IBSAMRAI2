//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxTreeNode.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2172 $
// Modified:    $LastChangedDate: 2008-05-02 11:02:08 -0700 (Fri, 02 May 2008) $
// Description: Basic building block used by BoxTree class.
//

#ifndef included_hier_BoxTreeNode
#define included_hier_BoxTreeNode

#include "SAMRAI_config.h"

#include <stdlib.h>
#include <iostream>
#include <vector>

#include "Box.h"
#include "BoxList.h"
#include "BoxArray.h"
#include "tbox/Array.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
   namespace hier {


/*!
 * @brief Building block used by BoxTree<DIM>.
 * 
 * This is a private class that is used by the BoxTree<DIM> class;
 * there is no reason that anyone should make direct use of this class.
 */

template<int DIM> class BoxTreeNode : public tbox::DescribedClass
{
public:

   struct Triple {
      Box<DIM> box;
      int       idx;
      int       owner;
   };

   /*!
    * @brief Constructs a BoxTreeNode that represents the physical
    * domain specified by the bounding box that includes all boxes
    * in \b box_triples.
    * 
    * Constructs a BoxTreeNode that represents the physical
    * domain specified by \b box.  The \b dim parameter is
    * the dimension along which the domain will be cut when
    * constructing child nodes.
    * 
    * @param box_triples input
    * @param have_mapping true is the calling BoxTree was passed a
    *                     processor mapping.
    * @param min_length if \b box_triples contains less than \b min_length
    *                   elements, than the list will never be partitioned
    *                   amongst child nodes.  Setting to a larger value
    *                   tends to decrease the total number of nodes in the tree
    *                   (and hence reduces memory requirements).
    * @param dim the dimension along which the domain will be cut
    *            when constructing child nodes.
    * @param recurse_level is a counter that is used internally to
    *        keep track of the recurse level.  Do not change the
    *        default value.
    */
   BoxTreeNode(
      const tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples,
      bool have_mapping,
      int min_length = 10,
      int dim = DIM-1,
      int recurse_level = 0);

   /*!
    * The dtor does nothing even slightly interesting.
    */
   ~BoxTreeNode();

   /*!
    * @brief Append the triples of this node and its subtrees to \b triples:
    * those of d_tree, or this node's own, then those of the left and right
    * children.  BoxTree lists overlapping boxes in this order.
    */
   void appendTriplesInQueryOrder(std::vector<Triple>& triples) const;

   /*!
    * Counts the total number of nodes in all trees rooted at this
    * node, and add this to the totals storeed in BoxTree<DIM>::s_node_count
    * and BoxTree<DIM>::s_box_count.
    * This is a recursive call, i.e, all nodes in all subtrees
    * rooted at this node are also examined.
    * 
    * This call is primarily of interest for gathering information
    * on performance, during development and testing.
    */
   void count();

private:

   /*!
    * The copy ctor is not implemented.
    */
   BoxTreeNode(const BoxTreeNode<DIM>& box);

   /*!
    * The assignemnt operator is not implemented.
    */
   BoxTreeNode<DIM>& operator=(const BoxTreeNode<DIM>& box);

   /*!
    * The physical domain that this node represents.
    */
   Box<DIM> d_domain;

   /*!
    * Pointers to familial nodes.
    */
   tbox::Pointer< BoxTreeNode<DIM> > d_left_child;
   tbox::Pointer< BoxTreeNode<DIM> > d_right_child;

   /*!
    * The list of boxes, and their associated indices and processor
    * mappings, that are contained within the physical domain that 
    * this node represents.
    */
   tbox::Array<Triple> d_box_triples;

   /*!
    * The dimension along which the input box triples will be partitioned.
    */
   int d_dim;

   /*!
    * A tree rooted at this node.  The tree is constructed
    * from d_box_triples; these boxes are partitioned along dimension
    * d_dim - 1.
    */
   tbox::Pointer< BoxTreeNode<DIM> > d_tree;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxTreeNode.C"
#endif

#endif

