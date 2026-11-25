//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxTreeNode.C $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Basic building block used by BoxTree class.
//

#ifndef included_hier_BoxTreeNode_C
#define included_hier_BoxTreeNode_C

#include "BoxTreeNode.h"


namespace SAMRAI {
   namespace hier {

/*
 * ************************************************************************
 *                                                                       *
 * Constructs a BoxTreeNodeX that represents the physical
 * domain specified by box.
 *                                                                       *
 * ************************************************************************
 */


template<int DIM>  BoxTreeNode<DIM>::BoxTreeNode(
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples,
   bool have_mapping,
   int min_length,
   int dim,
   int recurse_level)
{
   d_dim = dim;

   /*
    * Compute this node's domain, which is the bounding box
    * for the list of boxes.
    */
   int len = box_triples.size();
   for (int j=0; j<len; ++j) {
      d_domain += box_triples[j].box;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_tree.
    */
   min_length = (min_length < 1) ? 1 : min_length;
   if (len <= min_length || dim == 0) {
      d_box_triples = box_triples;
      return;
   }

   /*
    * Partition the boxes into three sets:
    *    - those that belong to me
    *    - those that belong to my left child
    *    - those that belong to my right child
    */
   tbox::Array< typename BoxTreeNode<DIM>::Triple> left, right;
   int mid = (d_domain.lower(d_dim) + d_domain.upper(d_dim))/2;

   //first, compute the number of items in each list
   //(this loop would be eliminated if we had an Array::appendItem())
   int my_ct, left_ct, right_ct;
   my_ct = left_ct = right_ct = 0;
   for (int i=0; i<len; ++i) {
      if (box_triples[i].box.upper(d_dim)  <= mid) {
         ++left_ct;
      } else if (box_triples[i].box.lower(d_dim) > mid) {
         ++right_ct;
      } else {
         ++my_ct;
      }
   }

   //second, resize the arrays and insert the elements
   d_box_triples.resizeArray(my_ct); 
   left.resizeArray(left_ct);
   right.resizeArray(right_ct);
   my_ct = left_ct = right_ct = 0;

   for (int i=0; i<len; ++i) {
      if (box_triples[i].box.upper(d_dim)  <= mid) {
         left[left_ct++] = box_triples[i];
      } else if (box_triples[i].box.lower(d_dim) > mid) {
         right[right_ct++] = box_triples[i];
      } else {
         d_box_triples[my_ct++] = box_triples[i];
      }
   }

   /*
    * Recurse to build a private tree on this node;
    * this tree contains the boxes assigned to this node,
    * but they will be partitioned with referenc to
    * the next lower dimension.
    */
   if (d_dim > 0) {
      d_tree = new BoxTreeNode<DIM>(d_box_triples, have_mapping, 
                                    min_length, d_dim-1, recurse_level+1);
      d_box_triples.setNull();
   }

   /*
    * Recurse to build this node's left and right children.
    */
   if (!left.isNull()) {
      d_left_child = new BoxTreeNode<DIM>(left, have_mapping, 
                                          min_length, d_dim-1, recurse_level+1);
   }

   if (!right.isNull()) {
      d_right_child = new BoxTreeNode<DIM>(right, have_mapping, 
                                           min_length, d_dim-1, recurse_level+1);
   }
}


/*
 * ************************************************************************
 *                                                                       *
 * dtor
 *                                                                       *
 * ************************************************************************
 */

template<int DIM>  BoxTreeNode<DIM>::~BoxTreeNode()
{
}

template<int DIM> void BoxTreeNode<DIM>::appendTriplesInQueryOrder(
   std::vector<Triple>& triples) const
{
   if (d_tree) {
      d_tree->appendTriplesInQueryOrder(triples);
   } else {
      for (int i = 0; i < d_box_triples.size(); ++i) {
         triples.push_back(d_box_triples[i]);
      }
   }
   if (d_left_child) {
      d_left_child->appendTriplesInQueryOrder(triples);
   }
   if (d_right_child) {
      d_right_child->appendTriplesInQueryOrder(triples);
   }
}

}
}

#endif
