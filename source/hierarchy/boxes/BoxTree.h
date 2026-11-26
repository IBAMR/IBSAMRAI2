//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxTree.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2172 $
// Modified:    $LastChangedDate: 2008-05-02 11:02:08 -0700 (Fri, 02 May 2008) $
// Description: Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTree
#define included_hier_BoxTree

#include "SAMRAI_config.h"

#include <stdlib.h>
#include <utility>
#include <vector>

#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "BoxTreeNode.h"
#include "ProcessorMapping.h"
#include "tbox/DescribedClass.h"

namespace SAMRAI {
   namespace hier {


/*!
 * Class BoxTree is a utility class that provides functionality
 * that can be used to reduce the runtime complexity of certain box
 * calculus operations.
 * 
 * A BoxTree object is constructed by passing the ctor a BoxArray;
 * this array is used internally to setup private data structures that
 * will be used in subsequent methods, described below. If the
 * BoxArray contains \b n boxes, then the setup phase (which is
 * invoked by the constructor) has an expected runtime complexity 
 * of O(n log(n)).
 *
 * The constructor builds the tree only to fix the order in which queries
 * list boxes.  The queries use a BinIndex of the boxes and their periodic
 * images, so their cost usually grows with the number of boxes near the
 * query box rather than with n.  The queries return the same results, in
 * the same order, as a search of the tree.
 * 
 * Following construction, two types of operations are supported by the 
 * member functions findOverlapIndices() and findLocalOverlapIndices().  
 * Each function takes an integer array and a box as input, and computes 
 * the indices of the subset of the boxes (in the array that was passed 
 * to the constructor) that overlap with the argument box.  On return from 
 * the function, the integer array holds the indices of the boxes in the 
 * array that overlap with the argument box.  The function findOverlapIndices()
 * returns indices of all boxes that overlap the argument box.  The function
 * findLocalOverlapIndices() returns only the indices of the boxes that are
 * assigned to this processor.
 * 
 * Second, the removeIntersections call takes a list of boxes as
 * input, and removes those portions that intersect with the boxes
 * in the BoxArray that was passed to the constructor.
 *
 * This class provides functionality that is similar to that provided
 * by the BoxTop and BoxGraph classes.  However, it employs different
 * algorithms to do so, and is expected to have different runtime
 * characteristics.  Questions such as "which is faster for
 * removeIntersections(): BoxTop, BoxGraph, or BoxTree?" depends on many
 * factors for which analysis is difficult.  The answer is largely a
 * function of a particular problem and domain topology, and is best 
 * answered experimentally.
 *
 * 
 * @see hier::BoxTop
 * @see hier::BoxGraph
 * @see hier::Box
 * @see hier::BoxArray
 */
template<int DIM> class BoxTree : public tbox::DescribedClass
{
public:

   //! @{
   /*!
    * @name BoxTree constructors
    */

   /*!
    * @brief Constructor for BoxTree.
    *
    * The primary difference between the constructors is that the first takes
    * processor mapping information while the others do not.  This is because
    * processor mapping information is required in the method
    * findLocalOverlapIndices(), but is not required for the other 
    * methods.
    *
    * If you do not pass processor mappying information, and subsequently call
    * findLocalOverlapIndices(), an unrecoverable assertion will be thrown.
    *
    * @param boxes input array of boxes.
    * @param box_shifts the amount by which each box is shifted when there
    *               are periodic boundary conditions.  If there are no
    *               periodic boundary conditions you can pass an array
    *               of length zero; otherwise the array must contain
    *               an entry for each box in \b boxes. 
    * @param box_mapping maps each box in \b boxes to a processor.
    * @param min_length controls the partitioning of \b boxes amongst
    *                   child nodes in the tree, and so the order in which
    *                   findOverlapBoxes() lists boxes.  It does not affect
    *                   the cost of the queries.  Please use the default
    *                   value.
    */
   BoxTree(const BoxArray<DIM>& boxes,
                 const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
                 const ProcessorMapping& box_mapping,
                 int min_length = 10);

   BoxTree(const BoxArray<DIM>&  boxes,
                 const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
                 int min_length = 10);

   BoxTree(const BoxArray<DIM>& boxes,
                 int min_length = 10);
   //@}

   
   /*!
    * The dtor does nothing interesting.
    */
   ~BoxTree();

   //! @{
   /*!
    * @name BoxTree findIndices
    */

   /*!
    * @brief Compute the indices of all boxes that overlap the specified \b box.
    *
    * The integer array will contain the indices of all boxes (in the BoxArray 
    * that was passed to the constructor) and that overlap the specified \b box.  
    *
    * @param indices           integer array to hold the array indices of 
    *                          the overlapping boxes in the box array.
    * @param box               box whose overlaps are requested.
    */
   void findOverlapIndices(tbox::Array<int>& indices,
                           const Box<DIM>& box) const;

   /*!
    * @brief Compute the indices of the boxes local to this processor that 
    * overlap the specified \b box.
    *
    * The integer array will contain the indices of the boxes (in the BoxArray
    * that was passed to the constructor) and that overlap the specified \b box
    * and which are assigned to this processor.  The processor assignment is
    * determined from the mapping information passed to the constructor.  This
    * routine will result in an unrecoverable exception if a constructor that
    * does not accept the mapping information was used to create this object.
    * See the comments for the constructors. 
    *
    * @param indices           integer array to hold the array indices of
    *                          the overlapping boxes in the box array.
    * @param box               box whose overlaps are requested.
    */
   void findLocalOverlapIndices(tbox::Array<int>& indices,
                                const Box<DIM>& box) const;

   //! @}

   /*!
    * @brief Compute list of boxes that overlap the specified \b box.
    *
    * @param overlap_boxes boxes in BoxArray that was passed to the constructor
    *                      that overlap with argument box.
    * @param box           box whose overlaps are requested.
    */
   void findOverlapBoxes(BoxList<DIM>& overlap_boxes, 
                         const Box<DIM>& box) const;

  /*!
   * @brief Remove from \b list the portions that intersect the boxes
   * in the BoxArray that was passed to the constructor.
   *
   * CAUTION: the semantics of this call differ from that of
   * BoxList::removeIntersections(const BoxList<DIM> takeaway).
   * Here, the list that is being modified is the list that is
   * passed as an argument; the "takeaway" list is the list that
   * was passed when the BoxTop object was constructed.
   *
   * @param list the list of boxes from which intersections are to be removed.
   */
  void removeIntersections(BoxList<DIM>& list) const;

   /*!
    * @brief Bin index of a sequence of boxes, for finding the boxes that
    * intersect a given box.
    *
    * The boxes are numbered in the order in which they are given, and a
    * query returns the numbers of the boxes it finds in increasing order,
    * so the caller fixes the order of the results through the order of the
    * boxes.  The boxes are sorted by size into classes whose bins grow by
    * factors of two from the smallest box extent in each direction, so each
    * box meets at most 2^DIM bins of its class.  A class stores every bin
    * of its bounding box if there are at most twice as many bins as
    * entries, and otherwise only its B occupied bins, in a sorted list.
    * In each class a query visits the bins its range covers, at most
    * (query extent / bin size + 1)^DIM and never more than the occupied
    * bins, at a cost of log B per bin for a sorted list, and tests the
    * boxes of those bins.  Sequences of at most 32 boxes are scanned
    * instead.
    */
   class BinIndex
   {
   public:
      /*!
       * @brief Construct an index of no boxes.
       */
      BinIndex();

      /*!
       * @brief Construct an index of the boxes of \b boxes, numbered in
       * list order.
       */
      explicit BinIndex(const BoxList<DIM>& boxes);

      /*!
       * @brief Index \b boxes, numbered in order, in place of the boxes
       * indexed before.
       */
      void initialize(const std::vector< Box<DIM> >& boxes);

      /*!
       * @brief Return box number \b i.
       */
      const Box<DIM>& getBox(int i) const
      {
         return d_boxes[i];
      }

      /*!
       * @brief Set \b numbers to the numbers of the boxes that intersect
       * \b region, in increasing order.
       */
      void findIntersecting(const Box<DIM>& region,
                            std::vector<int>& numbers) const;

      /*!
       * @brief Return true if any box intersects \b region.
       */
      bool intersects(const Box<DIM>& region) const;

   private:
      // One size class.  Bin i holds bin_boxes[bin_start[i]..bin_start[i+1]);
      // dense grids have every bin, sparse grids only occupied ones, in bins.
      struct Grid
      {
         IntVector<DIM> bin_size;
         Box<DIM> bounding_box;
         IntVector<DIM> number_bins;
         bool dense;
         std::vector< IntVector<DIM> > bins;
         std::vector<int> bin_start;
         std::vector<int> bin_boxes;
      };

      // Sequences of at most this many boxes are scanned.
      static const int s_max_scanned_boxes = 32;

      /*
       * Append to candidates the boxes of grid's bins that meet region,
       * and return the number of bins that contributed.
       */
      int findCandidates(const Grid& grid, const Box<DIM>& region,
                         std::vector<int>& candidates) const;
      /*
       * Return the position of bin in grid's list of bins, or -1 if it is
       * not occupied.
       */
      static int findBin(const Grid& grid, const IntVector<DIM>& bin);
      static bool getBinRange(const Grid& grid, const Box<DIM>& region,
                              IntVector<DIM>& lower, IntVector<DIM>& upper);
      static int denseBinNumber(const Grid& grid, const IntVector<DIM>& bin);
      static bool lessBin(const IntVector<DIM>& a, const IntVector<DIM>& b);
      static bool lessEntry(const std::pair< IntVector<DIM>, int >& a,
                            const std::pair< IntVector<DIM>, int >& b);
      static bool nextBin(IntVector<DIM>& bin, const IntVector<DIM>& lower,
                          const IntVector<DIM>& upper);

      std::vector< Box<DIM> > d_boxes;
      bool d_use_bins;
      std::vector<Grid> d_grids;
   };

private:
   /*!
    * The copy ctor is not implemented.
    */
   BoxTree(const BoxTree<DIM> &box);

   /*!
    * The assignment operator is not implemented.
    */
   BoxTree<DIM>& operator=(const BoxTree<DIM>& box);

   /*!
    * Private routine to construct overlap indices.
    */
   void privateFindOverlapIndices(tbox::Array<int>& indices,
                                  const Box<DIM>& box,
                                  bool find_local_overlaps) const;

   /*
    * Build a tree of box_triples to fix their order, and index them in
    * that order.
    */
   void buildIndex(
      const tbox::Array<typename BoxTreeNode<DIM>::Triple>& box_triples,
      int min_length);

   /*
    * The boxes and periodic images of the tree, in the order in which the
    * tree lists them, and for each its index in d_boxes and its owner.
    */
   BinIndex d_index;
   std::vector<int> d_entry_indices;
   std::vector<int> d_entry_owners;

   BoxArray<DIM> d_boxes;

   bool d_have_mapping;

   /*!
    * This processor's rank.
    */
   int d_rank;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxTree.C"
#endif

#endif

