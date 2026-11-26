//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxTree.C $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTree_C
#define included_hier_BoxTree_C

#include "BoxTree.h"
#include "BoxGraphUtilities.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <iostream>
#include <iomanip>

// Ignore incorrect Intel warning
#ifdef __INTEL_COMPILER
#pragma warning (disable:177)
#endif

namespace SAMRAI {
   namespace hier {

/*
 * ************************************************************************
 *
 * ctors and dtor
 *
 * ************************************************************************
 */

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
   const ProcessorMapping & box_mapping,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(true)
{
   tbox::Array<int> indices;
   BoxArray<DIM>  boxes;
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      boxes, indices, box_array, box_shifts);

   int len = boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      int idx = indices[i];
      box_triples[i].box = boxes[i];
      box_triples[i].idx = idx;
      box_triples[i].owner = box_mapping.getProcessorAssignment(idx);
   }

   buildIndex(box_triples, min_length);
}

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(false)
{
   tbox::Array<int> indices;
   BoxArray<DIM>  boxes;
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      boxes, indices, box_array, box_shifts);

   int len = boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      int idx = indices[i];
      box_triples[i].box = boxes[i];
      box_triples[i].idx = idx;
      box_triples[i].owner = -1;
   }

   buildIndex(box_triples, min_length);
}

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(false)
{
   int len = d_boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      box_triples[i].box = d_boxes[i];
      box_triples[i].idx = i;
      box_triples[i].owner = -1;
   }

   buildIndex(box_triples, min_length);
}


template<int DIM>  BoxTree<DIM>::~BoxTree()
{
}

/*
 *************************************************************************
 *
 * Routines to return an array of integer indices indicating elements
 * of the box array passed to the ctor that overlap with the specified box.  
 * Depending on which routine is called, the indices will either correspond 
 * to all boxes in the original box array or only those that are local to 
 * this processor. 
 *
 *************************************************************************
 */

template<int DIM> void BoxTree<DIM>::findOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box) const
{
   bool find_local_overlap_indices = false;
   privateFindOverlapIndices(indices, box, find_local_overlap_indices);
}

template<int DIM> void BoxTree<DIM>::findLocalOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box) const
{
   if (!d_have_mapping) {
      TBOX_ERROR("BoxTree<DIM>::findLocalOverlapIndices() error!!"
                 << "\n Processor mapping must be passed to class constructor"
                 << "to use this method. None was passed." << std::endl);
   }

   bool find_local_overlap_indices = true;
   privateFindOverlapIndices(indices, box, find_local_overlap_indices);
}

template<int DIM> void BoxTree<DIM>::privateFindOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box,
   bool find_local_overlaps) const
{
   std::vector<int> entries;
   d_index.findIntersecting(box, entries);

   std::vector<int> found;
   found.reserve(entries.size());
   for (std::size_t e = 0; e < entries.size(); ++e) {
      if (!find_local_overlaps || d_entry_owners[entries[e]] == d_rank) {
         found.push_back(d_entry_indices[entries[e]]);
      }
   }
   std::sort(found.begin(), found.end());
   found.erase(std::unique(found.begin(), found.end()), found.end());

   indices.setNull();
   indices.resizeArray(static_cast<int>(found.size()));
   for (std::size_t i = 0; i < found.size(); ++i) {
      indices[static_cast<int>(i)] = found[i];
   }
}

/*
 *************************************************************************
 *
 * Routine that returns a list of boxes that overlap with the argument box.
 *
 *************************************************************************
 */

template<int DIM> void BoxTree<DIM>::findOverlapBoxes(
   BoxList<DIM>& overlap_boxes,
   const Box<DIM>& box) const
{
   overlap_boxes.clearItems();
   std::vector<int> entries;
   d_index.findIntersecting(box, entries);
   for (std::size_t e = 0; e < entries.size(); ++e) {
      overlap_boxes.appendItem(d_index.getBox(entries[e]));
   }
}

/*
 * **************************************************************
 * CAUTION: the semantics of this call differ from that of
 * BoxList::removeIntersections(const BoxList<DIM> takeaway);
 * here, the list that is being modified is the list that is
 * passed as an argument; the "takeaway" list is the list that
 * was passed when the BoxTop object was constructed.
 * **************************************************************
 */

template<int DIM> void BoxTree<DIM>::removeIntersections(BoxList<DIM>& fragments) const
{
   //This will hold boxes in "list" minus the intersection with "takeaway"
   BoxArray<DIM> boxes_in(fragments);
   fragments.clearItems();

   /*
    * Compare each box in "list" against all boxes in "take_away"
    * Note: this is opposite the approach of the current implementation;
    * here, the "take_away" list is the list with which the tree was
    * instantiated.
    */
   int len = boxes_in.getNumberOfBoxes();
   for (int k=0; k<len; ++k) {
      Box<DIM> tryme = boxes_in[k];

      /*
       * Find a shorter list of boxes from "takeaway" that might intersect
       * with "tryme."  Recall that "takeaway" was the BoxList with which
       * this object was instantiated; the next call returns a list of
       * boxes that are in "takeaway," and are also nabors (i.e., overlap with)
       * the "tryme" box.
       */
      BoxList<DIM> nabors;
      findOverlapBoxes(nabors, tryme);
      int nabor_count = nabors.size();

      /*
       * if "tryme" doesn't intersect any boxes on takeaway
       * i.e., if "nabors" is empty" then keep "tryme."
       */
      if (nabor_count == 0) {
        fragments.appendItem(tryme);
      } 

      /*
       * Else, burst up tryme against the nabors.
       * Here, we call the BoxList<DIM> version of
       * removeIntersections.
       */
      else {

         //possibly todo: hand-code a removeIntersections arrays ...
         BoxList<DIM> tryme_burst(tryme);
         tryme_burst.removeIntersections(nabors);
         fragments.copyItems(tryme_burst);
      }
   }
}



/*
 * Build the tree only to list the triples in the order of its queries,
 * and index them in that order.
 */
template<int DIM> void BoxTree<DIM>::buildIndex(
   const tbox::Array<typename BoxTreeNode<DIM>::Triple>& box_triples,
   int min_length)
{
   d_rank = tbox::SAMRAI_MPI::getRank();
   std::vector<typename BoxTreeNode<DIM>::Triple> triples;
   triples.reserve(box_triples.size());
   {
      const BoxTreeNode<DIM> tree(box_triples, d_have_mapping, min_length);
      tree.appendTriplesInQueryOrder(triples);
   }

   std::vector< Box<DIM> > boxes;
   boxes.reserve(triples.size());
   d_entry_indices.reserve(triples.size());
   d_entry_owners.reserve(triples.size());
   for (std::size_t i = 0; i < triples.size(); ++i) {
      boxes.push_back(triples[i].box);
      d_entry_indices.push_back(triples[i].idx);
      d_entry_owners.push_back(triples[i].owner);
   }
   d_index.initialize(boxes);
}

/*
 *************************************************************************
 *
 * BinIndex
 *
 *************************************************************************
 */

template<int DIM> BoxTree<DIM>::BinIndex::BinIndex()
:  d_use_bins(false)
{
}

template<int DIM> BoxTree<DIM>::BinIndex::BinIndex(
   const BoxList<DIM>& boxes)
:  d_use_bins(false)
{
   std::vector< Box<DIM> > box_vector;
   box_vector.reserve(boxes.getNumberOfItems());
   for (typename BoxList<DIM>::Iterator b(boxes); b; b++) {
      box_vector.push_back(b());
   }
   initialize(box_vector);
}

/*
 * Class l holds the boxes no larger than base * 2^l in every direction and
 * larger than base * 2^(l-1) in some; base(d) is the smallest extent in d.
 */
template<int DIM> void BoxTree<DIM>::BinIndex::initialize(
   const std::vector< Box<DIM> >& boxes)
{
   d_boxes = boxes;
   d_use_bins = false;
   d_grids.clear();

   const int num_boxes = static_cast<int>(d_boxes.size());
   if (num_boxes <= s_max_scanned_boxes) {
      return;
   }

   IntVector<DIM> base(0);
   bool have_box = false;
   for (int i = 0; i < num_boxes; ++i) {
      if (!d_boxes[i].empty()) {
         if (!have_box) {
            base = d_boxes[i].numberCells();
            have_box = true;
         } else {
            base.min(d_boxes[i].numberCells());
         }
      }
   }
   if (!have_box) {
      return;
   }

   std::vector<int> size_class(num_boxes, -1);
   int number_classes = 0;
   for (int i = 0; i < num_boxes; ++i) {
      if (d_boxes[i].empty()) {
         continue;
      }
      int l = 0;
      for (int d = 0; d < DIM; ++d) {
         while ((static_cast<long>(base(d)) << l) < d_boxes[i].numberCells(d)) {
            ++l;
         }
      }
      size_class[i] = l;
      number_classes = tbox::MathUtilities<int>::Max(number_classes, l + 1);
   }

   std::vector< std::pair< IntVector<DIM>, int > > entries;
   for (int l = 0; l < number_classes; ++l) {
      Grid grid;
      for (int i = 0; i < num_boxes; ++i) {
         if (size_class[i] == l) {
            grid.bounding_box += d_boxes[i];
         }
      }
      if (grid.bounding_box.empty()) {
         continue;
      }
      double total_bins = 1.0;
      for (int d = 0; d < DIM; ++d) {
         grid.bin_size(d) = base(d) << l;
         const int extent = grid.bounding_box.numberCells(d);
         grid.number_bins(d) =
            (extent + grid.bin_size(d) - 1) / grid.bin_size(d);
         total_bins *= grid.number_bins(d);
      }

      entries.clear();
      IntVector<DIM> lower(0), upper(0);
      for (int i = 0; i < num_boxes; ++i) {
         if (size_class[i] == l) {
            getBinRange(grid, d_boxes[i], lower, upper);
            IntVector<DIM> bin(lower);
            do {
               entries.push_back(std::make_pair(bin, i));
            } while (nextBin(bin, lower, upper));
         }
      }
      std::sort(entries.begin(), entries.end(), lessEntry);

      grid.dense = total_bins <= 2.0 * static_cast<double>(entries.size());
      grid.bin_boxes.resize(entries.size());
      if (grid.dense) {
         const int number_bins = static_cast<int>(total_bins);
         grid.bin_start.assign(number_bins + 1, 0);
         for (std::size_t e = 0; e < entries.size(); ++e) {
            ++grid.bin_start[denseBinNumber(grid, entries[e].first) + 1];
         }
         for (int b = 0; b < number_bins; ++b) {
            grid.bin_start[b + 1] += grid.bin_start[b];
         }
      } else {
         for (std::size_t e = 0; e < entries.size(); ++e) {
            if (e == 0 || lessBin(entries[e - 1].first, entries[e].first)) {
               grid.bins.push_back(entries[e].first);
               grid.bin_start.push_back(static_cast<int>(e));
            }
         }
         grid.bin_start.push_back(static_cast<int>(entries.size()));
      }
      for (std::size_t e = 0; e < entries.size(); ++e) {
         grid.bin_boxes[e] = entries[e].second;
      }
      d_grids.push_back(grid);
   }
   d_use_bins = true;
}

template<int DIM> void BoxTree<DIM>::BinIndex::findIntersecting(
   const Box<DIM>& region,
   std::vector<int>& numbers) const
{
   numbers.clear();
   if (!d_use_bins) {
      const int num_boxes = static_cast<int>(d_boxes.size());
      for (int i = 0; i < num_boxes; ++i) {
         if (d_boxes[i].intersects(region)) {
            numbers.push_back(i);
         }
      }
      return;
   }

   std::vector<int> candidates;
   int contributing_bins = 0;
   for (std::size_t g = 0; g < d_grids.size(); ++g) {
      contributing_bins += findCandidates(d_grids[g], region, candidates);
   }
   numbers.reserve(candidates.size());
   for (std::size_t c = 0; c < candidates.size(); ++c) {
      if (d_boxes[candidates[c]].intersects(region)) {
         numbers.push_back(candidates[c]);
      }
   }
   if (contributing_bins > 1) {
      // A box may meet several bins, and grids list their boxes separately.
      std::sort(numbers.begin(), numbers.end());
      numbers.erase(std::unique(numbers.begin(), numbers.end()),
                    numbers.end());
   }
}

template<int DIM> bool BoxTree<DIM>::BinIndex::intersects(
   const Box<DIM>& region) const
{
   if (!d_use_bins) {
      for (std::size_t i = 0; i < d_boxes.size(); ++i) {
         if (d_boxes[i].intersects(region)) {
            return true;
         }
      }
      return false;
   }

   std::vector<int> candidates;
   for (std::size_t g = 0; g < d_grids.size(); ++g) {
      candidates.clear();
      findCandidates(d_grids[g], region, candidates);
      for (std::size_t c = 0; c < candidates.size(); ++c) {
         if (d_boxes[candidates[c]].intersects(region)) {
            return true;
         }
      }
   }
   return false;
}

template<int DIM> int BoxTree<DIM>::BinIndex::findCandidates(
   const Grid& grid,
   const Box<DIM>& region,
   std::vector<int>& candidates) const
{
   IntVector<DIM> lower(0), upper(0);
   if (!getBinRange(grid, region, lower, upper)) {
      return 0;
   }
   int contributing_bins = 0;
   double range_bins = 1.0;
   for (int d = 0; d < DIM; ++d) {
      range_bins *= upper(d) - lower(d) + 1;
   }
   if (grid.dense ||
       range_bins <= static_cast<double>(grid.bins.size())) {
      IntVector<DIM> bin(lower);
      do {
         const int b = findBin(grid, bin);
         if (b >= 0 && grid.bin_start[b] < grid.bin_start[b + 1]) {
            candidates.insert(candidates.end(),
                              grid.bin_boxes.begin() + grid.bin_start[b],
                              grid.bin_boxes.begin() + grid.bin_start[b + 1]);
            ++contributing_bins;
         }
      } while (nextBin(bin, lower, upper));
   } else {
      // The region spans more bins than are occupied: test each occupied bin.
      for (std::size_t b = 0; b < grid.bins.size(); ++b) {
         bool inside = true;
         for (int d = 0; d < DIM; ++d) {
            inside = inside && lower(d) <= grid.bins[b](d) &&
                     grid.bins[b](d) <= upper(d);
         }
         if (inside) {
            candidates.insert(candidates.end(),
                              grid.bin_boxes.begin() + grid.bin_start[b],
                              grid.bin_boxes.begin() + grid.bin_start[b + 1]);
            ++contributing_bins;
         }
      }
   }
   return contributing_bins;
}

template<int DIM> int BoxTree<DIM>::BinIndex::findBin(
   const Grid& grid,
   const IntVector<DIM>& bin)
{
   if (grid.dense) {
      return denseBinNumber(grid, bin);
   }
   typename std::vector< IntVector<DIM> >::const_iterator it =
      std::lower_bound(grid.bins.begin(), grid.bins.end(), bin, lessBin);
   if (it == grid.bins.end() || lessBin(bin, *it)) {
      return -1;
   }
   return static_cast<int>(it - grid.bins.begin());
}

template<int DIM> int BoxTree<DIM>::BinIndex::denseBinNumber(
   const Grid& grid,
   const IntVector<DIM>& bin)
{
   int b = bin(DIM - 1);
   for (int d = DIM - 2; d >= 0; --d) {
      b = b * grid.number_bins(d) + bin(d);
   }
   return b;
}

/*
 * Set lower and upper to the range of bins of grid that meet region, and
 * return false if region misses the grid's bounding box.
 */
template<int DIM> bool BoxTree<DIM>::BinIndex::getBinRange(
   const Grid& grid,
   const Box<DIM>& region,
   IntVector<DIM>& lower,
   IntVector<DIM>& upper)
{
   const Box<DIM> clipped = region * grid.bounding_box;
   if (clipped.empty()) {
      return false;
   }
   for (int d = 0; d < DIM; ++d) {
      lower(d) = (clipped.lower(d) - grid.bounding_box.lower(d)) /
                 grid.bin_size(d);
      upper(d) = (clipped.upper(d) - grid.bounding_box.lower(d)) /
                 grid.bin_size(d);
   }
   return true;
}

/*
 * Order bins by their last index first, as dense bin numbers are ordered.
 */
template<int DIM> bool BoxTree<DIM>::BinIndex::lessBin(
   const IntVector<DIM>& a,
   const IntVector<DIM>& b)
{
   for (int d = DIM - 1; d >= 0; --d) {
      if (a(d) != b(d)) {
         return a(d) < b(d);
      }
   }
   return false;
}

template<int DIM> bool BoxTree<DIM>::BinIndex::lessEntry(
   const std::pair< IntVector<DIM>, int >& a,
   const std::pair< IntVector<DIM>, int >& b)
{
   if (lessBin(a.first, b.first)) {
      return true;
   }
   if (lessBin(b.first, a.first)) {
      return false;
   }
   return a.second < b.second;
}

/*
 * Advance bin through the range lower to upper, first index fastest, and
 * return false after the last bin.
 */
template<int DIM> bool BoxTree<DIM>::BinIndex::nextBin(
   IntVector<DIM>& bin,
   const IntVector<DIM>& lower,
   const IntVector<DIM>& upper)
{
   for (int d = 0; d < DIM; ++d) {
      if (bin(d) < upper(d)) {
         ++bin(d);
         return true;
      }
      bin(d) = lower(d);
   }
   return false;
}

}
}

#endif
