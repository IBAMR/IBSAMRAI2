//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/hierarchy/boxes/BoxUtilities.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Routines for processing boxes within a domain of index space.
//

#ifndef included_hier_BoxUtilities
#define included_hier_BoxUtilities

#include "SAMRAI_config.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "IntVector.h"
#include "tbox/Array.h"
#include "tbox/List.h"


namespace SAMRAI {
   namespace hier {

/**
 * Class BoxUtilities<DIM> provides several utility routines for processing
 * boxes or collections of boxes.  Many of these operations require
 * information about the location of the input boxes within some region of
 * index space (domain) or are used to compute this sort of information.
 * Often these routines are used in load balancing and communication routines
 * to determine the relationship of a box or set of boxes to some domain 
 * boundary, where the domain is specified by a list of box regions.
 *
 * The following provides an explanation some of the concepts
 * common to many of the functions in this class.
 *
 *     - \b cut point
 *
 *        A cut point for a coordinate direction is an integer that
 *        specifies the index value of the cell immediately to the
 *        right of the boundary along which a box will be cut.  For
 *        instance, if we have a cut point value of 12 in the y-coordinate
 *        direction for box [(0,0,0), (5,15,10)], then the box
 *        will be cut into [(0,0,0),(5,11,10)] and [(0,12,0),(5,15,10)].
 *        That is, the box has been cut along the boundary between
 *        the y = 11 and y = 12 cells.
 *
 *     - \b bad cut point
 *        If a box is cut at a "bad cut point", the resulting boxes
 *        will violate some of the box constraints. For example, a
 *        cut point that does not satify the cut factor restriction
 *        or violates the bad_interval constraint (see below)
 *        is a bad cut point.
 *
 *     - \b irregular boundary
 *        An irregular boundary results when the domain cannot be
 *        described by a single box and the domain is non-convex.
 *
 *
 * The following provides an explanation some of the arguments
 * common to many of the functions in this class.
 *
 *    - \b min_size
 *       min_size is a IntVector that specifies the minimum
 *       allowable box length in each dimension.  For example, if
 *       min_size = (10,4,15), then the minimum box length in the
 *       x, y, and z dimensions are 10, 4, and 15 respectively.
 *
 *    - \b max_size
 *       max_size is a IntVector that specifies the maximum
 *       allowable box length in each dimension.  For example, if
 *       max_size = (10,40,50), then the maximum box length in the
 *       x, y, and z dimensions are 10, 40, and 50 respectively.
 *
 *       It should be noted that the max_size constraint has lower
 *       priority than the other constraints.  In instances where
 *       all constraints cannot be simultaneously satisfied, the max_size
 *       constraint is sacrificed.
 *
 *    - \b cut_factor
 *       cut_factor is a IntVector that constrains the
 *       dimensions of a box to be multiples of the components of the
 *       cut_factor.  For instance, if cut_factor = (2,4,5), then the
 *       x, y, and z dimensions of a 8 box that satisfies the cut_factor
 *       constraint would be multiples of 2, 4, and 5 respectively.
 *
 *       This constraint is usually enforced with the cut_factor equal
 *       to a multiple of the refinement ratio between levels to ensure
 *       that if the box is coarsened, the resulting box is aligned with
 *       the coarse grid (it is assumed that the boundary of the fine box
 *       is aligned with a coarse cell boundary).
 *
 *    - \b bad_interval
 *       bad_interval is a IntVector that limits the distance
 *       a box can be from the boundary of the box list domain so that
 *       the outer boundaries of the box and the box list domain which are
 *       perpendicular to the i-th direction are no closer than the i-th
 *       component of bad_interval.  Another way to think of this is that
 *       the boundary of the box list domain in the i-th direction is not
 *       allowed to lie strictly within the interior of the box after
 *       it has been grown in the i-th direction by the i-th component
 *       of bad_interval.  For example, if bad_interval = (2,3,4), then
 *       the x, y, and z boundaries of the box must be at least 2, 3, and
 *       4 cells away from the x, y, and z boundaries of the box list
 *       domain.
 *
 *       The bad_interval constraint is enforced to avoid the situation
 *       where the ghost region for a box resides partially inside and
 *       partially outside the box list domain which complicates ghost
 *       cell filling.  In addition, this constraint avoids complicated
 *       issues with respect to the numerical accuracy of the solution.
 *
 *       Typically, bad_interval is based on the maximum ghost cell
 *       width over all patch data objects and some coarsen ratio.
 *
 *
 * Note that all member functions of this class are static.  The main intent
 * of the class is to group the functions into one name space.  Thus, you
 * should never attempt to instantiate a class of type BoxUtilities<DIM>;
 * simply call the functions as static functions; e.g.,
 * BoxUtilities<DIM>::function(...).  These routines are placed here rather
 * than in the box, box list, box array classes to avoid circular dependencies
 * among these classes.
 *
 * @see hier::Box
 * @see hier::BoxArray
 * @see hier::BoxList
 */

template<int DIM> struct BoxUtilities 
{
   /**
    * Check the given box for violation of minimum size, cut factor,
    * and box list domain constraints.  If a patch is generated from a box 
    * that violates any of these constraints, then some other routine 
    * (e.g., ghost cell filling, or inter-patch communication) may fail.  
    * Thus, this routine prints an error message describing the violation and 
    * executes a program abort.  
    *
    * Arguments:
    *
    *    - \b box (input) 
    *       box whose constraints are to be checked
    *
    *    - \b min_size (input) 
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    * 
    *       If there is no constraint on the box location within
    *       the box list domain, pass in an empty box array for the
    *       physical_boxes argument.
    *
    *    - \b physical_boxes (input) 
    *       box array representing the index space of box list domain
    *
    *
    * Assertion checks:
    * 
    *    - all components of min_size and bad_interval must be nonnegative
    *
    *    - all components of cut_factor must be positive 
    *
    */
   static void checkBoxConstraints(const Box<DIM>& box,
                                   const IntVector<DIM>& min_size,
                                   const IntVector<DIM>& cut_factor,
                                   const IntVector<DIM>& bad_interval,
                                   const BoxArray<DIM>& physical_boxes);
 
   /**
    * Replace each box in the list that is too large with a list of non-
    * overlapping smaller boxes whose union covers the same region of 
    * index space as the original box.  
    * 
    * Arguments:
    * 
    *    - \b boxes (input) 
    *       list of boxes to be chopped
    * 
    *    - \b max_size (input) 
    *       maximum allowed box size.  See class header for further
    *       description.
    * 
    *    - \b min_size (input) 
    *       minimum allowed box size.  See class header for further
    *       description.
    * 
    *    - \b cut_factor (input) 
    *       See class header for description.
    * 
    *    - \b bad_interval (input) 
    *       See class header for description.
    * 
    *    - \b physical_boxes (input) 
    *       box array representing the index space of box list domain
    *
    *
    * Notes:
    * 
    *    - The resulting boxes will obey the minimum size and cut factor 
    *      restrictions if the each of the original boxes does.  
    *
    *    - Any box with side length not equal to a multiple of the
    *      cut factor for that direction, will not be chopped along that
    *      direction.
    * 
    *    - The maximum size restriction may be sacrificed if the box 
    *      cannot be chopped at appropriate points.  However, this is 
    *      generally the case only when the box is adjacent to the 
    *      box list domain boundary and an irregular boundary configuration 
    *      restricts the cut locations or if the maximum size is not a 
    *      multiple of the cut factor.  
    * 
    */
   static void chopBoxes(BoxList<DIM>& boxes,
                         const IntVector<DIM>& max_size,
                         const IntVector<DIM>& min_size,
                         const IntVector<DIM>& cut_factor,
                         const IntVector<DIM>& bad_interval,
                         const BoxArray<DIM>& physical_boxes);

   /**
    * Chop the box into a collection of boxes according to the collection
    * of cut points specified along each coordinate direction.  Cut points 
    * that do not reside within the range of box indices are ignored.  
    *
    * Arguments:
    * 
    *    - \b boxes (output) 
    *       list of boxes into which the "box" argument was chopped
    * 
    *    - \b box (input) 
    *       box which is to be chopped
    * 
    *    - \b cut_points (input)
    *       cut_points is an array of integer lists, each of which 
    *       indicates the indices where the box will be cut in one of
    *       the coordinate directions
    * 
    *
    * Assertion checks:
    * 
    *    - The cut point array must have size equal to the number of 
    *      spatial dimensions for the box.
    *
    *    - The cut points for each direction must be on the list in 
    *      increasing order.
    *
    *
    * Notes:
    * 
    *    - The "boxes" BoxList is cleared before any box 
    *      operations are performed.  Thus, any boxes on the list when
    *      the function is called will be lost.  
    * 
    */
   static void chopBox(BoxList<DIM>& boxes,
                       const Box<DIM>& box,
                       const tbox::Array< tbox::List<int> > cut_points);

   /**
    * Extend the box in the list to domain boundary as needed so that
    * the domain boundary does not intersect the ghost cell region around
    * the box in an inappropriate manner.  Intersections that are 
    * disallowed are those in which a portion of the domain boundary is 
    * parallel to a box face and lies strictly in the interior of the ghost 
    * cell layer adjacent to that face.  In other words, we eliminate
    * ghost cell regions residing outside of a given domain and which are
    * narrower than the specified ghost width.   The boolean return value
    * is true if the input box was extended to the boundary and thus
    * is changed by the routine.  Otherwise, the return value is false.
    * 
    * See description of bad_interval in the class header comments for 
    * more details.
    * 
    * Arguments:
    * 
    *    - \b box (input/ouput)
    *       box to be extended 
    *
    *    - \b domain (input)
    *       some domain whose interior is the union of boxes in a list
    *
    *    - \b ext_ghosts (input)
    *       IntVector that specifies the size of the desired 
    *       ghost cell region 
    *
    *  
    * Assertion checks:
    * 
    *    - domain must not be an empty boxlist.
    *
    *    - All components of ext_ghosts must be nonnegative.
    *
    *
    * Notes:
    * 
    *    - The ext_ghosts argument often corresponds to the bad_interval 
    *      argument in many of the other functions in class.
    *    
    *    - This operation may produce overlap regions among boxes on the 
    *      list.  
    * 
    *    - There exist some bizarre domain configurations for which it is 
    *      impossible to grow a box to the boundary and eliminate bad 
    *      ghost region intersections.  This routine will extend each box as 
    *      far as it can, but will not remedy these degenerate situations 
    *      in general.  
    * 
    */
   static bool extendBoxToDomainBoundary(Box<DIM>& box,
                                         const BoxList<DIM>& domain,
                                         const IntVector<DIM>& ext_ghosts);

   /**
    * Same function as extendBoxToDomainBoundary() above except that it
    * extends each box in a list of boxes to the domain boundary specified
    * by the box list argument as needed.  The boolean return value 
    * is true if any box in the input box list was extended to the boundary 
    * and thus is changed by the routine.  Otherwise, the return value 
    * is false.
    */
   static bool extendBoxesToDomainBoundary(BoxList<DIM>& boxes,
                                           const BoxList<DIM>& domain,
                                           const IntVector<DIM>& ext_ghosts);

   /**
    * Grow each box in the list that is smaller than the specified minimum
    * size.  
    * 
    * Arguments:
    * 
    *    - \b boxes (input/output)
    *       list of boxes to be grown to satisfy the min_size constraint 
    * 
    *    - \b domain (input)
    *       list of boxes whose union is some domain
    * 
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       description.
    * 
    *
    * Assertion checks:
    * 
    *    - The minimum box size argument must have all nonnegative entries. 
    *
    *
    * Notes:
    * 
    *    - Each box that is grown must remain within the union of the
    *      boxes of the given domain.  
    *
    *    - If the specified domain is an empty box list, then each box 
    *      will be grown to be as large as the minimum size with no 
    *      particular restrictions applied.  
    *
    *    - This operation may produce overlap regions among boxes on list 
    * 
    *    - There exist some bizarre domain configurations for which it is 
    *      impossible to grow a box sufficiently within the domain.  
    *
    *      For instance if the domain is given by 
    *      [(0,0),(2,10)], [(0,3),(1,4)], [(0,5),(10,10)] 
    *      and the box is given by [(4,1),(6,2)] with a minimum size
    *      of (4,4), there is no way the box can be grown to the minimum
    *      size without have to "cross" the gap in the box list domain.
    *
    *      This routine will grow each box as far as it can, but will not 
    *      remedy these situations, generally. 
    */
   static void growBoxesWithinDomain(BoxList<DIM>& boxes,
                                     const BoxList<DIM>& domain,
                                     const IntVector<DIM>& min_size);

   /**
    * Determine whether the box may be chopped according to specified
    * min_size, max_size and cut_factor constraints.  For those 
    * directions along which the box may be chopped, the cut points are 
    * computed.  The cut points for the j-th coordinate direction are 
    * placed into a list of integers corresponding to the j-th component 
    * of the cut_point array.
    *
    * Return value:  
    * 
    *    - true is returned if the box may be chopped along any 
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    * 
    *    - \b cut_points (output)
    *       array of list of cut points for the box
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b max_size(input)
    *       minimum allowed box size.  See class header for further 
    *       description.
    *
    *    - \b min_size(input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b cut_factor(input)
    *       See class header for description.
    *
    *
    * Assertion checks:
    * 
    *    - all components of min_size must be nonnegative
    *
    *    - all components of max_size and cut_factor must be positive 
    * 
    *    - for each i between 0 and (DIM-1), the i-th component of 
    *      min_size must be less than or equal to the i-th component 
    *      of max_size
    * 
    */
   static bool findBestCutPointsGivenMax(
      tbox::Array< tbox::List<int> >& cut_points, 
      const Box<DIM>& box, 
      const IntVector<DIM>& max_size, 
      const IntVector<DIM>& min_size, 
      const IntVector<DIM>& cut_factor);

   /**
    * Determine whether the box may be chopped according to specified
    * min_size, max_size and cut_factor constraints along given
    * coordinate direction.  If the box may be chopped, the cut points 
    * are computed and placed into a list of integers.
    *
    * Return value:  
    * 
    *    - true is returned if the box may be chopped along the specified
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    * 
    *    - \b idir (input)
    *       coordinate direction along which cut points will be computed
    *
    *    - \b cut_points (output)
    *       list of cut points for the box along the idir coordinate 
    *       direction
    *
    *    - \b box (input)
    *       box to be chopped
    *
    *    - \b max_size (input)
    *       maximum allowed box size in idir coordinate direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size in idir coordinate direction.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * 
    * Assertion checks:
    * 
    *    - min_size must be nonnegative
    * 
    *    - max_size and cut_factor must be positive
    * 
    *    - min_size must be less than or equal to max_size 
    * 
    */
   static bool findBestCutPointsForDirectionGivenMax(
      const int idir,
      tbox::List<int>& cut_points,
      const Box<DIM>& box, 
      const int max_size,
      const int min_size,
      const int cut_factor);

   /**
    * Determine whether the box may be chopped into the specified
    * number of boxes along each coordinate direction.  For those 
    * directions along which the box may be chopped, the cut points are
    * computed.  The cut points for the j-th coordinate direction are 
    * placed into a list of integers corresponding to the j-th component 
    * of the cut_point array.
    *
    * Return value:  
    * 
    *    - true is returned if the box may be chopped along any
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    *
    *    - \b cut_points (output)
    *       array of list of cut points for the box
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b number_boxes (input)
    *       the i-th component of number_boxes specifies the desired 
    *       number of cuts to be made along the i-th coordinate 
    *       direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size. See class header for further
    *       description.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *
    * Important note: By convention, each integer cut point that is computed
    *                 corresponds to the cell index to the right of cut point.
    *
    * Assertion checks:
    * 
    *    - all components of min_size must be nonnegative
    *
    *    - all components of cut_factor and number_boxes must be positive
    *
    */
   static bool findBestCutPointsGivenNumber(
      tbox::Array< tbox::List<int> >& cut_points,
      const Box<DIM>& box,
      const IntVector<DIM>& number_boxes,
      const IntVector<DIM>& min_size,
      const IntVector<DIM>& cut_factor);

   /**
    * Determine whether the box may be chopped into the specified
    * number of boxes along along given coordinate direction.  If the 
    * box may be chopped, the cut points are computed and placed 
    * into a list of integers.
    * 
    * Return value:  
    * 
    *    - true is returned if the box may be chopped along the specified
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    * 
    *    - \b idir (input)
    *       coordinate direction along which cut points will be computed
    *
    *    - \b cut_points (output)
    *       list of cut points for the box along the idir coordinate 
    *       direction
    *
    *    - \b box (input)
    *       box to be chopped
    *
    *    - \b num_boxes (input)
    *       num_boxes specifies the desired number of cuts to be made 
    *       along the idir coordinate direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size in idir coordinate direction.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * 
    * Assertion checks:
    * 
    *    - min_size must be nonnegative
    * 
    *    - cut_factor and num_boxes must be positive
    * 
    */
   static bool findBestCutPointsForDirectionGivenNumber(
      const int idir,
      tbox::List<int>& cut_points,
      const Box<DIM>& box,
      const int num_boxes,
      const int min_size,
      const int cut_factor);

   /**
    * Determine whether box has any bad cut points based on its
    * position within the box list domain.  Information about the potentially
    * bad directions is returned in the IntVector 
    * bad_cut_information.  An entry of zero indicates that there are no 
    * bad cut points for the box along that coordinate direction.  An entry 
    * of one indicates that there may be a bad cut point along that direction.
    *  
    * Return value:
    * 
    *    - true is returned if the box may potentially have a bad point 
    *      along some coordinate direction. Otherwise false is returned.  
    * 
    * Arguments:
    * 
    *    - \b bad_cut_information (output)
    *       A value of 0 in the i-th component of bad_cut_information 
    *       indicates that there are no bad cut points in the i-th
    *       coordinate direction.
    *
    *    - \b box (input)
    *       box to be cut.
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    *
    * Assertion checks:
    * 
    *    - all components of bad_interval must be nonnegative
    * 
    */
   static bool checkBoxForBadCutPoints(
      IntVector<DIM>& bad_cut_information,
      const Box<DIM>& box,
      const BoxArray<DIM>& physical_boxes,
      const IntVector<DIM>& bad_interval);

   /**
    * Determine whether box may have any bad cut points along the specified
    * coordinate direction based on its position within the box array domain.  
    * 
    * Return value:
    * 
    *    - true is returned if the box may potentially have a bad point;
    *      otherwise false is returned.
    *
    * Arguments:
    * 
    *    - \b dir (input)
    *       coordinate direction to be checked for bad cut points
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    *
    * Assertion checks:
    * 
    *    - all components of bad_interval must be nonnegative.
    *
    */
   static bool checkBoxForBadCutPointsInDirection(
      const int dir,
      const Box<DIM>& box,
      const BoxArray<DIM>& physical_boxes,
      const IntVector<DIM>& bad_interval);

   /**
    * Determine bad cut points for box based on the specified box array domain
    * and bad interval.  
    *
    * The cut information is returned as an array (size = DIM) of arrays
    * (size = number of cells along edge of box) of boolean values.  A
    * false value indicates a good cut point, a true value indicates that 
    * the box should not be cut at that point. 
    * 
    * Arguments:
    * 
    *    - \b bad_cuts (output)
    *        stores an array of boolean arrays that indicates whether
    *        a potential cut point is bad.  A value of false indicates 
    *        a good cut point, and a true value indicates a bad cut point.
    *
    *    - \b box (input)
    *        box to be cut
    *
    *    - \b physical_boxes (input)
    *        box array that represents some domain
    *    
    *    - \b bad_interval (input)
    *        See class header for description.
    *
    *
    * Assertion checks:
    * 
    *    - all components of bad_interval must be nonnegative
    * 
    *    - bad_cuts must have size equal to DIM
    *
    */
   static void findBadCutPoints(
      tbox::Array< tbox::Array<bool> >& bad_cuts, 
      const Box<DIM>& box,
      const BoxArray<DIM>& physical_boxes, 
      const IntVector<DIM>& bad_interval);

   /**
    * Find bad cut points for a box given a single coordinate direction.
    * The cut information is returned as an array of boolean values
    * (size = number of cells along specified edge of box).  A false value 
    * indicates a good cut point, a true value indicates that the box should 
    * not be cut at that point.
    *
    * Arguments:
    * 
    *    - \b dir (input)
    *       coordinate direction to be checked for bad cut points
    *
    *    - \b bad_cuts (output)
    *        boolean arrays whose entries indicates whether
    *        a potential cut point is bad.  
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    *
    * Assertion checks:
    * 
    *    - all components of bad_interval must be nonnegative
    *
    */
   static void findBadCutPointsForDirection(
      const int dir,
      tbox::Array<bool>& bad_cuts,
      const Box<DIM>& box,
      const BoxArray<DIM>& physical_boxes, 
      const IntVector<DIM>& bad_interval);

   /**
    * Given a set of potential cut points and a set of bad cut points for
    * a box, adjust the cut points so that they do not coincide with 
    * bad cut points.  Typically, the cuts are generated using either of 
    * the findBestCutPoints...() functions, and the bad cut points are 
    * generated using the findBadCutPoints() function.
    *
    * Arguments:
    * 
    *    - \b cuts (input/output)
    *       array of integer lists each of which holds a list of 
    *       cut points for the box.  Each list is adjusted so that
    *       no cut points coincide with bad cut points
    *
    *    - \b bad_cuts (input)
    *       array of boolean arrays each of which stores information
    *       about which offsets from the lower corner of the box 
    *       are bad cut points
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       details.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * 
    * Assertion checks:
    * 
    *    - All components of min_size must be nonnegative.
    *
    *    - All components of cut_factor must be positive.
    * 
    *    - cuts and bad_cuts must have size equal to DIM
    * 
    *    - Each array of integers in bad_cuts must have length
    *      equal to the number of cells for the box in that dimension.
    *    
    *    - The cut points for each direction must be strictly increasing
    *      and all satisfy the cut_factor restriction.
    *
    */
   static void fixBadCutPoints(
      tbox::Array< tbox::List<int> >& cuts, 
      const tbox::Array< tbox::Array<bool> >& bad_cuts, 
      const Box<DIM>& box,
      const IntVector<DIM>& min_size, 
      const IntVector<DIM>& cut_factor);

   /**
    * Given a set of potential cut points and a set of bad cut points for
    * a box, adjust the cut points in the specified coordinate direction
    * so that they do not coincide with bad cut points.  Typically, the 
    * cuts are generated using either of the findBestCutPoints...() 
    * functions, and the bad cut points are generated using the 
    * findBadCutPoints() function.
    *
    * Arguments:
    * 
    *    - \b dir (input)
    *       coordinate direction along which to fix cut points
    *
    *    - \b cuts (input/output)
    *       list of integers which holds a list of cut points for the box.  
    *       This list is adjusted so that no cut points coincide with bad 
    *       cut points.
    *
    *    - \b bad_cuts (input)
    *       array of booleans which stores information about which 
    *       offsets from the lower corner of the box are bad cut points
    *       
    *    - \b box (input) 
    *       box to be cut
    * 
    *    - \b min_size (input) 
    *       minimum allowed box size along specified coordinate direction.
    *
    *    - \b cut_factor (input) 
    *       See class header for description.
    *
    * 
    * Assertion checks:
    * 
    *    - min_size must be nonnegative. 
    *
    *    - cut_factor must be positive. 
    *
    *    - bad_cut_points must have size equal to the number of
    *      cells in the box along the specified coordinate direction.
    *
    *    - The cut points must be strictly increasing and all
    *      satisfy the cut_factor constraint.
    *
    */
   static void fixBadCutPointsForDirection(
      const int dir,
      tbox::List<int>& cuts, 
      const tbox::Array<bool>& bad_cuts, 
      const Box<DIM>& box,
      const int min_size, 
      const int cut_factor);

   /**
     *                                                                       
     * This static private member function is called by findBadCutPoints(),  
     * and the findBadCutPointsForDirection() member functions.  It sets bad 
     * cut points near the lower and upper ends of the border box in the     
     * given coordinate direction.                                           
     *                                                                       
     */
   static void findBadCutPointsForBorderAndDirection(
      const int id,
      tbox::Array<bool>& bad_cuts,
      const Box<DIM>& box,
      const Box<DIM>& border,
      const int bad_interval);


   /**
    * Construct an array of box lists so that each list contains
    * a non-overlapping set of boxes covering some portion of the
    * box at the same array location in the box array.  The regions 
    * of index space formed by composing the union of boxes on each 
    * box list are mutually disjoint and the union of all boxes in 
    * the box lists exactly covers the union of boxes in the original
    * box array.  In other words, this routine partitions the boxes 
    * in the "boxes" argument into a set of non-overlapping box collections.  
    * If none of the boxes in this box array overlap then each box list 
    * in the resulting array has a single box equal to the corresponding 
    * box in the box array.  This routine is especially useful for 
    * determining a unique set of index points given an array of boxes 
    * in some index space.
    * 
    * Arguments:
    * 
    *    - \b box_list_array (output) 
    *       array of box lists which cover mutually exclusive portions
    *       of the index space covered by the "boxes" argument
    *
    *    - \b boxes (input) 
    *       an arbitrary box array
    */
   static void makeNonOverlappingBoxLists(
      tbox::Array< BoxList<DIM> >& box_list_array,
      const BoxArray<DIM>& boxes);

};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxUtilities.C"
#endif

#endif

