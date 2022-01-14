//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/apputils/embedded_boundary/EmbeddedBoundaryDefines.h $
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Definitions of cells and nodes on embedded boundaries.
// 

#ifndef included_appu_EmbeddedBoundaryDefines
#define included_appu_EmbeddedBoundaryDefines

#include "SAMRAI_config.h"

namespace SAMRAI {
   namespace appu {

/*!
 * @brief Class EmbeddedBoundaryDefines sets the enumerated types used
 * to define cells and nodes on the embedded boundary level.
 * 
 * @see appu::EmbeddedBoundaryGeometry 
 * @see appu::CutCell 
 * @see appu::BoundaryNode 
 */

class EmbeddedBoundaryDefines 
{
public:

   /*!
    * Enumerated type for the different cell classifications.
    *
    * - \b SOLID          {Cell is located in the "solid" region.}
    * - \b CUT            {Cell is cut, meaning a CutCell data
    *                      structure will be maintained at this cell.}
    * - \b BORDER         {Cell neighbors a cut cell, in the "flow" region.}
    * - \b FLOW           {Cell is located in the "flow" region.}
    */
   enum CELL_TYPE{ SOLID  = 0,
                   CUT    = 1,
                   BORDER = 2,
                   FLOW   = 3 };

  /*!
    * Enumerated type for inside/outside node classification.
    *
    * - \b INSIDE         {Node is located "inside" the prescribed geometry.}
    * - \b OUTSIDE        {Node is outside the geometry.}
    * - \b BOUNDARY       {Node is on the boundary of the geometry.  That is
    *                      it is the first one "inside" the geometry.}
    * - \b ONBOUNDARY     {Node is located exactly on the boundary of the
    *                      geometry (used to avoid divide-by-zero problems)
    *                      in numerical operations at embedded boundary.}
    */
   enum NODE_TYPE{ OUTSIDE    = 0,
                   INSIDE     = 1,
                   BOUNDARY   = 2,
                   ONBOUNDARY = -9};
};

}
}
#endif


