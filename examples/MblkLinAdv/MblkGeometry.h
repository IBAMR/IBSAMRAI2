//
// File:        MblkGeometry.h
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: set geometry for multiblock domain
//              
// 

#ifndef included_MblkGeometryXD
#define included_MblkGeometryXD

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "Box.h"
#include "BoxArray.h"
#include "CellData.h"
#include "tbox/Database.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include "BlockGridGeometry.h"

using namespace std;
using namespace SAMRAI;

/*!
 * This class creates the mapped multi-block grid geometry used
 * for calculations in the MblkLinAdv code.  The supported grid types
 * include Cartesian, Wedge, and Spherical shell.  The spherical shell
 * case is a full multi-block grid with 3 blocks.
 */
class MblkGeometry
{
public:

   /*!
    * Reads geometry information from the "MblkGeometry" input file
    * entry.                                                                 
    */
   MblkGeometry(const string& object_name,
                tbox::Pointer<tbox::Database> input_db,
                const int nblocks);
   
   ~MblkGeometry();
   
   /*!
    * Return the geometry type (CARTESIAN, WEDGE, or SPHERICAL_SHELL)
    */
   string getGeometryType();

   /*!
    * Return the user-specified refine boxes, given a block and 
    * level number
    */
   bool getRefineBoxes(hier::BoxArray<NDIM>& refine_boxes,
                       const int block_number,
                       const int level_number);
   
   /*!
    * Build mapped grid on patch.  The method defers the actual grid 
    * construction to private members, depending on the geometry 
    * choice in input.
    */
   void buildGridOnPatch(const hier::Patch<NDIM>& patch,
                         const hier::Box<NDIM>& domain,
                         const int xyz_id,
                         const int level_number,
                         const int block_number);

   /*!
    * Access the stored dx
    */

   void getDx(const hier::Box<NDIM>& domain,
              const int level_number,
              double* dx);

   void getDx(const int level_number,
              double* dx);

   /*!
    * Access the block rotation
    */
   int getBlockRotation(const int block_number);

   /*!
    * Tag cells for the octant problem.
    */
   void tagOctantCells(hier::Patch<NDIM>& patch,
                       const int xyz_id,
                       tbox::Pointer<pdat::CellData<NDIM,int> >& temp_tags,
                       const double regrid_time,
                       const int refine_tag_val);
   

private:

   /*
    * Read data members from input.  
    */
   void getFromInput(tbox::Pointer<tbox::Database> input_db,
                     bool is_from_restart);
   
   /*
    * Cartesian grid construction.
    */
   void setCartesianMetrics(const hier::Box<NDIM>& domain,
                            const int level_number);

   void buildCartesianGridOnPatch(const hier::Patch<NDIM>& patch,
                                  const int xyz_id,
                                  const int level_number,
                                  const int block_number);

   
   /*
    * Wedge grid construction.
    */
   void setWedgeMetrics(const hier::Box<NDIM>& domain,
                        const int level_number);

   void buildWedgeGridOnPatch(const hier::Patch<NDIM>& patch,
                              const int xyz_id,
                              const int level_number,
                              const int block_number);


   /*
    * Spherical shell grid construction
    */
   void setSShellMetrics(const hier::Box<NDIM>& domain,
                         const int level_number);

   void buildSShellGridOnPatch(const hier::Patch<NDIM>& patch,
                               const hier::Box<NDIM>& domain,
                               const int xyz_id,
                               const int level_number,
                               const int block_number);

   /*
    * For the spherical shell construction, i always points in the r direction
    * and j,k are points on the shell.  Given a certain j,k compute the 
    * unit sphere locations for block face (actual xyz is computed
    * by x = r*xface, y = r*yface, z = r*zface.  Note that the dimension
    * in the theta direction (nth) should be the same for each block.
    */
   void computeUnitSphereOctant(int nblock, 
                                int nth, 
                                int j, int k, 
                                double *xface, 
                                double *yface, 
                                double *zface );

   /*
    * Geometry type.  Choices are CARTESIAN, WEDGE, SPHERICAL_SHELL
    */
   string d_geom_problem;
   string d_object_name;
   

   /*
    * The number of blocks and the set of skelton grid geometries that make
    * up a multiblock mesh.
    */
   int d_nblocks;
   tbox::Array<bool> d_metrics_set;

   
   /*
    * The grid spacing.  For cartesian, d_dx = (dx,dy,dz).  For wedge,
    * d_dx = (dr, dth, dz). For spherical shell, d_dx = (dr, dth, dphi)
    */
   tbox::Array<tbox::Array<double> > d_dx;

   /*
    * Cartesian inputs
    */
   tbox::Array<tbox::Array<double> > d_cart_xlo;
   tbox::Array<tbox::Array<double> > d_cart_xhi;

   /*
    * Wedge inputs
    */
   tbox::Array<double> d_wedge_rmin;
   tbox::Array<double> d_wedge_rmax;
   double d_wedge_thmin;
   double d_wedge_thmax;
   double d_wedge_zmin;
   double d_wedge_zmax;
   
   /*
    * Shell inputs
    */
   double d_sshell_rmin;
   double d_sshell_rmax;

   // options are SOLID, OCTANT
   string d_sshell_type;

   /*
    * For tagging in the spherical octant case
    */
   double d_tag_velocity;
   double d_tag_width;


   // if SOLID, need to read in these...
   double d_sangle_degrees;
   double d_sangle_thmin;
   

   /*
    * Specify block rotation.
    */
   tbox::Array<int> d_block_rotation;

   /*
    * Refine boxes for different blocks/levels
    */
   tbox::Array< tbox::Array< hier::BoxArray<NDIM> > > d_refine_boxes;
   
};


#endif // included_MblkGeometry
