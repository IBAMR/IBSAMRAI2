/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/test/FAC/PoissonPolynomialSolution.h $
  Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1993 $
  Modified:	$LastChangedDate: 2008-02-19 08:24:52 -0800 (Tue, 19 Feb 2008) $
  Description:	PoissonPolynomialSolution class declaration
*/

#ifndef included_PoissonPolynomialSolution
#define included_PoissonPolynomialSolution


#include "QuarticFcn.h"


#include <string>
#include "tbox/Database.h"


/*
  SAMRAI classes
*/
#include "CellData.h"
#include "SideData.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"


using namespace SAMRAI;


/*!
  @brief Specialized class to provide polynomial-solution-specific stuff.
*/
class PoissonPolynomialSolution :
  public solv::RobinBcCoefStrategy<NDIM>
{

public:

  PoissonPolynomialSolution();

  PoissonPolynomialSolution(
    /*! Ojbect name */
    const string &object_name
  , /*! Input database */
    tbox::Database &database
  , /*! Standard output stream */ ostream *out_stream=NULL
  , /*! Log output stream */ ostream *log_stream=NULL
  );

  virtual ~PoissonPolynomialSolution();

  void setFromDatabase( tbox::Database &database );

  void setPoissonSpecifications(
    /*! Object to set */ solv::PoissonSpecifications &sps,
    /*! C id, if used */ int C_patch_data_id,
    /*! D id, if used */ int D_patch_data_id  ) const;

  /*!
    @brief Set parameters living on grid.

    Ignored data are: diffcoef_data and ccoef_data
    because they are constant.
  */
  void setGridData( hier::Patch<NDIM> &patch ,
		    pdat::SideData<NDIM,double> &diffcoef_data ,
		    pdat::CellData<NDIM,double> &ccoef_data ,
		    pdat::CellData<NDIM,double> &exact_data ,
		    pdat::CellData<NDIM,double> &source_data );

  virtual void setBcCoefs (
      tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
      tbox::Pointer<pdat::ArrayData<NDIM,double> > &bcoef_data ,
      tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
      const tbox::Pointer< hier::Variable<NDIM> > &variable ,
      const hier::Patch<NDIM> &patch ,
      const hier::BoundaryBox<NDIM> &bdry_box,
      const double fill_time=0.0) const;

  hier::IntVector<NDIM> numberOfExtensionsFillable() const;


  friend ostream &operator<<( ostream &os,
			      const PoissonPolynomialSolution &r );

private:

  QuarticFcn d_exact;
  QuarticFcn d_source;

};


#endif	// included_PoissonPolynomialSolution
