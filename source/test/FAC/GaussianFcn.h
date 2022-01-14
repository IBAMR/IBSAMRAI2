/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/test/FAC/GaussianFcn.h $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Gaussian function support for FAC solver tests.
 */

#ifndef included_GaussianFcn_h
#define included_GaussianFcn_h

#include <iostream>

using namespace std;

/*!
  @brief Gaussian function functor.

  Computes the function
  @f[ e^{ \lambda |r-r_0|^2 } @f]

  lambda is generally, but not necessarily, negative.

  This functor exists for convenience, not efficiency.
  If you use it heavily, you can improve your efficiency
  by replacing its usage with lower level codes.
*/
class GaussianFcn {

public:

  GaussianFcn();


  /*!
    @brief Return the function value.
  */
#if NDIM == 1
  double operator()( double x ) const;
#endif
#if NDIM == 2
  double operator()( double x, double y ) const;
#endif
#if NDIM == 3
  double operator()( double x, double y, double z ) const;
#endif


  /*!
    @brief Set amplitude.
  */
  int setAmplitude( const double amp );

  /*!
    @brief Set all wave numbers.

    Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int setLambda( const double lambda );

  /*!
    @brief Set all phase angles.

    Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
  */
  int setCenter( const double *center );


  /*!
    @brief Get amplitude.
  */
  double getAmplitude() const;

  /*!
    @brief Get lambda.
  */
  double getLambda() const;

  /*!
    @brief Get center coordinates.
  */
  int getCenter( double *center ) const;


  //@{
  /*!
    @name IO operators.
  */
  /*!
    @brief Input from stream.

    Stream extract function reads input in the format used
    by the stream insert function (see
    opertor<<(ostream&, GaussianFcn &).
    Except for allowing for missing centers (set to zero)
    and lambda (set to 1), this function requires the correct
    syntax or the result is undefined.

    @see GaussianFcn::ctype
  */
  friend istream &operator>>( istream &ci, GaussianFcn &cf );
  /*!
    @brief Output to stream.

    Outputs sets of name=double values where the name is one
    of nx, px, ny, nz or pz, and the double is the value of
    the coefficient.

    Format of output with example values is
    @verbatim
    { lambda=1.0 cx=0 cy=0.0 cz=0.0 }
    @endverbatim
    where cx, cy and cz are the center of the Gaussian function.

    @see GaussianFcn::ctype
  */
  friend ostream &operator<<( ostream &co, const GaussianFcn &cf );
  //@}


private:

  //! Amplitude
  double d_amp;
  //! Center
  double d_center[NDIM];
  //! Lambda
  double d_lambda;

};




#endif
