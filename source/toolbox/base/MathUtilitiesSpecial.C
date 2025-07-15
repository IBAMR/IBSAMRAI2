//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/toolbox/templates/special/MathUtilitiesSpecial.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	MathUtilities routines to set up handlers and get signaling NaNs
//

#include "tbox/MathUtilities.h"

#include "tbox/Complex.h"

#include <cmath>
#include <limits>

namespace SAMRAI {
   namespace tbox {

#define SAMRAI_FLT_SNAN std::numeric_limits<float>::signaling_NaN()
#define SAMRAI_DBL_SNAN std::numeric_limits<double>::signaling_NaN()

template<> bool   MathUtilities<bool>::s_zero           = false;
template<> bool   MathUtilities<bool>::s_one            = true;
template<> bool   MathUtilities<bool>::s_signaling_nan  = false;
template<> bool   MathUtilities<bool>::s_max            = true;
template<> bool   MathUtilities<bool>::s_min            = false;
template<> bool   MathUtilities<bool>::s_epsilon        = true;

template<> char   MathUtilities<char>::s_zero           = 0;
template<> char   MathUtilities<char>::s_one            = 1;
template<> char   MathUtilities<char>::s_signaling_nan  = std::numeric_limits<char>::max();
template<> char   MathUtilities<char>::s_max            = std::numeric_limits<char>::max();
template<> char   MathUtilities<char>::s_min            = std::numeric_limits<char>::min();
template<> char   MathUtilities<char>::s_epsilon        = 1;

template<> int    MathUtilities<int>::s_zero           = 0;
template<> int    MathUtilities<int>::s_one            = 1;
template<> int    MathUtilities<int>::s_signaling_nan   = std::numeric_limits<int>::max();
template<> int    MathUtilities<int>::s_max             = std::numeric_limits<int>::max();
template<> int    MathUtilities<int>::s_min             = std::numeric_limits<int>::min();
template<> int    MathUtilities<int>::s_epsilon         = 1;

template<> float  MathUtilities<float>::s_zero          = 0.0;
template<> float  MathUtilities<float>::s_one           = 1.0;
template<> float  MathUtilities<float>::s_signaling_nan = SAMRAI_FLT_SNAN;
template<> float  MathUtilities<float>::s_max           = std::numeric_limits<float>::max();
template<> float  MathUtilities<float>::s_min           = std::numeric_limits<float>::min();
template<> float  MathUtilities<float>::s_epsilon       = std::numeric_limits<float>::epsilon();

template<> double MathUtilities<double>::s_zero          = 0.0;
template<> double MathUtilities<double>::s_one           = 1.0;
template<> double MathUtilities<double>::s_signaling_nan = SAMRAI_DBL_SNAN;
template<> double MathUtilities<double>::s_max           = std::numeric_limits<double>::max();
template<> double MathUtilities<double>::s_min           = std::numeric_limits<double>::min();
template<> double MathUtilities<double>::s_epsilon       = std::numeric_limits<double>::epsilon();

template<> dcomplex   MathUtilities<dcomplex>::s_zero             = dcomplex(0.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_one              = dcomplex(1.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_signaling_nan  = dcomplex(SAMRAI_DBL_SNAN,SAMRAI_DBL_SNAN);
template<> dcomplex   MathUtilities<dcomplex>::s_max            = dcomplex(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
template<> dcomplex   MathUtilities<dcomplex>::s_min            = dcomplex(std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
template<> dcomplex   MathUtilities<dcomplex>::s_epsilon        = dcomplex(std::numeric_limits<double>::min(), 0.0);

template<>
bool MathUtilities<float>::isNaN(const float& value)
{
    return std::isnan(value);
}

template<> 
bool MathUtilities<double>::isNaN(const double& value)
{
    return std::isnan(value);

}

template<>
bool MathUtilities<dcomplex>::isNaN(const dcomplex& value)
{
  int i_re;
  int i_im;
  i_re = std::isnan( std::real(value) );
  i_im = std::isnan( std::imag(value) );

  return( ( (i_re != 0) || (i_im !=0) ) ? true : false );
}

template<>
bool MathUtilities<float>::equalEps(const float& a, const float& b)
{
   float absmax = MathUtilities<float>::Max(
                     MathUtilities<float>::Abs(a),
                     MathUtilities<float>::Abs(b) );
   float numerator = MathUtilities<float>::Abs(a-b);
   float denomenator =
      MathUtilities<float>::Max(absmax,
           MathUtilities<float>::s_epsilon);

   return( numerator/denomenator < std::sqrt(MathUtilities<float>::s_epsilon) );
}

template<>
bool MathUtilities<double>::equalEps(const double& a, const double& b)
{
   double absmax = MathUtilities<double>::Max(
                      MathUtilities<double>::Abs(a),
                      MathUtilities<double>::Abs(b) );
   double numerator = MathUtilities<double>::Abs(a-b);
   double denomenator =
      MathUtilities<double>::Max(absmax,
           MathUtilities<double>::s_epsilon);

   return( numerator/denomenator < std::sqrt(MathUtilities<double>::s_epsilon) );
}

template<>
bool MathUtilities<dcomplex>::equalEps(const dcomplex& a, const dcomplex& b)
{
   double a_re = std::real(a);
   double a_im = std::imag(a);
   double b_re = std::real(b);
   double b_im = std::imag(b);

   return( MathUtilities<double>::equalEps(a_re,b_re) && 
           MathUtilities<double>::equalEps(a_im,b_im) );
}

template<> 
dcomplex MathUtilities<dcomplex>::Min(dcomplex a, dcomplex b)
{
   return(norm(a) < norm(b) ? a : b);
}

template<> 
dcomplex MathUtilities<dcomplex>::Max(dcomplex a, dcomplex b)
{
   return(norm(a) > norm(b) ? a : b);
}

template<>
int MathUtilities<int>::Abs(int a)
{
   return(a > 0 ? a : -a);
}

template<>
float MathUtilities<float>::Abs(float a)
{
   return(a > 0.0 ? a : -a);
}

template<>
double MathUtilities<double>::Abs(double a)
{
   return(a > 0.0 ? a : -a);
}

template<> 
bool MathUtilities<bool>::Rand(const bool& low, const bool& width)
{
   (void) low;
   (void) width;
   return( mrand48() > 0 ? true : false );
}

template<> 
char MathUtilities<char>::Rand(const char& low, const char& width)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   return static_cast<char>( static_cast<double>(width) * drand48() ) + low;
}

template<> 
int MathUtilities<int>::Rand(const int& low, const int& width)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   return static_cast<int>( static_cast<double>(width) * drand48() ) + low;
}

template<> 
float MathUtilities<float>::Rand(const float& low, const float& width)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   return width * drand48() + low;
}

template<> 
double MathUtilities<double>::Rand(const double& low, const double& width)
{
   return width * drand48() + low;
}

template<> 
dcomplex MathUtilities<dcomplex>::Rand(const dcomplex& low, const dcomplex& width)
{
   double real_part = std::real(width) * drand48() + std::real(low);
   double imag_part = std::imag(width) * drand48() + std::imag(low);
   return dcomplex(real_part, imag_part);
}

template<class TYPE>  
TYPE round_internal(TYPE x)
{
   /* algorithm used from Steven G. Kargl */
   double t;
   if (x >= 0.0) {
      t = ceil(x);
      if (t - x > 0.5)
        t -= 1.0;
      return (t);
    } else {
      t = ceil(-x);
      if (t + x > 0.5)
        t -= 1.0;
      return (-t);
    }
}

template float round_internal< float >(float x);
template double round_internal< double >(double x);

template<>  
float MathUtilities<float>::round(float x) {
   return round_internal<float>(x);
}

template<>  
double MathUtilities<double>::round(double x) {
   return round_internal<double>(x);
}

}
}

