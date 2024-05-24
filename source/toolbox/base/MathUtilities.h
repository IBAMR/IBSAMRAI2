//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/base/MathUtilities.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Utilities class to access common POSIX constants and math ops
//

#ifndef included_tbox_MathUtilities
#define included_tbox_MathUtilities

#include "SAMRAI_config.h"

#include "tbox/Array.h"

#include "tbox/Complex.h"

namespace SAMRAI {
   namespace tbox {

/*!
 * Class MathUtilities is a utility that provides some basic math-related 
 * functions and routines for initializing data to signaling NaNs and
 * POSIX constants like INT_MAX, FLT_MAX, DBL_MAX, etc.  Signaling
 * NaNs force a trap if they are used in a numerical operation, so
 * they are a useful way to track uninitialized floating point data.
 * For example, setting integer values to INT_MAX is a useful way to track
 * uninitialized integer values.
 *
 * IMPORTANT: To properly trap operations based on signaling NaN values,
 *            the routine IEEE::setupFloatingPointExceptionHandlers()
 *            must be called.  This is normally done in the
 *            SAMRAIManager::startup() routine.
 *
 * The implementation of this class depends heavily on the particular
 * computer architecture and how it implements floating point arithmetic
 * and hardware traps.  
 * 
 * Note that the class @see tbox::IEEE also provides operations for 
 * dealing with signaling NaNs.  This class provides the actual 
 * implementation for such operations with functions in tbox::IEEE
 * calling the approriate operations provided here.  The class
 * tbox::IEEE is not templated on the data type and so calling the
 * operations provided there may be easier in some cases, such as in
 * codes built based on earlier versions of SAMRAI.
 */

template<class TYPE> class MathUtilities
{
public:
   /*!
    * @brief Return the value 0 for the template type.
    */
   static TYPE getZero();

   /*!
    * @brief Return the value 1 for the template type.
    */
   static TYPE getOne();

   /*!
    * @brief Return the IEEE signaling NaN for the template type on 
    * architectures that support it.  
    *
    * Using this value in a numerical expression will
    * cause a program abort.
    *
    * For float, double, and dcomplex types, the usual signaling
    * Nan value will be returned if defined for the architecture
    * and compiler.  For char and int types, the POSIX max value for
    * the type is returned.  For other template types, zero is returned.
    */
   static TYPE getSignalingNaN();

   /*!
    * @brief Return true if the supplied value is NaN; else, false.
    *
    * For float and double, will check value against signaling NaN.
    * For dcomplex will return true if either real and imaginary part 
    * is signaling NaN and false otherwise.  
    * For other types will return false always.
    *
    * @param value Value to test
    */
   static bool isNaN(const TYPE& value);

   /*!
    * @brief Set array entries to value given by getSignalingNaN().
    *
    * @param array SAMRAI array to set
    */
   static void setArrayToSignalingNaN(Array<TYPE>& array);

   /*!
    * @brief Set array entries to value given by getSignalingNaN().
    *
    * @param array Pointer to first array entry 
    * @param n     Integer to number of array entries.
    */
   static void setArrayToSignalingNaN(TYPE* array, int n = 1);

   /*!
    * @brief Return true if given values have a relative difference
    * smaller than sqrt(mach_eps) for the template type.
    *
    * Valid for float/double/dcomplex only.  For dcomplex will 
    * return true if both real and imaginary part values 
    * have a relative difference smaller than sqrt(mach_eps). 
    *
    * For non float/double/dcomplex will just return ( a == b ).  
    * 
    * @param a
    * @param b
    */
   static bool equalEps(const TYPE& a, const TYPE& b);

   /*!
    * @brief Return max value for the template type.
    * 
    * For boolean type, will return "true".  For dcomplex type,
    * will return a dcomplex value with both real and imaginary
    * parts set to the POSIX max value for type double.  For 
    * other types, will return the POSIX max value for the type.
    */
   static TYPE getMax();

   /*!
    * @brief Set array entries to value given by getMax().
    *
    * @param array SAMRAI array to set
    */
   static void setArrayToMax(Array<TYPE>& array);
 
   /*!
    * @brief Set array entries to value given by getMax().
    *
    * @param array Pointer to first array entry
    * @param n     Integer to number of array entries.
    */
   static void setArrayToMax(TYPE* array, int n = 1);

   /*!
    * @brief Return min value for the template type.
    * 
    * For boolean type, will return "false".  For dcomplex type,
    * will return a dcomplex value with both real and imaginary
    * parts set to the POSIX min value for type double.  For 
    * other types, will return the POSIX min value for the type.
    */
   static TYPE getMin();

   /*!
    * @brief Set array entries to value given by getMin().
    *
    * @param array SAMRAI array to set
    */
   static void setArrayToMin(Array<TYPE>& array);
 
   /*!
    * @brief Set array entries to value given by getMin().
    *
    * @param array Pointer to first array entry
    * @param n     Integer to number of array entries.
    */
   static void setArrayToMin(TYPE* array, int n = 1);

   /*!
    * @brief Return epsilon value for the template type.
    * 
    * For boolean type, will return "true".  For integer type, will
    * return "1".  For dcomplex type, will return a dcomplex value with 
    * both real and imaginary parts set to the POSIX epsilon value for 
    * type double.  For other types, will return the POSIX epsilon value 
    * for the type.
    */
   static TYPE getEpsilon();

   /*!
    * @brief Set array entries to value given by getEpsilon().
    *
    * @param array SAMRAI array to set
    */
   static void setArrayToEpsilon(Array<TYPE>& array);
 
   /*!
    * @brief Set array entries to value given by getEpsilon().
    *
    * @param array Pointer to first array entry
    * @param n     Integer to number of array entries.
    */
   static void setArrayToEpsilon(TYPE* array, int n = 1);

   /*!
    * @brief Return the minimum value of a and b.
    *
    * For dcomplex type will return the value with the minimum norm.
    * 
    * @param a
    * @param b
    */
   static TYPE Min(TYPE a, TYPE b);

   /*!
    * @brief Return the maximum value of a and b.
    *
    * For dcomplex type will return the value with the maximum norm.
    * 
    * @param a
    * @param b
    */
   static TYPE Max(TYPE a, TYPE b);

   /*!
    * @brief Return absolute value of a.
    *
    * For int, float, and double types will return the absolute
    * numerical value.  For other types, will return the input value.
    *
    * @param a
    */
   static TYPE Abs(TYPE a);


   /*!
    * @brief Return nearest integral value.
    *
    * If the argument is halfway between two integral
    * values then round away from zero for float and double types.
    * For other types, will return the input value.
    *
    * @param a
    */
   static TYPE round(TYPE a);

   /*!
    * @brief Generate and return a random value from low to low+width.
    *
    * For boolean type, value is either "true" or "false" based on
    * value generated by mrand48().  For other types, returned value is 
    * computed as width * drand48() + low.  When type is char, this value
    * is cast to a char. When type is dcomplex, real and imaginary
    * parts are computed separately this way. 
    * 
    * @param low   Starting value for range
    * @param width Width of the range.
    */
   static TYPE Rand(const TYPE& low, const TYPE& width);

private:
   static TYPE  s_zero;
   static TYPE  s_one;
   static TYPE  s_signaling_nan;
   static TYPE  s_max;
   static TYPE  s_min;
   static TYPE  s_epsilon;
};

/*
 * Template specializations.
 */
template<> bool MathUtilities<float>::isNaN(const float& value);
template<> bool MathUtilities<double>::isNaN(const double& value);
template<> bool MathUtilities<dcomplex>::isNaN(const dcomplex& value);

template<> bool MathUtilities<float>::equalEps(const float& a, const float& b);
template<> bool MathUtilities<double>::equalEps(const double& a, const double& b);
template<> bool MathUtilities<dcomplex>::equalEps(const dcomplex& a, const dcomplex& b);

template<> bool     MathUtilities<bool>::Rand(const bool& low, const bool& width);
template<> char     MathUtilities<char>::Rand(const char& low, const char& width);
template<> int      MathUtilities<int>::Rand(const int& low, const int& width);
template<> float    MathUtilities<float>::Rand(const float& low, const float& width);
template<> double   MathUtilities<double>::Rand(const double& low, const double& width);
template<> dcomplex MathUtilities<dcomplex>::Rand(const dcomplex& low, const dcomplex& width);

template<> dcomplex MathUtilities<dcomplex>::Min(dcomplex a, dcomplex b);
template<> dcomplex MathUtilities<dcomplex>::Max(dcomplex a, dcomplex b);

template<> int      MathUtilities<int>::Abs(int a);
template<> float    MathUtilities<float>::Abs(float a);
template<> double   MathUtilities<double>::Abs(double a);


template<> float    MathUtilities<float>::round(float a);
template<> double   MathUtilities<double>::round(double a);

template<> bool   MathUtilities<bool>::s_zero;
template<> bool   MathUtilities<bool>::s_one;
template<> bool   MathUtilities<bool>::s_signaling_nan;
template<> bool   MathUtilities<bool>::s_max;
template<> bool   MathUtilities<bool>::s_min;
template<> bool   MathUtilities<bool>::s_epsilon;

template<> char   MathUtilities<char>::s_zero;
template<> char   MathUtilities<char>::s_one;
template<> char   MathUtilities<char>::s_signaling_nan;
template<> char   MathUtilities<char>::s_max;
template<> char   MathUtilities<char>::s_min;
template<> char   MathUtilities<char>::s_epsilon;

template<> int    MathUtilities<int>::s_zero;
template<> int    MathUtilities<int>::s_one;
template<> int    MathUtilities<int>::s_signaling_nan;
template<> int    MathUtilities<int>::s_max;
template<> int    MathUtilities<int>::s_min;
template<> int    MathUtilities<int>::s_epsilon;

template<> float  MathUtilities<float>::s_zero;
template<> float  MathUtilities<float>::s_one;
template<> float  MathUtilities<float>::s_signaling_nan;
template<> float  MathUtilities<float>::s_max;
template<> float  MathUtilities<float>::s_min;
template<> float  MathUtilities<float>::s_epsilon;

template<> double MathUtilities<double>::s_zero;
template<> double MathUtilities<double>::s_one;
template<> double MathUtilities<double>::s_signaling_nan;
template<> double MathUtilities<double>::s_max;
template<> double MathUtilities<double>::s_min;
template<> double MathUtilities<double>::s_epsilon;

template<> dcomplex   MathUtilities<dcomplex>::s_zero;
template<> dcomplex   MathUtilities<dcomplex>::s_one;
template<> dcomplex   MathUtilities<dcomplex>::s_signaling_nan;
template<> dcomplex   MathUtilities<dcomplex>::s_max;
template<> dcomplex   MathUtilities<dcomplex>::s_min;
template<> dcomplex   MathUtilities<dcomplex>::s_epsilon;

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/MathUtilities.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "tbox/MathUtilities.C"
#endif

#endif
