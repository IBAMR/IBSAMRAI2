//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/vectors/PETSc_SAMRAIVectorReal.h $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 3281 $
// Modified:    $LastChangedDate: 2009-06-17 19:10:03 -0700 (Wed, 17 Jun 2009) $
// Description: "Glue code" between PETSc vector interface and SAMRAI vectors. 
//

#ifndef included_solv_PETSc_SAMRAIVectorReal
#define included_solv_PETSc_SAMRAIVectorReal

#include "SAMRAI_config.h"

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
************************************************************************
*/
#ifdef HAVE_PETSC

#include "tbox/SAMRAI_MPI.h"


/* petsc defines skip so avoid redefining */
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "petscvec.h"

#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#include "PETScAbstractVectorReal.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#ifdef REQUIRES_CMATH
#include <cmath>
#endif 


namespace SAMRAI {
    namespace solv {

/**
 * Class PETSc_SAMRAIVectorReal<DIM> wraps a real-valued SAMRAI vector 
 * (see SAMRAIVectorReal class) object so that it may be used with 
 * the PETSc solver package.  This class is derived from the abstract base 
 * class PETScAbstractVectorReal, the <TT>C++</TT> interface for 
 * PETSc vectors where the underlying data is float or double.  It also 
 * maintains a pointer to a SAMRAI vector object.  A SAMRAI vector is defined 
 * as a collection of patch data components and associated operations living 
 * on some subset of levels in a structured AMR mesh hierarchy.  
 *
 * Observe that there are only three public member functions in this class
 * They are used to create and destroy PETSc vector objects (i.e., "Vec"s),
 * and to obtain the SAMRAI vector associated with the PETSc vector.
 * In particular, note that the constructor and destructor of this class
 * are protected members.  The construction and destruction of instances of
 * this class may occur only through the static member functions that 
 * create and destroy PETSc vector objects.
 *
 * Finally, we remark that PETSc allows vectors with complex-valued entries. 
 * This class and the class VectorReal assume real-values vectors 
 * (i.e., data of type <TT>double</TT> or <TT>float</TT>.  The class 
 * PETSc_SAMRAIVectorComplex must be used for complex data.
 *
 * @see solv::PETScAbstractVectorReal
 * @see solv::SAMRAIVectorReal 
 */

template<int DIM, class TYPE>
class PETSc_SAMRAIVectorReal : public PETScAbstractVectorReal<TYPE>
{
public:
   /**
    * Create and return a new PETSc vector object.  The SAMRAI vector
    * object is wrapped so that it may be manipulated through the PETSc 
    * vector structure in the PETScAbstractVectorReal base class.  It is 
    * important to note that this function does not allocate storage for 
    * the vector data.  Data must be allocated through the SAMRAI vector 
    * object directly.  For output of the data through PETSc "ViewVec" calls,
    * the output stream to which the SAMRAI vector object writes will be used.
    */
   static Vec createPETScVector(
      tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > samrai_vec,
      MPI_Comm comm=PETSC_COMM_SELF);

   /**
    * Destroy a given PETSc vector object. It is important to note that
    * this function does not deallocate storage for the vector data.
    * Vector data must be deallocated through the SAMRAI vector object.
    */
   static void destroyPETScVector(Vec petsc_vec);

   /**
    * Return pointer to the SAMRAI vector object associated with the  
    * given PETSc vector.
    */
   static tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > 
      getSAMRAIVector(Vec petsc_vec);

protected:
   /*
    * Constructor for PETSc_SAMRAIVectorReal<DIM> is protected so that an
    * object of this class cannot be constructed without calling the static
    * member function createPETScVector, which is used to construct a PETSc
    * vector and associated "wrapper" that allows the PETSc vector to 
    * manipulate the SAMRAIVectorReal data.  
    *
    * The boolean argument is used to control whether the SAMRAI vector is 
    * destroyed when the associated PETSc vector is destroyed.  This should 
    * happen if the PETSc vector is created within PETSc via a duplicate 
    * (i.e., clone) operation, but not otherwise.
    */
   PETSc_SAMRAIVectorReal(
      tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > samrai_vector,
      bool vector_created_via_duplicate,
      MPI_Comm comm);

   /*
    * Virtual destructor for PETSc_SAMRAIVectorReal<DIM> is protected so that 
    * an object of this class cannot be destroyed without calling the static
    * member function destroyPETScVector, which is used to destroy the PETSc
    * vector and associated "wrapper".
    */
   virtual ~PETSc_SAMRAIVectorReal<DIM,TYPE>();

private:
   /*
    * Return SAMRAI vector owned by this SAMRAI_PETScVector object.
    */
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > getSAMRAIVector();

   /*
    * Clone this vector structure and allocate storage for the vector
    * data.  Then, return a pointer to the new vector instance.  This
    * function is distinct from the vector constructor since it is called
    * from within PETSc to allocate new vectors.
    */
   PETScAbstractVectorReal<TYPE>* makeNewVector();


   /*
    * Destroy vector structure and its storage. This function is
    * distinct from the destructor since it will be called from
    * within PETSc to deallocate vectors.
    */
   virtual void freeVector();

   /*
    * Print all vector data as implemented in SAMRAI vector class.
    */
   virtual void viewVector() const;
   
   /*
    * Return \f$(x,y) = \sum_i (x_i * conj(y_i)), where \f$x\f$ is
    * this vector.  Note that we assume real vectors; so this is the
    * same as TdotWith().  If local_only is true, the operation is
    * performed only on the local parts.
    */
   virtual double dotWith(
      const PETScAbstractVectorReal<TYPE>* y,
        bool local_only=false) const;

    /*
     * Return \f$(x,y) = \sum_i (x_i * (y_i)), where \f$x\f$ is this
     * vector.  Note that we assume real vectors; so this is the same
     * as dotWith().  If local_only is true, the operation is
     * performed only on the local parts.
     */
    virtual double TdotWith(
        const PETScAbstractVectorReal<TYPE>* y,
        bool local_only=false) const;

    /*
     * Return \f$L_1\f$-norm of this vector.
     *
     * \param local_only Flag to get result for local data only.
     */
    virtual double L1Norm(
        bool local_only=false) const;

    /*
     * Return \f$L_2\f$-norm of this vector.
     *
     * \param local_only Flag to get result for local data only.
     */
    virtual double L2Norm(
        bool local_only=false) const;

    /*
     * Return \f$L_{\infty}\f$-norm of this vector.
     *
     * \param local_only Flag to get result for local data only.
     */
    virtual double maxNorm(
        bool local_only=false) const;

    /*
     * Multiply each entry of this vector by given scalar.
     */
    virtual void scaleVector(
        const TYPE alpha);

    /*
     * Copy source vector data to this vector.
     */
    virtual void copyVector(
        const PETScAbstractVectorReal<TYPE>* v_src);

    /*
     * Set each entry of this vector to given scalar.
     */
    virtual void setToScalar(
        const TYPE alpha);

    /*
     * Swap data between this vector and argument vector.
     */
    virtual void swapWith(
        PETScAbstractVectorReal<TYPE>* v_other);

    /*
     * Set \f$y = \alpha x + y\f$, where \f$y\f$ is this vector.
     */
    virtual void setAXPY(
        const TYPE alpha,
        const PETScAbstractVectorReal<TYPE>* x);

    /*
     * Set \f$y = \alpha x + \beta y\f$, where \f$y\f$ is this vector.
     */
    virtual void setAXPBY(
        const TYPE alpha,
        const PETScAbstractVectorReal<TYPE>* x,
        const TYPE beta);

    /*
     * Set \f$w = \alpha x + y\f$, where \f$w\f$ is this vector.
     */
    virtual void setWAXPY(
        const TYPE alpha,
        const PETScAbstractVectorReal<TYPE>* x,
        const PETScAbstractVectorReal<TYPE>* y);

    /*
     * Set \f$w_i = x_i y_i\f$, where \f$w_i\f$ is \f$i\f$-th entry of
     * this vector.
     */
    virtual void pointwiseMultiply(
        const PETScAbstractVectorReal<TYPE>* x,
        const PETScAbstractVectorReal<TYPE>* y);

    /*
     * Set \f$w_i = x_i / y_i\f$, where \f$w_i\f$ is \f$i\f$-th entry
     * of this vector.
     */
    virtual void pointwiseDivide(
        const PETScAbstractVectorReal<TYPE>* x,
        const PETScAbstractVectorReal<TYPE>* y);

    /*
     * Compute \f$ max_i = abs(w_i / y_i) \f$ where \f$ y_i \neq 0
     * \f$, and \f$ max_i = abs(w_i) \f$ where \f$ y_i = 0 \f$, where
     * \f$ w_i \f$ is \f$ i \f$ -th entry of this vector.
     */
    virtual double maxPointwiseDivide(
        const PETScAbstractVectorReal<TYPE>* y);

    /*
     * Find maximum vector entry and vector index at which maximum
     * occurs.  Note that this function sets the index to a bogus
     * value since the it is not clear how to define the vector index
     * for the the SAMRAI vector class.
     */
    virtual void vecMax(
        int& i,
        TYPE& max) const;

    /*
     * Find minimum vector entry and vector index at which minimum
     * occurs.  Note that this function sets the index to a bogus
     * value since the it is not clear how to define the vector index
     * for the the SAMRAI vector class.
     */
    virtual void vecMin(
        int& i,
        TYPE& min) const;


    /*
     * Set each entry of this vector to random value.  To be
     * consistent with PETSc, we define the random value to be
     * width*drand48()+low.
     */
    virtual void setRandomValues(
        const TYPE width,
        const TYPE low);

    /*
     * Set argument to pointer to vector data in contiguous array
     * (local to processor).
     *
     * \note This function returns a null pointer.  Hopefully, this
     * will not cause problems.
     */
    virtual void getDataArray(
        TYPE** array);

    /*
     * Restore pointer to vector data in contiguous array (local to
     * processor).
     */
    virtual void restoreDataArray(
        TYPE** array);

    /*
     * Return total length of vector.
     *
     * \note This function returns zero since PETSc requires a valid
     * integer return value.  However, the manner in which we are
     * using PETSc currently does not require the value to be used
     * within PETSc.  Hopefully, this will not cause problems.
     */
    virtual int getDataSize() const;

    /*
     * Return length of vector (local to processor).
     *
     * \note This function returns zero since PETSc requires a valid
     * integer return value.  However, the manner in which we are
     * using PETSc currently does not require the value to be used
     * within PETSc.  Hopefully, this will not cause problems.
     */
    virtual int getLocalDataSize() const;

   /*
    * Vector data is maintained in SAMRAI vector structure.
    */
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > d_samrai_vector; 

   /*
    * Boolean flag to control whether SAMRAI vector is destroyed when
    * the associated PETSc vector is destroyed. 
    */
   bool d_vector_created_via_duplicate;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "PETSc_SAMRAIVectorReal.I"
#endif
#endif

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PETSc_SAMRAIVectorReal.C"
#endif
