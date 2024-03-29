//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/solvers/packages/petsc/PETScAbstractVectorReal.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2195 $
// Modified:    $LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description: Interface to C++ vector implementation for PETSc package.
//

#ifndef included_solv_PETScAbstractVectorReal_C
#define included_solv_PETScAbstractVectorReal_C

#include "PETScAbstractVectorReal.h"

#ifdef HAVE_PETSC

#include <string>

#include "tbox/IOStream.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
namespace solv {


#define PABSVEC_CAST(v) (static_cast<PETScAbstractVectorReal<TYPE>*>(v->data))


/*
*************************************************************************
*                                                                       *
* PETScAbstractVectorReal constructor and destructor.              *
*                                                                       *
*************************************************************************
*/

template <class TYPE>
PETScAbstractVectorReal<TYPE>::PETScAbstractVectorReal(
   bool vector_created_via_duplicate,
   MPI_Comm comm)
   : d_petsc_vector(static_cast<Vec>(NULL)),
     d_vector_created_via_duplicate(vector_created_via_duplicate),
     d_comm(comm)
{
   
   int ierr = 0;

   ierr = VecCreate(d_comm, &d_petsc_vector); PETSC_SAMRAI_ERROR(ierr);
   
   // Set PETSc vector data to this abstract vector object
   d_petsc_vector->precision                    = PETSC_SCALAR;
   d_petsc_vector->data                         = this;
   d_petsc_vector->petscnative                  = PETSC_FALSE;
   d_petsc_vector->map.n                        = 0;
   d_petsc_vector->map.N                        = 0;
   d_petsc_vector->map.bs                       = 1;
   
   // Assign vector operations to PETSc vector object.
   d_petsc_vector->ops->duplicate               = PETScAbstractVectorReal<TYPE>::vecDuplicate;
   d_petsc_vector->ops->duplicatevecs           = PETScAbstractVectorReal<TYPE>::vecDuplicateVecs;
   d_petsc_vector->ops->destroyvecs             = PETScAbstractVectorReal<TYPE>::vecDestroyVecs;
   d_petsc_vector->ops->dot                     = PETScAbstractVectorReal<TYPE>::vecDot;
   d_petsc_vector->ops->mdot                    = PETScAbstractVectorReal<TYPE>::vecMDot;
   d_petsc_vector->ops->norm                    = PETScAbstractVectorReal<TYPE>::vecNorm;
   d_petsc_vector->ops->tdot                    = PETScAbstractVectorReal<TYPE>::vecTDot;
   d_petsc_vector->ops->mtdot                   = PETScAbstractVectorReal<TYPE>::vecMTDot;
   d_petsc_vector->ops->scale                   = PETScAbstractVectorReal<TYPE>::vecScale;
   d_petsc_vector->ops->copy                    = PETScAbstractVectorReal<TYPE>::vecCopy;
   d_petsc_vector->ops->set                     = PETScAbstractVectorReal<TYPE>::vecSet;
   d_petsc_vector->ops->swap                    = PETScAbstractVectorReal<TYPE>::vecSwap;
   d_petsc_vector->ops->axpy                    = PETScAbstractVectorReal<TYPE>::vecAXPY;
   d_petsc_vector->ops->axpby                   = PETScAbstractVectorReal<TYPE>::vecAXPBY;
   d_petsc_vector->ops->maxpy                   = PETScAbstractVectorReal<TYPE>::vecMAXPY;
   d_petsc_vector->ops->aypx                    = PETScAbstractVectorReal<TYPE>::vecAYPX;
   d_petsc_vector->ops->waxpy                   = PETScAbstractVectorReal<TYPE>::vecWAXPY;
   d_petsc_vector->ops->pointwisemult           = PETScAbstractVectorReal<TYPE>::vecPointwiseMult;
   d_petsc_vector->ops->pointwisedivide         = PETScAbstractVectorReal<TYPE>::vecPointwiseDivide;
   d_petsc_vector->ops->getarray                = PETScAbstractVectorReal<TYPE>::vecGetArray;
   d_petsc_vector->ops->getsize                 = PETScAbstractVectorReal<TYPE>::vecGetSize;
   d_petsc_vector->ops->getlocalsize            = PETScAbstractVectorReal<TYPE>::vecGetLocalSize;
   d_petsc_vector->ops->restorearray            = PETScAbstractVectorReal<TYPE>::vecRestoreArray;
   d_petsc_vector->ops->max                     = PETScAbstractVectorReal<TYPE>::vecMax;
   d_petsc_vector->ops->min                     = PETScAbstractVectorReal<TYPE>::vecMin;
   d_petsc_vector->ops->setrandom               = PETScAbstractVectorReal<TYPE>::vecSetRandom;
   d_petsc_vector->ops->destroy                 = PETScAbstractVectorReal<TYPE>::vecDestroy;
   d_petsc_vector->ops->view                    = PETScAbstractVectorReal<TYPE>::vecView;
   d_petsc_vector->ops->dot_local               = PETScAbstractVectorReal<TYPE>::vecDot_local;
   d_petsc_vector->ops->tdot_local              = PETScAbstractVectorReal<TYPE>::vecTDot_local;
   d_petsc_vector->ops->norm_local              = PETScAbstractVectorReal<TYPE>::vecNorm_local;
   d_petsc_vector->ops->mdot_local              = PETScAbstractVectorReal<TYPE>::vecMDot_local;
   d_petsc_vector->ops->mtdot_local             = PETScAbstractVectorReal<TYPE>::vecMTDot_local;
   d_petsc_vector->ops->loadintovector          = VecLoadIntoVector_Default;
   d_petsc_vector->ops->maxpointwisedivide      = PETScAbstractVectorReal<TYPE>::vecMaxPointwiseDivide;

   // The remaining functions will result in program abort.
   d_petsc_vector->ops->setvalues               = PETScAbstractVectorReal<TYPE>::vecSetValues;
   d_petsc_vector->ops->assemblybegin           = PETScAbstractVectorReal<TYPE>::vecAssemblyBegin;
   d_petsc_vector->ops->assemblyend             = PETScAbstractVectorReal<TYPE>::vecAssemblyEnd;
   d_petsc_vector->ops->setoption               = PETScAbstractVectorReal<TYPE>::vecSetOption;
   d_petsc_vector->ops->setvaluesblocked        = PETScAbstractVectorReal<TYPE>::vecSetValuesBlocked;
   d_petsc_vector->ops->placearray              = PETScAbstractVectorReal<TYPE>::vecPlaceArray;
   d_petsc_vector->ops->replacearray            = PETScAbstractVectorReal<TYPE>::vecReplaceArray;
   d_petsc_vector->ops->reciprocal              = PETScAbstractVectorReal<TYPE>::vecReciprocal;
   d_petsc_vector->ops->viewnative              = PETScAbstractVectorReal<TYPE>::vecViewNative;
   d_petsc_vector->ops->conjugate               = PETScAbstractVectorReal<TYPE>::vecConjugate;
   d_petsc_vector->ops->setlocaltoglobalmapping = PETScAbstractVectorReal<TYPE>::vecSetLocalToGlobalMapping;
   d_petsc_vector->ops->setvalueslocal          = PETScAbstractVectorReal<TYPE>::vecSetValuesLocal;
   d_petsc_vector->ops->resetarray              = PETScAbstractVectorReal<TYPE>::vecResetArray;
   d_petsc_vector->ops->setfromoptions          = PETScAbstractVectorReal<TYPE>::vecSetFromOptions;
   d_petsc_vector->ops->load                    = PETScAbstractVectorReal<TYPE>::vecLoad;
   d_petsc_vector->ops->pointwisemax            = PETScAbstractVectorReal<TYPE>::vecPointwiseMax;
   d_petsc_vector->ops->pointwisemaxabs         = PETScAbstractVectorReal<TYPE>::vecPointwiseMaxAbs;
   d_petsc_vector->ops->pointwisemin            = PETScAbstractVectorReal<TYPE>::vecPointwiseMin;
   d_petsc_vector->ops->getvalues               = PETScAbstractVectorReal<TYPE>::vecGetValues;

   ierr = PetscMapInitialize(d_comm, &d_petsc_vector->map); PETSC_SAMRAI_ERROR(ierr);
   ierr = PetscMapSetBlockSize(&d_petsc_vector->map, 1); PETSC_SAMRAI_ERROR(ierr);
   ierr = PetscMapSetSize(&d_petsc_vector->map, 0); PETSC_SAMRAI_ERROR(ierr);
   ierr = PetscMapSetLocalSize(&d_petsc_vector->map, 0); PETSC_SAMRAI_ERROR(ierr);

   const std::string my_name = "PETScAbstractVectorReal";

   if (d_petsc_vector->type_name) {
      ierr = PetscFree(d_petsc_vector->type_name);  PETSC_SAMRAI_ERROR(ierr);
   }

   ierr = PetscObjectChangeTypeName(reinterpret_cast<PetscObject>(d_petsc_vector), my_name.c_str()); PETSC_SAMRAI_ERROR(ierr);

   return;
}


template <class TYPE>
PETScAbstractVectorReal<TYPE>::~PETScAbstractVectorReal()
{
   int ierr = 0;
   
   if (!d_vector_created_via_duplicate)
   {
      d_petsc_vector->ops->destroy = 0;
      ierr = VecDestroy(d_petsc_vector); PETSC_SAMRAI_ERROR(ierr);
   }
   
   d_petsc_vector = 0;
   
   return;
}

template <class TYPE>
PetscErrorCode PETScAbstractVectorReal<TYPE>::vecDuplicateVecs(Vec v_in, 
							       int n, Vec** varr_new)
{
   int ierr = 0;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(v_in == (Vec)NULL));
#endif
   ierr = PetscMalloc( n * sizeof(Vec *), varr_new);
   PETSC_SAMRAI_ERROR(ierr);
   
   for (int i = 0; i < n; i++) {
      PETScAbstractVectorReal<TYPE>::vecDuplicate(v_in, *varr_new+i);
   }
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode PETScAbstractVectorReal<TYPE>::vecDestroyVecs(Vec *v_arr, PetscInt n)
{
   int i;
   int ierr = 0;

#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < n; i++) {
      TBOX_ASSERT(!(v_arr[i] == (Vec)NULL));
   }
#endif
   for (i = 0; i < n; i++) {
      Vec petsc_vec = v_arr[i];
      vecDestroy(petsc_vec);

      // There is some assymetry in the allocation/deallocation 
      // The PETSc Vec structures are created in the constructor
      // and normally deallocated by the PETSc VecDestroy call.
      // This holds for VecCreated and VecDuplicate Vec.
      // However in the case of DuplicateVecs the VecDestroy 
      // needs to be called on the PETSc Vec structure.
      petsc_vec->ops->destroy = 0;
      ierr = VecDestroy(petsc_vec); PETSC_SAMRAI_ERROR(ierr);
   }
	 
   // SGS foobar
   ierr = PetscFree(v_arr);  PETSC_SAMRAI_ERROR(ierr);

   PetscFunctionReturn(0);
}


template <class TYPE>
Vec
PETScAbstractVectorReal<TYPE>::getPETScVector()
{
   return d_petsc_vector;
}// getPETScVector

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecDuplicate(
   Vec v,
   Vec* newv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(v == (Vec)NULL));
#endif

   PETScAbstractVectorReal<TYPE>* new_pav = PABSVEC_CAST(v)->makeNewVector();
   *newv = new_pav->getPETScVector();
   
   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*newv)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


   template <class TYPE>
   PetscErrorCode
   PETScAbstractVectorReal<TYPE>::vecDot(
      Vec x,
      Vec y,
      TYPE* val)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif

   *val = PABSVEC_CAST(x)->dotWith(PABSVEC_CAST(y));

   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMDot(
   Vec x,
   PetscInt nv,
   const Vec* y,
   TYPE* val)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(x == (Vec)NULL));
   for (int i = 0; i < nv; i++) {
      TBOX_ASSERT(!(y[i] == static_cast<Vec>(NULL)));
   }
#endif

   for (PetscInt i = 0; i < nv; ++i)
   {
      val[i] = PABSVEC_CAST(x)->dotWith(PABSVEC_CAST(y[i]));
   }

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecNorm(
   Vec x,
   NormType type,
   TYPE* val)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   if (type == NORM_1)
   {
      *val = PABSVEC_CAST(x)->L1Norm();
   }
   else if (type == NORM_2)
   {
      *val = PABSVEC_CAST(x)->L2Norm();
   }
   else if (type == NORM_INFINITY)
   {
      *val = PABSVEC_CAST(x)->maxNorm();
   }
   else if (type == NORM_1_AND_2)
   {
      val[0] = PABSVEC_CAST(x)->L1Norm();
      val[1] = PABSVEC_CAST(x)->L2Norm();
   }
   else
   {
      TBOX_ERROR("PETScAbstractVectorReal<TYPE>::norm()\n" <<
		 "  vector norm type " << static_cast<int>(type) << " unsupported" << std::endl);
   }

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecTDot(
   Vec x,
   Vec y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   *val = PABSVEC_CAST(x)->TdotWith(PABSVEC_CAST(y));

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMTDot(
   Vec x,
   PetscInt nv,
   const Vec* y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   for (PetscInt i = 0; i < nv; ++i)
   {
      TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
   }
#endif
   for (PetscInt i = 0; i < nv; ++i)
   {
      val[i] = PABSVEC_CAST(x)->TdotWith(PABSVEC_CAST(y[i]));
   }

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecScale(
   Vec x,
   TYPE alpha)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(x)->scaleVector(alpha);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecCopy(
   Vec x,
   Vec y)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(y)->copyVector(PABSVEC_CAST(x));

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSet(
   Vec x,
   TYPE alpha)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(x)->setToScalar(alpha);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSwap(
   Vec x,
   Vec y)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(x)->swapWith(PABSVEC_CAST(y));

   int ierr;
   ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecAXPY(
   Vec y,
   TYPE alpha,
   Vec x)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(y)->setAXPY(alpha, PABSVEC_CAST(x));

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecAXPBY(
   Vec y,
   TYPE alpha,
   TYPE beta,
   Vec x)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(y)->setAXPBY(alpha, PABSVEC_CAST(x), beta);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMAXPY(
   Vec y,
   PetscInt nv,
   const TYPE* alpha,
   Vec* x)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
   for (PetscInt i = 0; i < nv; ++i)
   {
      TBOX_ASSERT(x[i] != static_cast<Vec>(NULL));
   }
#endif
   for (PetscInt i = 0; i < nv; ++i)
   {
      PABSVEC_CAST(y)->setAXPY(alpha[i], PABSVEC_CAST(x[i]));
   }

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecAYPX(
   Vec y,
   const TYPE alpha,
   Vec x)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(y)->setAXPBY(1.0, PABSVEC_CAST(x), alpha);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecWAXPY(
   Vec w,
   TYPE alpha,
   Vec x,
   Vec y)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
   TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(w)->setWAXPY(alpha, PABSVEC_CAST(x), PABSVEC_CAST(y));

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPointwiseMult(
   Vec w,
   Vec x,
   Vec y)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
   TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(w)->pointwiseMultiply(PABSVEC_CAST(x), PABSVEC_CAST(y));
   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPointwiseDivide(
   Vec w,
   Vec x,
   Vec y)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
   TBOX_ASSERT(w != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(w)->pointwiseDivide(PABSVEC_CAST(x), PABSVEC_CAST(y));

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecGetArray(
   Vec x,
   TYPE** a)
{
   NULL_USE(x);
   *a = NULL;

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecGetSize(
   Vec x,
   PetscInt* size)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   *size = PABSVEC_CAST(x)->getDataSize();

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecGetLocalSize(
   Vec x,
   PetscInt* size)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   *size = PABSVEC_CAST(x)->getLocalDataSize();

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


// SGS this looks odd
template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecRestoreArray(
   Vec x,
   TYPE** a)
{
   NULL_USE(x);
   *a = NULL;
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMax(
   Vec x,
   PetscInt* p,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(x)->vecMax(*p, *val);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMin(
   Vec x,
   PetscInt* p,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   PABSVEC_CAST(x)->vecMin(*p, *val);

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetRandom(
   Vec x,
   PetscRandom rctx)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif
   TYPE lo, hi;
   int ierr;
   ierr = PetscRandomGetInterval(rctx, &lo, &hi); PETSC_SAMRAI_ERROR(ierr);
   PABSVEC_CAST(x)->setRandomValues(hi-lo, lo);

   ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecDestroy(
   Vec v)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif

   PABSVEC_CAST(v)->freeVector();

   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecView(
   Vec v,
   PetscViewer viewer)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(v != static_cast<Vec>(NULL));
#endif

   NULL_USE(viewer);
   PABSVEC_CAST(v)->viewVector();

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecDot_local(
   Vec x,
   Vec y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif

   *val = PABSVEC_CAST(x)->dotWith(PABSVEC_CAST(y),true);

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecTDot_local(
   Vec x,
   Vec y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif

   *val = PABSVEC_CAST(x)->TdotWith(PABSVEC_CAST(y),true);

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecNorm_local(
   Vec x,
   NormType type,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
#endif

   if (type == NORM_1)
   {
      *val = PABSVEC_CAST(x)->L1Norm(true);
   } else if (type == NORM_2)   {
      *val = PABSVEC_CAST(x)->L2Norm(true);
   } else if (type == NORM_INFINITY) {
      *val = PABSVEC_CAST(x)->maxNorm(true);
   }else if (type == NORM_1_AND_2) {
      val[0] = PABSVEC_CAST(x)->L1Norm(true);
      val[1] = PABSVEC_CAST(x)->L2Norm(true);
   } else {
      TBOX_ERROR("PETScAbstractVectorReal<TYPE>::norm()\n" <<
		 "  vector norm type " << static_cast<int>(type) << " unsupported" << std::endl);
   }

   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x)); PETSC_SAMRAI_ERROR(ierr);
   PetscFunctionReturn(0);
}// VecNorm_local

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMDot_local(
   Vec x,
   PetscInt nv,
   const Vec* y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   for (PetscInt i = 0; i < nv; ++i)
   {
      TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
   }
#endif
   for (PetscInt i = 0; i < nv; ++i)
   {
      val[i] = PABSVEC_CAST(x)->dotWith(PABSVEC_CAST(y[i]),true);
   }

   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMTDot_local(
   Vec x,
   PetscInt nv,
   const Vec* y,
   TYPE* val)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   for (PetscInt i = 0; i < nv; ++i)
   {
      TBOX_ASSERT(y[i] != static_cast<Vec>(NULL));
   }
#endif

   for (PetscInt i = 0; i < nv; ++i)
   {
      val[i] = PABSVEC_CAST(x)->TdotWith(PABSVEC_CAST(y[i]),true);
   }

   PetscFunctionReturn(0);
}


template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecMaxPointwiseDivide(
   Vec x,
   Vec y,
   TYPE* max)
{
#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(x != static_cast<Vec>(NULL));
   TBOX_ASSERT(y != static_cast<Vec>(NULL));
#endif
   *max = PABSVEC_CAST(x)->maxPointwiseDivide(PABSVEC_CAST(y));

   PetscFunctionReturn(0);
}

///
/// The remaining functions are not implemented and will result in an
/// unrecoverable exception being thrown and program abort if called.
///

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetValues(
   Vec x,
   PetscInt ni,
   const PetscInt* ix,
   const TYPE* y,
   InsertMode iora)
{
   NULL_USE(x);
   NULL_USE(ni);
   NULL_USE(ix);
   NULL_USE(y);
   NULL_USE(iora);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetValues() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecAssemblyBegin(
   Vec vec)
{
   NULL_USE(vec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecAssemblyBegin() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecAssemblyEnd(
   Vec vec)
{
   NULL_USE(vec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecAssemblyEnd() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetOption(
   Vec x,
   VecOption op)
{
   NULL_USE(x);
   NULL_USE(op);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetOption() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetValuesBlocked(
   Vec x,
   PetscInt ni,
   const PetscInt* ix,
   const TYPE* y,
   InsertMode iora)
{
   NULL_USE(x);
   NULL_USE(ni);
   NULL_USE(ix);
   NULL_USE(y);
   NULL_USE(iora);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetValuesBlocked() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPlaceArray(
   Vec vec,
   const TYPE* array)
{
   NULL_USE(vec);
   NULL_USE(array);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecPlaceArray() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecReplaceArray(
   Vec vec,
   const TYPE* array)
{
   NULL_USE(vec);
   NULL_USE(array);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecReplaceArray() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecReciprocal(
   Vec vec)
{
   NULL_USE(vec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecReciprocal() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecViewNative(
   Vec v,
   PetscViewer viewer)
{
   NULL_USE(v);
   NULL_USE(viewer);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecViewNative() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecConjugate(
   Vec x)
{
   NULL_USE(x);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecConjugate() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetLocalToGlobalMapping(
   Vec x,
   ISLocalToGlobalMapping mapping)
{
   NULL_USE(x);
   NULL_USE(mapping);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetLocalToGlobalMapping() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetValuesLocal(
   Vec x,
   PetscInt ni,
   const PetscInt* ix,
   const TYPE* y,
   InsertMode iora)
{
   NULL_USE(x);
   NULL_USE(ni);
   NULL_USE(ix);
   NULL_USE(y);
   NULL_USE(iora);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetValuesLocal() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecResetArray(
   Vec vec)
{
   NULL_USE(vec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecResetArray() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecSetFromOptions(
   Vec vec)
{
   NULL_USE(vec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecSetFromOptions() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecLoad(
   PetscViewer viewer,
   VecType outtype,
   Vec* newvec)
{
   NULL_USE(viewer);
   NULL_USE(outtype);
   NULL_USE(newvec);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecLoad() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPointwiseMax(
   Vec w,
   Vec x,
   Vec y)
{
   NULL_USE(w);
   NULL_USE(x);
   NULL_USE(y);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecPointwiseMax() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPointwiseMaxAbs(
   Vec w,
   Vec x,
   Vec y)
{
   NULL_USE(w);
   NULL_USE(x);
   NULL_USE(y);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecPointwiseMaxAbs() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecPointwiseMin(
   Vec w,
   Vec x,
   Vec y)
{
   NULL_USE(w);
   NULL_USE(x);
   NULL_USE(y);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecPointwiseMin() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

template <class TYPE>
PetscErrorCode
PETScAbstractVectorReal<TYPE>::vecGetValues(
   Vec x,
   PetscInt ni,
   const PetscInt* ix,
   PetscScalar* y)
{
   NULL_USE(x);
   NULL_USE(ni);
   NULL_USE(ix);
   NULL_USE(y);
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::vecGetValues() unimplemented" <<std::endl);
   PetscFunctionReturn(0);
}

}
}
#endif

#endif
