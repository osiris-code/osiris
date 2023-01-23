/*****************************************************************************************

Relativistic particle pusher, SIMD optimized version

*****************************************************************************************/

#ifdef SIMD

/**************************************** SSE *******************************************/

#ifdef SIMD_SSE

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-push-sse.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-push-ssed.c"

#endif 

#endif

/**************************************** AVX *******************************************/

#ifdef SIMD_AVX

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-push-avx.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-push-avxd.c"

#endif 

#endif


/************************************ BlueGene/Q ****************************************/

#ifdef SIMD_BGQ

#define HAS_SIMD_CODE 1
#include "os-spec-push-bgq.c"

#endif

/************************************* Intel MIC ****************************************/

#ifdef SIMD_MIC

#define HAS_SIMD_CODE 1
#include "os-spec-push-mic.c"

#endif 

/*********************************** Intel AVX512 ***************************************/

#ifdef SIMD_AVX512

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-push-avx512.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-push-avx512d.c"

#endif 

#endif 

/*************************************** AVX2 ******************************************/

#ifdef SIMD_AVX2

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-push-avx2.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-push-avx2d.c"

#endif 

#endif 


/*********************************** Sanity Check ***************************************/

#ifndef HAS_SIMD_CODE

#error SIMD code is requested, but no specific hardware architecture has been defined
#error When using SIMD code define SIMD = SSE | AVX | AVX2 | BGQ | MIC | KNL in the config file

#endif

#endif
