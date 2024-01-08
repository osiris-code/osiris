/*****************************************************************************************

Current deposition, SIMD optimized version

*****************************************************************************************/

#ifdef SIMD



/**************************************** SSE *******************************************/

#ifdef SIMD_SSE

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-sse.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-ssed.c"

#endif 

#endif

/**************************************** AVX *******************************************/

#ifdef SIMD_AVX

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-avx.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-avxd.c"

#endif 

#endif

/************************************ BlueGene/Q ****************************************/

#ifdef SIMD_BGQ
  #define HAS_SIMD_CODE 1
  #include "os-spec-current-bgq.c"
#endif 

/************************************ Intel MIC ****************************************/

#ifdef SIMD_MIC
#define HAS_SIMD_CODE 1

#include "os-spec-current-mic.c"

#endif 

/************************************ Intel KNL ****************************************/

#ifdef SIMD_AVX512

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-avx512.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-avx512d.c"

#endif 

#endif

/**************************************** AVX2 *****************************************/

#ifdef SIMD_AVX2

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-avx2.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-avx2d.c"

#endif 

#endif

/*********************************** Sanity Check ***************************************/

#ifndef HAS_SIMD_CODE
  #error When using SIMD code define SIMD = SSE | AVX | AVX2 | BGQ | MIC | KNL in the config file
#endif

#endif
