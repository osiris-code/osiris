/****************************************************************************************/
/* SIMD utilites.                                                                       */
/* Definitions for BlueGene/Q                                                           */
/****************************************************************************************/

#ifndef _VECTOR_BGQ_H
#define _VECTOR_BGQ_H

#include <stdio.h>

#define VEC_WIDTH 4

typedef vector4double vector;

#define DECLARE_ALIGNED_32(v)       v __attribute__ ((aligned (32)))

/*
Depending on the precision chosen t_real will be set to float or double
*/

#ifdef PRECISION_SINGLE
typedef float t_real ;
#else
typedef double t_real ;
#endif

/*
BG/Q has no (single) instructions for unaligned memory access, you need to use these
macros instead. These work for both double and single precision data, if you need integers
just replace vec_ld with vec_ldi
*/

#define VEC_LDU( addr, v ) { \
  vector4double v1, v2, vp; \
  /* Load the left part of the vector */ \
  v1 = vec_ld( 0x00, addr ); \
  /* Load the right part of the vector */ \
  v2 = vec_ld( 0x00, add + 4 ); \
  /* Generate control value */ \
  vp = vec_lvsl( 0x00, addr ); \
  /* Generate the aligned vector */ \
  v = vec_perm( v1, v2, vp ); \
}

#define VEC_STU( v, addr ) { \
   vector4double v1, v2, v3, p, m1, m2, m3; \
   /* load existing data */ \
   v1 = vec_ld( 0x00, addr     ); \
   v2 = vec_ld( 0x00, addr + 4 ); \
   /* generate permutation mask */ \
   p = vec_lvsr( 0x00, addr ); \
   /* generate vector of all FALSE */ \
   m1 = vec_cmplt(p,p); \
   /* generate vector of all TRUE */ \
   m2 = vec_neg(m1); \
   /* Generate insert mask */ \
   m3 = vec_perm(m1,m2,p); \
   /* permute and insert */ \
   v3 = vec_perm(v,v,p); \
   v1 = vec_sel(v1,v3,m3); \
   v2 = vec_sel(v3,v2,m3); \
   /* store data back */ \
   vec_st( 0x00, addr    , v1); \
   vec_st( 0x00, addr + 4, v2); \
}


/* 
Number of Newton-Raphson refinements for the full precision operations. The maximum
is 2 (default) that for a = M_PI gives exactly the same results as the normal scalar op. */

#ifndef __MATH_PRECISION__
#define __MATH_PRECISION__ 2
#endif

/* Full precision reciprocal, reciprocal square root                                    */

/*   vec_r - Parallel Reciprocal  */
/*

In terms of performance this needs to be compared with

vec_r( a )  =  vec_swdiv_nochk( vec_splats( 1.0 ), a );

*/

inline vector vec_r( const vector a )
{
   
   // Get the reciprocal estimate
   vector estimate = vec_re( a );

#if __MATH_PRECISION__ > 0   
   vector one     = vec_splats( 1.0 );
   estimate = vec_madd( vec_nmsub( estimate, a, one ), estimate, estimate );
#endif

#if __MATH_PRECISION__ > 1   
   estimate = vec_madd( vec_nmsub( estimate, a, one ), estimate, estimate );
#endif

   return estimate;
}

/* vec_rsqrt - Parallel Reciprocal square root */

/*

In terms of performance this needs to be compared with

vec_rsqrt( a )  =  vec_swdiv_nochk( vec_splats( 1.0 ), vec_swsqrt_nochk( a ) );

*/


inline vector vec_rsqrt( const vector a )
{
   // Get the reciprocal square root estimate
   vector estimate = vec_rsqrte( a );

#if __MATH_PRECISION__ > 0   
   vector one     = vec_splats(1.0);
   vector oneHalf = vec_splats(0.5);
  
   // One round of Newton-Raphson refinement
   vector estimateSquared = vec_mul( estimate, estimate );
   vector halfEstimate    = vec_mul( estimate, oneHalf );
   estimate = vec_madd( vec_nmsub( a, estimateSquared, one ), halfEstimate, estimate );
#endif

#if __MATH_PRECISION__ > 1   
   estimateSquared = vec_mul( estimate, estimate );
   halfEstimate    = vec_mul( estimate, oneHalf );
   estimate = vec_madd( vec_nmsub( a, estimateSquared, one ), halfEstimate, estimate );
#endif
   
   return estimate;
}


/* vec_reduce_add - add all vector elements */

/*

I could not find a better way of doing this with QPX intrinsics.

*/

inline double vec_reduce_add( const vector a )
{
   double sum;
   
   sum = vec_extract( a, 0 ) + vec_extract( a, 1 ) +
	     vec_extract( a, 2 ) + vec_extract( a, 3 ) ;

   return sum;
}


/***********************************************************************
VEC_LD4v3

Loads 4 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4] -> VEC_LD4v3(v1,v2,v3,b)

-> v1 = [x1,x2,x3,x4]
-> v2 = [y1,y2,y3,y4]
-> v3 = [z1,z2,z3,z4]

************************************************************************/

#ifdef PRECISION_SINGLE

#define VEC_LDA4v3(v1, v2, v3, buf) {               \
                                                   \
  register vector a0,a1,a2, t0, t1;                \
                                                   \
  a0 = vec_lda( 0x00, (float*) buf );              \
  a1 = vec_lda( 0x10, (float*) buf );              \
  a2 = vec_lda( 0x20, (float*) buf );              \
                                                   \
  t0 = vec_perm( a1, a2, vec_gpci( 02536 ) );      \
  t1 = vec_perm( a0, a1, vec_gpci( 01425 ) );      \
                                                   \
  v1 = vec_perm( a0, t0, vec_gpci( 00345 ) );      \
  v2 = vec_perm( t0, t1, vec_gpci( 04523 ) );      \
  v3 = vec_perm( t1, a2, vec_gpci( 02347 ) );      \
}


#else

#define VEC_LDA4v3(v1, v2, v3, buf) {               \
                                                   \
  register vector a0,a1,a2, t0, t1;                \
                                                   \
  a0 = vec_lda( 0x00, (double*) buf );             \
  a1 = vec_lda( 0x20, (double*) buf );             \
  a2 = vec_lda( 0x40, (double*) buf );             \
                                                   \
  t0 = vec_perm( a1, a2, vec_gpci( 02536 ) );      \
  t1 = vec_perm( a0, a1, vec_gpci( 01425 ) );      \
                                                   \
  v1 = vec_perm( a0, t0, vec_gpci( 00345 ) );      \
  v2 = vec_perm( t0, t1, vec_gpci( 04523 ) );      \
  v3 = vec_perm( t1, a2, vec_gpci( 02347 ) );      \
}

#endif

/***********************************************************************
VEC_ST4v3

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,x3,x4], v2 = [y1,y2,y3,y4], v3 = [z1,z2,z3,z4]

-> VEC_ST4v3(b, v1, v2, v3)

-> b = [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4]

************************************************************************/

#ifdef PRECISION_SINGLE

#define VEC_STA4v3(buf, v1, v2, v3) {                                      \
                                                                          \
  register vector t0, t1;                                                 \
                                                                          \
  t0 = vec_perm( v2, v3, vec_gpci( 00415 ) );                             \
  t1 = vec_perm( v1, v2, vec_gpci( 03726 ) );                             \
                                                                          \
  vec_sta( vec_perm( v1, t0, vec_gpci( 00451 ) ), 0x00, (float*) buf );   \
  vec_sta( vec_perm( t0, t1, vec_gpci( 02367 ) ), 0x10, (float*) buf );   \
  vec_sta( vec_perm( t1, v3, vec_gpci( 06017 ) ), 0x20, (float*) buf );   \
}


#else

#define VEC_STA4v3(buf, v1, v2, v3) {                                      \
                                                                          \
  register vector t0, t1;                                                 \
                                                                          \
  t0 = vec_perm( v2, v3, vec_gpci( 00415 ) );                             \
  t1 = vec_perm( v1, v2, vec_gpci( 03726 ) );                             \
                                                                          \
  vec_sta( vec_perm( v1, t0, vec_gpci( 00451 ) ), 0x00, (double*) buf );  \
  vec_sta( vec_perm( t0, t1, vec_gpci( 02367 ) ), 0x20, (double*) buf );  \
  vec_sta( vec_perm( t1, v3, vec_gpci( 06017 ) ), 0x40, (double*) buf );  \
}

#endif



/***********************************************************************
VEC_LDA4v2

Loads 2 x 4 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,x2,y2,x3,y3,x4,y4] -> VEC_LD4v2(v1,v2,b)

-> v1 = [x1,x2,x3,x4]
-> v2 = [y1,y2,y3,y4]

************************************************************************/

#ifdef PRECISION_SINGLE

#define VEC_LDA4v2(v1, v2, buf) {                       \
                                                       \
  register vector a0,a1;                               \
                                                       \
  a0 = vec_lda( 0x00, (float *) buf );                 \
  a1 = vec_lda( 0x10, (float *) buf );                 \
                                                       \
  (v1) = vec_perm( a0, a1, vec_gpci( 00246 ) );        \
  (v2) = vec_perm( a0, a1, vec_gpci( 01357 ) );        \
}


#else

#define VEC_LDA4v2(v1, v2, buf) {                       \
                                                       \
  register vector a0,a1;                               \
                                                       \
  a0 = vec_lda( 0x00, (double *) buf );                \
  a1 = vec_lda( 0x20, (double *) buf );                \
                                                       \
  (v1) = vec_perm( a0, a1, vec_gpci( 00246 ) );        \
  (v2) = vec_perm( a0, a1, vec_gpci( 01357 ) );        \
}

#endif

/***********************************************************************
VEC_STA4v2

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,x3,x4], v2 = [y1,y2,y3,y4]

-> VEC_ST4v2(b, v1, v2)

-> b = [x1,y1,x2,y2,x3,y3,x4,y4] 

*********************************************************************/

#ifdef PRECISION_SINGLE

#define VEC_STA4v2(buf, v1, v2) {                                          \
                                                                          \
  vec_sta( vec_perm( v1, v2, vec_gpci( 00415 ) ), 0x00, (float *) buf );  \
  vec_sta( vec_perm( v1, v2, vec_gpci( 02637 ) ), 0x10, (float *) buf );  \
}


#else

#define VEC_STA4v2(buf, v1, v2) {                                          \
                                                                          \
  vec_sta( vec_perm( v1, v2, vec_gpci( 00415 ) ), 0x00, (double *) buf ); \
  vec_sta( vec_perm( v1, v2, vec_gpci( 02637 ) ), 0x20, (double *) buf ); \
}

#endif

/***********************************************************************
VEC_ZERO

Set vector to 0.
***********************************************************************/

//#define VEC_ZERO( v ) ( v = vec_xor(v,v))

#define VEC_ZERO( v ) ( v = vec_splats(0.0) )


/***********************************************************************
printfv

Print vector components
***********************************************************************/

static void printfv(char *tag, vector v)
{
  printf("%s = %f, %f, %f, %f \n", tag, vec_extract(v,0), vec_extract(v,1), 
                                        vec_extract(v,2), vec_extract(v,3) );
}

/***********************************************************************
test vector components

***********************************************************************/
#define TEST_VECTOR( a ) {							             \
if ( isinf((a)[0]) || isinf((a)[1]) || isinf((a)[2]) || isinf((a)[3]) || \
     isnan((a)[0]) || isnan((a)[1]) || isnan((a)[2]) || isnan((a)[3]) ) {		\
     printf( "%s:%d %s = %g, %g, %g, %g\n", __FILE__, __LINE__, #a, (a)[0], (a)[1], (a)[2], (a)[3] );	\
     exit(-1);										\
   }												\
}

#endif
