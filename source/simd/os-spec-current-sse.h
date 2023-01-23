#ifndef _OS_SPEC_CURRENT_H
#define _OS_SPEC_CURRENT_H

#include "vector-sse.h"
#include "os-spec-push-sse.h"

void vwl_s1( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] );
void vwl_s2( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] );
void vwl_s3( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] );
void vwl_s4( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] );

/* Current deposition routines for x86 SSE */

void vdepcurrent_1d_s1(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s2(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s3(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s4(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf1D * const part);


void vdepcurrent_2d_s1(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s2(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s3(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s4(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_split_buf2D * const part);

void vdepcurrent_3d_s1(float * const current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s2(float * const current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s3(float * const current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s4(float * const current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part);


#endif
