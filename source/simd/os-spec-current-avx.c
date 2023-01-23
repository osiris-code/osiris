/*****************************************************************************************

Charge conserving current deposition, AVX optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avx.h"

#include "vector-avx.h"
#include "splines-avx.h"


/***********************************************************************
vwl_s1
***********************************************************************/

inline void vwl_s1( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{

  vwl[0] = _mm256_mul_ps( vqn, _mm256_sub_ps( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2
***********************************************************************/

inline void vwl_s2( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{

  __m256 const c1_2 = _mm256_set1_ps( 0.5f );
  
  __m256 d    = _mm256_sub_ps( vx1, vx0 );
  __m256 s1_2 = _mm256_sub_ps( c1_2, vx0 );
  __m256 p1_2 = _mm256_add_ps( c1_2, vx0 );
  
  __m256 n = _mm256_mul_ps( vqn, d );
  
  vwl[0] = _mm256_mul_ps( n, _mm256_sub_ps( s1_2, _mm256_mul_ps( c1_2, d ) ) );
  vwl[1] = _mm256_mul_ps( n, _mm256_add_ps( p1_2, _mm256_mul_ps( c1_2, d ) ) );  

}



/***********************************************************************
vwl_s3
***********************************************************************/

inline void vwl_s3( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{
   __m256 const c1_2 = _mm256_set1_ps( 0.5f );
   __m256 const c3_4 = _mm256_set1_ps( 0.75f );
   __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );
   
   __m256   d = _mm256_sub_ps(  vx1, vx0 );
   __m256   s = _mm256_sub_ps( c1_2, vx0 );
   __m256   p = _mm256_add_ps( c1_2, vx0 );
   __m256 d_3 = _mm256_mul_ps( c1_3,   d );

   vwl[0] = _mm256_mul_ps( c1_2, _mm256_sub_ps( _mm256_mul_ps( s, s ),       _mm256_mul_ps( d, _mm256_sub_ps( s, d_3 ) ) ) );
   vwl[1] = _mm256_sub_ps( _mm256_sub_ps( c3_4, _mm256_mul_ps( vx0, vx0 ) ), _mm256_mul_ps( d, _mm256_add_ps( vx0, d_3 ) ) );
   vwl[2] = _mm256_mul_ps( c1_2, _mm256_add_ps( _mm256_mul_ps( p, p ),       _mm256_mul_ps( d, _mm256_add_ps( p, d_3 ) ) ) );

   __m256 n = _mm256_mul_ps( vqn, d );
   vwl[0] = _mm256_mul_ps( n, vwl[0] );
   vwl[1] = _mm256_mul_ps( n, vwl[1] );
   vwl[2] = _mm256_mul_ps( n, vwl[2] );


}

/***********************************************************************
vwl_s4
***********************************************************************/
inline void vwl_s4( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{
   __m256 const c1_2 = _mm256_set1_ps( 0.5f );
   __m256 const c1_4 = _mm256_set1_ps( 0.25f );
   __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
   __m256 const c3_2 = _mm256_set1_ps( 1.5f );
   __m256 const c3_4 = _mm256_set1_ps( 0.75f );
   __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );
   __m256 const c2_3 = _mm256_set1_ps( 2.0f/3.0f );
   
   __m256 d = _mm256_sub_ps( vx1, vx0 );
   __m256 s = _mm256_sub_ps( c1_2, vx0 );
   __m256 p = _mm256_add_ps( c1_2, vx0 );

   __m256 t = _mm256_sub_ps( vx0, c1_6 );
   __m256 u = _mm256_add_ps( vx0, c1_6 );

   __m256 s2 = _mm256_mul_ps( s, s );
   __m256 s3 = _mm256_mul_ps( s2, s );

   __m256 p2 = _mm256_mul_ps( p, p );
   __m256 p3 = _mm256_mul_ps( p2, p );

   __m256 d_2 = _mm256_mul_ps( d, c1_2 );
   __m256 d_4 = _mm256_mul_ps( d, c1_4 );
   

   vwl[0] = _mm256_mul_ps( c1_6, _mm256_add_ps( s3, _mm256_mul_ps( d, _mm256_sub_ps( _mm256_mul_ps( d, _mm256_sub_ps( s, d_4 ) ),  _mm256_mul_ps( c3_2, s2 ) ) ) ) );
   vwl[1] = _mm256_add_ps( _mm256_add_ps( _mm256_sub_ps( c2_3, p2 ), _mm256_mul_ps( c1_2, p3 ) ), _mm256_mul_ps( d, _mm256_add_ps( _mm256_sub_ps( _mm256_mul_ps( c3_4, _mm256_mul_ps( t, t ) ), c1_3 ), _mm256_mul_ps( d_2, _mm256_add_ps( t, d_4 ) ) ) ) );
   vwl[2] = _mm256_sub_ps( _mm256_add_ps( _mm256_sub_ps( c2_3, s2 ), _mm256_mul_ps( c1_2, s3 ) ), _mm256_mul_ps( d, _mm256_add_ps( _mm256_sub_ps( _mm256_mul_ps( c3_4, _mm256_mul_ps( u, u ) ), c1_3 ), _mm256_mul_ps( d_2, _mm256_add_ps( u, d_4 ) ) ) ) );
   vwl[3] = _mm256_mul_ps( c1_6, _mm256_add_ps( p3, _mm256_mul_ps( d, _mm256_add_ps( _mm256_mul_ps( d, _mm256_add_ps( p, d_4 ) ),  _mm256_mul_ps( c3_2, p2 ) ) ) ) );

   __m256 n = _mm256_mul_ps( vqn, d );
   vwl[0] = _mm256_mul_ps( n, vwl[0] );
   vwl[1] = _mm256_mul_ps( n, vwl[1] );
   vwl[2] = _mm256_mul_ps( n, vwl[2] );
   vwl[3] = _mm256_mul_ps( n, vwl[3] );

}

/****************************************************************************************

  Generate specific current deposition functions for 1D, 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o

/********************************** Linear interpolation ********************************/

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Current grid offset (1 + ORDER/2)
#define OFFSET 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ORDER 4
#define NP 5
#define OFFSET 3

#include __FILE__

#else

/****************************************************************************************

  Template function definitions for 2D and 3D current deposition

****************************************************************************************/

/*

Note: In 2D there was a division by 2 removed from the calculation of the perpendicular
      weights, so jnorm must be dx / dt / (2*NORM). In 3D it was a division by 3, so 
      jnorm must be dx / dt / (3*NORM). To avoid repeating this every time the 
      DEP_CURRENT functions are called, this is done in the ADVANCE_DEPOSIT functions.

*/

/***************   Generate Function names based on interpolation order  ****************/

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

#define SPLINE  ONAME( vspline, ORDER )
#define WL      ONAME( vwl, ORDER )

/*
  Reference implementation, serial deposition
*/

inline void DEP_CURRENT_1D
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_split_buf1D * const part) {
        
  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *p0; 
    
  fvec j2[NP], j3[NP];
  
  fvec vwl1[ORDER];
  
  __m256 vx0, vx1, vq, vvy, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
      
  __m256 vqnx, vqvy, vqvz;
  __m256 vs0x[NP], vs1x[NP]; 

  unsigned int i, np;

  __m256 const c1_2  = _mm256_set1_ps( 0.5f );

  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 8 add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;
       
       part->x0[k] = part->x1[k] = 0.0f;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0f;
       
       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    unsigned k, k1;
    __m256i vix;

    // load 8 particles
    LOAD8P1D( part, i, vx0, vx1, vq, vvy, vvz, vix );
    
    // Store cell index
    _mm256_store_si256( (__m256i *) idx, vix ); 
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    
    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqvz = _mm256_mul_ps( vq, vvz );
    vqvy = _mm256_mul_ps( vq, vvy );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
     
    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m256 vwp1 = _mm256_mul_ps( c1_2, _mm256_add_ps(vs0x[k1], vs1x[k1]) );
      j2[k1].v8 = _mm256_mul_ps( vqvy, vwp1 );
      j3[k1].v8 = _mm256_mul_ps( vqvz, vwp1 );
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      p0 = pj + idx[k];

      // accumulate j1
      for ( k1 = 0; k1 < ORDER; k1++ ) {
        p0[k1].j1 += vwl1[k1].v[k] ;
      }

      // accumulate j2, j3
      for ( k1 = 0; k1 < NP; k1++ ) {
        p0[k1].j2 += j2[k1].v[k];
        p0[k1].j3 += j3[k1].v[k];
      }

    }

  }
}


#if 1

#warning Using mk. I AVX current deposition

inline void DEP_CURRENT_2D
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_split_buf2D * const part)
{
        
  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  
  fvec j3[NP][NP];
  fvec vwp1[NP], vwp2[NP]; 
  fvec vwl1[ORDER], vwl2[ORDER];
  
  __m256 vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
		  
  __m256 vqnx, vqny, vqvz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 
 
  unsigned int np;
  unsigned int i, k, k1, k2;

  __m256 const oneThird = _mm256_set1_ps( 1.0f/3.0f );

  unsigned int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       k = np + i;
       
       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;
       
       part->ix[k] = part->iy[k] = 1;
    }
  }
  
  // Each vp is 8 floats
  for( i = 0; i < np; i+=VEC_WIDTH ) {
    
    __m256i vix, viy;
    
    // load 8 particles
    LOAD8P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
    
    // Store cell index
    viadd_store( idx, vix, vimul_s( viy, Dy ) ) ;
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqvz = _mm256_mul_ps( vq, vvz);
    vqvz = _mm256_mul_ps( vqvz, oneThird );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v8 = _mm256_add_ps(vs0y[k], vs1y[k]);
      vwp2[k].v8 = _mm256_add_ps(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k1 = 0; k1 < NP; k1++ ) {
      for ( k2 = 0; k2 < NP; k2++ ) {
  		
  		__m256 const oneHalf  = _mm256_set1_ps( 0.5f );
  		__m256  tmp1, tmp2;
  		
		tmp1 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs1y[k2] ) );
	    
		tmp2 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs0y[k2] ) );

		j3[k1][k2].v8 = _mm256_mul_ps( vqvz, _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf ) ) );
      }
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      // This is not vectorially because it's a 64bit op
      pjpart = pj + idx[k];

	  // accumulate j1
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy;
    #pragma unroll(ORDER)
		for ( k1 = 0; k1 < ORDER; k1++ ) {
		  p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
		}
	  }
  
	  // accumulate j2 - making k2 the outside loop gives marginal perf. gain
	  for( k2 = 0; k2 < ORDER; k2++ ) {
		p0 = pjpart + k2*Dy;
    #pragma unroll(NP)
		for ( k1 = 0; k1 < NP; k1++ ) {
		   p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
		}
	  }
  
	  // accumulate j3
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy;
    #pragma unroll(NP)
		for ( k1 = 0; k1 < NP; k1++ ) {
		   p0[k1].j3 += (j3[k1][k2]).v[k];
		}
	  }
	  
    }
	  
  }
}


inline void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part)
{
        
  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
      
  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
  
  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  fvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  fvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  unsigned int np;
  unsigned int i, k, k1, k2, k3;
     
  unsigned int const Dy = size[0];    
  unsigned int const Dz = Dy * size[1]; 
    
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 + 
                                                  ( offset[1] - OFFSET  ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET  ) * 3 * Dz );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);
  
  np = part -> np;
    
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       k = np + i;
       
       part -> x0[k] = part -> x1[k] = 0.;
       part -> y0[k] = part -> y1[k] = 0.;
       part -> z0[k] = part -> z1[k] = 0.;
       part -> q[k] = 0.;
       part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }


  for( i = 0; i < np; i+=VEC_WIDTH ) {

    __m256i vix, viy, viz;
      
	// load 8 particles
    LOAD8P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
    
    // Do idx[k]= ix[k] + iy[k]*Dy + iz[k]*Dz vectorially 
    {
      __m128i vDy, vDz, a, b;
      vDy = _mm_set1_epi32( Dy );
      vDz = _mm_set1_epi32( Dz );
      
      a = _mm_mullo_epi32( _mm256_castsi256_si128(viz), vDz );
      b = _mm_mullo_epi32( _mm256_extractf128_si256(viz,1), vDz );

      a = _mm_add_epi32( a, _mm_mullo_epi32( _mm256_castsi256_si128(viy), vDy ));
      b = _mm_add_epi32( b, _mm_mullo_epi32( _mm256_extractf128_si256(viy,1), vDy ));
      
      a = _mm_add_epi32( a, _mm256_castsi256_si128(vix) );
      b = _mm_add_epi32( b, _mm256_extractf128_si256(vix,1) );
      
      _mm_store_si128((__m128i *)  &idx[0], a );
      _mm_store_si128((__m128i *)  &idx[4], b );
    }
    
    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqnz = _mm256_mul_ps( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );
    WL( vqnz, vz0, vz1, (__m256 *) vwl3 );
    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
         
         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         const __m256 oneHalf = _mm256_set1_ps(0.5f);
         __m256 tmp1, tmp2;        
        
         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs0z[k2]), _mm256_mul_ps( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs1z[k2]), _mm256_mul_ps( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0z[k2]), _mm256_mul_ps( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1z[k2]), _mm256_mul_ps( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2]), _mm256_mul_ps( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2]), _mm256_mul_ps( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

      }
    }

    // looping by particle on the outside loop yields the best performance
     
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

       pjpart = pj + idx[k];
          
	   // accumulate j1
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
       #pragma unroll(ORDER)
		   for ( k1 = 0; k1 < ORDER; k1++ ) {
			  p0[k1].j1 += vwl1[k1].v[k] * vwp1[k3][k2].v[k];
		   }
		 }
	   }
   
	   // accumulate j2
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < ORDER; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
       #pragma unroll(NP)
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j2 += vwl2[k2].v[k] * vwp2[k3][k1].v[k];
		   }
		 }
	   }
	   
	   // accumulate j3
	   for( k3 = 0; k3 < ORDER; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
       #pragma unroll(NP)
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j3  += vwl3[k3].v[k] * vwp3[k2][k1].v[k];
		   }
		 }
	   }

    }

  }
}

#endif

#if 0

#warning Using mk. II AVX current deposition

#if ( ORDER == 1 ) 

#define ACC_CURRENT2D( k, v, p ) {                                            \
    t_current *p0, *p1;                                                       \
    p0 = pj + idx[k];                                                         \
	p1 = p0 + Dy;                                                             \
	__m256 v0, v1, j0, j1;                                                    \
    j0 = _mm256_castps128_ps256(   _mm256_extractf128_ps( v[0][0], p ));      \
	j0 = _mm256_insertf128_ps( j0, _mm256_extractf128_ps( v[0][1], p ), 1);   \
	j1 = _mm256_castps128_ps256(   _mm256_extractf128_ps( v[1][0], p ));      \
	j1 = _mm256_insertf128_ps( j1, _mm256_extractf128_ps( v[1][1], p ), 1);   \
	v0 = _mm256_castps128_ps256(   _mm_loadu_ps( (float *) (p0+0) ) );        \
	v0 = _mm256_insertf128_ps( v0, _mm_loadu_ps( (float *) (p0+1) ), 1 );     \
	v1 = _mm256_castps128_ps256(   _mm_loadu_ps( (float *) (p1+0) ) );        \
	v1 = _mm256_insertf128_ps( v1, _mm_loadu_ps( (float *) (p1+1) ), 1 );     \
	v0 = _mm256_add_ps( v0, j0 );                                             \
	v1 = _mm256_add_ps( v1, j1 );                                             \
	_mm_storeu_ps( (float *) (p0+0),  _mm256_castps256_ps128( v0 )  );        \
	_mm_storeu_ps( (float *) (p0+1),   _mm256_extractf128_ps( v0, 1 ) );      \
	_mm_storeu_ps( (float *) (p1+0),  _mm256_castps256_ps128( v1 )  );        \
	_mm_storeu_ps( (float *) (p1+1),   _mm256_extractf128_ps( v1, 1 ) );      \
}

#elif ( ORDER == 2 ) 

#define ACC_CURRENT2D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2;                                                  \
    p0 = pj + idx[k];                                                         \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	__m256 v0a, v1a, v2a;                                                     \
	__m128 v0b, v1b, v2b;                                                     \
	__m256 j0a, j1a, j2a;                                                     \
	__m128 j0b, j1b, j2b;                                                     \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b =                            _mm256_extractf128_ps( v[0][2], p );     \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
    j1b =                            _mm256_extractf128_ps( v[1][2], p );     \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
    j2b =                            _mm256_extractf128_ps( v[2][2], p );     \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ), 1);    \
	v0b =                            _mm_loadu_ps( (float *) (p0+2) );        \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ), 1);    \
	v1b =                            _mm_loadu_ps( (float *) (p1+2) );        \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ), 1);    \
	v2b =                            _mm_loadu_ps( (float *) (p2+2) );        \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b =    _mm_add_ps( v0b, j0b );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b =    _mm_add_ps( v1b, j1b );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b =    _mm_add_ps( v2b, j2b );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a ) ) ;        \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ) );      \
	_mm_storeu_ps( (float *) (p0+2),                         v0b ) ;          \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a ) );         \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ) );      \
	_mm_storeu_ps( (float *) (p1+2),                         v1b ) ;          \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a ) );         \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ) );      \
	_mm_storeu_ps( (float *) (p2+2),                         v2b ) ;          \
}

#elif ( ORDER == 3 ) 

#define ACC_CURRENT2D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2, *p3;                                             \
    p0 = pj + idx[k];                                                         \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	p3 = p0 + 3*Dy;                                                           \
	__m256 v0a, v1a, v2a, v3a;                                                \
	__m256 v0b, v1b, v2b, v3b;                                                \
	__m256 j0a, j1a, j2a, j3a;                                                \
	__m256 j0b, j1b, j2b, j3b;                                                \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][2], p ));    \
	j0b = _mm256_insertf128_ps( j0b, _mm256_extractf128_ps( v[0][3], p ), 1); \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
	j1b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][2], p ));    \
	j1b = _mm256_insertf128_ps( j1b, _mm256_extractf128_ps( v[1][3], p ), 1); \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
	j2b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][2], p ));    \
	j2b = _mm256_insertf128_ps( j2b, _mm256_extractf128_ps( v[2][3], p ), 1); \
	j3a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][0], p ));    \
	j3a = _mm256_insertf128_ps( j3a, _mm256_extractf128_ps( v[3][1], p ), 1); \
	j3b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][2], p ));    \
	j3b = _mm256_insertf128_ps( j3b, _mm256_extractf128_ps( v[3][3], p ), 1); \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ), 1);    \
	v0b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+2) ) );      \
	v0b = _mm256_insertf128_ps( v0b, _mm_loadu_ps( (float *) (p0+3) ), 1);    \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ), 1);    \
	v1b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+2) ) );      \
	v1b = _mm256_insertf128_ps( v1b, _mm_loadu_ps( (float *) (p1+3) ), 1);    \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ), 1);    \
	v2b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+2) ) );      \
	v2b = _mm256_insertf128_ps( v2b, _mm_loadu_ps( (float *) (p2+3) ), 1);    \
	v3a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+0) ) );      \
	v3a = _mm256_insertf128_ps( v3a, _mm_loadu_ps( (float *) (p3+1) ), 1);    \
	v3b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+2) ) );      \
	v3b = _mm256_insertf128_ps( v3b, _mm_loadu_ps( (float *) (p3+3) ), 1);    \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b = _mm256_add_ps( v0b, j0b );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b = _mm256_add_ps( v1b, j1b );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b = _mm256_add_ps( v2b, j2b );                                          \
	v3a = _mm256_add_ps( v3a, j3a );                                          \
	v3b = _mm256_add_ps( v3b, j3b );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a )   );       \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ));       \
	_mm_storeu_ps( (float *) (p0+2), _mm256_castps256_ps128( v0b )   );       \
	_mm_storeu_ps( (float *) (p0+3),  _mm256_extractf128_ps( v0b, 1 ));       \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a )   );       \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ));       \
	_mm_storeu_ps( (float *) (p1+2), _mm256_castps256_ps128( v1b )   );       \
	_mm_storeu_ps( (float *) (p1+3),  _mm256_extractf128_ps( v1b, 1 ));       \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a )   );       \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ));       \
	_mm_storeu_ps( (float *) (p2+2), _mm256_castps256_ps128( v2b )   );       \
	_mm_storeu_ps( (float *) (p2+3),  _mm256_extractf128_ps( v2b, 1 ));       \
	_mm_storeu_ps( (float *) (p3+0), _mm256_castps256_ps128( v3a )   );       \
	_mm_storeu_ps( (float *) (p3+1),  _mm256_extractf128_ps( v3a, 1 ));       \
	_mm_storeu_ps( (float *) (p3+2), _mm256_castps256_ps128( v3b )   );       \
	_mm_storeu_ps( (float *) (p3+3),  _mm256_extractf128_ps( v3b, 1 ));       \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT2D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2, *p3, *p4;                                        \
    p0 = pj + idx[k];                                                         \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	p3 = p0 + 3*Dy;                                                           \
	p4 = p0 + 4*Dy;                                                           \
	__m256 v0a, v1a, v2a, v3a, v4a;                                           \
	__m256 v0b, v1b, v2b, v3b, v4b;                                           \
	__m128 v0c, v1c, v2c, v3c, v4c;                                           \
	__m256 j0a, j1a, j2a, j3a, j4a;                                           \
	__m256 j0b, j1b, j2b, j3b, j4b;                                           \
	__m128 j0c, j1c, j2c, j3c, j4c;                                           \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][2], p ));    \
	j0b = _mm256_insertf128_ps( j0b, _mm256_extractf128_ps( v[0][3], p ), 1); \
    j0c =                            _mm256_extractf128_ps( v[0][4], p );     \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
	j1b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][2], p ));    \
	j1b = _mm256_insertf128_ps( j1b, _mm256_extractf128_ps( v[1][3], p ), 1); \
    j1c =                            _mm256_extractf128_ps( v[1][4], p );     \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
	j2b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][2], p ));    \
	j2b = _mm256_insertf128_ps( j2b, _mm256_extractf128_ps( v[2][3], p ), 1); \
    j2c =                            _mm256_extractf128_ps( v[2][4], p );     \
	j3a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][0], p ));    \
	j3a = _mm256_insertf128_ps( j3a, _mm256_extractf128_ps( v[3][1], p ), 1); \
	j3b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][2], p ));    \
	j3b = _mm256_insertf128_ps( j3b, _mm256_extractf128_ps( v[3][3], p ), 1); \
    j3c =                            _mm256_extractf128_ps( v[3][4], p );     \
	j4a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[4][0], p ));    \
	j4a = _mm256_insertf128_ps( j4a, _mm256_extractf128_ps( v[4][1], p ), 1); \
	j4b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[4][2], p ));    \
	j4b = _mm256_insertf128_ps( j4b, _mm256_extractf128_ps( v[4][3], p ), 1); \
    j4c =                            _mm256_extractf128_ps( v[4][4], p );     \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ),1);     \
	v0b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+2) ) );      \
	v0b = _mm256_insertf128_ps( v0b, _mm_loadu_ps( (float *) (p0+3) ),1);     \
	v0c =                            _mm_loadu_ps( (float *) (p0+4) );        \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ),1);     \
	v1b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+2) ) );      \
	v1b = _mm256_insertf128_ps( v1b, _mm_loadu_ps( (float *) (p1+3) ),1);     \
	v1c =                            _mm_loadu_ps( (float *) (p1+4) );        \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ),1);     \
	v2b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+2) ) );      \
	v2b = _mm256_insertf128_ps( v2b, _mm_loadu_ps( (float *) (p2+3) ),1);     \
	v2c =                            _mm_loadu_ps( (float *) (p2+4) );        \
	v3a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+0) ) );      \
	v3a = _mm256_insertf128_ps( v3a, _mm_loadu_ps( (float *) (p3+1) ),1);     \
	v3b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+2) ) );      \
	v3b = _mm256_insertf128_ps( v3b, _mm_loadu_ps( (float *) (p3+3) ),1);     \
	v3c =                            _mm_loadu_ps( (float *) (p3+4) );        \
	v4a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p4+0) ) );      \
	v4a = _mm256_insertf128_ps( v4a, _mm_loadu_ps( (float *) (p4+1) ),1);     \
	v4b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p4+2) ) );      \
	v4b = _mm256_insertf128_ps( v4b, _mm_loadu_ps( (float *) (p4+3) ),1);     \
	v4c =                            _mm_loadu_ps( (float *) (p4+4) );        \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b = _mm256_add_ps( v0b, j0b );                                          \
	v0c =    _mm_add_ps( v0c, j0c );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b = _mm256_add_ps( v1b, j1b );                                          \
	v1c =    _mm_add_ps( v1c, j1c );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b = _mm256_add_ps( v2b, j2b );                                          \
	v2c =    _mm_add_ps( v2c, j2c );                                          \
	v3a = _mm256_add_ps( v3a, j3a );                                          \
	v3b = _mm256_add_ps( v3b, j3b );                                          \
	v3c =    _mm_add_ps( v3c, j3c );                                          \
	v4a = _mm256_add_ps( v4a, j4a );                                          \
	v4b = _mm256_add_ps( v4b, j4b );                                          \
	v4c =    _mm_add_ps( v4c, j4c );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a )   );       \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ));       \
	_mm_storeu_ps( (float *) (p0+2), _mm256_castps256_ps128( v0b )   );       \
	_mm_storeu_ps( (float *) (p0+3),  _mm256_extractf128_ps( v0b, 1 ));       \
	_mm_storeu_ps( (float *) (p0+4),                         v0c );           \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a )   );       \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ));       \
	_mm_storeu_ps( (float *) (p1+2), _mm256_castps256_ps128( v1b )   );       \
	_mm_storeu_ps( (float *) (p1+3),  _mm256_extractf128_ps( v1b, 1 ));       \
	_mm_storeu_ps( (float *) (p1+4),                         v1c );           \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a )   );       \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ));       \
	_mm_storeu_ps( (float *) (p2+2), _mm256_castps256_ps128( v2b )   );       \
	_mm_storeu_ps( (float *) (p2+3),  _mm256_extractf128_ps( v2b, 1 ));       \
	_mm_storeu_ps( (float *) (p2+4),                         v2c );           \
	_mm_storeu_ps( (float *) (p3+0), _mm256_castps256_ps128( v3a )   );       \
	_mm_storeu_ps( (float *) (p3+1),  _mm256_extractf128_ps( v3a, 1 ));       \
	_mm_storeu_ps( (float *) (p3+2), _mm256_castps256_ps128( v3b )   );       \
	_mm_storeu_ps( (float *) (p3+3),  _mm256_extractf128_ps( v3b, 1 ));       \
	_mm_storeu_ps( (float *) (p3+4),                         v3c );           \
	_mm_storeu_ps( (float *) (p4+0), _mm256_castps256_ps128( v4a )   );       \
	_mm_storeu_ps( (float *) (p4+1),  _mm256_extractf128_ps( v4a, 1 ));       \
	_mm_storeu_ps( (float *) (p4+2), _mm256_castps256_ps128( v4b )   );       \
	_mm_storeu_ps( (float *) (p4+3),  _mm256_extractf128_ps( v4b, 1 ));       \
	_mm_storeu_ps( (float *) (p4+4),                         v4c );           \
}
#else
#error Unsupported interpolation order
#endif

inline void DEP_CURRENT_2D
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_split_buf2D * const part)
{
        
  typedef struct Current { float j1, j2, j3; } t_current;
  
  __m256 vwp1[NP], vwp2[NP]; 
  __m256 vwl1[NP], vwl2[NP];
  
  __m256 vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
		  
  __m256 vqnx, vqny, vqvz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 
 
  unsigned int np;
  unsigned int i, k1, k2;

  __m256 const oneThird = _mm256_set1_ps( 1.0f/3.0f );

  unsigned int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;
       
       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;
       
       part->ix[k] = part->iy[k] = 1;
    }
  }
 
  // Each vp is 8 floats
  for( i = 0; i < np; i+=VEC_WIDTH ) {
    
    __m256i vix, viy;
    
    // load 8 particles
    LOAD8P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
    
    // Store cell index
    viadd_store( idx, vix, vimul_s( viy, Dy ) ) ;
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqvz = _mm256_mul_ps( vq, vvz);
    vqvz = _mm256_mul_ps( vqvz, oneThird );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();
    WL( vqny, vy0, vy1, (__m256 *) vwl2 ); vwl2[NP-1] = _mm256_setzero_ps();

    // get perpendicular weights
    for( k1 = 0; k1 < NP; k1++ ) {
      vwp1[k1] = _mm256_add_ps(vs0y[k1], vs1y[k1]);
      vwp2[k1] = _mm256_add_ps(vs0x[k1], vs1x[k1]);
    }

    __m256 a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];
     
    // get {j1,j2,j3} currents
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
  		
  		__m256 j1, j2, j3;
 
        j1 = _mm256_mul_ps( vwl1[k1] , vwp1[k2] );
        j2 = _mm256_mul_ps( vwl2[k2] , vwp2[k1] );
 		
  		__m256 const oneHalf  = _mm256_set1_ps( 0.5f );
  		__m256  tmp1, tmp2;
  		
		tmp1 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs1y[k2] ) );
	    
		tmp2 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs0y[k2] ) );

		j3 = _mm256_mul_ps( vqvz, _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf ) ) );
		
		// Do a 16x4 transpose, ignoring 4th current component
        __m256 t0 = _mm256_unpacklo_ps( j1, j3 );
        __m256 t1 = _mm256_unpackhi_ps( j1, j3 );
        __m256 t2 = _mm256_unpacklo_ps( j2, _mm256_setzero_ps() );  
        __m256 t3 = _mm256_unpackhi_ps( j2, _mm256_setzero_ps() );  
    
        a[k2][k1] = _mm256_unpacklo_ps( t0, t2 );
        b[k2][k1] = _mm256_unpackhi_ps( t0, t2 );
        c[k2][k1] = _mm256_unpacklo_ps( t1, t3 );
        d[k2][k1] = _mm256_unpackhi_ps( t1, t3 );
		
      }
    }

	ACC_CURRENT2D(  0, a, 0 )
	ACC_CURRENT2D(  1, b, 0 )
	ACC_CURRENT2D(  2, c, 0 )
	ACC_CURRENT2D(  3, d, 0 )
	ACC_CURRENT2D(  4, a, 1 )
	ACC_CURRENT2D(  5, b, 1 )
	ACC_CURRENT2D(  6, c, 1 )
	ACC_CURRENT2D(  7, d, 1 )
	  
  }
}

/********************************* 3D current deposition ********************************/

#if ( ORDER == 1 ) 

#define ACC_CURRENT3D( k, v, p ) {                                            \
    t_current *p0, *p1;                                                       \
    p0 = pj + idx[k] + k3*Dz;                                                 \
	p1 = p0 + Dy;                                                             \
	__m256 v0, v1, j0, j1;                                                    \
    j0 = _mm256_castps128_ps256(   _mm256_extractf128_ps( v[0][0], p ));      \
	j0 = _mm256_insertf128_ps( j0, _mm256_extractf128_ps( v[0][1], p ), 1);   \
	j1 = _mm256_castps128_ps256(   _mm256_extractf128_ps( v[1][0], p ));      \
	j1 = _mm256_insertf128_ps( j1, _mm256_extractf128_ps( v[1][1], p ), 1);   \
	v0 = _mm256_castps128_ps256(   _mm_loadu_ps( (float *) (p0+0) ) );        \
	v0 = _mm256_insertf128_ps( v0, _mm_loadu_ps( (float *) (p0+1) ),1);       \
	v1 = _mm256_castps128_ps256(   _mm_loadu_ps( (float *) (p1+0) ) );        \
	v1 = _mm256_insertf128_ps( v1, _mm_loadu_ps( (float *) (p1+1) ),1 );      \
	v0 = _mm256_add_ps( v0, j0 );                                             \
	v1 = _mm256_add_ps( v1, j1 );                                             \
	_mm_storeu_ps( (float *) p0,      _mm256_castps256_ps128( v0 ) );         \
	_mm_storeu_ps( (float *) (p0+1),   _mm256_extractf128_ps( v0, 1 ) );      \
	_mm_storeu_ps( (float *) p1,      _mm256_castps256_ps128( v1 ) );         \
	_mm_storeu_ps( (float *) (p1+1),   _mm256_extractf128_ps( v1, 1 ) );      \
}

#elif ( ORDER == 2 ) 

#define ACC_CURRENT3D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2;                                                  \
    p0 = pj + idx[k] + k3*Dz;                                                 \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	__m256 v0a, v1a, v2a;                                                     \
	__m128 v0b, v1b, v2b;                                                     \
	__m256 j0a, j1a, j2a;                                                     \
	__m128 j0b, j1b, j2b;                                                     \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b =                            _mm256_extractf128_ps( v[0][2], p );     \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
    j1b =                            _mm256_extractf128_ps( v[1][2], p );     \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
    j2b =                            _mm256_extractf128_ps( v[2][2], p );     \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ) , 1);   \
	v0b =                            _mm_loadu_ps( (float *) (p0+2) );        \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ) , 1);   \
	v1b =                            _mm_loadu_ps( (float *) (p1+2) );        \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ) , 1);   \
	v2b =                            _mm_loadu_ps( (float *) (p2+2) );        \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b =    _mm_add_ps( v0b, j0b );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b =    _mm_add_ps( v1b, j1b );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b =    _mm_add_ps( v2b, j2b );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a ) );         \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ) );      \
	_mm_storeu_ps( (float *) (p0+2),                         v0b );           \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a ) ) ;        \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ) );      \
	_mm_storeu_ps( (float *) (p1+2),                         v1b );           \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a ) );         \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ) );      \
	_mm_storeu_ps( (float *) (p2+2),                         v2b );           \
}

#elif ( ORDER == 3 ) 

#define ACC_CURRENT3D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2, *p3;                                             \
    p0 = pj + idx[k] + k3*Dz;                                                 \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	p3 = p0 + 3*Dy;                                                           \
	__m256 v0a, v1a, v2a, v3a;                                                \
	__m256 v0b, v1b, v2b, v3b;                                                \
	__m256 j0a, j1a, j2a, j3a;                                                \
	__m256 j0b, j1b, j2b, j3b;                                                \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][2], p ));    \
	j0b = _mm256_insertf128_ps( j0b, _mm256_extractf128_ps( v[0][3], p ), 1); \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
	j1b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][2], p ));    \
	j1b = _mm256_insertf128_ps( j1b, _mm256_extractf128_ps( v[1][3], p ), 1); \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
	j2b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][2], p ));    \
	j2b = _mm256_insertf128_ps( j2b, _mm256_extractf128_ps( v[2][3], p ), 1); \
	j3a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][0], p ));    \
	j3a = _mm256_insertf128_ps( j3a, _mm256_extractf128_ps( v[3][1], p ), 1); \
	j3b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][2], p ));    \
	j3b = _mm256_insertf128_ps( j3b, _mm256_extractf128_ps( v[3][3], p ), 1); \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ),1);     \
	v0b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+2) ) );      \
	v0b = _mm256_insertf128_ps( v0b, _mm_loadu_ps( (float *) (p0+3) ),1);     \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ),1);     \
	v1b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+2) ) );      \
	v1b = _mm256_insertf128_ps( v1b, _mm_loadu_ps( (float *) (p1+3) ),1);     \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ),1);     \
	v2b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+2) ) );      \
	v2b = _mm256_insertf128_ps( v2b, _mm_loadu_ps( (float *) (p2+3) ),1);     \
	v3a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+0) ) );      \
	v3a = _mm256_insertf128_ps( v3a, _mm_loadu_ps( (float *) (p3+1) ),1);     \
	v3b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+2) ) );      \
	v3b = _mm256_insertf128_ps( v3b, _mm_loadu_ps( (float *) (p3+3) ),1);     \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b = _mm256_add_ps( v0b, j0b );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b = _mm256_add_ps( v1b, j1b );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b = _mm256_add_ps( v2b, j2b );                                          \
	v3a = _mm256_add_ps( v3a, j3a );                                          \
	v3b = _mm256_add_ps( v3b, j3b );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a )   );       \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ));       \
	_mm_storeu_ps( (float *) (p0+2), _mm256_castps256_ps128( v0b )   );       \
	_mm_storeu_ps( (float *) (p0+3),  _mm256_extractf128_ps( v0b, 1 ));       \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a )   );       \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ));       \
	_mm_storeu_ps( (float *) (p1+2), _mm256_castps256_ps128( v1b )   );       \
	_mm_storeu_ps( (float *) (p1+3),  _mm256_extractf128_ps( v1b, 1 ));       \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a )   );       \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ));       \
	_mm_storeu_ps( (float *) (p2+2), _mm256_castps256_ps128( v2b )   );       \
	_mm_storeu_ps( (float *) (p2+3),  _mm256_extractf128_ps( v2b, 1 ));       \
	_mm_storeu_ps( (float *) (p3+0), _mm256_castps256_ps128( v3a )   );       \
	_mm_storeu_ps( (float *) (p3+1),  _mm256_extractf128_ps( v3a, 1 ));       \
	_mm_storeu_ps( (float *) (p3+2), _mm256_castps256_ps128( v3b )   );       \
	_mm_storeu_ps( (float *) (p3+3),  _mm256_extractf128_ps( v3b, 1 ));       \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v, p ) {                                            \
    t_current *p0, *p1, *p2, *p3, *p4;                                        \
    p0 = pj + idx[k] + k3*Dz;                                                 \
	p1 = p0 + Dy;                                                             \
	p2 = p0 + 2*Dy;                                                           \
	p3 = p0 + 3*Dy;                                                           \
	p4 = p0 + 4*Dy;                                                           \
	__m256 v0a, v1a, v2a, v3a, v4a;                                           \
	__m256 v0b, v1b, v2b, v3b, v4b;                                           \
	__m128 v0c, v1c, v2c, v3c, v4c;                                           \
	__m256 j0a, j1a, j2a, j3a, j4a;                                           \
	__m256 j0b, j1b, j2b, j3b, j4b;                                           \
	__m128 j0c, j1c, j2c, j3c, j4c;                                           \
    j0a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][0], p ));    \
	j0a = _mm256_insertf128_ps( j0a, _mm256_extractf128_ps( v[0][1], p ), 1); \
    j0b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[0][2], p ));    \
	j0b = _mm256_insertf128_ps( j0b, _mm256_extractf128_ps( v[0][3], p ), 1); \
    j0c =                            _mm256_extractf128_ps( v[0][4], p );     \
	j1a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][0], p ));    \
	j1a = _mm256_insertf128_ps( j1a, _mm256_extractf128_ps( v[1][1], p ), 1); \
	j1b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[1][2], p ));    \
	j1b = _mm256_insertf128_ps( j1b, _mm256_extractf128_ps( v[1][3], p ), 1); \
    j1c =                            _mm256_extractf128_ps( v[1][4], p );     \
	j2a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][0], p ));    \
	j2a = _mm256_insertf128_ps( j2a, _mm256_extractf128_ps( v[2][1], p ), 1); \
	j2b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[2][2], p ));    \
	j2b = _mm256_insertf128_ps( j2b, _mm256_extractf128_ps( v[2][3], p ), 1); \
    j2c =                            _mm256_extractf128_ps( v[2][4], p );     \
	j3a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][0], p ));    \
	j3a = _mm256_insertf128_ps( j3a, _mm256_extractf128_ps( v[3][1], p ), 1); \
	j3b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[3][2], p ));    \
	j3b = _mm256_insertf128_ps( j3b, _mm256_extractf128_ps( v[3][3], p ), 1); \
    j3c =                            _mm256_extractf128_ps( v[3][4], p );     \
	j4a = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[4][0], p ));    \
	j4a = _mm256_insertf128_ps( j4a, _mm256_extractf128_ps( v[4][1], p ), 1); \
	j4b = _mm256_castps128_ps256(    _mm256_extractf128_ps( v[4][2], p ));    \
	j4b = _mm256_insertf128_ps( j4b, _mm256_extractf128_ps( v[4][3], p ), 1); \
    j4c =                            _mm256_extractf128_ps( v[4][4], p );     \
	v0a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+0) ) );      \
	v0a = _mm256_insertf128_ps( v0a, _mm_loadu_ps( (float *) (p0+1) ),1);     \
	v0b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p0+2) ) );      \
	v0b = _mm256_insertf128_ps( v0b, _mm_loadu_ps( (float *) (p0+3) ),1);     \
	v0c =                            _mm_loadu_ps( (float *) (p0+4) );        \
	v1a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+0) ) );      \
	v1a = _mm256_insertf128_ps( v1a, _mm_loadu_ps( (float *) (p1+1) ),1);     \
	v1b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p1+2) ) );      \
	v1b = _mm256_insertf128_ps( v1b, _mm_loadu_ps( (float *) (p1+3) ),1);     \
	v1c =                            _mm_loadu_ps( (float *) (p1+4) );        \
	v2a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+0) ) );      \
	v2a = _mm256_insertf128_ps( v2a, _mm_loadu_ps( (float *) (p2+1) ),1);     \
	v2b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p2+2) ) );      \
	v2b = _mm256_insertf128_ps( v2b, _mm_loadu_ps( (float *) (p2+3) ),1);     \
	v2c =                            _mm_loadu_ps( (float *) (p2+4) );        \
	v3a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+0) ) );      \
	v3a = _mm256_insertf128_ps( v3a, _mm_loadu_ps( (float *) (p3+1) ),1);     \
	v3b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p3+2) ) );      \
	v3b = _mm256_insertf128_ps( v3b, _mm_loadu_ps( (float *) (p3+3) ),1);     \
	v3c =                            _mm_loadu_ps( (float *) (p3+4) );        \
	v4a = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p4+0) ) );      \
	v4a = _mm256_insertf128_ps( v4a, _mm_loadu_ps( (float *) (p4+1) ),1);     \
	v4b = _mm256_castps128_ps256(    _mm_loadu_ps( (float *) (p4+2) ) );      \
	v4b = _mm256_insertf128_ps( v4b, _mm_loadu_ps( (float *) (p4+3) ),1);     \
	v4c =                            _mm_loadu_ps( (float *) (p4+4) );        \
	v0a = _mm256_add_ps( v0a, j0a );                                          \
	v0b = _mm256_add_ps( v0b, j0b );                                          \
	v0c =    _mm_add_ps( v0c, j0c );                                          \
	v1a = _mm256_add_ps( v1a, j1a );                                          \
	v1b = _mm256_add_ps( v1b, j1b );                                          \
	v1c =    _mm_add_ps( v1c, j1c );                                          \
	v2a = _mm256_add_ps( v2a, j2a );                                          \
	v2b = _mm256_add_ps( v2b, j2b );                                          \
	v2c =    _mm_add_ps( v2c, j2c );                                          \
	v3a = _mm256_add_ps( v3a, j3a );                                          \
	v3b = _mm256_add_ps( v3b, j3b );                                          \
	v3c =    _mm_add_ps( v3c, j3c );                                          \
	v4a = _mm256_add_ps( v4a, j4a );                                          \
	v4b = _mm256_add_ps( v4b, j4b );                                          \
	v4c =    _mm_add_ps( v4c, j4c );                                          \
	_mm_storeu_ps( (float *) (p0+0), _mm256_castps256_ps128( v0a )   );       \
	_mm_storeu_ps( (float *) (p0+1),  _mm256_extractf128_ps( v0a, 1 ));       \
	_mm_storeu_ps( (float *) (p0+2), _mm256_castps256_ps128( v0b )   );       \
	_mm_storeu_ps( (float *) (p0+3),  _mm256_extractf128_ps( v0b, 1 ));       \
	_mm_storeu_ps( (float *) (p0+4),                         v0c );           \
	_mm_storeu_ps( (float *) (p1+0), _mm256_castps256_ps128( v1a )   );       \
	_mm_storeu_ps( (float *) (p1+1),  _mm256_extractf128_ps( v1a, 1 ));       \
	_mm_storeu_ps( (float *) (p1+2), _mm256_castps256_ps128( v1b )   );       \
	_mm_storeu_ps( (float *) (p1+3),  _mm256_extractf128_ps( v1b, 1 ));       \
	_mm_storeu_ps( (float *) (p1+4),                         v1c );           \
	_mm_storeu_ps( (float *) (p2+0), _mm256_castps256_ps128( v2a )   );       \
	_mm_storeu_ps( (float *) (p2+1),  _mm256_extractf128_ps( v2a, 1 ));       \
	_mm_storeu_ps( (float *) (p2+2), _mm256_castps256_ps128( v2b )   );       \
	_mm_storeu_ps( (float *) (p2+3),  _mm256_extractf128_ps( v2b, 1 ));       \
	_mm_storeu_ps( (float *) (p2+4),                         v2c );           \
	_mm_storeu_ps( (float *) (p3+0), _mm256_castps256_ps128( v3a )   );       \
	_mm_storeu_ps( (float *) (p3+1),  _mm256_extractf128_ps( v3a, 1 ));       \
	_mm_storeu_ps( (float *) (p3+2), _mm256_castps256_ps128( v3b )   );       \
	_mm_storeu_ps( (float *) (p3+3),  _mm256_extractf128_ps( v3b, 1 ));       \
	_mm_storeu_ps( (float *) (p3+4),                         v3c );           \
	_mm_storeu_ps( (float *) (p4+0), _mm256_castps256_ps128( v4a )   );       \
	_mm_storeu_ps( (float *) (p4+1),  _mm256_extractf128_ps( v4a, 1 ));       \
	_mm_storeu_ps( (float *) (p4+2), _mm256_castps256_ps128( v4b )   );       \
	_mm_storeu_ps( (float *) (p4+3),  _mm256_extractf128_ps( v4b, 1 ));       \
	_mm_storeu_ps( (float *) (p4+4),                         v4c );           \
}

#else
#error Unsupported interpolation order
#endif



inline void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part)
{
        
  typedef struct Current { float j1, j2, j3; } t_current;
      
  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
  
  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  __m256 vwl1[NP], vwl2[NP], vwl3[NP];
  __m256 vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  unsigned int np;
  unsigned int i, k1, k2, k3;
     
  unsigned int const Dy = size[0];    
  unsigned int const Dz = Dy * size[1]; 
    
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 + 
                                                  ( offset[1] - OFFSET  ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET  ) * 3 * Dz );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);
  
  np = part -> np;
    
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;
       
       part -> x0[k] = part -> x1[k] = 0.;
       part -> y0[k] = part -> y1[k] = 0.;
       part -> z0[k] = part -> z1[k] = 0.;
       part -> q[k] = 0.;
       part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( i = 0; i < np; i+=VEC_WIDTH ) {

    __m256i vix, viy, viz;
      
	// load 8 particles
    LOAD8P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
    
    // Do idx[k]= ix[k] + iy[k]*Dy + iz[k]*Dz vectorially 
    {
      __m128i vDy, vDz, a, b;
      vDy = _mm_set1_epi32( Dy );
      vDz = _mm_set1_epi32( Dz );
      
      a = _mm_mullo_epi32( _mm256_castsi256_si128(viz), vDz );
      b = _mm_mullo_epi32( _mm256_extractf128_si256(viz,1), vDz );

      a = _mm_add_epi32( a, _mm_mullo_epi32( _mm256_castsi256_si128(viy), vDy ));
      b = _mm_add_epi32( b, _mm_mullo_epi32( _mm256_extractf128_si256(viy,1), vDy ));
      
      a = _mm_add_epi32( a, _mm256_castsi256_si128(vix) );
      b = _mm_add_epi32( b, _mm256_extractf128_si256(vix,1) );
      
      _mm_store_si128((__m128i *)  &idx[0], a );
      _mm_store_si128((__m128i *)  &idx[4], b );
    }
    
    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqnz = _mm256_mul_ps( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();
    WL( vqny, vy0, vy1, (__m256 *) vwl2 ); vwl2[NP-1] = _mm256_setzero_ps();
    WL( vqnz, vz0, vz1, (__m256 *) vwl3 ); vwl3[NP-1] = _mm256_setzero_ps();
    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
         
         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         const __m256 oneHalf = _mm256_set1_ps(0.5f);
         __m256 tmp1, tmp2;        
        
         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs0z[k2]), _mm256_mul_ps( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs1z[k2]), _mm256_mul_ps( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1] = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0z[k2]), _mm256_mul_ps( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1z[k2]), _mm256_mul_ps( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1] = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2]), _mm256_mul_ps( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2]), _mm256_mul_ps( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1] = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

      }
    }

    // Accumulate current 1 plane at a time using the 2D algorithm
    for( k3 = 0; k3 < NP; k3++ ) {
      __m256 a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];
      
	  for( k2 = 0; k2 < NP; k2++ ) {
		
		for ( k1 = 0; k1 < NP; k1++ ) {
		   __m256 j1, j2, j3;
		 
		   // Calculate current components
		   j1 = _mm256_mul_ps( vwl1[k1] , vwp1[k3][k2] );
		   j2 = _mm256_mul_ps( vwl2[k2] , vwp2[k3][k1] );
		   j3 = _mm256_mul_ps( vwl3[k3] , vwp3[k2][k1] );
    
		   // Do a 8x4 transpose, ignoring 4th current component
		   __m256 t0 = _mm256_unpacklo_ps( j1, j3 );
		   __m256 t1 = _mm256_unpackhi_ps( j1, j3 );
		   __m256 t2 = _mm256_unpacklo_ps( j2, _mm256_setzero_ps() );  
		   __m256 t3 = _mm256_unpackhi_ps( j2, _mm256_setzero_ps() );  
	
		   a[k2][k1] = _mm256_unpacklo_ps( t0, t2 );
		   b[k2][k1] = _mm256_unpackhi_ps( t0, t2 );
		   c[k2][k1] = _mm256_unpacklo_ps( t1, t3 );
		   d[k2][k1] = _mm256_unpackhi_ps( t1, t3 );
        }
      }

	  // Accumulate electric current in k3 plane
	  ACC_CURRENT3D(  0, a, 0 )
	  ACC_CURRENT3D(  1, b, 0 )
	  ACC_CURRENT3D(  2, c, 0 )
	  ACC_CURRENT3D(  3, d, 0 )
	  ACC_CURRENT3D(  4, a, 1 )
	  ACC_CURRENT3D(  5, b, 1 )
	  ACC_CURRENT3D(  6, c, 1 )
	  ACC_CURRENT3D(  7, d, 1 )

    }
     

  }
}

#endif

#undef ACC_CURRENT2D
#undef ACC_CURRENT3D

#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D
#undef SPLINE
#undef WL

#undef OFFSET
#undef ORDER
#undef NP

#endif
