/*****************************************************************************************

Current deposition, AVX optimized version (double precision)
 
*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avxd.h"

#include <stdio.h>

#include "fortran.h"

#include "vector-avx.h"
#include "splines-avx.h"


/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] ) {

  vwl[0] = _mm256_mul_pd( vqn, _mm256_sub_pd( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition
***********************************************************************/

inline void vwl_s2( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] ) {

  __m256d const c1_2 = _mm256_set1_pd( 0.5 );
  
  __m256d d    = _mm256_sub_pd( vx1, vx0 );
  __m256d s1_2 = _mm256_sub_pd( c1_2, vx0 );
  __m256d p1_2 = _mm256_add_pd( c1_2, vx0 );
  
  __m256d n = _mm256_mul_pd( vqn, d );
  
  // qn * d * ( s1_2 - 0.5 * d ) 
  vwl[0] = _mm256_mul_pd( n, _mm256_sub_pd( s1_2, _mm256_mul_pd( c1_2, d ) ) );
  // qn * d * ( p1_2 + 0.5 * d ) 
  vwl[1] = _mm256_mul_pd( n, _mm256_add_pd( p1_2, _mm256_mul_pd( c1_2, d ) ) );  

}

/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition
***********************************************************************/

inline void vwl_s3( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] ) {

   __m256d const c1_2 = _mm256_set1_pd( 0.5 );
   __m256d const c3_4 = _mm256_set1_pd( 0.75 );
   __m256d const c1_3 = _mm256_set1_pd( 1.0/3.0 );
   
   __m256d d = _mm256_sub_pd( vx1, vx0 );
   __m256d s = _mm256_sub_pd( c1_2, vx0 );
   __m256d p = _mm256_add_pd( c1_2, vx0 );
   __m256d d_3 = _mm256_mul_pd( c1_3, d );

   vwl[0] = _mm256_mul_pd( c1_2, _mm256_sub_pd( _mm256_mul_pd( s, s ),       _mm256_mul_pd( d, _mm256_sub_pd( s, d_3 ) ) ) );
   vwl[1] = _mm256_sub_pd( _mm256_sub_pd( c3_4, _mm256_mul_pd( vx0, vx0 ) ), _mm256_mul_pd( d, _mm256_add_pd( vx0, d_3 ) ) );
   vwl[2] = _mm256_mul_pd( c1_2, _mm256_add_pd( _mm256_mul_pd( p, p ),       _mm256_mul_pd( d, _mm256_add_pd( p, d_3 ) ) ) );

   __m256d n = _mm256_mul_pd( vqn, d );
   vwl[0] = _mm256_mul_pd( n, vwl[0] );
   vwl[1] = _mm256_mul_pd( n, vwl[1] );
   vwl[2] = _mm256_mul_pd( n, vwl[2] );

}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition
***********************************************************************/
inline void vwl_s4( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] ) {

   __m256d const c1_2 = _mm256_set1_pd( 0.5 );
   __m256d const c1_4 = _mm256_set1_pd( 0.25 );
   __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
   __m256d const c3_2 = _mm256_set1_pd( 1.5 );
   __m256d const c3_4 = _mm256_set1_pd( 0.75 );
   __m256d const c1_3 = _mm256_set1_pd( 1.0/3.0 );
   __m256d const c2_3 = _mm256_set1_pd( 2.0/3.0 );
   
   __m256d d = _mm256_sub_pd( vx1, vx0 );
   __m256d s = _mm256_sub_pd( c1_2, vx0 );
   __m256d p = _mm256_add_pd( c1_2, vx0 );

   __m256d t = _mm256_sub_pd( vx0, c1_6 );
   __m256d u = _mm256_add_pd( vx0, c1_6 );

   __m256d s2 = _mm256_mul_pd( s, s );
   __m256d s3 = _mm256_mul_pd( s2, s );

   __m256d p2 = _mm256_mul_pd( p, p );
   __m256d p3 = _mm256_mul_pd( p2, p );

   __m256d d_2 = _mm256_mul_pd( d, c1_2 );
   __m256d d_4 = _mm256_mul_pd( d, c1_4 );
   

   vwl[0] = _mm256_mul_pd( c1_6, _mm256_add_pd( s3, _mm256_mul_pd( d, _mm256_sub_pd( _mm256_mul_pd( d, _mm256_sub_pd( s, d_4 ) ),  _mm256_mul_pd( c3_2, s2 ) ) ) ) );
   vwl[1] = _mm256_add_pd( _mm256_add_pd( _mm256_sub_pd( c2_3, p2 ), _mm256_mul_pd( c1_2, p3 ) ), _mm256_mul_pd( d, _mm256_add_pd( _mm256_sub_pd( _mm256_mul_pd( c3_4, _mm256_mul_pd( t, t ) ), c1_3 ), _mm256_mul_pd( d_2, _mm256_add_pd( t, d_4 ) ) ) ) );
   vwl[2] = _mm256_sub_pd( _mm256_add_pd( _mm256_sub_pd( c2_3, s2 ), _mm256_mul_pd( c1_2, s3 ) ), _mm256_mul_pd( d, _mm256_add_pd( _mm256_sub_pd( _mm256_mul_pd( c3_4, _mm256_mul_pd( u, u ) ), c1_3 ), _mm256_mul_pd( d_2, _mm256_add_pd( u, d_4 ) ) ) ) );
   vwl[3] = _mm256_mul_pd( c1_6, _mm256_add_pd( p3, _mm256_mul_pd( d, _mm256_add_pd( _mm256_mul_pd( d, _mm256_add_pd( p, d_4 ) ),  _mm256_mul_pd( c3_2, p2 ) ) ) ) );

   __m256d n = _mm256_mul_pd( vqn, d );
   vwl[0] = _mm256_mul_pd( n, vwl[0] );
   vwl[1] = _mm256_mul_pd( n, vwl[1] );
   vwl[2] = _mm256_mul_pd( n, vwl[2] );
   vwl[3] = _mm256_mul_pd( n, vwl[3] );
      
}



/****************************************************************************************

  Generate specific current deposition functions for 2D and 3D, 1st to 4th order

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

#define SPLINE  ONAME( vsplined, ORDER )
#define WL      ONAME( vwl, ORDER )

inline void DEP_CURRENT_1D
(double * const current, int const * const size, int const * const offset, 
 double * const norm, t_split_buf1D * const part) {
        
  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *p0; 
    
  dvec j2[NP], j3[NP];
  
  dvec vwl1[ORDER];
  
  __m256d vx0, vx1, vq, vvy, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

      
  __m256d vqnx, vqvz, vqvy;
  __m256d vs0x[NP], vs1x[NP]; 

  unsigned int np;
  unsigned int i, k, k1;

  __m256d const c1_2  = _mm256_set1_pd( 0.5 );

  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       k = np + i;
       
       part->x0[k] = part->x1[k] = 0.0;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0;
       
       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    __m128i vix;

    // load 4 particles
    LOAD4P1D( part, i, vx0, vx1, vq, vvy, vvz, vix );
    
    // Store cell index
    _mm_store_si128( (__m128i *) idx, vix ); 
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    
    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqvz = _mm256_mul_pd( vq, vvz );
    vqvy = _mm256_mul_pd( vq, vvy );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );
     
    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m256d vwp1 = _mm256_mul_pd( c1_2, _mm256_add_pd(vs0x[k1], vs1x[k1]) );
      j2[k1].v4 = _mm256_mul_pd( vqvy, vwp1 );
      j3[k1].v4 = _mm256_mul_pd( vqvz, vwp1 );
    }

    // New version loop by particle on the outside loop
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


inline void DEP_CURRENT_2D
(double * const current, int const * const size, int const * const offset, 
 double * const norm, t_split_buf2D * const part)
{
        
  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  dvec j3[NP][NP];
  

  dvec vwp1[NP], vwp2[NP]; 
  dvec vwl1[ORDER], vwl2[ORDER];
  
  __m256d vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

		  
  __m256d vqnx, vqny, vqvz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  unsigned int np;
  unsigned int i, k, k1, k2;

  __m256d const oneThird = _mm256_set1_pd( 1.0/3.0 );

  const unsigned int Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );


  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);
  
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

  for( i = 0; i < np; i+=VEC_WIDTH ) {
    
     __m128i vix, viy;
    
    // load 4 particles
    LOAD4P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
    
    // Store cell index
    _mm_store_si128( (__m128i *) idx , 
                     _mm_add_epi32( vix, _mm_mullo_epi32( viy, _mm_set1_epi32(Dy) ) ) ); 
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm256_mul_pd( vq, vnorm1);
    vqny = _mm256_mul_pd( vq, vnorm2);
    vqvz = _mm256_mul_pd( vq, vvz);
    vqvz = _mm256_mul_pd( vqvz, oneThird );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );
    WL( vqny, vy0, vy1, (__m256d *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v4 = _mm256_add_pd(vs0y[k], vs1y[k]);
      vwp2[k].v4 = _mm256_add_pd(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        __m256d const oneHalf  = _mm256_set1_pd( 0.5 );
        __m256d tmp1, tmp2;
  		
		tmp1 = _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs0y[k2] ), _mm256_mul_pd( vs1x[k1], vs1y[k2] ) );
		tmp2 = _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs1y[k2] ), _mm256_mul_pd( vs1x[k1], vs0y[k2] ) );

		j3[k2][k1].v4 = _mm256_mul_pd( vqvz, _mm256_add_pd( tmp1, _mm256_mul_pd( tmp2, oneHalf ) ) );
      }
    }

    // New version loop by particle on the outside loop
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
		   p0[k1].j3 += (j3[k2][k1]).v[k];
		}
	  }
	  
    }
	  
  }
}


inline void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{
        
  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  dvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  dvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  __m256d vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
  
  __m256d vqnx, vqny, vqnz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  
  unsigned int np;
  unsigned int i, k, k1, k2, k3;
  
  unsigned int const Dy = size[0];    
  unsigned int const Dz = Dy * size[1]; 
  
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                                  ( offset[1] - OFFSET ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET ) * 3 * Dz );

  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);
  __m256d const vnorm3 = _mm256_set1_pd(norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
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

  for( i = 0; i < np; i += VEC_WIDTH ) {

     __m128i vix, viy, viz;
      
	// load 4 particles
    LOAD4P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
         
    // Do idx[k] = ix[k] + iy[k]*Dy + iz[k]*Dz vectorially 
    _mm_store_si128( (__m128i *)  &idx[0], _mm_add_epi32(  
                       _mm_add_epi32( _mm_mullo_epi32( viz, _mm_set1_epi32( Dz ) ), 
                                      _mm_mullo_epi32( viy, _mm_set1_epi32( Dy ) ) ), 
                       vix ) );
   
    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm256_mul_pd( vq, vnorm1);
    vqny = _mm256_mul_pd( vq, vnorm2);
    vqnz = _mm256_mul_pd( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );
    WL( vqny, vy0, vy1, (__m256d *) vwl2 );
    WL( vqnz, vz0, vz1, (__m256d *) vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
		 __m256d tmp1, tmp2;        
		 const __m256d oneHalf = _mm256_set1_pd(0.5);
         
         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 =  _mm256_add_pd( _mm256_mul_pd( vs0y[k1], vs0z[k2]), _mm256_mul_pd( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm256_add_pd( _mm256_mul_pd( vs0y[k1], vs1z[k2]), _mm256_mul_pd( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v4 = _mm256_add_pd( tmp1, _mm256_mul_pd( tmp2, oneHalf));

         tmp1 =  _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs0z[k2]), _mm256_mul_pd( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs1z[k2]), _mm256_mul_pd( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v4 = _mm256_add_pd( tmp1, _mm256_mul_pd( tmp2, oneHalf));

         tmp1 =  _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs0y[k2]), _mm256_mul_pd( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm256_add_pd( _mm256_mul_pd( vs0x[k1], vs1y[k2]), _mm256_mul_pd( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v4 = _mm256_add_pd( tmp1, _mm256_mul_pd( tmp2, oneHalf));

      }
    }
    
    // loop by particle on the outside loop
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

#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D
#undef SPLINE
#undef WL

#undef OFFSET
#undef ORDER
#undef NP

#endif
