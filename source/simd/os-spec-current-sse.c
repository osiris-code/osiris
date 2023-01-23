/*****************************************************************************************

Current deposition, SSE optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-sse.h"

#include <stdio.h>

#include "fortran.h"

#include "vector-sse.h"
#include "splines-sse.h"


/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

  vwl[0] = _mm_mul_ps( vqn, _mm_sub_ps( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition
***********************************************************************/

inline void vwl_s2( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

  __m128 const c1_2 = _mm_set1_ps( 0.5f );
  
  __m128 d    = _mm_sub_ps( vx1, vx0 );
  __m128 s1_2 = _mm_sub_ps( c1_2, vx0 );
  __m128 p1_2 = _mm_add_ps( c1_2, vx0 );
  
  __m128 n = _mm_mul_ps( vqn, d );
  
  // qn * d * ( s1_2 - 0.5 * d ) 
  vwl[0] = _mm_mul_ps( n, _mm_sub_ps( s1_2, _mm_mul_ps( c1_2, d ) ) );
  // qn * d * ( p1_2 + 0.5 * d ) 
  vwl[1] = _mm_mul_ps( n, _mm_add_ps( p1_2, _mm_mul_ps( c1_2, d ) ) );  

}



/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition
***********************************************************************/

inline void vwl_s3( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {
   
   __m128 const c1_2 = _mm_set1_ps( 0.5f );
   __m128 const c3_4 = _mm_set1_ps( 0.75f );
   __m128 const c1_3 = _mm_set1_ps( 1.0f/3.0f );
   
   __m128 d = _mm_sub_ps( vx1, vx0 );
   __m128 s = _mm_sub_ps( c1_2, vx0 );
   __m128 p = _mm_add_ps( c1_2, vx0 );
   __m128 d_3 = _mm_mul_ps( c1_3, d );

   vwl[0] = _mm_mul_ps( c1_2, _mm_sub_ps( _mm_mul_ps( s, s ),       _mm_mul_ps( d, _mm_sub_ps( s, d_3 ) ) ) );
   vwl[1] = _mm_sub_ps( _mm_sub_ps( c3_4, _mm_mul_ps( vx0, vx0 ) ), _mm_mul_ps( d, _mm_add_ps( vx0, d_3 ) ) );
   vwl[2] = _mm_mul_ps( c1_2, _mm_add_ps( _mm_mul_ps( p, p ),       _mm_mul_ps( d, _mm_add_ps( p, d_3 ) ) ) );

   __m128 n = _mm_mul_ps( vqn, d );
   vwl[0] = _mm_mul_ps( n, vwl[0] );
   vwl[1] = _mm_mul_ps( n, vwl[1] );
   vwl[2] = _mm_mul_ps( n, vwl[2] );


}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition
***********************************************************************/
inline void vwl_s4( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

   __m128 const c1_2 = _mm_set1_ps( 0.5f );
   __m128 const c1_4 = _mm_set1_ps( 0.25f );
   __m128 const c1_6 = _mm_set1_ps( 1.0f/6.0f );
   __m128 const c3_2 = _mm_set1_ps( 1.5f );
   __m128 const c3_4 = _mm_set1_ps( 0.75f );
   __m128 const c1_3 = _mm_set1_ps( 1.0f/3.0f );
   __m128 const c2_3 = _mm_set1_ps( 2.0f/3.0f );
   
   __m128 d = _mm_sub_ps( vx1, vx0 );
   __m128 s = _mm_sub_ps( c1_2, vx0 );
   __m128 p = _mm_add_ps( c1_2, vx0 );

   __m128 t = _mm_sub_ps( vx0, c1_6 );
   __m128 u = _mm_add_ps( vx0, c1_6 );

   __m128 s2 = _mm_mul_ps( s, s );
   __m128 s3 = _mm_mul_ps( s2, s );

   __m128 p2 = _mm_mul_ps( p, p );
   __m128 p3 = _mm_mul_ps( p2, p );

   __m128 d_2 = _mm_mul_ps( d, c1_2 );
   __m128 d_4 = _mm_mul_ps( d, c1_4 );
   

   vwl[0] = _mm_mul_ps( c1_6, _mm_add_ps( s3, _mm_mul_ps( d, _mm_sub_ps( _mm_mul_ps( d, _mm_sub_ps( s, d_4 ) ),  _mm_mul_ps( c3_2, s2 ) ) ) ) );
   vwl[1] = _mm_add_ps( _mm_add_ps( _mm_sub_ps( c2_3, p2 ), _mm_mul_ps( c1_2, p3 ) ), _mm_mul_ps( d, _mm_add_ps( _mm_sub_ps( _mm_mul_ps( c3_4, _mm_mul_ps( t, t ) ), c1_3 ), _mm_mul_ps( d_2, _mm_add_ps( t, d_4 ) ) ) ) );
   vwl[2] = _mm_sub_ps( _mm_add_ps( _mm_sub_ps( c2_3, s2 ), _mm_mul_ps( c1_2, s3 ) ), _mm_mul_ps( d, _mm_add_ps( _mm_sub_ps( _mm_mul_ps( c3_4, _mm_mul_ps( u, u ) ), c1_3 ), _mm_mul_ps( d_2, _mm_add_ps( u, d_4 ) ) ) ) );
   vwl[3] = _mm_mul_ps( c1_6, _mm_add_ps( p3, _mm_mul_ps( d, _mm_add_ps( _mm_mul_ps( d, _mm_add_ps( p, d_4 ) ),  _mm_mul_ps( c3_2, p2 ) ) ) ) );

   __m128 n = _mm_mul_ps( vqn, d );
   vwl[0] = _mm_mul_ps( n, vwl[0] );
   vwl[1] = _mm_mul_ps( n, vwl[1] );
   vwl[2] = _mm_mul_ps( n, vwl[2] );
   vwl[3] = _mm_mul_ps( n, vwl[3] );
      
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

#define SPLINE  ONAME( vspline, ORDER )
#define WL      ONAME( vwl, ORDER )


inline void DEP_CURRENT_1D
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_split_buf1D * const part) {
        
  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *p0; 
    
  fvec j2[NP], j3[NP];
  
  fvec vwl1[ORDER];
  
  __m128 vx0, vx1, vq, vvy, vvz;
  
  DECLARE_ALIGNED_16( int idx[VEC_WIDTH] );

		  
  __m128 vqnx, vqvz, vqvy;
  __m128 vs0x[NP], vs1x[NP]; 

  unsigned int np;
  unsigned int i, k, k1;

  __m128 const c1_2  = _mm_set1_ps( 0.5f );

  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  __m128 const vnorm1 = _mm_set1_ps(norm[0]);
  
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
    _mm_store_si128( (void *) idx, vix ); 
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    
    vqnx = _mm_mul_ps( vq, vnorm1 );
    vqvy = _mm_mul_ps( vq, vvy );
    vqvz = _mm_mul_ps( vq, vvz );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128 *) vwl1 );
     
    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m128 vwp1 = _mm_mul_ps( c1_2, _mm_add_ps(vs0x[k1], vs1x[k1]) );
      j2[k1].v4 = _mm_mul_ps( vqvy, vwp1 );
      j3[k1].v4 = _mm_mul_ps( vqvz, vwp1 );
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
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_split_buf2D * const part)
{
        
  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
    
  fvec j3[NP][NP];
  
  fvec vwp1[NP], vwp2[NP]; 
  fvec vwl1[ORDER], vwl2[ORDER];
  
  __m128 vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_16( int idx[VEC_WIDTH] );

      
  __m128 vqnx, vqny, vqvz;
  __m128 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  unsigned int np;
  unsigned int i, k, k1, k2;

  __m128 const c1_3 = _mm_set_ps1( 1.0f/3.0f );

  int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  __m128 const vnorm1 = _mm_set_ps1(norm[0]);
  __m128 const vnorm2 = _mm_set_ps1(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       k = np + i;
       
       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;
       
       part->ix[k] = part->iy[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    __m128i vix, viy;

    // load 4 particles
    LOAD4P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
    
    // Store cell index
    _mm_store_si128( (void *) idx, _mm_add_epi32( vix, vsmul( viy, Dy ) ) ); 
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm_mul_ps( vq, vnorm1);
    vqny = _mm_mul_ps( vq, vnorm2);
    vqvz = _mm_mul_ps( vq, vvz);
    vqvz = _mm_mul_ps( vqvz, c1_3 );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128 *) vwl1 );
    WL( vqny, vy0, vy1, (__m128 *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v4 = _mm_add_ps(vs0y[k], vs1y[k]);
      vwp2[k].v4 = _mm_add_ps(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
    __m128 const c1_2  = _mm_set_ps1( 0.5f );
    __m128 tmp1, tmp2;

    tmp1 = _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0y[k2] ), 
                       _mm_mul_ps( vs1x[k1], vs1y[k2] ) );
     
    tmp2 = _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1y[k2] ), 
           _mm_mul_ps( vs1x[k1], vs0y[k2] ) );
    
    j3[k1][k2].v4 = _mm_mul_ps( vqvz, _mm_add_ps( tmp1, _mm_mul_ps( tmp2, c1_2 ) ) );
      }
    }

    // New version loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      pjpart = pj + idx[k];

    // accumulate j1
    for( k2 = 0; k2 < NP; k2++ ) {
    p0 = pjpart + k2*Dy;
    for ( k1 = 0; k1 < ORDER; k1++ ) {
      p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
    }
    }
  
    // accumulate j2 - making k2 the outside loop gives marginal perf. gain
    for( k2 = 0; k2 < ORDER; k2++ ) {
    p0 = pjpart + k2*Dy;
    for ( k1 = 0; k1 < NP; k1++ ) {
       p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
    }
    }
  
    // accumulate j3
    for( k2 = 0; k2 < NP; k2++ ) {
    p0 = pjpart + k2*Dy;
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
    
  fvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  fvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  __m128 vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_16( int idx[4] );
  
  __m128 vqnx, vqny, vqnz;
  __m128 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  
  unsigned int np;
  unsigned int i, k, k1, k2, k3;
  
  int const Dy = size[0];    
  int const Dz = Dy * size[1]; 
  
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                                  ( offset[1] - OFFSET ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET ) * 3 * Dz );

  __m128 const vnorm1 = _mm_set_ps1(norm[0]);
  __m128 const vnorm2 = _mm_set_ps1(norm[1]);
  __m128 const vnorm3 = _mm_set_ps1(norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
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


  for( i = 0; i < np; i+= VEC_WIDTH ) {
      
    __m128i vix, viy, viz;
      
	// load 4 particles
    LOAD4P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
    
    // Store cell index
    _mm_store_si128( (void *) idx, _mm_add_epi32( vix, 
                                                  _mm_add_epi32( vsmul( viy, Dy ), 
                                                                 vsmul( viz, Dz ) ) ) ); 

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm_mul_ps( vq, vnorm1);
    vqny = _mm_mul_ps( vq, vnorm2);
    vqnz = _mm_mul_ps( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128 *) vwl1 );
    WL( vqny, vy0, vy1, (__m128 *) vwl2 );
    WL( vqnz, vz0, vz1, (__m128 *) vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
         
         __m128 tmp1, tmp2;        
         const __m128 c1_2 = _mm_set_ps1(0.5f);

         // wp1[k1][k2] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0y[k1], vs0z[k2]), _mm_mul_ps( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0y[k1], vs1z[k2]), _mm_mul_ps( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, c1_2));

         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0z[k2]), _mm_mul_ps( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1z[k2]), _mm_mul_ps( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, c1_2));

         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0y[k2]), _mm_mul_ps( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1y[k2]), _mm_mul_ps( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, c1_2));

      }
    }

    // looping by particle on the outside loop yields the best performance
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

       pjpart = pj + idx[k];
          
	   // accumulate j1
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < ORDER; k1++ ) {
			  p0[k1].j1 += vwl1[k1].v[k] * vwp1[k3][k2].v[k];
		   }
		 }
	   }
   
	   // accumulate j2
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < ORDER; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j2 += vwl2[k2].v[k] * vwp2[k3][k1].v[k];
		   }
		 }
	   }
   
	   
	   // accumulate j3
	   for( k3 = 0; k3 < ORDER; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
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
