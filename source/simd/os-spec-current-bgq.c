/*****************************************************************************************

Charge conserving current deposition, BG/Q optimized version

*****************************************************************************************/


/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-bgq.h"

#include "vector-bgq.h"
#include "splines-bgq.h"
#include "os-spec-push-bgq.h"

/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
   
  vwl[0] = vec_mul(vqn, vec_sub(vx1, vx0));
    
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition.
This expects that the normalization factor has been multiplied by 1/2 i.e.
vqn = q dx / dt / 2
***********************************************************************/

inline void vwl_s2( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {

  vector const c1_2 = vec_splats( 0.5 );
  
  vector d    = vec_sub( vx1, vx0 );
  vector s1_2 = vec_sub( c1_2, vx0 );
  vector p1_2 = vec_add( c1_2, vx0 );
  
  vector n     = vec_mul( vqn, d );

  // qn * d * ( s1_2 - 0.5 * d ) 
  vwl[0] = vec_mul( n, vec_nmsub( c1_2, d, s1_2 ));
  // qn * d * ( p1_2 + 0.5 * d ) 
  vwl[1] = vec_mul( n, vec_madd(  c1_2, d, p1_2 ));
    
}


/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition.
This expects that the normalization factor has been multiplied by 1/8 i.e.
vqn = q dx / dt / 8
***********************************************************************/

inline void vwl_s3( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
  
   vector const c1_2 = vec_splats( 0.5 );
   vector const c3_4 = vec_splats( 0.75 );
   vector const c1_3 = vec_splats( 1.0/3.0 );
   
   vector   d = vec_sub(  vx1, vx0 );
   vector   s = vec_sub( c1_2, vx0 );
   vector   p = vec_add( c1_2, vx0 );
   vector d_3 = vec_mul( c1_3,   d );
   vector   n = vec_mul( vqn, d );

   vwl[0] = vec_mul( c1_2, vec_msub( s, s, vec_mul( d, vec_sub( s, d_3 ) ) ) );
   vwl[1] = vec_nmsub( d, vec_add( vx0, d_3 ), vec_nmsub( vx0, vx0, c3_4 ) );
   vwl[2] = vec_mul( c1_2, vec_madd( p, p, vec_mul( d, vec_add( p, d_3 ) ) ) );

   vwl[0] = vec_mul( n, vwl[0] );
   vwl[1] = vec_mul( n, vwl[1] );
   vwl[2] = vec_mul( n, vwl[2] );
    
}

/***********************************************************************
vwl_s4
***********************************************************************/

inline void vwl_s4( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {

   vector const c1_2 = vec_splats( 0.5 );
   vector const c1_4 = vec_splats( 0.25 );
   vector const c1_6 = vec_splats( 1.0/6.0 );
   vector const c3_2 = vec_splats( 1.5 );
   vector const c3_4 = vec_splats( 0.75 );
   vector const c1_3 = vec_splats( 1.0/3.0 );
   vector const c2_3 = vec_splats( 2.0/3.0 );
   
   vector d = vec_sub( vx1, vx0 );
   vector s = vec_sub( c1_2, vx0 );
   vector p = vec_add( c1_2, vx0 );

   vector t = vec_sub( vx0, c1_6 );
   vector u = vec_add( vx0, c1_6 );

   vector s2 = vec_mul( s, s );
   vector s3 = vec_mul( s2, s );

   vector p2 = vec_mul( p, p );
   vector p3 = vec_mul( p2, p );

   vector d_2 = vec_mul( d, c1_2 );
   vector d_4 = vec_mul( d, c1_4 );

   vector n = vec_mul( vqn, d );
   
   vwl[0] = vec_mul( vec_madd( vec_nmsub( c3_2, s2, vec_mul( d, vec_sub( s, d_4 ) ) ), d, s3 ), c1_6 );
   vwl[1] = vec_madd ( vec_madd( d_2, vec_add( t, d_4 ),  vec_msub( vec_mul( c3_4, t ), t, c1_3 ) ), d, vec_madd( c1_2, p3, vec_sub( c2_3, p2 ) ) );
   vwl[2] = vec_nmsub( vec_madd( d_2, vec_add( u, d_4 ),  vec_msub( vec_mul( c3_4, u ), u, c1_3 ) ), d, vec_madd( c1_2, s3, vec_sub( c2_3, s2 ) ) );
   vwl[3] = vec_mul( vec_madd( vec_madd ( c3_2, p2, vec_mul( d, vec_add( p, d_4 ) ) ), d, p3 ), c1_6 );

   vwl[0] = vec_mul( n, vwl[0] );
   vwl[1] = vec_mul( n, vwl[1] );
   vwl[2] = vec_mul( n, vwl[2] );
   vwl[3] = vec_mul( n, vwl[3] );
    
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
      weights, so jnorm must be dx / dt / 2. In 3D it was a division by 3, so  jnorm must
      be dx / dt / 3. To avoid repeating this every time the DEP_CURRENT functions are
      called, this is done in the ADVANCE_DEPOSIT functions.

Optimizations

Several optimizations were considered for the current deposition, but the standard
serialized current accumulation strategy gave the best results. The algorithms tested
were:

mk. I    - Standard Algorithm, baseline performance (1.0)
mk. II   - Add one cell at a time, directly on main memory (0.904)
mk. IIa  - Copy one cell to local variable, add components, copy back to main memory (0.804)
mk. IIb  - Same as mk. IIa using memcpy
mk. III  - Copy one line to local buffer, add components, copy back to main memory (0.246)
mk. IIIa - Same as mk. III using memcpy (0.246)
mk. IV   - Copy plane to local buffer, add components, copy back to main memory (0.242)
mk. IVa  - Same as mk. IV using memcpy (0.609)

I believe that, since memory access is not aligned, scalar memory access is the most 
efficient, and mk. I minimizes the number of scalar accesses.

The current deposition may benefit from aligning the current buffer memory by adding a 
dummy 4th component. We could then read 1 aligned vector (cell) from memory, accumulate
current, and store it.

*/

/***************   Generate Function names based on interpolation order  ****************/

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

#define SPLINE  ONAME( vspline, ORDER )
#define WL      ONAME( vwl, ORDER )

void DEP_CURRENT_1D 
(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf1D * restrict const part)
{
 
  typedef struct { t_real j1, j2, j3; } t_current;

  vector j2[NP], j3[NP];

  vector vwl1[ORDER];
  
  vector vx0, vx1, vq, vvy, vvz;
  
  DECLARE_ALIGNED_32( int ix[VEC_WIDTH] );
      
  vector vqnx, vqvy, vqvz;
  vector vs0x[NP], vs1x[NP]; 

  
  unsigned int np;
  unsigned int i, k, k1;
  
  vector const c1_2 = vec_splats( 0.5 );

  unsigned int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  vector const vnorm1 = vec_splats(norm[0]);
  
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
  
    // load 4 particles
    // BGQ does not have integer vector registers
    LOAD4P1D( part, i, vx0, vx1, vq, vvy, vvz, ix );
 
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    
    vqnx = vec_mul( vq, vnorm1 );
    vqvy = vec_mul( vq, vvy );
    vqvz = vec_mul( vq, vvz );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );

    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      vector vwp1 = vec_mul( c1_2, vec_add(vs0x[k1], vs1x[k1]) );
      j2[k1] = vec_mul( vqvy, vwp1 );
      j3[k1] = vec_mul( vqvz, vwp1 );
    }

    // New version loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      t_current * restrict p0 = pj + ix[k];

      // accumulate j1
      for ( k1 = 0; k1 < ORDER; k1++ ) {
        p0[k1].j1 += vec_extract(vwl1[k1],k);
      }

      // accumulate j2, j3
      for ( k1 = 0; k1 < NP; k1++ ) {
        p0[k1].j2 += vec_extract(j2[k1],k);
        p0[k1].j3 += vec_extract(j3[k1],k);
      }

    }

  }
  
}

void DEP_CURRENT_2D 
(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf2D * restrict const part)
{
 
  typedef struct { t_real j1, j2, j3; } t_current;

  vector j3[NP][NP];

  vector vwp1[NP], vwp2[NP]; 
  vector vwl1[ORDER], vwl2[ORDER];
  
  vector vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int ix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int iy[VEC_WIDTH] );
		  
  vector vqnx, vqny, vqvz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  
  unsigned int np;
  unsigned int i, k, k1, k2;
  

  unsigned int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  vector const vnorm1 = vec_splats(norm[0]);
  vector const vnorm2 = vec_splats(norm[1]);
  
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
  
    vector const c1_3 = vec_splats( 1.0/3.0 );
    
    // load 4 particles
    // BGQ does not have integer vector registers
    LOAD4P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, ix, iy );
 
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = vec_mul( vq, vnorm1 );
    vqny = vec_mul( vq, vnorm2 );
    vqvz = vec_mul( vec_mul( vq, vvz), c1_3 );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );
    WL( vqny, vy0, vy1, vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k] = vec_add(vs0y[k], vs1y[k]);
      vwp2[k] = vec_add(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        vector const c1_2 = vec_splats( 0.5 );
        vector tmp1, tmp2;        
		
		tmp1 = vec_madd( vs1y[k2], vs1x[k1], vec_mul( vs0x[k1], vs0y[k2] ) );  // tmp1 = s0x*s0y + s1x*s1y
		tmp2 = vec_madd( vs0y[k2], vs1x[k1], vec_mul( vs0x[k1], vs1y[k2] ) );  // tmp2 = s0x*s1y + s1x*s0y
		
		j3[k1][k2] = vec_mul( vqvz, vec_madd( tmp2, c1_2, tmp1 ) );          // j3 = vqvz * (tmp1 + 0.5*tmp2)
      }
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {
       t_current* restrict pjpart = pj + ix[k] + iy[k] * Dy;
   
	   // accumulate j1
	   for( k2 = 0; k2 < NP; k2++ ) {
		 for ( k1 = 0; k1 < ORDER; k1++ ) {
		   pjpart[k2*Dy + k1].j1 += vec_extract(vec_mul(vwl1[k1],vwp1[k2]),k);
		 }
	   }
   
	   // accumulate j2
	   for( k2 = 0; k2 < ORDER; k2++ ) {
		 for ( k1 = 0; k1 < NP; k1++ ) {
			pjpart[k2*Dy + k1].j2 += vec_extract(vec_mul(vwl2[k2],vwp2[k1]),k);
		 }
	   }
   
	   // accumulate j3
	   for( k2 = 0; k2 < NP; k2++ ) {
		 for ( k1 = 0; k1 < NP; k1++ ) {
			pjpart[k2*Dy + k1].j3 += vec_extract(j3[k1][k2],k);
		 }
	   }
	}
  }
  
}


void DEP_CURRENT_3D
(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                              double * restrict const norm, t_split_buf3D * restrict const part)
{
        
  typedef struct { t_real j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  vector vnorm1, vnorm2, vnorm3;
  
  vector vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  vector vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  vector vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int ix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int iy[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int iz[VEC_WIDTH] );
  
  vector vqnx, vqny, vqnz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 

  const vector c1_2  = vec_splats( 0.5);
  
  
  unsigned int np;
  unsigned int i, k, k1, k2, k3;
  unsigned int Dx,Dy,Dz;
  
  Dy = size[0];    
  Dz = Dy * size[1]; 
  

  // The particle cell indexes in t_vpbuf3D are indexed to 1 	(subtract 1)
  // The deposition starts at current[*][ix-1][ij-1][ik-1] 	    (subtract 1)
  // 												Total :  subtract 2
  

   t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                                   ( offset[1] - OFFSET ) * 3 * Dy + 
                                                   ( offset[2] - OFFSET ) * 3 * Dz );


  // norm[i] = dx[i] / dt / 6
  vnorm1 = vec_splats( norm[0]);
  vnorm2 = vec_splats( norm[1]);
  vnorm3 = vec_splats( norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np % 4 != 0 ) {
    for( i = 0; i < 4 - np%4; i ++ ) {
       k = np + i;
       
       part -> x0[k] = part -> x1[k] = 0.;
       part -> y0[k] = part -> y1[k] = 0.;
       part -> z0[k] = part -> z1[k] = 0.;
       part -> q[k] = 0.;
       part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }
  
  // Each virtual particle uses 12 doubles, so 4 vp = 48 doubles
  for( i = 0; i < np; i += 4 ) {
      
	// load 2 particles
    LOAD4P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz );
        
    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = vec_mul( vq, vnorm1 );
    vqny = vec_mul( vq, vnorm2 );
    vqnz = vec_mul( vq, vnorm3 );
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );
    WL( vqny, vy0, vy1, vwl2 );
    WL( vqnz, vz0, vz1, vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
	  for( k1 = 0; k1 < NP; k1++ ) {
 
         vector tmp1, tmp2;        
        
         // wp1[k1][k2] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 = vec_mul( vs0y[k1], vs0z[k2]);
         tmp2 = vec_mul( vs0y[k1], vs1z[k2] );
         tmp1 = vec_madd( vs1z[k2], vs1y[k1], tmp1 );
         tmp2 = vec_madd( vs0z[k2], vs1y[k1], tmp2 );

         vwp1[k2][k1] = vec_madd( tmp2, c1_2, tmp1 );

         // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] + 
         //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

         tmp1 = vec_mul( vs0x[k1], vs0z[k2] );
         tmp2 = vec_mul( vs0x[k1], vs1z[k2] );
         tmp1 = vec_madd( vs1z[k2], vs1x[k1], tmp1 );
         tmp2 = vec_madd( vs0z[k2], vs1x[k1], tmp2 );
         vwp2[k2][k1] = vec_madd( tmp2, c1_2, tmp1 );

         // wp2[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] + 
         //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

         tmp1 = vec_mul( vs0x[k1], vs0y[k2] );
         tmp2 = vec_mul( vs0x[k1], vs1y[k2] );
         tmp1 = vec_madd( vs1y[k2], vs1x[k1], tmp1 );
         tmp2 = vec_madd(vs0y[k2],  vs1x[k1], tmp2 );
         vwp3[k2][k1] = vec_madd( tmp2, c1_2, tmp1 );

      }
    }
    
    // loop by particle on the outside loop
   for ( k = 0; k < 4; k ++ ) {
    
	   pjpart = pj + ix[k] + iy[k] * Dy + iz[k] * Dz;
		  
	   // accumulate j1
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < ORDER; k1++ ) {
			  p0[k1].j1 += vec_extract( vec_mul(vwl1[k1],vwp1[k3][k2]),k);
		   }
		 }
	   }
   
	   // accumulate j2
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < ORDER; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j2 += vec_extract( vec_mul(vwl2[k2],vwp2[k3][k1]),k);
		   }
		 }
	   }
   
	   
	   // accumulate j3
	   for( k3 = 0; k3 < ORDER; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j3  += vec_extract( vec_mul(vwl3[k3],vwp3[k2][k1]),k);
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







