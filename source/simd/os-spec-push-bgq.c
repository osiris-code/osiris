/*****************************************************************************************

Relativistic particle pusher, BG/Q optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-push-bgq.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-bgq.h"
#include "splines-bgq.h"

#include "os-spec-current-bgq.h"

// Internal calculations on BG/Q QPX code are always done in double precision, even when
// global data is store in single precision

#define REAL double
#include "split-vec.h"
#undef REAL

/*****************************************************************************************
vhp_shift

Gets the nearest half points of the 4 particles loaded as a vector for interpolating
scattered grids. This routine is for positions defined as distance to the nearest grid
points.
*****************************************************************************************/

#define vhp_shift(dx,ix, h, ixh, delta ) { \
  h = vec_sel( vec_splats( 1.0 ), vec_splats( 0.0 ), dx ); \
  vec_sta( vec_ctiw( h ), 0x00, ixh ); \
  for( unsigned i = 0; i < 4; i++ ) \
    ixh[i] = ix[i] - (ixh[i] * delta); \
}						

/*****************************************************************************************
vntrim

Returns the required cell shift (trim) for particles that must remain in [ -0.5, 0.5 [
to the nearest cell point.

dx = ( vx < -0.5 ) ? -1 : 0 + ( vx >= +0.5 ) ? +1 : 0

x'  = x  - dx
ix' = ix + dx

*****************************************************************************************/

inline vector vntrim( const vector vx )
{
  const vector oneHalf = vec_splats( 0.5 );
  const vector one     = vec_splats( 1.0 );
  const vector zero    = vec_splats( 0.0 );
  
  vector va = vec_sel( one, zero, vec_add( oneHalf, vx ) );
  vector vb = vec_sel( zero, one, vec_sub( vx, oneHalf ) );
  
  return vec_sub( vb, va );
}


/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene ) \
{																			\
   vector const one = vec_splats( 1.0 ); \
   register vector vut1, vut2, vut3; \
   register vector u2; \
   register vector vgamma_tem, votsq; \
   \
   /* Perform first half of electric field acceleration.*/\
   ve1 = vec_mul( ve1, vtem ); \
   ve2 = vec_mul( ve2, vtem ); \
   ve3 = vec_mul( ve3, vtem ); \
   \
   vut1 = vec_add( vu1, ve1 ); \
   vut2 = vec_add( vu2, ve2 ); \
   vut3 = vec_add( vu3, ve3 ); \
   \
   /* Perform first half of the rotation and store in u. */ \
   u2 = vec_madd( vut1, vut1, vec_madd( vut2, vut2, vec_mul( vut3, vut3 ) ) ); \
   \
   vgamma_tem = vec_swsqrt_nochk( vec_add( u2, one ) ); \
   \
   vene = vec_swdiv_nochk( u2, vec_add( vgamma_tem, one ) ); \
   \
   vgamma_tem = vec_swdiv_nochk( vtem, vgamma_tem ); \
   \
   vb1 = vec_mul( vb1, vgamma_tem ); \
   vb2 = vec_mul( vb2, vgamma_tem ); \
   vb3 = vec_mul( vb3, vgamma_tem ); \
   \
   vu1 = vec_madd( vb3, vut2, vut1 ); \
   vu2 = vec_madd( vb1, vut3, vut2 ); \
   vu3 = vec_madd( vb2, vut1, vut3 ); \
   \
   vu1 = vec_nmsub( vb2, vut3, vu1 ); \
   vu2 = vec_nmsub( vb3, vut1, vu2 ); \
   vu3 = vec_nmsub( vb1, vut2, vu3 ); \
   \
   /* Perform second half of the rotation. */\
   votsq = vec_madd( vb1,   vb1, one ); \
   votsq = vec_madd( vb2, vb2, votsq ); \
   votsq = vec_madd( vb3, vb3, votsq ); \
   \
   votsq = vec_r(votsq); \
   votsq = vec_add( votsq, votsq );	\
   \
   vb1 = vec_mul( vb1, votsq );	\
   vb2 = vec_mul( vb2, votsq );	\
   vb3 = vec_mul( vb3, votsq );	\
   \
   vut1 = vec_madd( vb3, vu2, vut1 ); \
   vut2 = vec_madd( vb1, vu3, vut2 ); \
   vut3 = vec_madd( vb2, vu1, vut3 ); \
   \
   vut1 = vec_nmsub( vb2, vu3, vut1 ); \
   vut2 = vec_nmsub( vb3, vu1, vut2 ); \
   vut3 = vec_nmsub( vb1, vu2, vut3 ); \
   \
   /* Perform second half of electric field acceleration.*/ \
   vu1 = vec_add( vut1, ve1 ); \
   vu2 = vec_add( vut2, ve2 ); \
   vu3 = vec_add( vut3, ve3 ); \
}


/****************************************************************************************

  Generate specific particle push functions for 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
// These cannot be used to define the fortran function names because of recursive macro
// expansion issues.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o


/********************************** Linear interpolation ********************************/

#define ADVANCE_DEPOSIT_1D FNAME( vadvance_deposit_1d_s1 )
#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s1 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s1 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s1 )

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Grid offset ( 1 + ORDER/2 )
#define OFFSET 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ADVANCE_DEPOSIT_1D FNAME( vadvance_deposit_1d_s2 )
#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s2 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s2 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s2 )

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ADVANCE_DEPOSIT_1D FNAME( vadvance_deposit_1d_s3 )
#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s3 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s3 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s3 )

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ADVANCE_DEPOSIT_1D FNAME( vadvance_deposit_1d_s4 )
#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s4 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s4 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s4 )

#define ORDER 4
#define NP 5
#define OFFSET 3

#include __FILE__

/****************************************************************************************/

#else

/****************************************************************************************

  Template function definitions for 2D and 3D particle push

****************************************************************************************/

/***************   Generate Function names based on interpolation order  ****************/

#define SPLINE  ONAME( vspline, ORDER )
#define SPLINEH  ONAME( vsplineh, ORDER )

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

/********************************** 1D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_1D
( int* restrict const ix, t_real* restrict const x, t_real* restrict const u, t_real* restrict const q, 
  int const* restrict const i0, int const * restrict const i1, double const * restrict const rqm,
  t_real * restrict const e, t_real * restrict const b, 
  int const * restrict const emf_size, int const * restrict const emf_offset,  
  t_real * restrict const j, int const * restrict const j_size, int const * restrict const j_offset, 
  double const * restrict const dx, double const * restrict const dt, double *ene )
{  
  int *pix;
  t_real *px, *pu, *pq;
  
  unsigned k, i, np, np_total;
    
  DECLARE_ALIGNED_32( double jnorm[1] );
  
  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf1D vpbuf;

   // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components 

  t_real* const e1 = e   + (emf_offset[0]-OFFSET)*deltaX;
  t_real* const e2 = e1 + 1;
  t_real* const e3 = e1 + 2;

  t_real* const b1 = b   + (emf_offset[0]-OFFSET)*deltaX;
  t_real* const b2 = b1  + 1;
  t_real* const b3 = b1  + 2;
  
  vector const vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vector const vdt_dx = vec_splats( (*dt) / dx[0] );
  vector const vdt_dy = vec_splats( (*dt) / dx[1] );

  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt);
  
  
  // jump to 1st particle
  px  = x  + (*i0-1);
  pix = ix + (*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  + (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;
  
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	for( i = 0; i < VEC_WIDTH - np_total % VEC_WIDTH; i ++ ) {
	  pix[ (np_total + i)     ] = 1;
	  px[  (np_total + i)     ] = 0.;
	  
	  pu[ 3 * (np_total + i)     ] = 0.;
	  pu[ 3 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i) + 2 ] = 0.;
	  
	  pq[ (np_total + i ) ] = 0.;
	}
	
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }
  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) p_cache_size_1D is not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

	// Initialize energy
	vector vene_group = vec_splats( 0.0 );
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_1D ) ? np_total - k : p_cache_size_1D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i += VEC_WIDTH ) {
	  vector vene;
	
	  vector vx0;
	  DECLARE_ALIGNED_32( int vix[VEC_WIDTH] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;	  
	  vector vq;
	  
	  unsigned int l;
  
	  // Load particle positions 
	  vx0 = vec_lda( 0x00, px );
	  
	  for( l = 0; l < VEC_WIDTH; l++ ) {
	    vix[l] = pix[l];
	  }
 
	  // Interpolate fields
	  {

		 vector hx;
		 int idxi[VEC_WIDTH], idxih[VEC_WIDTH];
		 vector vwx[NP], vwxh[NP];
	     
	     unsigned int l;
	     
	     for( l=0; l<VEC_WIDTH; l++) {
	       idxi[l] =      3 * vix[l];
	     }

		 SPLINE( vx0, vwx );

		 vhp_shift( vx0, idxi, hx, idxih, 3 ); 
	 
		 SPLINEH( vx0, hx, vwxh );

		 // Interpolate E field
         {
			unsigned int k1;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = e1 + idxih[k1];
			   fld2[k1] = e2 +  idxi[k1];
			   fld3[k1] = e3 +  idxi[k1];
			}
						 
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		    		    
			register vector f1point[NP], f2point[NP], f3point[NP];

			for ( k1 = 0; k1 < NP; k1++ ) {
				unsigned int shift = k1*3;

				f1point[k1] = LOADFLD4( fld1, shift ); 
				f2point[k1] = LOADFLD4( fld2, shift ); 
				f3point[k1] = LOADFLD4( fld3, shift ); 
			}
			       
			for ( k1 = 0; k1 < NP; k1 ++ ) {
				ve1 = vec_madd( vwxh[k1], f1point[k1], ve1 );
				ve2 = vec_madd( vwx[k1],  f2point[k1], ve2 );
				ve3 = vec_madd( vwx[k1],  f3point[k1], ve2 );
			}
			        
         }		 

		 
		 // Interpolate B field
         {
			unsigned int k1;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
				fld1[k1] = b1 +  idxi[k1];
				fld2[k1] = b2 + idxih[k1];
				fld3[k1] = b3 + idxih[k1];
			}

			vb1 = zero;
			vb2 = zero;
			vb3 = zero;

			register vector f1point[NP], f2point[NP], f3point[NP];
			for ( k1 = 0; k1 < NP; k1++ ) {
				unsigned int shift = k1*3;

				f1point[k1] = LOADFLD4( fld1, shift ); 
				f2point[k1] = LOADFLD4( fld2, shift ); 
				f3point[k1] = LOADFLD4( fld3, shift ); 
			}

			for ( k1 = 0; k1 < NP; k1 ++ ) {
				vb1 = vec_madd( vwx[k1],  f1point[k1], vb1 );
				vb2 = vec_madd( vwxh[k1], f2point[k1], vb2 );
				vb3 = vec_madd( vwxh[k1], f3point[k1], vb3 );
			}
         
         }

	  }

	  // ---------- advance momenta
	  
	  // Load momenta
	  VEC_LDA4v3( vu1, vu2, vu3, pu );
	  
	  // results are stored in vu1, vu2, vu3 and vene
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
      
      // Store momenta
      VEC_STA4v3(pu, vu1, vu2, vu3);
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg;
		register vector vx1, vtrx;
		unsigned int l;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
		
        vq = vec_lda( 0x00, pq );

		// Accumulate energy
		vene_group = vec_madd( vq, vene, vene_group );
		
		// get velocities
		vu1 = vec_mul( vu1, vrg );
		vu2 = vec_mul( vu2, vrg ); // this is required for current deposition
		vu3 = vec_mul( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = vec_madd( vdt_dx, vu1, vx0 );

        // Store virtual particles with positions still indexed to the original cell
        STOREU4P1D( vpbuf, vpbuf.np, vx0, vx1, vq, vu2, vu3, vix );
	  
	    // Trim positions and store results
		vtrx = vntrim(vx1);

		vx1  = vec_sub( vx1, vtrx );
		vec_sta( vx1, 0x00, px );
			
		vtrx = vec_ctiw( vtrx ) ;
		vec_sta( vtrx, 0x00, dix );
		
        for( l=0; l < VEC_WIDTH; l++) {
          // store new particle cell indexes
          pix[l] = vix[l] + dix[l];
        }

	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit1D( &vpbuf, dix );
	  
	  // advance pointers
	  px  += VEC_WIDTH;  
	  pix += VEC_WIDTH;
	  pu  += 3 * VEC_WIDTH;
	  pq  += VEC_WIDTH;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_1D( j, j_size, j_offset, jnorm, &vpbuf ); 

	// Accumulate energy from group
	*ene += vec_reduce_add( vene_group );

  }
    
}


/********************************** 2D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_2D
( int* restrict const ix, t_real* restrict const x, t_real* restrict const u, t_real* restrict const q, 
  int const* restrict const i0, int const * restrict const i1, double const * restrict const rqm,
  t_real * restrict const e, t_real * restrict const b, 
  int const * restrict const emf_size, int const * restrict const emf_offset,  
  t_real * restrict const j, int const * restrict const j_size, int const * restrict const j_offset, 
  double const * restrict const dx, double const * restrict const dt, double *ene )
{  
  int *pix;
  t_real *px, *pu, *pq;
  
  unsigned k, i, np, np_total;
    
  DECLARE_ALIGNED_32( double jnorm[2] );
  
  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

   // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components 
  unsigned int const deltaY = emf_size[0] * deltaX;

  t_real* const e1 = e   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  t_real* const e2 = e1 + 1;
  t_real* const e3 = e1 + 2;

  t_real* const b1 = b   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  t_real* const b2 = b1  + 1;
  t_real* const b3 = b1  + 2;
  
  vector const vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vector const vdt_dx = vec_splats( (*dt) / dx[0] );
  vector const vdt_dy = vec_splats( (*dt) / dx[1] );

  
  // Normalization for currents
  // The factor 2 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;
  
  
  // jump to 1st particle
  px  = x  + 2*(*i0-1);
  pix = ix + 2*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  + (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;
  
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	for( i = 0; i < VEC_WIDTH - np_total % VEC_WIDTH; i ++ ) {
	  pix[ 2 * (np_total + i)     ] = 1;
	  pix[ 2 * (np_total + i) + 1 ] = 1;
	  
	  px[ 2 * (np_total + i)     ] = 0.;
	  px[ 2 * (np_total + i) + 1 ] = 0.;
	  
	  pu[ 3 * (np_total + i)     ] = 0.;
	  pu[ 3 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i) + 2 ] = 0.;
	  
	  pq[ (np_total + i ) ] = 0.;
	}
	
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }
  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) p_cache_size_2D is not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

	// Initialize energy
	vector vene_group = vec_splats( 0.0 );
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i += VEC_WIDTH ) {
	  vector vene;
	
	  vector vx0, vy0;
	  DECLARE_ALIGNED_32( unsigned int vix[VEC_WIDTH] );
	  DECLARE_ALIGNED_32( unsigned int viy[VEC_WIDTH] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;	  
	  vector vq;
	  
	  unsigned int l;
  
	  // Load particle positions 
	  VEC_LDA4v2( vx0, vy0, px );
	  
	  for( l = 0; l < VEC_WIDTH; l++ ) {
	    vix[l] = pix[2*l    ];
	    viy[l] = pix[2*l + 1];
	  }
 
	  // Interpolate fields
	  {

		 vector hx, hy;
		 vector  vx0h, vy0h;
		 int idxi[VEC_WIDTH], idxj[VEC_WIDTH], idxih[VEC_WIDTH], idxjh[VEC_WIDTH];
		 register vector vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	     
	     unsigned int l;
	     
	     for( l=0; l<VEC_WIDTH; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	     }

		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );

		 vhp_shift( vx0, idxi, hx, idxih, 3 ); 
		 vhp_shift( vy0, idxj, hy, idxjh, deltaY ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );

		 // Interpolate E field
         {
			unsigned int k1, k2;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = e1 + idxih[k1] +  idxj[k1];
			   fld2[k1] = e2 +  idxi[k1] + idxjh[k1];
			   fld3[k1] = e3 +  idxi[k1] +  idxj[k1];
			}
						 
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwxh[k1], f1point[k1], f1line );
				 f2line = vec_madd( vwx[k1],  f2point[k1], f2line );
				 f3line = vec_madd( vwx[k1],  f3point[k1], f3line );
			  }

			  ve1 = vec_madd( vwy[k2],  f1line, ve1 );
			  ve2 = vec_madd( vwyh[k2], f2line, ve2 );
			  ve3 = vec_madd( vwy[k2],  f3line, ve3 );
			}
			        
         }		 

		 
		 // Interpolate B field
         {
			unsigned int k1, k2;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = b1 +  idxi[k1] + idxjh[k1];
			   fld2[k1] = b2 + idxih[k1] +  idxj[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];
			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
              
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwx[k1],  f1point[k1], f1line );
				 f2line = vec_madd( vwxh[k1], f2point[k1], f2line );
				 f3line = vec_madd( vwxh[k1], f3point[k1], f3line );
			  }

			  vb1 = vec_madd( vwyh[k2], f1line, vb1 );
			  vb2 = vec_madd( vwy[k2],  f2line, vb2 );
			  vb3 = vec_madd( vwyh[k2], f3line, vb3 );
			}
         
         }

	  }

	  // ---------- advance momenta
	  
	  // Load momenta
	  VEC_LDA4v3( vu1, vu2, vu3, pu );
	  
	  // results are stored in vu1, vu2, vu3 and vene
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
      
      // Store momenta
      VEC_STA4v3(pu, vu1, vu2, vu3);
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg;
		register vector vx1, vy1, vtrx, vtry;
		unsigned int l;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
		
        vq = vec_lda( 0x00, pq );

		// Accumulate energy
		vene_group = vec_madd( vq, vene, vene_group );
		
		// get velocities
		vu1 = vec_mul( vu1, vrg );
		vu2 = vec_mul( vu2, vrg );
		vu3 = vec_mul( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = vec_madd( vdt_dx, vu1, vx0 );
		vy1 = vec_madd( vdt_dy, vu2, vy0 );

        // Store virtual particles with positions still indexed to the original cell
        STOREU4P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vu3, vix, viy );
	  
	    // Trim positions and store results
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);

		vx1  = vec_sub( vx1, vtrx );
		vy1  = vec_sub( vy1, vtry );
		VEC_STA4v2( px, vx1, vy1 );
			
		vtrx = vec_ctiw( vtrx ) ;
		vtry = vec_ctiw( vtry ) ;
		vec_sta( vtrx, 0x00, dix );
		vec_sta( vtry, 0x00, diy );
		
        for( l=0; l < VEC_WIDTH; l++) {
          // store new particle cell indexes
          pix[2*l     ] = vix[l] + dix[l];
          pix[2*l + 1 ] = viy[l] + diy[l];
          
          // calculate crossings - this can be done using vtr1 / vtr2
          cross[l] = (( dix[l] ) ? 1 : 0 ) + (( diy[l] ) ? 2 : 0 );
        }

	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );
	  
	  // advance pointers
	  px  += 2 * VEC_WIDTH;  
	  pix += 2 * VEC_WIDTH;
	  pu  += 3 * VEC_WIDTH;
	  pq  += VEC_WIDTH;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf ); 

	// Accumulate energy from group
	*ene += vec_reduce_add( vene_group );

  }
    
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int* restrict const ix, t_real* restrict const x, t_real* restrict const u, t_real* restrict const q, 
  int const* restrict const i0, int const * restrict const i1, double const * restrict const rqm,
  int const *ilb2 , 
  t_real * restrict const e, t_real * restrict const b, 
  int const * restrict const emf_size, int const * restrict const emf_offset,  
  t_real * restrict const j, int const * restrict const j_size, int const * restrict const j_offset, 
  double const * restrict const dx, double const * restrict const dt, double *ene )
{
  
  int *pix;
  t_real *px, *pu, *pq;
  
  int k, i, np, np_total;
    
  DECLARE_ALIGNED_32( double jnorm[2] );
  
  DECLARE_ALIGNED_32( unsigned int cross[4] );
  DECLARE_ALIGNED_32( int dix[4] );
  DECLARE_ALIGNED_32( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

   // get pointers to position 0,0 of each field component
  int const deltaX = 3; // 3 field components 
  int const deltaY = emf_size[0] * deltaX;

  t_real * const e1 = e   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  t_real * const e2 = e1 + 1;
  t_real * const e3 = e1 + 2;

  t_real * const b1 = b   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  t_real * const b2 = b1  + 1;
  t_real * const b3 = b1  + 2;
  
  vector const vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vector const vdt_dx = vec_splats( (*dt) / dx[0] );
  vector const vdt_dy = vec_splats( (*dt) / dx[1] );
  
  // Normalization for currents
  // The factor 2 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;

  #if ( ORDER == 1 ) || ( ORDER == 3 )
  vector const gshift2 = vec_splats( (double) (*ilb2 - 2 ) );
  #else  
  vector const gshift2 = vec_splats( (double) (*ilb2 - 2 ) - 0.5 );
  #endif
  
  vector const dr      = vec_splats( dx[1] );
  vector const rdr     = vec_splats( 1.0/dx[1] );
  vector const vdt     = vec_splats( *dt );
  
  // jump to 1st particle
  px  = x  + 2*(*i0-1);
  pix = ix + 2*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  + (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	for( i = 0; i < VEC_WIDTH - np_total % VEC_WIDTH; i ++ ) {
	  pix[ 2 * (np_total + i)     ] = 1;
	  pix[ 2 * (np_total + i) + 1 ] = 1;
	  
	  px[ 2 * (np_total + i)     ] = 0.;
	  px[ 2 * (np_total + i) + 1 ] = 0.25;
	  
	  pu[ 3 * (np_total + i)     ] = 0.;
	  pu[ 3 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i) + 2 ] = 0.;
	  
	  pq[ (np_total + i ) ] = 0.;
	}
	
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) Number of particles to push is still not a multiple of VEC_WIDTH!\n");
       exit(1);
    }
  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

	// Initialize energy
	vector vene_group = vec_splats(0.0);
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i+=VEC_WIDTH ) {
	  vector vene;
	
	  vector vx0, vy0;
	  DECLARE_ALIGNED_32( unsigned int vix[VEC_WIDTH] );
	  DECLARE_ALIGNED_32( unsigned int viy[VEC_WIDTH] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;	  
	  vector vq;

	  int l;
  
	  // Load particle positions 
	  VEC_LDA4v2( vx0, vy0, px );
	  
	  for( l = 0; l < VEC_WIDTH; l++ ) {
	    vix[l] = pix[2*l    ];
	    viy[l] = pix[2*l + 1];
	  }
 
	  // Interpolate fields
	  {

		 vector  hx, hy;
		 int idxi[VEC_WIDTH], idxj[VEC_WIDTH], idxih[VEC_WIDTH], idxjh[VEC_WIDTH];
		 register vector vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     int l;
	     
	     for( l=0; l<VEC_WIDTH; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	     }

		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );

		 vhp_shift( vx0, idxi, hx, idxih, 3 ); 
		 vhp_shift( vy0, idxj, hy, idxjh, deltaY ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = e1 + idxih[k1] + idxj[k1];
			   fld2[k1] = e2 + idxi[k1] + idxjh[k1];
			   fld3[k1] = e3 + idxi[k1] + idxj[k1];
			}
						 
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwxh[k1], f1point[k1], f1line );
				 f2line = vec_madd( vwx[k1],  f2point[k1], f2line );
				 f3line = vec_madd( vwx[k1],  f3point[k1], f3line );
			  }

			  ve1 = vec_madd( vwy[k2],  f1line, ve1 );
			  ve2 = vec_madd( vwyh[k2], f2line, ve2 );
			  ve3 = vec_madd( vwy[k2],  f3line, ve3 );
			}
			        
         }		 

		 
		 // Interpolate B field
        {
			int k1, k2;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = b1 + idxi[k1]  + idxjh[k1];
			   fld2[k1] = b2 + idxih[k1] + idxj[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];
			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
              
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwx[k1],  f1point[k1], f1line );
				 f2line = vec_madd( vwxh[k1], f2point[k1], f2line );
				 f3line = vec_madd( vwxh[k1], f3point[k1], f3line );
			  }

			  vb1 = vec_madd( vwyh[k2], f1line, vb1 );
			  vb2 = vec_madd( vwy[k2],  f2line, vb2 );
			  vb3 = vec_madd( vwyh[k2], f3line, vb3 );
			}
         
         }

	  }

	  // ---------- advance momenta

	  // Load momenta
	  VEC_LDA4v3( vu1, vu2, vu3, pu );
	  
	  // results are stored in vu1, vu2, vu3 and vene
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
	  
	  // The position advance will further modify u2 and u3
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg;
		register vector vx1, vy1, vtrx, vtry;
		register vector vv1, vv2, vv3;

		int l;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
		
        vq = vec_lda( 0, pq );

		// Accumulate energy
		vene_group = vec_madd( vq, vene, vene_group );
	
		// get velocities
		vv1 = vec_mul( vu1, vrg );
		vv2 = vec_mul( vu2, vrg );
		vv3 = vec_mul( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = vec_madd( vdt_dx, vv1, vx0 );

        // x2 - 3D like push
        // This further changes u2 and u3
        {  
		   
		   vector gix2;
		   vector r_old, r_new;
		   vector x2_new, x3_new; 
		   vector tmp;
		   
		   
		   // gshift2 = shift_ix2 - 0.5d0
   
		   // gix2   = viy0 + gshift2;          // global cell
	       gix2   = vec_add( vec_cfid( vec_ldiza( 0x00, viy ) ), gshift2 );
		   r_old  = vec_mul( vec_add( vy0, gix2 ), dr );

		   x2_new = vec_madd( vdt, vv2, r_old );
		   x3_new = vec_mul( vv3, vdt ) ;

		   r_new  = vec_swsqrt_nochk( vec_madd( x3_new, x3_new, vec_mul( x2_new, x2_new ) ) ); 
		   
		   // This is a protection against roundoff for cold plasmas
		   // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
		   vy1 = vec_sel(  vec_sub( vec_mul( r_new, rdr ), gix2 ), vy0, 
		                   vec_nabs( vec_sub( r_old, r_new ) ) );
		   		   		   
		   // Correct p_r and p_\theta to conserve angular momentum
		   tmp  = vec_r( r_new );
		   vu2  = vec_mul( vec_add( vec_mul( vu2, x2_new ), vec_mul( vu3, x3_new ) ), tmp );
		   vu3  = vec_mul( vu3 , vec_mul( r_old, tmp ) );
   
		}        

        // Store Momenta
        VEC_STA4v3(pu, vu1, vu2, vu3);

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU4P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix, viy );

	    // Trim positions and store results
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);

		vx1  = vec_sub( vx1, vtrx );
		vy1  = vec_sub( vy1, vtry );
		VEC_STA4v2( px, vx1, vy1 );
			
		vtrx = vec_ctiw( vtrx ) ;
		vtry = vec_ctiw( vtry ) ;
		vec_sta( vtrx, 0, dix );
		vec_sta( vtry, 0, diy );
		
        for( l=0; l<4; l++) {
          // store new particle cell indexes
          pix[2*l     ] = vix[l] + dix[l];
          pix[2*l + 1 ] = viy[l] + diy[l];
          
          // calculate crossings - this can be done using vtr1 / vtr2
          cross[l] = (( dix[l] ) ? 1 : 0 ) + (( diy[l] ) ? 2 : 0 );
        }

	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );
      
	  // advance x, ix 4 particles * 2 components
	  px  += 2 * VEC_WIDTH;  
	  pix += 2 * VEC_WIDTH;
	  pu  += 3 * VEC_WIDTH; 
	  pq  += VEC_WIDTH;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

	// Accumulate energy from group
	*ene += vec_reduce_add( vene_group );

  }
    
}


/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int* restrict const ix, t_real * restrict const x, t_real * restrict const u, t_real * restrict const q, 
  int const* restrict const i0, int const* restrict const i1, double const * restrict const rqm,
  t_real  * restrict const e, t_real  * restrict const b, 
  int const * restrict const emf_size, int const * restrict const emf_offset,  
  t_real  * restrict const j, int const * restrict const j_size, int const * restrict const j_offset, 
  double const * restrict const dx, double const * restrict const dt, double *ene )
{  
  int *pix;
  t_real  *px, *pu, *pq;
  
  int k, i, np, np_total;
      
  DECLARE_ALIGNED_32( double jnorm[3] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf3D vpbuf;

  // get pointers to position -1,-1,-1 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;
  int const deltaZ = emf_size[1] * deltaY;


  t_real * const e1 = e  + (emf_offset[0]-OFFSET)*deltaX +  
                           (emf_offset[1]-OFFSET)*deltaY + 
                           (emf_offset[2]-OFFSET)*deltaZ;
  t_real * const e2 = e1 + 1;
  t_real * const e3 = e1 + 2;

  t_real * const b1 = b   + (emf_offset[0]-OFFSET)*deltaX + 
                            (emf_offset[1]-OFFSET)*deltaY + 
                            (emf_offset[2]-OFFSET)*deltaZ;
  t_real * const b2 = b1  + 1;
  t_real * const b3 = b1  + 2;
  
  //
  
  vector const vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vector const vdt_dx = vec_splats( (*dt) / dx[0] );
  vector const vdt_dy = vec_splats( (*dt) / dx[1] );
  vector const vdt_dz = vec_splats( (*dt) / dx[2] );

  
  // Normalization for currents
  // The factor 3 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / 3;
  jnorm[1] = dx[1] / (*dt) / 3;
  jnorm[2] = dx[2] / (*dt) / 3;
  
  // jump to 1st particle
  px  = x  + 3*(*i0-1);
  pix = ix + 3*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  +   (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total % VEC_WIDTH; i ++ ) {
	  
	  pix[ 3 * (np_total + i)     ] = 1;
	  pix[ 3 * (np_total + i) + 1 ] = 1;
	  pix[ 3 * (np_total + i) + 2 ] = 1;
	   
	   px[ 3 * (np_total + i)     ] = 0.;
	   px[ 3 * (np_total + i) + 1 ] = 0.;
	   px[ 3 * (np_total + i) + 2 ] = 0.;
	   
	   pu[ 3 * (np_total + i)     ] = 0.;
	   pu[ 3 * (np_total + i) + 1 ] = 0.;
	   pu[ 3 * (np_total + i) + 2 ] = 0.;
	   
	   pq[     (np_total + i)     ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total% VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       fprintf(stderr, "(*error*) Number of particles to push is still not a multiple of VEC_WIDTH!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_3D % 4 ) {
       fprintf(stderr, "(*error*) p_cache_size_3D is not a multiple of 4!\n");
       exit(1);
  }
  
  for ( k = 0; k < np_total; k += p_cache_size_3D ) {

	// Initialize energy
	vector vene_group = vec_splats( 0.0 );
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i+= VEC_WIDTH ) {

  	  vector vene;

	  vector vx0, vy0, vz0;
	  DECLARE_ALIGNED_32( unsigned int vix[VEC_WIDTH] );
	  DECLARE_ALIGNED_32( unsigned int viy[VEC_WIDTH] );
	  DECLARE_ALIGNED_32( unsigned int viz[VEC_WIDTH] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;

	  // Load particle positions 
	  VEC_LDA4v3( vx0, vy0, vz0, px );

	  for( int l = 0; l < VEC_WIDTH; l++ ) {
	    vix[l] = pix[3*l    ];
	    viy[l] = pix[3*l + 1];
	    viz[l] = pix[3*l + 2];
	  }

	  // Interpolate fields
	  {
		 vector  hx, hy, hz;
		 int  idxi[VEC_WIDTH],  idxj[VEC_WIDTH],  idxk[VEC_WIDTH], 
		              idxih[VEC_WIDTH], idxjh[VEC_WIDTH], idxkh[VEC_WIDTH];
		 vector vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];
	     
	     for( unsigned l = 0; l < VEC_WIDTH; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	       idxk[l] = deltaZ * viz[l];
	     }
	     	 
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );
		 SPLINE( vz0, vwz );

		 vhp_shift( vx0, idxi, hx, idxih, 3 ); 
		 vhp_shift( vy0, idxj, hy, idxjh, deltaY ); 
		 vhp_shift( vz0, idxk, hz, idxkh, deltaZ ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );
		 SPLINEH( vz0, hz, vwzh );
	       
         // Interpolate E field
         {
			int k1, k2, k3;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			vector const zero = vec_splats( 0.0 );
						  
			// get pointers to fields in particle cells (i-1, j-1, k-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = e1 + idxih[k1] +  idxj[k1] +  idxk[k1]; 
			   fld2[k1] = e2 +  idxi[k1] + idxjh[k1] +  idxk[k1];
			   fld3[k1] = e3 +  idxi[k1] +  idxj[k1] + idxkh[k1];
            }
            			
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register vector f1plane, f2plane, f3plane;
			   f1plane = zero;
			   f2plane = zero;
			   f3plane = zero;
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register vector f1line, f2line, f3line;
				 f1line = zero;
				 f2line = zero;
				 f3line = zero;

                 #pragma unroll(NP)
                 for( k1 = 0; k1 < NP; k1++ ) {
                    register vector f1point, f2point, f3point;
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point = LOADFLD4( fld1, shift ); 
					f2point = LOADFLD4( fld2, shift ); 
					f3point = LOADFLD4( fld3, shift ); 

					f1line = vec_madd( vwxh[k1], f1point, f1line );
				    f2line = vec_madd(  vwx[k1], f2point, f2line );
				    f3line = vec_madd(  vwx[k1], f3point, f3line );
                 }
				 
				 f1plane = vec_madd(  vwy[k2], f1line, f1plane );
				 f2plane = vec_madd( vwyh[k2], f2line, f2plane );
				 f3plane = vec_madd(  vwy[k2], f3line, f3plane );
			   }
			   
			   ve1 = vec_madd(  vwz[k3], f1plane, ve1 );
			   ve2 = vec_madd(  vwz[k3], f2plane, ve2 );
			   ve3 = vec_madd( vwzh[k3], f3plane, ve3 );
			}
         
         }		 

         // Interpolate B field
         {
			int k1, k2, k3;
			t_real *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1, k-1)
			for( k1 = 0; k1 < VEC_WIDTH; k1++) {
			   fld1[k1] = b1 +  idxi[k1] + idxjh[k1] + idxkh[k1]; 
			   fld2[k1] = b2 + idxih[k1] +  idxj[k1] + idxkh[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1] +  idxk[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register vector f1plane, f2plane, f3plane;
			   f1plane = zero;
			   f2plane = zero;
			   f3plane = zero;
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register vector f1line, f2line, f3line;
				 f1line = zero;
				 f2line = zero;
				 f3line = zero;

				 #pragma unroll(NP)
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
                    register vector f1point, f2point, f3point;
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point = LOADFLD4( fld1, shift ); 
					f2point = LOADFLD4( fld2, shift ); 
					f3point = LOADFLD4( fld3, shift ); 

					f1line = vec_madd(  vwx[k1], f1point, f1line );
					f2line = vec_madd( vwxh[k1], f2point, f2line );
					f3line = vec_madd( vwxh[k1], f3point, f3line );
				 }
				 				 
				 f1plane = vec_madd( vwyh[k2], f1line, f1plane  );
				 f2plane = vec_madd(  vwy[k2], f2line, f2plane  );
				 f3plane = vec_madd( vwyh[k2], f3line, f3plane  );
			   }
			   
			   vb1 = vec_madd( vwzh[k3], f1plane, vb1 );
			   vb2 = vec_madd( vwzh[k3], f2plane, vb2 );
			   vb3 = vec_madd(  vwz[k3], f3plane, vb3 );
			}
         
         }		 
      
	  }	  

	  // ---------- advance momenta
	  // Load momenta
	  VEC_LDA4v3( vu1, vu2, vu3, pu );

	  // results are stored in vu1, vu2, vu3 and vene
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store momenta
      VEC_STA4v3(pu, vu1, vu2, vu3);
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg, vq;
		register vector vtrx, vtry, vtrz;
		register vector vx1, vy1, vz1;
		
	    DECLARE_ALIGNED_32( unsigned int cross[4] );
		DECLARE_ALIGNED_32( int dix[4] );
		DECLARE_ALIGNED_32( int diy[4] );
		DECLARE_ALIGNED_32( int diz[4] );

  
		VRGAMMA( vrg, vu1, vu2, vu3 );
	    
	    vq = vec_lda( 0, pq );

		// Accumulate energy
		vene_group = vec_madd( vq, vene, vene_group );
	    
		// get velocities
		vu1 = vec_mul( vu1, vrg );
		vu2 = vec_mul( vu2, vrg );
		vu3 = vec_mul( vu3, vrg );
			
		vx1 = vec_madd( vdt_dx, vu1, vx0 );
		vy1 = vec_madd( vdt_dy, vu2, vy0 );
		vz1 = vec_madd( vdt_dz, vu3, vz0 );
	  
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);
		vtrz = vntrim(vz1);
				
		STOREU4P3D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );

		vx1   = vec_sub( vx1, vtrx );
		vy1   = vec_sub( vy1, vtry );
		vz1   = vec_sub( vz1, vtrz );
        VEC_STA4v3( px, vx1, vy1, vz1 );

		vtrx = vec_ctiw( vtrx );
		vtry = vec_ctiw( vtry );
		vtrz = vec_ctiw( vtrz );
		vec_sta( vtrx, 0, dix );
		vec_sta( vtry, 0, diy );
		vec_sta( vtrz, 0, diz );
				
        for( int l = 0; l < 4; l++) {
          // store new particle cell indexes
		  pix[3*l    ] = vix[l] + dix[l];
		  pix[3*l + 1] = viy[l] + diy[l];
		  pix[3*l + 2] = viz[l] + diz[l];
        
		  // calculate crossings
		  cross[l] = (( dix[l] )?1:0) + (( diy[l] )?2:0) + (( diz[l] )?4:0);
        }
        
		// ---------- split trajectories for current deposition
  
		vsplit3D( &vpbuf, cross, dix, diy, diz );

	  }

	  // advance pointers
	  px  += 3 * VEC_WIDTH;  
	  pix += 3 * VEC_WIDTH;
	  pu  += 3 * VEC_WIDTH; 
	  pq  +=     VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );

	// Accumulate energy from group
	*ene += vec_reduce_add( vene_group ) ;
		
  }
 
}

/****************************************************************************************
   Clear template definitions
****************************************************************************************/

#undef ORDER
#undef ADVANCE_DEPOSIT_1D
#undef ADVANCE_DEPOSIT_2D
#undef ADVANCE_DEPOSIT_2D_CYL
#undef ADVANCE_DEPOSIT_3D
#undef OFFSET
#undef SPLINE
#undef NP

#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D


#endif

