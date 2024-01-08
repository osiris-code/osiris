/*****************************************************************************************

Relativistic particle pusher, AVX optimized version (single precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-avx.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-avx.h"
#include "splines-avx.h"
#include "os-spec-current-avx.h"

#include "split-vec.h"

/*****************************************************************************************
Single pass particle advance deposit

This version of the code does the full push in a single pass: field interpolation, du/dt and dx/dt. 
This allows for savings on memory accesses, since particle positions need only be read once and it 
avoids storing temporary field and momenta values. 

Each group of VEC_WIDTH particles is immediatly split to avoid storing old position values in a special
buffer optimized for the current deposition routines.

The routine then deposits the current for all virtual (split) particles.
*****************************************************************************************/


/*****************************************************************************************
vmntrim

Returns the required shift (trim) for particles that must remain in [ -0.5, 0.5 [
to the nearest cell point:
          x >= 0.5 : +1
  -0.5 <= x <  0.5 :  0
  -0.5 <  x        : -1

*****************************************************************************************/

// Reference implementation
/*
inline __m256 vmntrim( const __m256 vx )
{
   register __m256 va, vb;
 
   va = _mm256_cmp_ps( vx, _mm256_set1_ps( -0.5f ), _CMP_LT_OS );
   va = _mm256_and_ps( va, _mm256_set1_ps( -1.0f ) );
   
   vb = _mm256_cmp_ps( vx, _mm256_set1_ps( +0.5f ), _CMP_GE_OS );
   vb = _mm256_and_ps( vb, _mm256_set1_ps( +1.0f ) );
   
   return _mm256_add_ps( va, vb );
}
*/

// This version uses one less constant
static inline __m256 vmntrim( const __m256 vx )
{
  __m256 va, vb;

  va = _mm256_cmp_ps( vx, _mm256_set1_ps( -0.5f ), _CMP_LT_OS );
  va = _mm256_and_ps( va, _mm256_set1_ps( +1.0f ) );

  vb = _mm256_cmp_ps( vx, _mm256_set1_ps( +0.5f ), _CMP_GE_OS );
  vb = _mm256_and_ps( vb, _mm256_set1_ps( +1.0f ) );

  return  _mm256_sub_ps( vb, va );
}

/*****************************************************************************************
LOADFLD8

Loads 8 field values corresponding to 8 particle positions into a vector
variable. Using _mm_set_ps yields the most efficient code in all scenarios

*****************************************************************************************/


#define LOADFLD8( fp, shift )                                                            \
  _mm256_set_ps( (fp[7])[shift], (fp[6])[shift], (fp[5])[shift], (fp[4])[shift],         \
                 (fp[3])[shift], (fp[2])[shift], (fp[1])[shift], (fp[0])[shift] )

/*****************************************************************************************
vhp_shift
  
Gets the required h shift and cell index for interpolating quantities on a staggered grid
(quantities located at the half cell position).

  h  = ( dx < 0 ) ? 1.0 : 0.0
  ixh = ix - h * delta

*****************************************************************************************/

#define vhp_shift(dx,ix, h,ixh,delta ) {						              \
   register __m256 cmp;                                           \
                                                                  \
   cmp  = _mm256_cmp_ps( dx, _mm256_setzero_ps(), _CMP_LT_OS );   \
   h    =  _mm256_and_ps( cmp, _mm256_set1_ps(  1.0f ) );         \
   ixh = visub( ix, _mm256_castps_si256(                          \
					   _mm256_and_ps(                                       \
						  cmp,                                                \
						  _mm256_castsi256_ps( _mm256_set1_epi32( delta ) )   \
					   )                                                    \
					) );                                                    \
}


/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene )        \
{                                                                                        \
   __m256 const c1 = _mm256_set1_ps( 1.0f ); \
   __m256 vut1, vut2, vut3, vutsq;\
   __m256 vgamma, vtem_gamma, votsq;\
   \
   /* Perform first half of electric field acceleration.*/\
   ve1 = _mm256_mul_ps( ve1, vtem );\
   ve2 = _mm256_mul_ps( ve2, vtem );\
   ve3 = _mm256_mul_ps( ve3, vtem );\
   \
   vut1 = _mm256_add_ps( vu1, ve1 );\
   vut2 = _mm256_add_ps( vu2, ve2 );\
   vut3 = _mm256_add_ps( vu3, ve3 );\
   \
   /* Perform first half of the rotation and store in u. */\
   vutsq = _mm256_add_ps( _mm256_mul_ps( vut3, vut3 ),\
                          _mm256_add_ps( _mm256_mul_ps( vut2, vut2 ),\
                                         _mm256_mul_ps( vut1, vut1 ) ) );\
   \
   vgamma = _mm256_sqrt_ps( _mm256_add_ps( vutsq, c1 ) );\
   \
   vene = _mm256_div_ps( vutsq, _mm256_add_ps( vgamma, c1 ) );\
   \
   vtem_gamma = _mm256_div_ps( vtem, vgamma );\
   \
   vb1 = _mm256_mul_ps( vb1, vtem_gamma );\
   vb2 = _mm256_mul_ps( vb2, vtem_gamma );\
   vb3 = _mm256_mul_ps( vb3, vtem_gamma );\
   \
   vu1 = _mm256_add_ps( _mm256_mul_ps(vb3, vut2), vut1 );\
   vu2 = _mm256_add_ps( _mm256_mul_ps(vb1, vut3), vut2 );\
   vu3 = _mm256_add_ps( _mm256_mul_ps(vb2, vut1), vut3 );\
   \
   vu1 = _mm256_sub_ps( vu1, _mm256_mul_ps(vb2, vut3)  );\
   vu2 = _mm256_sub_ps( vu2, _mm256_mul_ps(vb3, vut1)  );\
   vu3 = _mm256_sub_ps( vu3, _mm256_mul_ps(vb1, vut2)  );\
   \
   /* Perform second half of the rotation. */\
  votsq = _mm256_add_ps( _mm256_mul_ps(vb1, vb1),    c1 );\
  votsq = _mm256_add_ps( _mm256_mul_ps(vb2, vb2), votsq );\
  votsq = _mm256_add_ps( _mm256_mul_ps(vb3, vb3), votsq );\
  \
  votsq = _mm256_div_ps( _mm256_set1_ps( 2.0f ), votsq ) ;\
   \
   vb1 = _mm256_mul_ps( vb1, votsq );\
   vb2 = _mm256_mul_ps( vb2, votsq );\
   vb3 = _mm256_mul_ps( vb3, votsq );\
   \
  vut1 = _mm256_add_ps( _mm256_mul_ps(vb3, vu2), vut1 );\
  vut2 = _mm256_add_ps( _mm256_mul_ps(vb1, vu3), vut2 );\
  vut3 = _mm256_add_ps( _mm256_mul_ps(vb2, vu1), vut3 );\
  \
  vut1 = _mm256_sub_ps( vut1, _mm256_mul_ps(vb2, vu3) );\
  vut2 = _mm256_sub_ps( vut2, _mm256_mul_ps(vb3, vu1) );\
  vut3 = _mm256_sub_ps( vut3, _mm256_mul_ps(vb1, vu2) );\
  \
   /* Perform second half of electric field acceleration.*/\
   vu1 = _mm256_add_ps( vut1, ve1 );\
   vu2 = _mm256_add_ps( vut2, ve2 );\
   vu3 = _mm256_add_ps( vut3, ve3 );\
\
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. Using a macro is beneficial.

Note: Need to check use of rsqrt_ps
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                                                  \
  __m256 const c1 = _mm256_set1_ps( 1.0f );                                              \
  (vrg) = _mm256_add_ps( _mm256_mul_ps((vu3),(vu3)),                                     \
                         _mm256_add_ps(_mm256_mul_ps((vu2),(vu2)),                       \
                                       _mm256_add_ps(_mm256_mul_ps((vu1),(vu1)),c1 )));  \
  (vrg) = _mm256_div_ps( c1, _mm256_sqrt_ps(vrg) );                                      \
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

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s1 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s1 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s1 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s1 )

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Grid offset ( 1 + ORDER/2 )
#define OFFSET 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s2 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s2 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s2 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s2 )

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s3 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s3 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s3 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s3 )

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s4 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s4 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s4 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s4 )

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
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;

  DECLARE_ALIGNED_32( float jnorm[1] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf1D vpbuf;

  // Simulation constants
  __m256 const vtem    = _mm256_set1_ps( (float) (0.5 * (*dt) / (*rqm)) );
  __m256 const vdt_dx1 = _mm256_set1_ps( (float) ((*dt) / dx[0]) );

  // get pointers to position 0 of each field component
  unsigned int const deltaX = 3; // 3 field components

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;


  // Normalization for currents
  jnorm[0] = (float) (dx[0] / (*dt));

  // jump to 1st particle
  x  += (*i0-1);
  ix += (*i0-1);
  u  += 3*(*i0-1);
  q  += (*i0-1);

  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 8 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

    for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
      ix[ (np_total + i)     ] = 1;
      x[  (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i) + 1 ] = 0.;
      u[ 3 * (np_total + i) + 2 ] = 0.;
      q[ (np_total + i ) ] = 0.;
    }

    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
      printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
      exit(1);
    }

  }

  if ( p_cache_size_1D % VEC_WIDTH ) {
    printf("(*error*) p_cache_size_1D is not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
    exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_1D ) {

    // Initialize energy
    __m256 vene_group = _mm256_setzero_ps();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_1D ) ? np_total - k : p_cache_size_1D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i += VEC_WIDTH ) {
    __m256 vene;
    __m256 vx0;
    __m256i vix0;

    __m256 ve1, ve2, ve3;
    __m256 vb1, vb2, vb3;

    __m256 vu1, vu2, vu3;

    // Load particle positions 
    vx0  = _mm256_load_ps( (float const*) x ); 
    vix0 = _mm256_load_si256( (__m256i const* ) ix );

    // Interpolate fields
    {
      __m256  hx;
      __m256i vidxi, vidxih;
      __m256 vwx[NP], vwxh[NP];

      // idxi = 3 * ix0
      vidxi = viadd3( vix0, vix0, vix0 );

      SPLINE( vx0, vwx );

      vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 

      SPLINEH( vx0, hx, vwxh );

      // Interpolate E field
      {
        unsigned int k1;

        float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
        DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
        DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
        DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

        // get pointers to fields 
        _mm256_store_si256( (__m256i *) idx1, vidxih );
        _mm256_store_si256( (__m256i *) idx2,  vidxi );
        _mm256_store_si256( (__m256i *) idx3,  vidxi );

        // This cannot be done vectorially since pointers are 64 bit
        for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
          fld1[k1] = e1 + idx1[k1];
          fld2[k1] = e2 + idx2[k1];
          fld3[k1] = e3 + idx3[k1];
        }

        ve1 = _mm256_setzero_ps();
        ve2 = _mm256_setzero_ps();
        ve3 = _mm256_setzero_ps();

        __m256 f1point[NP], f2point[NP], f3point[NP];

        for ( k1 = 0; k1 < NP; k1++ ) {
          unsigned int shift = k1*3;

          f1point[k1] = LOADFLD8( fld1, shift ); 
          f2point[k1] = LOADFLD8( fld2, shift ); 
          f3point[k1] = LOADFLD8( fld3, shift ); 
        }

        for ( k1 = 0; k1 < NP; k1 ++ ) {
          ve1 = _mm256_add_ps( ve1, _mm256_mul_ps( f1point[k1], vwxh[k1] ) );
          ve2 = _mm256_add_ps( ve2, _mm256_mul_ps( f2point[k1],  vwx[k1] ) );
          ve3 = _mm256_add_ps( ve3, _mm256_mul_ps( f3point[k1],  vwx[k1] ) );
        }

      }     

      // Interpolate B field
      {
        unsigned int k1;

        float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
        DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
        DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
        DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

        // get pointers to fields 
        _mm256_store_si256( (__m256i *) idx1,  vidxi );
        _mm256_store_si256( (__m256i *) idx2, vidxih );
        _mm256_store_si256( (__m256i *) idx3, vidxih );

        // This cannot be done vectorially since pointers are 64 bit
        for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
          fld1[k1] = b1 + idx1[k1];
          fld2[k1] = b2 + idx2[k1];
          fld3[k1] = b3 + idx3[k1];
        }

        vb1 = _mm256_setzero_ps();
        vb2 = _mm256_setzero_ps();
        vb3 = _mm256_setzero_ps();

        __m256 f1point[NP], f2point[NP], f3point[NP];

        for ( k1 = 0; k1 < NP; k1++ ) {
          unsigned int shift = k1*3;

          f1point[k1] = LOADFLD8( fld1, shift ); 
          f2point[k1] = LOADFLD8( fld2, shift ); 
          f3point[k1] = LOADFLD8( fld3, shift ); 
        }

        for ( k1 = 0; k1 < NP; k1 ++ ) {
          vb1 = _mm256_add_ps( vb1, _mm256_mul_ps( f1point[k1],  vwx[k1] ) );
          vb2 = _mm256_add_ps( vb2, _mm256_mul_ps( f2point[k1], vwxh[k1] ) );
          vb3 = _mm256_add_ps( vb3, _mm256_mul_ps( f3point[k1], vwxh[k1] ) );
        }

      }     
    }

    // load momenta
    _MM256_LOAD8v3_PS( vu1, vu2, vu3, u );

    // advance momenta
    VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

    // Store Results
    _MM256_STORE8v3_PS(u, vu1, vu2, vu3); 

    // ---------- advance position and get velocities
    {
      __m256 vrg;

      // Get 1 / \gamma
      VRGAMMA( vrg, vu1, vu2, vu3 );

      // Load charge
      __m256 vq = _mm256_load_ps( (float const*) q );

      // Accumulate energy
      // Energy for each group is accumulated in single precision
      vene_group = _mm256_add_ps(_mm256_mul_ps( vq, vene ),  vene_group  );

      // get velocities
      __m256 vv1 = _mm256_mul_ps( vu1, vrg );
      __m256 vv2 = _mm256_mul_ps( vu2, vrg ); // this is required for current deposition
      __m256 vv3 = _mm256_mul_ps( vu3, vrg ); // this is required for current deposition

      // Push positions
      __m256 vx1 = _mm256_add_ps( _mm256_mul_ps( vv1, vdt_dx1 ), vx0 );

      // Store virtual particles with positions still indexed to the original cell
      STOREU8P1D( vpbuf, vpbuf.np, vx0, vx1, vq, vv2, vv3, vix0 );

      // Trim positions and store results
      __m256 vtr1 = vmntrim(vx1);

      vx1  = _mm256_sub_ps( vx1, vtr1 );
      _mm256_store_ps( x, vx1 );

      __m256i vitr1 = _mm256_cvttps_epi32( vtr1 ) ;
      _mm256_store_si256( (__m256i *) dix, vitr1 );

      __m256i vix1 = viadd( vix0, vitr1 ); 
      _mm256_store_si256( (__m256i *) ix, vix1 );
    }

    // ---------- split trajectories for current deposition
    vsplit1D( &vpbuf, dix );

    // ---------- advance pointers
    x  += VEC_WIDTH;  
    ix += VEC_WIDTH;
    u  += 3 * VEC_WIDTH; 
    q  += VEC_WIDTH;
    }

    // Deposit current from all virtual particles
    DEP_CURRENT_1D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    // this is done in double precision
    *ene += _mm256_reduce_add_ps( vene_group );

  }

}

/********************************** 2D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_2D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;
    
  DECLARE_ALIGNED_32( float jnorm[2] );
  
  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m256 const vtem    = _mm256_set1_ps( (float) (0.5 * (*dt) / (*rqm)) );
  __m256 const vdt_dx1 = _mm256_set1_ps( (float) ((*dt) / dx[0]) );
  __m256 const vdt_dx2 = _mm256_set1_ps( (float) ((*dt) / dx[1]) );

  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;
  
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;
  
  
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
  	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH;  i ++ ) {
  	  ix[ 2 * (np_total + i)     ] = 1;
  	  ix[ 2 * (np_total + i) + 1 ] = 1;
  	  x[ 2 * (np_total + i)     ] = 0.;
  	  x[ 2 * (np_total + i) + 1 ] = 0.;
  	  u[ 3 * (np_total + i)     ] = 0.;
  	  u[ 3 * (np_total + i) + 1 ] = 0.;
  	  u[ 3 * (np_total + i) + 2 ] = 0.;
  	  q[ (np_total + i ) ] = 0.;
  	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D )  {

    // Initialize energy
    __m256 vene_group = _mm256_setzero_ps();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m256 vene;
      __m256 vx0, vy0;
      __m256i vix0, viy0;

      __m256 ve1, ve2, ve3;
      __m256 vb1, vb2, vb3;

      __m256 vu1, vu2, vu3;

      // Load particle positions 
      _MM256_LOAD8v2_PS( vx0, vy0, x );
      _MM256_LOAD8v2_EPI32( vix0, viy0, ix  );

      // Interpolate fields
      {
        __m256  hx, hy;
        __m256i vidxi, vidxj, vidxih, vidxjh;
        __m256 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = viadd3( vix0, vix0, vix0 );
        // idxj = deltaY * iy0
        vidxj = vimul_s( viy0, deltaY );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
        vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd_store( idx1, vidxih,  vidxj );
          viadd_store( idx2,  vidxi, vidxjh );
          viadd_store( idx3,  vidxi,  vidxj );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = e1 + idx1[k1];
            fld2[k1] = e2 + idx2[k1];
            fld3[k1] = e3 + idx3[k1];
          }

          ve1 = _mm256_setzero_ps();
          ve2 = _mm256_setzero_ps();
          ve3 = _mm256_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256 f1line, f2line, f3line;
            f1line = _mm256_setzero_ps();
            f2line = _mm256_setzero_ps();
            f3line = _mm256_setzero_ps();

            __m256 f1point[NP], f2point[NP], f3point[NP];

            for ( k1 = 0; k1 < NP; k1++ ) {
              unsigned int shift = k1*3 + k2*deltaY;

              f1point[k1] = LOADFLD8( fld1, shift ); 
              f2point[k1] = LOADFLD8( fld2, shift ); 
              f3point[k1] = LOADFLD8( fld3, shift ); 
            }

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1], vwxh[k1] ) );
              f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1],  vwx[k1] ) );
              f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1],  vwx[k1] ) );
            }

            ve1 = _mm256_add_ps( ve1, _mm256_mul_ps( f1line,  vwy[k2] ) );
            ve2 = _mm256_add_ps( ve2, _mm256_mul_ps( f2line, vwyh[k2] ) );
            ve3 = _mm256_add_ps( ve3, _mm256_mul_ps( f3line,  vwy[k2] ) );
          }

        }		 

        // Interpolate B field
        {
          unsigned int k1, k2;

          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd_store( idx1,  vidxi, vidxjh );
          viadd_store( idx2, vidxih,  vidxj );
          viadd_store( idx3, vidxih, vidxjh );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = b1 + idx1[k1];
            fld2[k1] = b2 + idx2[k1];
            fld3[k1] = b3 + idx3[k1];
          }

          vb1 = _mm256_setzero_ps();
          vb2 = _mm256_setzero_ps();
          vb3 = _mm256_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256 f1line, f2line, f3line;
            f1line = _mm256_setzero_ps();
            f2line = _mm256_setzero_ps();
            f3line = _mm256_setzero_ps();

            __m256 f1point[NP], f2point[NP], f3point[NP];

            for ( k1 = 0; k1 < NP; k1++ ) {
              unsigned int shift = k1*3 + k2*deltaY;

              f1point[k1] = LOADFLD8( fld1, shift ); 
              f2point[k1] = LOADFLD8( fld2, shift ); 
              f3point[k1] = LOADFLD8( fld3, shift ); 
            }

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1],  vwx[k1] ) );
              f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1], vwxh[k1] ) );
              f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1], vwxh[k1] ) );
            }

            vb1 = _mm256_add_ps( vb1, _mm256_mul_ps( f1line, vwyh[k2] ) );
            vb2 = _mm256_add_ps( vb2, _mm256_mul_ps( f2line,  vwy[k2] ) );
            vb3 = _mm256_add_ps( vb3, _mm256_mul_ps( f3line, vwyh[k2] ) );
          }

        }		 
      }

      // Load momenta
      _MM256_LOAD8v3_PS( vu1, vu2, vu3, u );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM256_STORE8v3_PS(u, vu1, vu2, vu3); 

      // ---------- advance position and get velocities
      {
        __m256 vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256 vq = _mm256_load_ps( q );

        // Accumulate energy
        vene_group = _mm256_add_ps( _mm256_mul_ps( vq, vene ), vene_group );

        // get velocities
        __m256 vv1 = _mm256_mul_ps( vu1, vrg );
        __m256 vv2 = _mm256_mul_ps( vu2, vrg );
        __m256 vv3 = _mm256_mul_ps( vu3, vrg ); // this is required for current deposition

        // Push positions
        __m256 vx1 = _mm256_add_ps( _mm256_mul_ps( vv1, vdt_dx1 ), vx0 );
        __m256 vy1 = _mm256_add_ps( _mm256_mul_ps( vv2, vdt_dx2 ), vy0 );

        // Store virtual particles with positions still indexed to the original cell
        STOREU8P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim
        __m256 vtr1 = vmntrim(vx1);
        __m256 vtr2 = vmntrim(vy1);

        // Calculate crossings and store result
        {
          __m256 const zero = _mm256_setzero_ps();
          __m256 vcross;

          // x crossings
          vcross = _mm256_and_ps( _mm256_cmp_ps(vtr1, zero, _CMP_NEQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(1) ));

          // y crossings
          vcross = _mm256_or_ps( vcross, 
          _mm256_and_ps( _mm256_cmp_ps(vtr2, zero, _CMP_NEQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(2) )));

          // Store vcross
          _mm256_store_si256( (__m256i *) cross, _mm256_castps_si256( vcross ) );
        }

        // Trim positions and store results
        vx1  = _mm256_sub_ps( vx1, vtr1 );
        vy1  = _mm256_sub_ps( vy1, vtr2 );
        _MM256_STORE8v2_PS( x, vx1, vy1 );

        // find cell crossings and store
        __m256i vitr1 = _mm256_cvtps_epi32( vtr1 ) ;
        __m256i vitr2 = _mm256_cvtps_epi32( vtr2 ) ;

        _mm256_store_si256( (__m256i *) dix, vitr1 );
        _mm256_store_si256( (__m256i *) diy, vitr2 );

        // Trim cell indexes and store
        __m256i vix1 = viadd( vix0, vitr1 ); 
        __m256i viy1 = viadd( viy0, vitr2 ); 
        _MM256_STORE8v2_EPI32( ix, vix1, viy1 );
      }


      // ---------- split trajectories for current deposition
      vsplit2D( &vpbuf, cross, dix, diy );


      // ---------- advance pointers
      x  += 2 * VEC_WIDTH;
      ix += 2 * VEC_WIDTH;
      u  += 3 * VEC_WIDTH;
      q  += VEC_WIDTH;

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_ps( vene_group );

  }
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  int *ilb2 , 
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;

  DECLARE_ALIGNED_32( float jnorm[2] );

  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m256 const vtem    = _mm256_set1_ps( (float) (0.5 * (*dt) / (*rqm)) );
  __m256 const vdt_dx1 = _mm256_set1_ps( (float) ((*dt) / dx[0]) );

  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
  (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
  (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;


  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;

  #if ( ORDER == 1 ) || ( ORDER == 3 )
  __m256d const gshift2 = _mm256_set1_pd( (*ilb2 - 2 ) );
  #else
  __m256d const gshift2 = _mm256_set1_pd( (*ilb2 - 2 ) - 0.5 );
  #endif

  __m256d const dr      = _mm256_set1_pd( dx[1] );
  __m256d const rdr     = _mm256_set1_pd( 1.0/dx[1] );
  __m256d const vdt     = _mm256_set1_pd( *dt );

  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;

  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

    for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH;  i ++ ) {
      ix[ 2 * (np_total + i)     ] = 1;
      ix[ 2 * (np_total + i) + 1 ] = 1;
      x[ 2 * (np_total + i)     ] = 0.;
      x[ 2 * (np_total + i) + 1 ] = 0.25f; // Avoid the axis
      u[ 3 * (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i) + 1 ] = 0.;
      u[ 3 * (np_total + i) + 2 ] = 0.;
      q[ (np_total + i ) ] = 0.;
    }

    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
      printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
      exit(1);
    }

  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
    printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
    exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

    // Initialize energy
    __m256 vene_group = _mm256_setzero_ps();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m256 vene;
      __m256 vx0, vy0;
      __m256i vix0, viy0;

      __m256 ve1, ve2, ve3;
      __m256 vb1, vb2, vb3;

      __m256 vu1, vu2, vu3;

      // Load particle positions 
      _MM256_LOAD8v2_PS( vx0, vy0, x );
      _MM256_LOAD8v2_EPI32( vix0, viy0, ix  );

      // Interpolate fields
      {
        __m256  hx, hy;
        __m256i vidxi, vidxj, vidxih, vidxjh;
        __m256 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = viadd3( vix0, vix0, vix0 );
        // idxj = deltaY * iy0
        vidxj = vimul_s( viy0, deltaY );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
        vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd_store( idx1, vidxih,  vidxj );
          viadd_store( idx2, vidxi, vidxjh );
          viadd_store( idx3, vidxi,  vidxj );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = e1 + idx1[k1];
            fld2[k1] = e2 + idx2[k1];
            fld3[k1] = e3 + idx3[k1];
          }

          ve1 = _mm256_setzero_ps();
          ve2 = _mm256_setzero_ps();
          ve3 = _mm256_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256 f1line, f2line, f3line;
            f1line = _mm256_setzero_ps();
            f2line = _mm256_setzero_ps();
            f3line = _mm256_setzero_ps();

            __m256 f1point[NP], f2point[NP], f3point[NP];

            for ( k1 = 0; k1 < NP; k1++ ) {
              unsigned int shift = k1*3 + k2*deltaY;

              f1point[k1] = LOADFLD8( fld1, shift ); 
              f2point[k1] = LOADFLD8( fld2, shift ); 
              f3point[k1] = LOADFLD8( fld3, shift ); 
            }

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1], vwxh[k1] ) );
              f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1],  vwx[k1] ) );
              f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1],  vwx[k1] ) );
            }

            ve1 = _mm256_add_ps( ve1, _mm256_mul_ps( f1line,  vwy[k2] ) );
            ve2 = _mm256_add_ps( ve2, _mm256_mul_ps( f2line, vwyh[k2] ) );
            ve3 = _mm256_add_ps( ve3, _mm256_mul_ps( f3line,  vwy[k2] ) );
          }

        }		 

        // Interpolate B field
        {
          unsigned int k1, k2;

          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd_store( idx1, vidxi,  vidxjh );
          viadd_store( idx2, vidxih,  vidxj );
          viadd_store( idx3, vidxih, vidxjh );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = b1 + idx1[k1];
            fld2[k1] = b2 + idx2[k1];
            fld3[k1] = b3 + idx3[k1];
          }

          vb1 = _mm256_setzero_ps();
          vb2 = _mm256_setzero_ps();
          vb3 = _mm256_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256 f1line, f2line, f3line;
            f1line = _mm256_setzero_ps();
            f2line = _mm256_setzero_ps();
            f3line = _mm256_setzero_ps();

            __m256 f1point[NP], f2point[NP], f3point[NP];

            for ( k1 = 0; k1 < NP; k1++ ) {
              unsigned int shift = k1*3 + k2*deltaY;

              f1point[k1] = LOADFLD8( fld1, shift ); 
              f2point[k1] = LOADFLD8( fld2, shift ); 
              f3point[k1] = LOADFLD8( fld3, shift ); 
            }

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1],  vwx[k1] ) );
              f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1], vwxh[k1] ) );
              f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1], vwxh[k1] ) );
            }

            vb1 = _mm256_add_ps( vb1, _mm256_mul_ps( f1line, vwyh[k2] ) );
            vb2 = _mm256_add_ps( vb2, _mm256_mul_ps( f2line,  vwy[k2] ) );
            vb3 = _mm256_add_ps( vb3, _mm256_mul_ps( f3line, vwyh[k2] ) );
          }

        }		 

      }

      // Load momenta
      _MM256_LOAD8v3_PS( vu1, vu2, vu3, u );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );


      // ---------- advance position in cylindrical geometry
      {
        __m256 vrg;
        __m256 vx1, vy1;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256 vq = _mm256_load_ps( q );

        // Accumulate energy
        vene_group = _mm256_add_ps( _mm256_mul_ps( vq, vene ), vene_group );

        // get velocities
        __m256 vv1 = _mm256_mul_ps( vu1, vrg );
        __m256 vv2 = _mm256_mul_ps( vu2, vrg );
        __m256 vv3 = _mm256_mul_ps( vu3, vrg ); // this is required for current deposition

        // Push positions
        vx1 = _mm256_add_ps( _mm256_mul_ps( vv1, vdt_dx1 ), vx0 );

        {  // This section needs to be performed in double precision

          __m256d gix2a, gix2b;
          __m256d vy0a, vy0b, vy1a, vy1b;
          __m256d r_olda, r_oldb, r_newa, r_newb;
          __m256d x2_newa, x2_newb, x3_newa, x3_newb; 

          __m256d vv2a, vv2b, vv3a, vv3b;
          __m256d vu2a, vu2b, vu3a, vu3b;

          __m256d tmpa, tmpb;

          // gix2   = viy0 + gshift2;          // global cell
          gix2a = _mm256_add_pd( _mm256_cvtepi32_pd( _mm256_castsi256_si128(viy0) ), gshift2 );
          gix2b = _mm256_add_pd( _mm256_cvtepi32_pd( _mm256_extractf128_si256( viy0, 1 ) ), gshift2 );

          // r_old = ( vy0 + gix2 ) * dr;      // global radial position
          _MM256_CVTPS_PD( vy0a, vy0b, vy0 );
          r_olda = _mm256_mul_pd( _mm256_add_pd( vy0a, gix2a ), dr );
          r_oldb = _mm256_mul_pd( _mm256_add_pd( vy0b, gix2b ), dr );

          // x2_new = r_old + vv2 * dt;
          // x3_new =         vv3 * dt;

          _MM256_CVTPS_PD( vv2a, vv2b, vv2 );
          _MM256_CVTPS_PD( vv3a, vv3b, vv3 );

          x2_newa = _mm256_add_pd( _mm256_mul_pd( vv2a, vdt ), r_olda  );
          x3_newa =                     _mm256_mul_pd( vv3a, vdt ) ;

          x2_newb = _mm256_add_pd( _mm256_mul_pd( vv2b, vdt ), r_oldb  );
          x3_newb =                     _mm256_mul_pd( vv3b, vdt );

          // r_new = sqrt( x2_new*x2_new + x3_new*x3_new );
          r_newa = _mm256_sqrt_pd( _mm256_add_pd( _mm256_mul_pd( x2_newa, x2_newa ),
          _mm256_mul_pd( x3_newa, x3_newa ) ) ); 

          r_newb = _mm256_sqrt_pd( _mm256_add_pd( _mm256_mul_pd( x2_newb, x2_newb ),
          _mm256_mul_pd( x3_newb, x3_newb ) ) ); 

          // This is a protection against roundoff for cold plasmas
          // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
          vy1a = _mm256_blendv_pd( _mm256_sub_pd( _mm256_mul_pd( r_newa, rdr ), gix2a ), vy0a, 
                                  _mm256_cmp_pd( r_olda, r_newa, _CMP_EQ_OQ ) );
          vy1b = _mm256_blendv_pd( _mm256_sub_pd( _mm256_mul_pd( r_newb, rdr ), gix2b ), vy0b, 
                                  _mm256_cmp_pd( r_oldb, r_newb, _CMP_EQ_OQ ) );

          // Convert vy1 to single precision
          _MM256_CVTPD_PS( vy1, vy1a, vy1b );

          // Correct p_r and p_\theta to conserve angular momentum
          // tmp = 1.0 / r_new;
          // vu2 = ( vu2 * x2_new + vu3 * x3_new ) * tmp;
          // vu3 = vu3 * r_old * tmp;

          // Convert vu2, vu3 to double precision
          _MM256_CVTPS_PD( vu2a, vu2b, vu2 );
          _MM256_CVTPS_PD( vu3a, vu3b, vu3 );

          tmpa  = _mm256_div_pd( _mm256_set1_pd(1.0), r_newa );
          tmpb  = _mm256_div_pd( _mm256_set1_pd(1.0), r_newb );

          vu2a  = _mm256_mul_pd( _mm256_add_pd( _mm256_mul_pd( vu2a, x2_newa ), 
                                                _mm256_mul_pd( vu3a, x3_newa ) ), 
                                                tmpa );
          vu2b = _mm256_mul_pd( _mm256_add_pd( _mm256_mul_pd( vu2b, x2_newb ), 
                                                _mm256_mul_pd( vu3b, x3_newb ) ), 
                                                tmpb );

          vu3a = _mm256_mul_pd( vu3a, _mm256_mul_pd( r_olda, tmpa ) );
          vu3b = _mm256_mul_pd( vu3b, _mm256_mul_pd( r_oldb, tmpb ) );

          // Convert vu2, vu3 to single precision
          _MM256_CVTPD_PS( vu2, vu2a, vu2b );
          _MM256_CVTPD_PS( vu3, vu3a, vu3b );
        }

        // Store Momenta
        _MM256_STORE8v3_PS(u, vu1, vu2, vu3); 

        // Store virtual particles with positions still indexed to the original cell
        STOREU8P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim
        __m256 vtr1 = vmntrim(vx1);
        __m256 vtr2 = vmntrim(vy1);

        // Calculate crossings and store result
        {
          const __m256 zero = _mm256_setzero_ps();
          __m256 vcross;

          // x crossings
          vcross = _mm256_andnot_ps( _mm256_cmp_ps(vtr1, zero, _CMP_EQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(1) ));

          // y crossings
          vcross = _mm256_or_ps( vcross, 
          _mm256_andnot_ps( _mm256_cmp_ps(vtr2, zero, _CMP_EQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(2) )));

          // Store vcross
          _mm256_store_si256( (__m256i *) cross, _mm256_castps_si256( vcross ) );
        }

        // Trim positions and store results
        vx1  = _mm256_sub_ps( vx1, vtr1 );
        vy1  = _mm256_sub_ps( vy1, vtr2 );
        _MM256_STORE8v2_PS( x, vx1, vy1 );

        // find cell crossings and store
        __m256i vitr1 = _mm256_cvtps_epi32( vtr1 ) ;
        __m256i vitr2 = _mm256_cvtps_epi32( vtr2 ) ;

        _mm256_store_si256( (__m256i *) dix, vitr1 );
        _mm256_store_si256( (__m256i *) diy, vitr2 );

        // Trim cell indexes and store
        __m256i vix1 = viadd( vix0, vitr1 ); 
        __m256i viy1 = viadd( viy0, vitr2 ); 
        _MM256_STORE8v2_EPI32( ix, vix1, viy1 );

      }

      // ---------- split trajectories for current deposition

      // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );


      // ---------- advance pointers
      x  += 2 * VEC_WIDTH;  
      u  += 3 * VEC_WIDTH; 
      ix += 2 * VEC_WIDTH;
      q  += VEC_WIDTH;

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_ps( vene_group );

  }
    
}


/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;

  DECLARE_ALIGNED_32( float jnorm[3] );

  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diz[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf3D vpbuf;

  // Simulation constants
  __m256 const vtem    = _mm256_set1_ps( (float) (0.5 * (*dt) / (*rqm)) );
  __m256 const vdt_dx1 = _mm256_set1_ps( (float) ((*dt) / dx[0]) );
  __m256 const vdt_dx2 = _mm256_set1_ps( (float) ((*dt) / dx[1]) );
  __m256 const vdt_dx3 = _mm256_set1_ps( (float) ((*dt) / dx[2]) );

  // get pointers to position 0,0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;
  unsigned int const deltaZ = emf_size[1] * deltaY;

  float* const e1 = efield + (emf_offset[0]-OFFSET)*deltaX + 
  (emf_offset[1]-OFFSET)*deltaY + 
  (emf_offset[2]-OFFSET)*deltaZ;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0]-OFFSET)*deltaX + 
  (emf_offset[1]-OFFSET)*deltaY + 
  (emf_offset[2]-OFFSET)*deltaZ;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;

  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 3;
  jnorm[1] = dx[1] / (*dt) / 3;
  jnorm[2] = dx[2] / (*dt) / 3;

  // jump to 1st particle
  x  += 3*(*i0-1);
  ix += 3*(*i0-1);
  u  += 3*(*i0-1);
  q  +=   (*i0-1);


  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

    for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
      ix[ 3 * (np_total + i)     ] = 1;
      ix[ 3 * (np_total + i) + 1 ] = 1;
      ix[ 3 * (np_total + i) + 2 ] = 1;
      x[ 3 * (np_total + i)     ] = 0.;
      x[ 3 * (np_total + i) + 1 ] = 0.;
      x[ 3 * (np_total + i) + 2 ] = 0.;
      u[ 3 * (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i) + 1 ] = 0.;
      u[ 3 * (np_total + i) + 2 ] = 0.;
      q[ (np_total + i ) ] = 0.;
    }

    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
      printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
      exit(1);
    }

  }

  if ( p_cache_size_3D % VEC_WIDTH ) {
    printf("(*error*) p_cache_size_3D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH );
    exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {

    // Initialize energy
    __m256 vene_group = _mm256_setzero_ps();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m256 vene;
      __m256 vx0, vy0, vz0;
      __m256i vix0, viy0, viz0;

      __m256 ve1, ve2, ve3;
      __m256 vb1, vb2, vb3;

      __m256 vu1, vu2, vu3;

      // Load particle positions 
      _MM256_LOAD8v3_PS( vx0, vy0, vz0, x );
      _MM256_LOAD8v3_EPI32( vix0, viy0, viz0, ix  );

      // Interpolate fields
      {
        __m256  hx, hy, hz;
        __m256i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
        __m256 vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

        // idxi = 3 * ix0
        vidxi = viadd3( vix0, vix0, vix0 );
        // idxj = deltaY * iy0
        vidxj = vimul_s( viy0, deltaY );
        // idxk = deltaZ * iz0
        vidxk = vimul_s( viz0, deltaZ );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy );
        SPLINE( vz0, vwz );

        vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
        vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 
        vhp_shift( vz0, vidxk, hz, vidxkh, deltaZ );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );
        SPLINEH( vz0, hz, vwzh );

        // Interpolate E field
        {
          unsigned int k1, k2, k3;
          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd3_store( idx1, vidxih, vidxj, vidxk  );
          viadd3_store( idx2, vidxi, vidxjh, vidxk  );
          viadd3_store( idx3, vidxi,  vidxj, vidxkh );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = e1 + idx1[k1];
            fld2[k1] = e2 + idx2[k1];
            fld3[k1] = e3 + idx3[k1];
          }

          ve1 = _mm256_setzero_ps();
          ve2 = _mm256_setzero_ps();
          ve3 = _mm256_setzero_ps();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m256 f1plane, f2plane, f3plane;
            f1plane = _mm256_setzero_ps();
            f2plane = _mm256_setzero_ps();
            f3plane = _mm256_setzero_ps();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m256 f1line, f2line, f3line;
              f1line = _mm256_setzero_ps();
              f2line = _mm256_setzero_ps();
              f3line = _mm256_setzero_ps();

              __m256 f1point[NP], f2point[NP], f3point[NP];
              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned int shift = k1*3 + k2*deltaY + k3*deltaZ;

                f1point[k1] = LOADFLD8( fld1, shift ); 
                f2point[k1] = LOADFLD8( fld2, shift ); 
                f3point[k1] = LOADFLD8( fld3, shift ); 
              }

              for ( k1 = 0; k1 < NP; k1 ++ ) {
                f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1], vwxh[k1] ) );
                f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1],  vwx[k1] ) );
                f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1],  vwx[k1] ) );
              }

              f1plane = _mm256_add_ps( f1plane, _mm256_mul_ps( f1line,  vwy[k2] ) );
              f2plane = _mm256_add_ps( f2plane, _mm256_mul_ps( f2line, vwyh[k2] ) );
              f3plane = _mm256_add_ps( f3plane, _mm256_mul_ps( f3line,  vwy[k2] ) );
            }

            ve1 = _mm256_add_ps( ve1, _mm256_mul_ps( f1plane,  vwz[k3] ) );
            ve2 = _mm256_add_ps( ve2, _mm256_mul_ps( f2plane,  vwz[k3] ) );
            ve3 = _mm256_add_ps( ve3, _mm256_mul_ps( f3plane, vwzh[k3] ) );
          }

        }		 

        // Interpolate B field
        {
          unsigned int k1, k2, k3;
          float *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_32( int idx1[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx2[VEC_WIDTH] );
          DECLARE_ALIGNED_32( int idx3[VEC_WIDTH] );

          // get pointers to fields 
          viadd3_store( idx1, vidxi,  vidxjh, vidxkh );
          viadd3_store( idx2, vidxih,  vidxj, vidxkh );
          viadd3_store( idx3, vidxih, vidxjh,  vidxk );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = b1 + idx1[k1];
            fld2[k1] = b2 + idx2[k1];
            fld3[k1] = b3 + idx3[k1];
          }

          vb1 = _mm256_setzero_ps();
          vb2 = _mm256_setzero_ps();
          vb3 = _mm256_setzero_ps();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m256 f1plane, f2plane, f3plane;
            f1plane = _mm256_setzero_ps();
            f2plane = _mm256_setzero_ps();
            f3plane = _mm256_setzero_ps();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m256 f1line, f2line, f3line;
              f1line = _mm256_setzero_ps();
              f2line = _mm256_setzero_ps();
              f3line = _mm256_setzero_ps();

              __m256 f1point[NP], f2point[NP], f3point[NP];
              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned int shift = k1*3 + k2*deltaY + k3*deltaZ;

                f1point[k1] = LOADFLD8( fld1, shift ); 
                f2point[k1] = LOADFLD8( fld2, shift ); 
                f3point[k1] = LOADFLD8( fld3, shift ); 
              }

              for ( k1 = 0; k1 < NP; k1 ++ ) {
                f1line = _mm256_add_ps( f1line, _mm256_mul_ps( f1point[k1],  vwx[k1] ) );
                f2line = _mm256_add_ps( f2line, _mm256_mul_ps( f2point[k1], vwxh[k1] ) );
                f3line = _mm256_add_ps( f3line, _mm256_mul_ps( f3point[k1], vwxh[k1] ) );
              }

              f1plane = _mm256_add_ps( f1plane, _mm256_mul_ps( f1line, vwyh[k2] ) );
              f2plane = _mm256_add_ps( f2plane, _mm256_mul_ps( f2line,  vwy[k2] ) );
              f3plane = _mm256_add_ps( f3plane, _mm256_mul_ps( f3line, vwyh[k2] ) );
            }

            vb1 = _mm256_add_ps( vb1, _mm256_mul_ps( f1plane, vwzh[k3] ) );
            vb2 = _mm256_add_ps( vb2, _mm256_mul_ps( f2plane, vwzh[k3] ) );
            vb3 = _mm256_add_ps( vb3, _mm256_mul_ps( f3plane,  vwz[k3] ) );
          }
        }
      }

      // Load momenta
      _MM256_LOAD8v3_PS( vu1, vu2, vu3, u );

      // ---------- advance momenta
      _MM256_LOAD8v3_PS( vu1, vu2, vu3, u );

      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM256_STORE8v3_PS(u, vu1, vu2, vu3); 

      // ---------- advance positions
      {
        __m256 vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256 vq = _mm256_load_ps( q );

        // Accumulate energy
        vene_group = _mm256_add_ps( _mm256_mul_ps( vq, vene ), vene_group );

        // Push positions
        __m256 vx1 = _mm256_add_ps( _mm256_mul_ps( _mm256_mul_ps( vu1, vrg ), vdt_dx1 ), vx0 );
        __m256 vy1 = _mm256_add_ps( _mm256_mul_ps( _mm256_mul_ps( vu2, vrg ), vdt_dx2 ), vy0 );
        __m256 vz1 = _mm256_add_ps( _mm256_mul_ps( _mm256_mul_ps( vu3, vrg ), vdt_dx3 ), vz0 );		

        // Store virtual particles with positions still indexed to the original cell
        STOREU8P3D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

        // Find position trim
        __m256 vtrx = vmntrim(vx1);
        __m256 vtry = vmntrim(vy1);
        __m256 vtrz = vmntrim(vz1);

        // Calculate crossings and store result
        {
          __m256 const zero = _mm256_setzero_ps();
          __m256 vcross;

          // x crossings
          vcross = _mm256_andnot_ps( _mm256_cmp_ps(vtrx, zero, _CMP_EQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(1) ));

          // y crossings
          vcross = _mm256_or_ps( vcross, 
          _mm256_andnot_ps( _mm256_cmp_ps(vtry, zero, _CMP_EQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(2) )));

          // z crossings
          vcross = _mm256_or_ps( vcross, 
          _mm256_andnot_ps( _mm256_cmp_ps(vtrz, zero, _CMP_EQ_OQ ),
          _mm256_castsi256_ps( _mm256_set1_epi32(4) )));

          // Store vcross
          _mm256_store_si256( (__m256i *) cross, _mm256_castps_si256( vcross ) );
        }

        // Trim positions and store results
        vx1   = _mm256_sub_ps( vx1, vtrx );
        vy1   = _mm256_sub_ps( vy1, vtry );
        vz1   = _mm256_sub_ps( vz1, vtrz );
        _MM256_STORE8v3_PS( x, vx1, vy1, vz1 );

        // find cell crossings and store
        __m256i vitrx = _mm256_cvtps_epi32( vtrx ) ;
        __m256i vitry = _mm256_cvtps_epi32( vtry ) ;
        __m256i vitrz = _mm256_cvtps_epi32( vtrz ) ;

        _mm256_store_si256( (__m256i *) dix, vitrx );
        _mm256_store_si256( (__m256i *) diy, vitry );
        _mm256_store_si256( (__m256i *) diz, vitrz );

        // Trim cell indexes and store
        vix0 = viadd( vix0, vitrx ); 
        viy0 = viadd( viy0, vitry ); 
        viz0 = viadd( viz0, vitrz ); 
        _MM256_STORE8v3_EPI32( ix, vix0, viy0, viz0 );

      }

      // ---------- split trajectories for current deposition
      vsplit3D( &vpbuf, cross, dix, diy, diz );

      // ---------- advance pointers
      x  += 3 * VEC_WIDTH;
      u  += 3 * VEC_WIDTH;
      ix += 3 * VEC_WIDTH;
      q  += VEC_WIDTH;

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_ps( vene_group );

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
#undef SPLINEH
#undef NP
#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D

#endif
