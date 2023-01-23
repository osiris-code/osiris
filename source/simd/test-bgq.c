#include <math.h>
#include <stdio.h>
#include <mpi.h>

#include <inttypes.h>

#include "vector-bgq.h"
#include "splines-bgq.h"
#include "os-spec-push-bgq.h"

/*****************************************************************************************
test_calc
Tests the vec_r and vec_rsqrt extensions (full precision reciprocal and reciprocal 
square root)
*****************************************************************************************/
void test_calc()
{
  double a, b;
  double ra, rb;
  vector v1, v2;
  
  a = M_PI; b = M_E;
  
  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Precision tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  //   vec_r - Parallel Reciprocal 
  
  ra = 1.0/a; rb = 1.0/b;
  
  printf( " Scalar:  1/a = %g, 1/b = %g \n", ra, rb );
  
  v1 = (vector) { a, b, 1.0, 2.0 };
  v2 = vec_r( v1 );

  printf( " Vector:  1/a = %g, 1/b = %g \n", v2[0], v2[1] );
  printf( " error :  eps = %g, eps = %g \n", (ra - v2[0])/ra, (rb - v2[1])/rb );

  //   vec_rsqrt
  printf( " ----------------------------------------------------------------------- \n" );
  
  ra = 1.0/sqrt(a); rb = 1.0/sqrt(b);

  printf( " Scalar:  1/sqrt(a) = %g, 1/sqrt(b) = %g \n", ra, rb );
  
  v1 = (vector) { a, b, 1.0, 2.0 };
  v2 = vec_rsqrt( v1 );

  printf( " Vector:  1/sqrt(a) = %g, 1/sqrt(b) = %g \n", v2[0], v2[1] );
  printf( " error :  eps       = %g, eps       = %g \n", (ra - v2[0])/ra, (rb - v2[0])/rb );

  printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_load_store
Tests the load / store macros
*****************************************************************************************/
void test_load_store()
{
  double b[12] __attribute__ ((aligned (32)));
  vector v1, v2, v3, s;

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Load/store macro tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  
  // Initialize buffer values
  for( int i = 0; i < 12 ; i ++ ) {
    b[i] = i;
  }
  
  VEC_LD4v3( v1, v2, v3, b );
  
  printf(" VEC_LD4v3 test \n" );
  printf(" v1 = %g, %g, %g, %g  | b[0,3,6,9] = %g, %g, %g, %g\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[3],  b[6],  b[9] );
  printf(" v2 = %g, %g, %g, %g  | b[1,4,7,10] = %g, %g, %g, %g\n", 
          v2[0], v2[1], v2[2], v2[3],
          b[1],  b[4],  b[7],  b[10] );
  printf(" v3 = %g, %g, %g, %g  | b[2,5,8,11] = %g, %g, %g, %g\n", 
          v3[0], v3[1], v3[2], v3[3],
          b[2],  b[5],  b[8],  b[11] );

  // Change the vector values and store them
  s  = vec_splats( 10.0 );
  v1 = vec_add( v1, s );
  v2 = vec_add( v2, s );
  v3 = vec_add( v3, s );
  
  VEC_ST4v3(b, v1, v2, v3);
  printf(" VEC_ST4v3 test \n" );
  printf(" v1 = %g, %g, %g, %g  | b[0,3,6,9] = %g, %g, %g, %g\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[3],  b[6],  b[9] );
  printf(" v2 = %g, %g, %g, %g  | b[1,4,7,10] = %g, %g, %g, %g\n", 
          v2[0], v2[1], v2[2], v2[3],
          b[1],  b[4],  b[7],  b[10] );
  printf(" v3 = %g, %g, %g, %g  | b[2,5,8,11] = %g, %g, %g, %g\n", 
          v3[0], v3[1], v3[2], v3[3],
          b[2],  b[5],  b[8],  b[11] );
  
  printf( " ----------------------------------------------------------------------- \n" );
 
  // Initialize buffer values
  for( int i = 0; i < 8 ; i ++ ) {
    b[i] = i;
  }

  VEC_LD4v2( v1, v2, b );
  
  printf(" VEC_LD4v2 test \n " );
  printf(" v1 = %g, %g, %g, %g  | b[0,2,4,6] = %g, %g, %g, %g\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[2],  b[4],  b[6] );
  printf(" v2 = %g, %g, %g, %g  | b[1,3,5,7] = %g, %g, %g, %g\n", 
          v2[0], v2[1], v2[2], v2[3],
          b[1],  b[3],  b[5],  b[7] );

  // Change the vector values and store them
  s = vec_splats( 10.0 );
  v1 = vec_add( v1, s );
  v2 = vec_add( v2, s );
  
  VEC_ST4v2(b, v1, v2);
  printf(" VEC_ST4v2 test \n " );
  printf(" v1 = %g, %g, %g, %g  | b[0,2,4,6] = %g, %g, %g, %g\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[2],  b[4],  b[6] );
  printf(" v2 = %g, %g, %g, %g  | b[1,3,5,7] = %g, %g, %g, %g\n", 
          v2[0], v2[1], v2[2], v2[3],
          b[1],  b[3],  b[5],  b[7] );

  printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_load_store_int
Tests the load / store of integer values
*****************************************************************************************/

void test_load_store_int()
{
  int b[4] __attribute__ ((aligned (32)));
  vector v1, v2;

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Load/store integer tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  
  // Initialize buffer values
  for( int i = 0; i < 4 ; i ++ ) {
    b[i] = i+1;
  }
  
  v1 = vec_ldiaa( 0, b );
  v1 = vec_cfid( v1 );
  
  
  printf(" Load integer test \n" );
  printf(" v1 = %g, %g, %g, %g  | b[0,1,2,3] = %d, %d, %d, %d\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[1],  b[2],  b[3] );

  // Change the vector values and store them
  v1 = vec_add( v1, vec_splats( 10.0 ));
  vec_sta( vec_ctiw( v1 ), 0, b );
  
  printf(" Store integer test \n" );
  printf(" v1 = %g, %g, %g, %g  | b[0,1,2,3] = %d, %d, %d, %d\n", 
          v1[0], v1[1], v1[2], v1[3],
          b[0],  b[1],  b[2],  b[3] );
  
  printf( " ----------------------------------------------------------------------- \n" );
 
}

/*****************************************************************************************
test_splines
Tests the spline calculation
*****************************************************************************************/
void test_splines()
{
  vector v, s[5];
  int i, j;

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Spline tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  
  v = (vector)( -0.25, -0.1, 0.1, 0.25 );
  
  for( i = 1; i <= 4; i++ ) {
     switch(i)
     {
       case 1:
          vspline_s1( v, s ); break;
       case 2:
          vspline_s2( v, s ); break;
       case 3:
          vspline_s3( v, s ); break;
       case 4:
          vspline_s4( v, s ); break;
     }
  
     printf(" vspline_s%1d test \n", i );
     printf(" v   = %g, %g, %g, %g\n", v[0], v[1], v[2], v[3] );
     for( j = 0; j < i+1; j++ )
       printf( " s[%2d] = %g, %g, %g, %g\n", j-i+1, (s[j])[0], (s[j])[1], (s[j])[2], (s[j])[3] );
     printf( "\n" );
     
  }  
  
  printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_trim
Tests the trimming calculations
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


void test_trim()
{
  vector v1, v2;
  

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Trim tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );

  v1 = (vector)( -0.1, 0.1, -0.5, 0.5 );
  v2 = vntrim( v1 );
  
  printf( " v       = %g, %g, %g, %g \n", v1[0], v1[1], v1[2], v1[3] );
  printf( " vntrim  = %g, %g, %g, %g \n", v2[0], v2[1], v2[2], v2[3] );
  
  printf( "\n" );
  
  v1 = (vector)( -0.6, 0.6, -1.0, 1.0 );
  v2 = vntrim( v1 );
  
  printf( " v       = %g, %g, %g, %g \n", v1[0], v1[1], v1[2], v1[3] );
  printf( " vntrim  = %g, %g, %g, %g \n", v2[0], v2[1], v2[2], v2[3] );
  
  printf( "\n" );
  
  printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_hp
Tests the half-point calculations
*****************************************************************************************/

#define vget_hp_near(dx,ix,dxh,ixh,delta ) 					\
{ vector s;													\
															\
  vector const oneHalf = vec_splats( 0.5 );				    \
  															\
  /* s = (dx < 0 ? -0.5 : 0.5) */							\
  s   = vec_sel( vec_splats( -0.5 ), oneHalf, dx );			\
  dxh = vec_sub( dx, s );									\
  															\
  s   = vec_sub( s, oneHalf );								\
  s   = vec_ctiw( s );										\
  vec_sta( s, 0, ixh );										\
  															\
  for( int i = 0; i < 4; i++ ) 								\
    ixh[i] = ix[i] - (ixh[i] & delta); 						\
}


void test_hp()
{
  vector dx, dxh;
  int ix[4], ixh[4];
  
  int delta;
  

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Half-Point tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  
  dx   = (vector)( 0.4, 0.1, -0.1, -0.4 );
  ix[0] = 10;
  ix[1] = 20;
  ix[2] = 30;
  ix[3] = 40;
  
  delta = 3;
  
  ixh[0] = -1;
  ixh[1] = -1;
  ixh[2] = -1;
  ixh[3] = -1;

  vget_hp_near( dx, ix, dxh, ixh, delta );
   															   
  printf( " dx      = %g, %g, %g, %g \n", dx[0], dx[1], dx[2], dx[3] );
  printf( " dxh     = %g, %g, %g, %g \n", dxh[0], dxh[1], dxh[2], dxh[3] );
  printf( " ix      = %d, %d, %d, %d \n", ix[0], ix[1], ix[2], ix[3] );
  printf( " delta   = %d \n", delta );
  printf( " ixh     = %d, %d, %d, %d \n", ixh[0], ixh[1], ixh[2], ixh[3] );
    
  printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_load_store_vp
Tests the load/store virtual particle macros
*****************************************************************************************/
void test_load_store_vp2d()
{
    
  t_vp2D vp[4];
    
  vector vx0, vx1;
  vector vy0, vy1;
  vector vq, vvz;
  DECLARE_ALIGNED_32( int vix[4] );
  DECLARE_ALIGNED_32( int viy[4] );

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Load/store 2D VP macro tests tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  printf( " sizeof( t_vp2D ) = %zu \n ", sizeof( t_vp2D ) );
  printf( " sizeof( t_vp2D ) / sizeof( double ) = %zu \n ", sizeof( t_vp2D ) / sizeof( double ) );
  printf( " ----------------------------------------------------------------------- \n" );

  vp[0].x0 = -0.11;  vp[1].x0 = -0.12;  vp[2].x0 = -0.13;  vp[3].x0 = -0.14;   
  vp[0].x1 =  0.11;  vp[1].x1 =  0.12;  vp[2].x1 =  0.13;  vp[3].x1 =  0.14;
    
  vp[0].y0 =  0.21;  vp[1].y0 =  0.22;  vp[2].y0 =  0.23;  vp[3].y0 =  0.24;
  vp[0].y1 = -0.21;  vp[1].y1 = -0.22;  vp[2].y1 = -0.23;  vp[3].y1 = -0.24;
  
  vp[0].q  = -0.01;  vp[1].q  = -0.02;  vp[2].q  = -0.03;  vp[3].q  = -0.04;
  vp[0].vz =  0.91;  vp[1].vz =  0.92;  vp[2].vz =  0.93;  vp[3].vz =  0.94;
  
  vp[0].ix = 10;     vp[1].ix =  11;    vp[2].ix = 12;     vp[3].ix =  13;
  vp[0].iy = 20;     vp[1].iy =  21;    vp[2].iy = 22;     vp[3].iy =  23;
  
  LOAD4VP2D( ((double *)vp), vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
  
  printf( "LOAD4VP2D test\n" );
  
  printf( " x0 ( %g, %g, %g, %g )\nvx0 ( %g, %g, %g, %g ) \n", 
          vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
          vx0[0], vx0[1], vx0[2], vx0[3] );

  printf( " x1 ( %g, %g, %g, %g )\nvx1 ( %g, %g, %g, %g ) \n", 
          vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
          vx1[0], vx1[1], vx1[2], vx1[3] );
  
  printf( " y0 ( %g, %g, %g, %g )\nvy0 ( %g, %g, %g, %g ) \n", 
          vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
          vy0[0], vy0[1], vy0[2], vy0[3] );

  printf( " y1 ( %g, %g, %g, %g )\nvy1 ( %g, %g, %g, %g ) \n", 
          vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
          vy1[0], vy1[1], vy1[2], vy1[3] );
  
  printf( " q ( %g, %g, %g, %g )\nvq ( %g, %g, %g, %g ) \n", 
          vp[0].q, vp[1].q, vp[2].q, vp[3].q,
          vq[0], vq[1], vq[2], vq[3] );

  printf( " vz ( %g, %g, %g, %g )\nvvz ( %g, %g, %g, %g ) \n", 
          vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz,
          vvz[0], vvz[1], vvz[2], vvz[3] );
  
  printf( " ix ( %d, %d, %d, %d )\nvix ( %d, %d, %d, %d ) \n", 
          vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
          vix[0], vix[1], vix[2], vix[3] );

  printf( " iy ( %d, %d, %d, %d )\nviy ( %d, %d, %d, %d ) \n", 
          vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
          viy[0], viy[1], viy[2], viy[3] );
  
  printf( "\n" );
  
  #define CHECK(A)\
    if ( vp[i]. A != v##A [i] ) printf("(*error*) Bad %s !\n", #A )
  
  for( int i = 0; i < 4; i ++ ) {
    CHECK(x0);
    CHECK(x1);
    CHECK(y0);
    CHECK(y1);
    CHECK(q);
    CHECK(vz);
    CHECK(ix);
    CHECK(iy);
  }

  printf( "STORE4VP2D test\n" );
    
  vector t = vec_splats(0.001);
  vx0 = vec_add( vx0, t );
  vx1 = vec_add( vx1, t );
  vy0 = vec_add( vy0, t );
  vy1 = vec_add( vy1, t );
  vq = vec_add( vq, t );
  vvz = vec_add( vvz, t );
      
  STORE4VP2D( (double *)vp, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
  
  printf( " x0 ( %g, %g, %g, %g )\nvx0 ( %g, %g, %g, %g ) \n", 
          vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
          vx0[0], vx0[1], vx0[2], vx0[3] );

  printf( " x1 ( %g, %g, %g, %g )\nvx1 ( %g, %g, %g, %g ) \n", 
          vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
          vx1[0], vx1[1], vx1[2], vx1[3] );
  
  printf( " y0 ( %g, %g, %g, %g )\nvy0 ( %g, %g, %g, %g ) \n", 
          vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
          vy0[0], vy0[1], vy0[2], vy0[3] );

  printf( " y1 ( %g, %g, %g, %g )\nvy1 ( %g, %g, %g, %g ) \n", 
          vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
          vy1[0], vy1[1], vy1[2], vy1[3] );
  
  printf( " q ( %g, %g, %g, %g )\nvq ( %g, %g, %g, %g ) \n", 
          vp[0].q, vp[1].q, vp[2].q, vp[3].q,
          vq[0], vq[1], vq[2], vq[3] );

  printf( " vz ( %g, %g, %g, %g )\nvvz ( %g, %g, %g, %g ) \n", 
          vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz,
          vvz[0], vvz[1], vvz[2], vvz[3] );
  
  printf( " ix ( %d, %d, %d, %d )\nvix ( %d, %d, %d, %d ) \n", 
          vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
          vix[0], vix[1], vix[2], vix[3] );

  printf( " iy ( %d, %d, %d, %d )\nviy ( %d, %d, %d, %d ) \n", 
          vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
          viy[0], viy[1], viy[2], viy[3] );
  

  for( int i = 0; i < 4; i ++ ) {
    CHECK(x0);
    CHECK(x1);
    CHECK(y0);
    CHECK(y1);
    CHECK(q);
    CHECK(vz);
    CHECK(ix);
    CHECK(iy);
  }
    
  printf( " ----------------------------------------------------------------------- \n" );

}


void test_load_store_vp3d()
{
    
  t_vp3D vp[4];
    
  vector vx0, vx1;
  vector vy0, vy1;
  vector vz0, vz1;
  vector vq;
  DECLARE_ALIGNED_32( int vix[4] );
  DECLARE_ALIGNED_32( int viy[4] );
  DECLARE_ALIGNED_32( int viz[4] );

  printf( " ----------------------------------------------------------------------- \n" );
  printf( " Load/store 3D VP macro tests tests \n" );
  printf( " ----------------------------------------------------------------------- \n" );
  printf( " sizeof( t_vp3D ) = %zu \n ", sizeof( t_vp3D ) );
  printf( " sizeof( t_vp3D ) / sizeof( double ) = %zu \n ", sizeof( t_vp3D ) / sizeof( double ) );
  printf( " ----------------------------------------------------------------------- \n" );

  vp[0].x0 = -0.11;  vp[1].x0 = -0.12;  vp[2].x0 = -0.13;  vp[3].x0 = -0.14;   
  vp[0].x1 =  0.11;  vp[1].x1 =  0.12;  vp[2].x1 =  0.13;  vp[3].x1 =  0.14;
    
  vp[0].y0 =  0.21;  vp[1].y0 =  0.22;  vp[2].y0 =  0.23;  vp[3].y0 =  0.24;
  vp[0].y1 = -0.21;  vp[1].y1 = -0.22;  vp[2].y1 = -0.23;  vp[3].y1 = -0.24;

  vp[0].z0 =  0.31;  vp[1].z0 =  0.32;  vp[2].z0 =  0.33;  vp[3].z0 =  0.34;
  vp[0].z1 = -0.31;  vp[1].z1 = -0.32;  vp[2].z1 = -0.33;  vp[3].z1 = -0.34;
  
  vp[0].q  = -0.01;  vp[1].q  = -0.02;  vp[2].q  = -0.03;  vp[3].q  = -0.04;
  
  vp[0].ix = 10;     vp[1].ix =  11;    vp[2].ix = 12;     vp[3].ix =  13;
  vp[0].iy = 20;     vp[1].iy =  21;    vp[2].iy = 22;     vp[3].iy =  23;
  vp[0].iz = 30;     vp[1].iz =  31;    vp[2].iz = 32;     vp[3].iz =  33;
  
  LOAD4VP3D( ((double *)vp), vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
  
  printf( "LOAD4VP3D test\n" );
  
  printf( " x0 ( %g, %g, %g, %g )\nvx0 ( %g, %g, %g, %g ) \n", 
          vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
          vx0[0], vx0[1], vx0[2], vx0[3] );

  printf( " x1 ( %g, %g, %g, %g )\nvx1 ( %g, %g, %g, %g ) \n", 
          vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
          vx1[0], vx1[1], vx1[2], vx1[3] );
  
  printf( " y0 ( %g, %g, %g, %g )\nvy0 ( %g, %g, %g, %g ) \n", 
          vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
          vy0[0], vy0[1], vy0[2], vy0[3] );

  printf( " y1 ( %g, %g, %g, %g )\nvy1 ( %g, %g, %g, %g ) \n", 
          vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
          vy1[0], vy1[1], vy1[2], vy1[3] );

  printf( " z0 ( %g, %g, %g, %g )\nvz0 ( %g, %g, %g, %g ) \n", 
          vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0,
          vz0[0], vz0[1], vz0[2], vz0[3] );

  printf( " z1 ( %g, %g, %g, %g )\nvz1 ( %g, %g, %g, %g ) \n", 
          vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1,
          vz1[0], vz1[1], vz1[2], vz1[3] );
  
  printf( " q ( %g, %g, %g, %g )\nvq ( %g, %g, %g, %g ) \n", 
          vp[0].q, vp[1].q, vp[2].q, vp[3].q,
          vq[0], vq[1], vq[2], vq[3] );
  
  printf( " ix ( %d, %d, %d, %d )\nvix ( %d, %d, %d, %d ) \n", 
          vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
          vix[0], vix[1], vix[2], vix[3] );

  printf( " iy ( %d, %d, %d, %d )\nviy ( %d, %d, %d, %d ) \n", 
          vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
          viy[0], viy[1], viy[2], viy[3] );

  printf( " iz ( %d, %d, %d, %d )\nviz ( %d, %d, %d, %d ) \n", 
          vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz,
          viz[0], viz[1], viz[2], viz[3] );
  
  printf( "\n" );
    
  for( int i = 0; i < 4; i ++ ) {
    CHECK(x0);  CHECK(x1);
    CHECK(y0);  CHECK(y1);
    CHECK(z0);  CHECK(z1);
    CHECK(q);
    CHECK(ix);  CHECK(iy);  CHECK(iz);
  }

  printf( "STORE4VP3D test\n" );
    
  vector t = vec_splats(0.001);
  vx0 = vec_add( vx0, t );
  vx1 = vec_add( vx1, t );
  vy0 = vec_add( vy0, t );
  vy1 = vec_add( vy1, t );
  vz0 = vec_add( vz0, t );
  vz1 = vec_add( vz1, t );
  vq = vec_add( vq, t );
      
  STORE4VP3D( (double *)vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
  
  printf( " x0 ( %g, %g, %g, %g )\nvx0 ( %g, %g, %g, %g ) \n", 
          vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
          vx0[0], vx0[1], vx0[2], vx0[3] );

  printf( " x1 ( %g, %g, %g, %g )\nvx1 ( %g, %g, %g, %g ) \n", 
          vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
          vx1[0], vx1[1], vx1[2], vx1[3] );
  
  printf( " y0 ( %g, %g, %g, %g )\nvy0 ( %g, %g, %g, %g ) \n", 
          vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
          vy0[0], vy0[1], vy0[2], vy0[3] );

  printf( " y1 ( %g, %g, %g, %g )\nvy1 ( %g, %g, %g, %g ) \n", 
          vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
          vy1[0], vy1[1], vy1[2], vy1[3] );

  printf( " z0 ( %g, %g, %g, %g )\nvz0 ( %g, %g, %g, %g ) \n", 
          vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0,
          vz0[0], vz0[1], vz0[2], vz0[3] );

  printf( " z1 ( %g, %g, %g, %g )\nvz1 ( %g, %g, %g, %g ) \n", 
          vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1,
          vz1[0], vz1[1], vz1[2], vz1[3] );
  
  printf( " q ( %g, %g, %g, %g )\nvq ( %g, %g, %g, %g ) \n", 
          vp[0].q, vp[1].q, vp[2].q, vp[3].q,
          vq[0], vq[1], vq[2], vq[3] );
  
  printf( " ix ( %d, %d, %d, %d )\nvix ( %d, %d, %d, %d ) \n", 
          vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
          vix[0], vix[1], vix[2], vix[3] );

  printf( " iy ( %d, %d, %d, %d )\nviy ( %d, %d, %d, %d ) \n", 
          vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
          viy[0], viy[1], viy[2], viy[3] );

  printf( " iz ( %d, %d, %d, %d )\nviz ( %d, %d, %d, %d ) \n", 
          vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz,
          viz[0], viz[1], viz[2], viz[3] );
  

 for( int i = 0; i < 4; i ++ ) {
    CHECK(x0);  CHECK(x1);
    CHECK(y0);  CHECK(y1);
    CHECK(z0);  CHECK(z1);
    CHECK(q);
    CHECK(ix);  CHECK(iy);  CHECK(iz);
  }
    
  printf( " ----------------------------------------------------------------------- \n" );

}


int main( int argc, char* argv[])
{
  int rank, size;
  
  MPI_Init (&argc, &argv);	
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
  MPI_Comm_size (MPI_COMM_WORLD, &size);	

  if ( rank == 0 ) {
     
     test_calc();
     
	 test_load_store();
  
	 test_splines();
	 
	 test_trim();
	 
	 test_hp();

	 test_load_store_int();
	 
	 test_load_store_vp2d();

	 test_load_store_vp3d();

  }
  
  MPI_Finalize();
  return 0;

}
