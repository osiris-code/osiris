#include <math.h>
#include <stdio.h>

#include "vector-sse.h"
#include "splines-sse.h"
#include "os-spec-push-ssed.h"


/*****************************************************************************************
test_load_store
Tests the load / store macros
*****************************************************************************************/
void test_load_store()
{
 DECLARE_ALIGNED_16( double b[6] );
 dvec v1, v2, v3;
 __m128d s;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store double tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 6 ; i ++ ) {
   b[i] = i;
 }

 _MM_LOAD2v3_PD( v1.v2, v2.v2, v3.v2, b );

 printf(" _MM_LOAD2v3_PD test \n" );
 printf(" v1 = %g, %g | b[0,3] = %g, %g\n",
         v1.v[0], v1.v[1], b[0],  b[3] );
 printf(" v2 = %g, %g | b[1,4] = %g, %g\n",
         v2.v[0], v2.v[1], b[1],  b[4] );
 printf(" v3 = %g, %g | b[2,5] = %g, %g\n",
         v3.v[0], v3.v[1], b[2],  b[5] );

 // Change the __m128d values and store them
 s  = _mm_set1_pd( 10.0 );
 v1.v2 = _mm_add_pd( v1.v2, s );
 v2.v2 = _mm_add_pd( v2.v2, s );
 v3.v2 = _mm_add_pd( v3.v2, s );

 _MM_STORE2v3_PD(b, v1.v2, v2.v2, v3.v2);
 printf(" _MM_STORE2v3_PD test \n" );
 printf(" v1 = %g, %g | b[0,3] = %g, %g\n",
         v1.v[0], v1.v[1], b[0],  b[3] );
 printf(" v2 = %g, %g | b[1,4] = %g, %g\n",
         v2.v[0], v2.v[1], b[1],  b[4] );
 printf(" v3 = %g, %g | b[2,5] = %g, %g\n",
         v3.v[0], v3.v[1], b[2],  b[5] );

 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 4 ; i ++ ) {
   b[i] = i;
 }

 _MM_LOAD2v2_PD( v1.v2, v2.v2, b );

 printf(" _MM_LOAD2v2_PD test \n" );
 printf(" v1 = %g, %g | b[0,2] = %g, %g\n",
         v1.v[0], v1.v[1], b[0],  b[2] );
 printf(" v2 = %g, %g | b[1,3] = %g, %g\n",
         v2.v[0], v2.v[1], b[1],  b[3] );

 // Change the __m128d values and store them
 s = _mm_set1_pd( 10.0 );
 v1.v2 = _mm_add_pd( v1.v2, s );
 v2.v2 = _mm_add_pd( v2.v2, s );

 _MM_STORE2v2_PD(b, v1.v2, v2.v2);
 printf(" _MM_STORE2v2_PD test \n" );
 printf(" v1 = %g, %g | b[0,2] = %g, %g\n",
         v1.v[0], v1.v[1], b[0],  b[2] );
 printf(" v2 = %g, %g | b[1,3] = %g, %g\n",
         v2.v[0], v2.v[1], b[1],  b[3] );

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_load_store_int
Tests the load / store of integer values
*****************************************************************************************/

void test_load_store_int()
{
 DECLARE_ALIGNED_16( int b[6] );
 ivec v1, v2, v3;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store integer tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 6 ; i ++ ) {
   b[i] = i+1;
 }

 _MM_LOAD2v2_EPI32( v1.v4, v2.v4, b )

 printf(" _MM_LOAD2v2_EPI32 test \n" );
 printf(" v1 = %d, %d | b[0,2] = %d, %d\n",
         v1.v[0], v1.v[1], b[0],  b[2] );
 printf(" v2 = %d, %d | b[1,3] = %d, %d\n",
         v2.v[0], v2.v[1], b[1],  b[3] );
 
 v1.v4 = _mm_add_epi32( v1.v4, _mm_set1_epi32(10) );
 v2.v4 = _mm_add_epi32( v2.v4, _mm_set1_epi32(20) );
 
 _MM_STORE2v2_EPI32( b, v1.v4, v2.v4 )

 printf(" _MM_STORE2v2_EPI32 \n" );
 printf(" v1 = %d, %d | b[0,2] = %d, %d\n",
         v1.v[0], v1.v[1], b[0],  b[2] );
 printf(" v2 = %d, %d | b[1,3] = %d, %d\n",
         v2.v[0], v2.v[1], b[1],  b[3] );

 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 6 ; i ++ ) {
   b[i] = i+1;
 }

 _MM_LOAD2v3_EPI32( v1.v4, v2.v4, v3.v4, b )

 printf(" _MM_LOAD2v3_EPI32 test \n" );
 printf(" v1 = %d, %d | b[0,3] = %d, %d\n",
         v1.v[0], v1.v[1], b[0],  b[3] );
 printf(" v2 = %d, %d | b[1,4] = %d, %d\n",
         v2.v[0], v2.v[1], b[1],  b[4] );
 printf(" v3 = %d, %d | b[2,5] = %d, %d\n",
         v3.v[0], v3.v[1], b[2],  b[5] );
 
 v1.v4 = _mm_add_epi32( v1.v4, _mm_set1_epi32(10) );
 v2.v4 = _mm_add_epi32( v2.v4, _mm_set1_epi32(20) );
 v3.v4 = _mm_add_epi32( v3.v4, _mm_set1_epi32(30) );
 
 _MM_STORE2v3_EPI32( b, v1.v4, v2.v4, v3.v4 )

 printf(" _MM_STORE2v3_EPI32 \n" );
 printf(" v1 = %d, %d | b[0,3] = %d, %d\n",
         v1.v[0], v1.v[1], b[0],  b[3] );
 printf(" v2 = %d, %d | b[1,4] = %d, %d\n",
         v2.v[0], v2.v[1], b[1],  b[4] );
 printf(" v3 = %d, %d | b[2,5] = %d, %d\n",
         v3.v[0], v3.v[1], b[2],  b[5] );

 printf( " ----------------------------------------------------------------------- \n" );


}


/*****************************************************************************************
test_splines
Tests the spline calculation
*****************************************************************************************/
void test_splines()
{
 __m128d v, s[5];
 int i, j;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Spline tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v = _mm_set_pd( -0.1, 0.1 );

 for( i = 1; i <= 4; i++ ) {
    dvec a,b ;
 
    switch(i)
    {
      case 1:
         vsplined_s1( v, s ); break;
      case 2:
         vsplined_s2( v, s ); break;
      case 3:
         vsplined_s3( v, s ); break;
      case 4:
         vsplined_s4( v, s ); break;
    }

    printf(" vsplined_s%1d test \n", i );
    a.v2 = v;
    printf(" v   = %g, %g\n", a.v[0], a.v[1] );
    for( j = 0; j < i+1; j++ ) {
      a.v2 = s[j];
      printf( " s[%2d] = %g, %g\n", j-i+1, a.v[0], a.v[1] );
    }
    printf( "\n" );

 }

 printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_trim
Tests the trimming calculations
*****************************************************************************************/

__m128d vmntrim( const __m128d vx )
{
  register __m128d va, vb;

  va = _mm_cmplt_pd( vx, _mm_set1_pd( -0.5 ) );
  va = _mm_and_pd(   va, _mm_set1_pd( -1.0 ) );
  
  vb = _mm_cmpge_pd( vx, _mm_set1_pd( +0.5 ) );
  vb = _mm_and_pd(   vb, _mm_set1_pd( +1.0 ) );
  
  return _mm_add_pd( va, vb );
}


void test_trim()
{
 __m128d v1, v2;
 dvec a, b;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Trim tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v1 = _mm_set_pd( -0.1, 0.1 );
 v2 = vmntrim( v1 );
 
 a.v2 = v1; b.v2 = v2;
 
 printf( " v       = %g, %g \n", a.v[0], a.v[1] );
 printf( " vntrim  = %g, %g \n", b.v[0], b.v[1] );

 printf( "\n" );

 v1 = _mm_set_pd( -0.6, 0.6 );
 v2 = vmntrim( v1 );
 
 a.v2 = v1; b.v2 = v2;
 
 printf( " v       = %g, %g \n", a.v[0], a.v[1] );
 printf( " vntrim  = %g, %g \n", b.v[0], b.v[1] );

 printf( "\n" );

 printf( " ----------------------------------------------------------------------- \n" );

}



/*****************************************************************************************
test_hp
Tests the half-point calculations
*****************************************************************************************/

#define vget_hp_near(dx,ix,dxh,ixh,delta ) {							                 \
  __m128d cmp  = _mm_cmplt_pd( dx, _mm_setzero_pd());                                    \
  dxh = _mm_add_pd( dx, _mm_sub_pd( _mm_and_pd( cmp, _mm_set1_pd( 1.0 ) ),               \
                                    _mm_set1_pd( 0.5 ) ) );                              \
  ixh  = _mm_sub_epi32( ix, _mm_and_si128(                                               \
               _mm_shuffle_epi32( _mm_castpd_si128( cmp ), _MM_SHUFFLE( 0, 0, 2, 0 ) ),  \
               _mm_set1_epi32( delta )) );                                               \
}


void test_hp()
{
 __m128d dx, dxh;
 __m128i ix, ixh;

 int delta;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Half-Point tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 dx   = _mm_set_pd( 0.1, -0.1 );
 ix   = _mm_set_epi32( 40, 30, 20, 10 );

 delta = 3;

 ixh   = _mm_set1_epi32( -1 );

 vget_hp_near( dx, ix, dxh, ixh, delta );

 dvec a;
 ivec b;

 a.v2 = dx;
 printf( " dx      = %g, %g \n", a.v[0], a.v[1] );
 a.v2 = dxh;
 printf( " dxh     = %g, %g \n", a.v[0], a.v[1] );
 
 b.v4 = ix;
 printf( " ix      = %d, %d \n", b.v[0], b.v[1] );
 
 printf( " delta   = %d \n", delta );
 
 b.v4 = ixh;
 printf( " ixh     = %d, %d \n", b.v[0], b.v[1] );

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_load_store_vp
Tests the load/store virtual particle macros
*****************************************************************************************/
void test_load_store_vp2d()
{

 t_vp2D vp[2];

 __m128d vx0, vx1;
 __m128d vy0, vy1;
 __m128d vq, vvz;
 DECLARE_ALIGNED_16( int ix[4] );
 DECLARE_ALIGNED_16( int iy[4] );
 
 int i;
  
 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 2D VP macro tests tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 vp[0].x0 = -0.11;  vp[1].x0 = -0.12; 
 vp[0].x1 =  0.11;  vp[1].x1 =  0.12; 

 vp[0].y0 =  0.21;  vp[1].y0 =  0.22;  
 vp[0].y1 = -0.21;  vp[1].y1 = -0.22;  

 vp[0].q  = -0.01;  vp[1].q  = -0.02; 
 vp[0].vz =  0.91;  vp[1].vz =  0.92; 

 vp[0].ix = 10;     vp[1].ix =  11;    
 vp[0].iy = 20;     vp[1].iy =  21;    

 LOAD2VP2D( ((double *)vp), vx0, vx1, vy0, vy1, vq, vvz, ix, iy );

 printf( "LOAD2VP2D test\n" );

 dvec a;
 a.v2 = vx0;
 printf( " x0 ( %g, %g )\nvx0 ( %g, %g ) \n",
         vp[0].x0, vp[1].x0, a.v[0], a.v[1] );

 a.v2 = vx1;
 printf( " x1 ( %g, %g )\nvx1 ( %g, %g ) \n",
         vp[0].x1, vp[1].x1,
         a.v[0], a.v[1] );

 a.v2 = vy0;
 printf( " y0 ( %g, %g )\nvy0 ( %g, %g ) \n",
         vp[0].y0, vp[1].y0, 
         a.v[0], a.v[1] );
 
 a.v2 = vy1;
 printf( " y1 ( %g, %g )\nvy1 ( %g, %g ) \n",
         vp[0].y1, vp[1].y1,
         a.v[0], a.v[1] );
 
 a.v2 = vq;
 printf( " q ( %g, %g )\nvq ( %g, %g ) \n",
         vp[0].q, vp[1].q, 
         a.v[0], a.v[1] );

 a.v2 = vvz;
 printf( " vz ( %g, %g )\nvvz ( %g, %g ) \n",
         vp[0].vz, vp[1].vz, 
         a.v[0], a.v[1] );
 
 printf( " ix ( %d, %d )\nvix ( %d, %d ) \n",
         vp[0].ix, vp[1].ix, 
         ix[0], ix[1] );
 
 printf( " iy ( %d, %d )\nviy ( %d, %d ) \n",
         vp[0].iy, vp[1].iy, 
         iy[0], iy[1] );

 printf( "\n" );

 #define CHECK(A) {\
   dvec a; \
   a.v2 = v##A; \
   if ( vp[i]. A !=  a.v[i] ) printf("(*error*) Bad %s !\n", #A ); }

 #define CHECKI(A) {\
   if ( vp[i]. A !=  A [i] ) printf("(*error*) Bad %s !\n", #A ); }

 for( i = 0; i < 2; i ++ ) {
   CHECK(x0);
   CHECK(x1);
   CHECK(y0);
   CHECK(y1);
   CHECK(q);
   CHECK(vz);
   CHECKI(ix);
   CHECKI(iy);
 }

 printf( "STORE2VP2D test\n" );

 __m128d t = _mm_set1_pd(0.001);
 vx0 = _mm_add_pd( vx0, t );
 vx1 = _mm_add_pd( vx1, t );
 vy0 = _mm_add_pd( vy0, t );
 vy1 = _mm_add_pd( vy1, t );
 vq  = _mm_add_pd( vq, t );
 vvz = _mm_add_pd( vvz, t );
 
 for(i=0; i<2; i++) { ix[i] +=100; iy[i]+=200; }
 
 __m128i vix, viy;
 vix = _mm_load_si128( (__m128i *)ix );
 viy = _mm_load_si128( (__m128i *)iy );
 
 
 STORE2VP2D( ((double *)vp), vx0, vx1, vy0, vy1, vq, vvz, vix, viy );

 a.v2 = vx0;
 printf( " x0 ( %g, %g )\nvx0 ( %g, %g ) \n",
         vp[0].x0, vp[1].x0, a.v[0], a.v[1] );

 a.v2 = vx1;
 printf( " x1 ( %g, %g )\nvx1 ( %g, %g ) \n",
         vp[0].x1, vp[1].x1,
         a.v[0], a.v[1] );

 a.v2 = vy0;
 printf( " y0 ( %g, %g )\nvy0 ( %g, %g ) \n",
         vp[0].y0, vp[1].y0, 
         a.v[0], a.v[1] );
 
 a.v2 = vy1;
 printf( " y1 ( %g, %g )\nvy1 ( %g, %g ) \n",
         vp[0].y1, vp[1].y1,
         a.v[0], a.v[1] );
 
 a.v2 = vq;
 printf( " q ( %g, %g )\nvq ( %g, %g ) \n",
         vp[0].q, vp[1].q, 
         a.v[0], a.v[1] );

 a.v2 = vvz;
 printf( " vz ( %g, %g )\nvvz ( %g, %g ) \n",
         vp[0].vz, vp[1].vz, 
         a.v[0], a.v[1] );
 
 printf( " ix ( %d, %d )\nvix ( %d, %d ) \n",
         vp[0].ix, vp[1].ix, 
         ix[0], ix[1] );
 
 printf( " iy ( %d, %d )\nviy ( %d, %d ) \n",
         vp[0].iy, vp[1].iy, 
         iy[0], iy[1] );

 printf( "\n" );


 for( i = 0; i < 2; i ++ ) {
   CHECK(x0);
   CHECK(x1);
   CHECK(y0);
   CHECK(y1);
   CHECK(q);
   CHECK(vz);
   CHECKI(ix);
   CHECKI(iy);
 }

 printf( " ----------------------------------------------------------------------- \n" );

}



void test_load_store_vp3d()
{

 t_vp3D vp[4];

 __m128d vx0, vx1;
 __m128d vy0, vy1;
 __m128d vz0, vz1;
 __m128d vq;
 DECLARE_ALIGNED_16( int ix[4] );
 DECLARE_ALIGNED_16( int iy[4] );
 DECLARE_ALIGNED_16( int iz[4] );
 
 int i;

 printf( " size( t_vp3D ) = %lu \n", sizeof( vp[0] ) );
 printf( " size( t_vp3D ) / sizeof( double )= %lu \n", sizeof( vp[0] ) / sizeof( double ) );

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 3D VP macro tests tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 vp[0].x0 = -0.11;  vp[1].x0 = -0.12; 
 vp[0].x1 =  0.11;  vp[1].x1 =  0.12; 

 vp[0].y0 =  0.21;  vp[1].y0 =  0.22;  
 vp[0].y1 = -0.21;  vp[1].y1 = -0.22;  

 vp[0].z0 =  0.31;  vp[1].z0 =  0.32;  
 vp[0].z1 = -0.31;  vp[1].z1 = -0.32;  

 vp[0].q  = -0.01;  vp[1].q  = -0.02; 
 
 vp[0].ix = 10;     vp[1].ix =  11;    
 vp[0].iy = 20;     vp[1].iy =  21;    
 vp[0].iz = 30;     vp[1].iz =  31;    

 LOAD2VP3D( ((double *)vp), vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz );

 printf( "LOAD2VP3D test\n" );

 dvec a;
 a.v2 = vx0;
 printf( " x0 ( %g, %g )\nvx0 ( %g, %g ) \n",
         vp[0].x0, vp[1].x0, a.v[0], a.v[1] );

 a.v2 = vx1;
 printf( " x1 ( %g, %g )\nvx1 ( %g, %g ) \n",
         vp[0].x1, vp[1].x1,
         a.v[0], a.v[1] );

 a.v2 = vy0;
 printf( " y0 ( %g, %g )\nvy0 ( %g, %g ) \n",
         vp[0].y0, vp[1].y0, 
         a.v[0], a.v[1] );
 
 a.v2 = vy1;
 printf( " y1 ( %g, %g )\nvy1 ( %g, %g ) \n",
         vp[0].y1, vp[1].y1,
         a.v[0], a.v[1] );

 a.v2 = vz0;
 printf( " z0 ( %g, %g )\nvz0 ( %g, %g ) \n",
         vp[0].z0, vp[1].z0, 
         a.v[0], a.v[1] );
 
 a.v2 = vz1;
 printf( " z1 ( %g, %g )\nvz1 ( %g, %g ) \n",
         vp[0].z1, vp[1].z1,
         a.v[0], a.v[1] );
 
 a.v2 = vq;
 printf( " q ( %g, %g )\nvq ( %g, %g ) \n",
         vp[0].q, vp[1].q, 
         a.v[0], a.v[1] );
 
 printf( " ix ( %d, %d )\nvix ( %d, %d ) \n",
         vp[0].ix, vp[1].ix, 
         ix[0], ix[1] );
 
 printf( " iy ( %d, %d )\nviy ( %d, %d ) \n",
         vp[0].iy, vp[1].iy, 
         iy[0], iy[1] );

 printf( " iz ( %d, %d )\nviz ( %d, %d ) \n",
         vp[0].iz, vp[1].iz, 
         iz[0], iz[1] );

 printf( "\n" );


 for( i = 0; i < 2; i ++ ) {
   CHECK(x0);   CHECK(x1);
   CHECK(y0);   CHECK(y1);
   CHECK(z0);   CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }
 
 printf( "STORE2VP3D test\n" );

 for(i=0; i<2; i++) { ix[i] +=100; iy[i]+=200; ; iz[i]+=300; }
 
 __m128i vix, viy, viz;
 vix = _mm_load_si128( (__m128i *)ix );
 viy = _mm_load_si128( (__m128i *)iy );
 viz = _mm_load_si128( (__m128i *)iz );
 
 
 STORE2VP3D( ((double *)vp), vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );

 a.v2 = vx0;
 printf( " x0 ( %g, %g )\nvx0 ( %g, %g ) \n",
         vp[0].x0, vp[1].x0, a.v[0], a.v[1] );

 a.v2 = vx1;
 printf( " x1 ( %g, %g )\nvx1 ( %g, %g ) \n",
         vp[0].x1, vp[1].x1,
         a.v[0], a.v[1] );

 a.v2 = vy0;
 printf( " y0 ( %g, %g )\nvy0 ( %g, %g ) \n",
         vp[0].y0, vp[1].y0, 
         a.v[0], a.v[1] );
 
 a.v2 = vy1;
 printf( " y1 ( %g, %g )\nvy1 ( %g, %g ) \n",
         vp[0].y1, vp[1].y1,
         a.v[0], a.v[1] );

 a.v2 = vz0;
 printf( " z0 ( %g, %g )\nvz0 ( %g, %g ) \n",
         vp[0].z0, vp[1].z0, 
         a.v[0], a.v[1] );
 
 a.v2 = vz1;
 printf( " z1 ( %g, %g )\nvz1 ( %g, %g ) \n",
         vp[0].z1, vp[1].z1,
         a.v[0], a.v[1] );

 
 a.v2 = vq;
 printf( " q ( %g, %g )\nvq ( %g, %g ) \n",
         vp[0].q, vp[1].q, 
         a.v[0], a.v[1] );

 
 printf( " ix ( %d, %d )\nvix ( %d, %d ) \n",
         vp[0].ix, vp[1].ix, 
         ix[0], ix[1] );
 
 printf( " iy ( %d, %d )\nviy ( %d, %d ) \n",
         vp[0].iy, vp[1].iy, 
         iy[0], iy[1] );

 printf( " iz ( %d, %d )\nviz ( %d, %d ) \n",
         vp[0].iz, vp[1].iz, 
         iz[0], iz[1] );

 printf( "\n" );


 for( i = 0; i < 2; i ++ ) {
   CHECK(x0);  CHECK(x1);
   CHECK(y0);  CHECK(y1);
   CHECK(z0);  CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }


 printf( " ----------------------------------------------------------------------- \n" );

}

#include "os-spec-current-ssed.c"

int main( int argc, char* argv[])
{
 
//  test_load_store();

//  test_load_store_int();

//  test_splines();

//  test_trim();
  
//  test_hp();

//  test_load_store_vp2d();
  
//  test_load_store_vp3d();

  __m128d q  = _mm_set_pd( 1.0, 1.0 );
     dvec a;
  
  __m128d s[5];
  __m128d wl[4];

  __m128d x0 = _mm_set_pd( 0.49999999999999994, 0.5 );
  __m128d x1 = _mm_set_pd( 0.49999999999999989, 0.49999999999999994 );

  vsplined_s4( x1, s );
  a.v2 = x1;
  printf( "x = %20.18g, %20.18g\n", a.v[0],a.v[1] );
  for( int i = 0; i < 5; i++ ) {
     a.v2 = s[i];
     printf("s[%d] = %20.18g, %20.18g\n",i, a.v[0], a.v[1]);
  }
  
  printf("-------------------------\n");
  
  a.v2 = x0;
  printf( "x0    = %18.18g, %18.18g\n", a.v[0],a.v[1] );
  a.v2 = x1;
  printf( "x1    = %18.18g, %18.18g\n", a.v[0],a.v[1] );
 
  vwl_s4( q, x0, x1, wl );
  for( int i = 0; i < 4; i++ ) {
     a.v2 = wl[i];
     printf("wl[%d] = %18.18g, %18.18g\n",i, a.v[0], a.v[1]);
  }

  printf("-------------------------\n");

  x0 = _mm_set_pd( 0., 0. );
  x1 = _mm_set_pd( .1, -.1 );

  vsplined_s4( x1, s );
  a.v2 = x1;
  printf( "x = %20.18g, %20.18g\n", a.v[0],a.v[1] );
  for( int i = 0; i < 5; i++ ) {
     a.v2 = s[i];
     printf("s[%d] = %20.18g, %20.18g\n",i, a.v[0], a.v[1]);
  }
  
  printf("-------------------------\n");
  
  a.v2 = x0;
  printf( "x0    = %18.18g, %18.18g\n", a.v[0],a.v[1] );
  a.v2 = x1;
  printf( "x1    = %18.18g, %18.18g\n", a.v[0],a.v[1] );
 
  vwl_s4( q, x0, x1, wl );
  for( int i = 0; i < 4; i++ ) {
     a.v2 = wl[i];
     printf("wl[%d] = %18.18g, %18.18g\n",i, a.v[0], a.v[1]);
  }

}
