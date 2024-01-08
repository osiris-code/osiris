/*****************************************************************************************
  Tests for AVX specific code (double precision)
*****************************************************************************************/

#include <math.h>
#include <stdio.h>

#include "vector-avx.h"
#include "splines-avx.h"

typedef union vecInt128 {int v[4]; __m128i v4;} ivec128;

#include "os-spec-push-avxd.h"

/*****************************************************************************************
test_load_store
Tests the load / store macros
*****************************************************************************************/

#define CHECKBUF_PD( vec, b, offset, stride ) { \
  int i; \
  for( i = 0; i < 4; i++ ) { \
     if ( vec.v[i] != b[ offset + i * stride ] ) { \
       printf( "\n(*error*) %g /= %g\n", vec.v[i], b[ offset + i * stride ] ); \
       exit(-1); \
     } \
  } \
}

void test_load_store()
{
 DECLARE_ALIGNED_32( double b[4*3] );
 dvec v1, v2, v3;
 __m256d s;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store double tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 4 ; i ++ ) {
   b[3*i  ] = 1.0 + 0.1 * i;
   b[3*i+1] = 2.0 + 0.1 * i;
   b[3*i+2] = 3.0 + 0.1 * i;
 }
 
 _MM256_LOAD4v3_PD(v1.v4, v2.v4, v3.v4, b)
 
 
 printf(" _MM256_LOAD4v3_PD test \n" );

 printf("             v1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf("b[ 0, 3, 6, 9 ] = %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9]  ); 

 CHECKBUF_PD( v1, b, 0, 3 );

 printf("             v2 = %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf("b[ 1, 4, 7,10 ] = %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10]  ); 

 CHECKBUF_PD( v2, b, 1, 3 );

 printf("             v3 = %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3]  );
 printf("b[ 2, 5, 8,11 ] = %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11] ); 

 CHECKBUF_PD( v3, b, 2, 3 );

 // Change the __m128d values and store them
 s  = _mm256_set1_pd( 10.0 );
 v1.v4 = _mm256_add_pd( v1.v4, s );
 v2.v4 = _mm256_add_pd( v2.v4, s );
 v3.v4 = _mm256_add_pd( v3.v4, s );

 printf(" _MM256_STORE4v3_PD test \n" );

 _MM256_STORE4v3_PD(b, v1.v4, v2.v4, v3.v4);

 printf("             v1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf("b[ 0, 3, 6, 9 ] = %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9] ); 

 CHECKBUF_PD( v1, b, 0, 3 );

 printf("             v2 = %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf("b[ 1, 4, 7,10 ] = %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10]  ); 

 CHECKBUF_PD( v2, b, 1, 3 );

 printf("             v3 = %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3]  );
 printf("b[ 2, 5, 8,11 ] = %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11]  ); 

 CHECKBUF_PD( v3, b, 2, 3 );
 
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 4 ; i ++ ) {
   b[2*i  ] = 1.0 + 0.1 * i;
   b[2*i+1] = 2.0 + 0.1 * i;
 }

 printf(" _MM256_LOAD4v2_PD test \n" );
 _MM256_LOAD4v2_PD( v1.v4, v2.v4, b );

 printf("             v1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf("b[ 0, 2, 4, 6 ] = %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6]  ); 

 CHECKBUF_PD( v1, b, 0, 2 );

 printf("             v2 = %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf("b[ 1, 3, 5, 7 ] = %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7]  ); 

 CHECKBUF_PD( v2, b, 1, 2 );

 // Change the __m256d values and store them
 s = _mm256_set1_pd( 10.0 );
 v1.v4 = _mm256_add_pd( v1.v4, s );
 v2.v4 = _mm256_add_pd( v2.v4, s );

 printf(" _MM256_STORE4v2_PD test \n" );
 _MM256_STORE4v2_PD(b, v1.v4, v2.v4);

 printf("             v1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf("b[ 0, 2, 4, 6 ] = %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6]  ); 

 CHECKBUF_PD( v1, b, 0, 2 );

 printf("             v2 = %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf("b[ 1, 3, 5, 7 ] = %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7]  ); 

 CHECKBUF_PD( v2, b, 1, 2 );

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_splines
Tests the spline calculation
*****************************************************************************************/
void test_splines()
{
 __m256d v, s[5];
 int i, j;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Spline tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v = _mm256_setr_pd( -0.5, -0.2, 0.0, 0.2 );

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
    a.v4 = v;
    printf(" v   = %g, %g, %g, %g\n", 
            a.v[0], a.v[1], a.v[2], a.v[3] );
    
    b.v4 = _mm256_setzero_pd();
    
    for( j = 0; j < i+1; j++ ) {
      a.v4 = s[j];
      b.v4 = _mm256_add_pd( b.v4, a.v4 );
      printf( " s[%2d] = %g, %g, %g, %g\n", j-i+1, 
               a.v[0], a.v[1], a.v[2], a.v[3] );
    }
    printf(" tot = %g, %g, %g, %g\n", 
            b.v[0], b.v[1], b.v[2], b.v[3] );
    printf( "\n" );

 }

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_trim
Tests the trimming calculations
*****************************************************************************************/

static inline __m256d vmntrim( const __m256d vx )
{
   register __m256d va, vb;
   
   va = _mm256_cmp_pd( vx, _mm256_set1_pd( -0.5 ), _CMP_LT_OS );
   va = _mm256_and_pd( va, _mm256_set1_pd( +1.0 ) );
   
   vb = _mm256_cmp_pd( vx, _mm256_set1_pd( +0.5 ), _CMP_GE_OS );
   vb = _mm256_and_pd( vb, _mm256_set1_pd( +1.0 ) );
   
   return  _mm256_sub_pd( vb, va );
}

void test_trim()
{
 __m256d v1, v2;
 dvec a, b;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Trim tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v1 = _mm256_setr_pd( -0.6, -0.4, 0.5, 0.6 );
 v2 = vmntrim( v1 );
 
 a.v4 = v1; b.v4 = v2;
 
 printf( " v       = %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3] );
 printf( " vntrim  = %g, %g, %g, %g\n", 
         b.v[0], b.v[1], b.v[2], b.v[3] );

 printf( "\n" );

 printf( " ----------------------------------------------------------------------- \n" );

}



/*****************************************************************************************
test_hp
Tests the half-point calculations
*****************************************************************************************/

#define vget_hp_near(dx, ix, dxh, ixh, delta ) 	{ \
   __m256d cmp  = _mm256_cmp_pd( (dx), _mm256_setzero_pd(), _CMP_LT_OS ); \
   (dxh) =  _mm256_add_pd( dx, _mm256_blendv_pd( _mm256_set1_pd( -0.5 ), \
											     _mm256_set1_pd(  0.5 ), cmp ) ); \
   cmp = _mm256_castps_pd( _mm256_permute_ps( _mm256_castpd_ps( cmp ), 0x08 ) ); \
   __m128i cmpi = _mm_castpd_si128( \
			  _mm_unpacklo_pd( _mm256_castpd256_pd128( _mm256_castpd_ps( cmp ) ), \
							   _mm256_extractf128_pd( _mm256_castpd_ps( cmp ), 1 ) ) ); \
   (ixh)  = _mm_sub_epi32( (ix), _mm_and_si128( cmpi,_mm_set1_epi32( delta )) ); \
}

void test_hp()
{
 __m256d dx, dxh;
 __m128i ix, ixh;

 DECLARE_ALIGNED_32( int b[4] );

 int delta;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Half-Point tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 dx   = _mm256_set_pd( 0.4, -0.2, 0.1, -0.3 );
 ix   = _mm_set_epi32( 10, 30, 50, 70 );

 delta = 3;

 ixh   = _mm_set1_epi32( -1 );

 vget_hp_near( dx, ix, dxh, ixh, delta );

 dvec a;

 a.v4 = dx;
 printf( " dx      = %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3] );
 a.v4 = dxh;
 printf( " dxh     = %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3] );
 
 _mm_store_si128( (__m128i * ) b, ix );
 printf( " ix      = %d, %d, %d, %d\n", 
         b[0], b[1], b[2], b[3] );
 
 printf( " delta   = %d \n", delta );
 
 _mm_store_si128( (__m128i * ) b, ixh );
 printf( " ixh     = %d, %d, %d, %d\n", 
         b[0], b[1], b[2], b[3] );

 printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_load_store_vp
Tests the load/store virtual particle macros
*****************************************************************************************/



void test_load_store_vp2d()
{

 DECLARE_ALIGNED_32( t_vp2D vp[4] );

 __m256d vx0, vx1;
 __m256d vy0, vy1;
 __m256d vq, vvz;
 __m128i vix, viy;
 
 int i;
  
 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 2D VP macro tests tests \n" );
 printf( " size( t_vp2D ) = %zu \n", sizeof( t_vp2D ) );
 printf( " size( t_vp2D ) / sizeof( double ) = %zu \n", sizeof( t_vp2D ) / sizeof( double ) );
 printf( " ----------------------------------------------------------------------- \n" );
 
 for( i = 0; i < 4; i ++ ) {
   vp[i].x0 = -0.11 - 0.01*i; 
   vp[i].x1 =  0.11 + 0.01*i; 

   vp[i].y0 =  0.21 + 0.01*i;  
   vp[i].y1 = -0.21 - 0.01*i;  

   vp[i].q  = -0.01 - 0.01*i; 
   vp[i].vz =  0.91 + 0.01*i; 

   vp[i].ix = 10 + i;    
   vp[i].iy = 20 + i;
 }   


 double* b = ((double *) vp);

 LOAD4VP2D( b, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

 
 printf( "LOAD4VP2D test\n" );

 dvec v1;
 ivec128 v2;

 v1.v4 = vx0;
 printf(" x0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x0 = %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0  );


 v1.v4 = vx1;
 printf(" x1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x1 = %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1  );

 v1.v4 = vy0;
 printf(" y0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y0 = %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0  );
 
 v1.v4 = vy1;
 printf(" y1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y1 = %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1  );
 
 v1.v4 = vq;
 printf(" q = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" q = %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q  );

v1.v4 = vvz;
 printf(" vvz = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" vvz = %g, %g, %g, %g\n",
         vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz  );
 
 v2.v4 = vix;
 printf(" ix = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf(" ix = %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix  );
         
 v2.v4 = viy;
 printf(" iy = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iy = %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy  );

 printf( "\n" );

 #define CHECK(A) {\
   dvec a; \
   a.v4 = v##A; \
   if ( vp[i]. A !=  a.v[i] ) printf("(*error*) Bad %s !\n", #A ); }

 #define CHECKI(A)  {\
   ivec128 a; \
   a.v4 = v##A; \
   if ( vp[i]. A !=  a.v[i] ) printf("(*error*) Bad %s !\n", #A ); }


 for( i = 0; i < 4; i ++ ) {
   CHECK(x0);
   CHECK(x1);
   CHECK(y0);
   CHECK(y1);
   CHECK(q);
   CHECK(vz);
   CHECKI(ix);
   CHECKI(iy);
 }


 printf( "STORE4VP2D test\n" );

 __m256d t = _mm256_set1_pd(0.001);
 vx0 = _mm256_add_pd( vx0, t );
 vx1 = _mm256_add_pd( vx1, t );
 vy0 = _mm256_add_pd( vy0, t );
 vy1 = _mm256_add_pd( vy1, t );
 vq  = _mm256_add_pd( vq, t );
 vvz = _mm256_add_pd( vvz, t );
 
 vix = _mm_add_epi32( vix, _mm_set1_epi32( 100 ) );
 viy = _mm_add_epi32( viy, _mm_set1_epi32( 200 ) );
 
 STORE4VP2D( b, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );

 v1.v4 = vx0;
 printf(" x0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x0 = %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0  );

 v1.v4 = vx1;
 printf(" x1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x1 = %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1  );

 v1.v4 = vy0;
 printf(" y0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y0 = %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0  );
 
 v1.v4 = vy1;
 printf(" y1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y1 = %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1  );
 
 v1.v4 = vq;
 printf(" q = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" q = %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q  );

v1.v4 = vvz;
 printf(" vvz = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" vvz = %g, %g, %g, %g\n",
         vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz  );
 
 v2.v4 = vix;
 printf(" ix = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf(" ix = %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix  );
         
 v2.v4 = viy;
 printf(" iy = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iy = %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy  );

 printf( "\n" );

 for( i = 0; i < 4; i ++ ) {
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

 __m256d vx0, vx1;
 __m256d vy0, vy1;
 __m256d vz0, vz1;
 __m256d vq;
 __m128i vix, viy, viz;
 
 int i;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 3D VP macro tests tests \n" );
 printf( " size( t_vp3D ) = %zu \n", sizeof( t_vp3D ) );
 printf( " size( t_vp3D ) / sizeof( float ) = %zu \n", sizeof( t_vp3D ) / sizeof( double ) );
 printf( " ----------------------------------------------------------------------- \n" );

  for( i = 0; i < 4; i ++ ) {
   vp[i].x0 = -0.11 - 0.01*i; 
   vp[i].x1 =  0.11 + 0.01*i; 

   vp[i].y0 =  0.21 + 0.01*i;  
   vp[i].y1 = -0.21 - 0.01*i;  

   vp[i].z0 =  0.31 + 0.01*i;  
   vp[i].z1 = -0.31 - 0.01*i;  

   vp[i].q  = -0.01 - 0.01*i; 

   vp[i].ix = 10 + i;    
   vp[i].iy = 20 + i;
   vp[i].iz = 30 + i;
 }   
 
 
 double* b = (double *) vp;
 LOAD4VP3D( b, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
 
 dvec v1;
 ivec128 v2;

 printf( "LOAD4VP3D test\n" );

 v1.v4 = vx0;
 printf(" x0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x0 = %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0  );

 v1.v4 = vx1;
 printf(" x1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x1 = %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1  );

 v1.v4 = vy0;
 printf(" y0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y0 = %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0  );
 
 v1.v4 = vy1;
 printf(" y1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y1 = %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1  );

 v1.v4 = vz0;
 printf(" z0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" z0 = %g, %g, %g, %g\n",
         vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0  );
 
 v1.v4 = vz1;
 printf(" z1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" z1 = %g, %g, %g, %g\n",
         vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1  );

 
 v1.v4 = vq;
 printf(" q = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" q = %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q  );
 
 v2.v4 = vix;
 printf(" ix = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf(" ix = %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix );
         
 v2.v4 = viy;
 printf(" iy = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iy = %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy  );

 v2.v4 = viz;
 printf(" iz = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iz = %d, %d, %d, %d\n",
         vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz );


 printf( "\n" );


 for( i = 0; i < 4; i ++ ) {
   CHECK(x0);   CHECK(x1);
   CHECK(y0);   CHECK(y1);
   CHECK(z0);   CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }

 printf( "STORE4VP3D test\n" );

 __m256d t = _mm256_set1_pd(0.001);
 vx0 = _mm256_add_pd( vx0, t );
 vx1 = _mm256_add_pd( vx1, t );
 vy0 = _mm256_add_pd( vy0, t );
 vy1 = _mm256_add_pd( vy1, t );
 vz0 = _mm256_add_pd( vz0, t );
 vz1 = _mm256_add_pd( vz1, t );
 vq  = _mm256_add_pd( vq, t );
 
 vix = _mm_add_epi32( vix, _mm_set1_epi32( 100 ) );
 viy = _mm_add_epi32( viy, _mm_set1_epi32( 200 ) );
 viz = _mm_add_epi32( viz, _mm_set1_epi32( 300 ) );
 
 
 STORE4VP3D( b, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );
 
 v1.v4 = vx0;
 printf(" x0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x0 = %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0  );

 v1.v4 = vx1;
 printf(" x1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" x1 = %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1  );

 v1.v4 = vy0;
 printf(" y0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y0 = %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0  );
 
 v1.v4 = vy1;
 printf(" y1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" y1 = %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1  );

 v1.v4 = vz0;
 printf(" z0 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" z0 = %g, %g, %g, %g\n",
         vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0  );
 
 v1.v4 = vz1;
 printf(" z1 = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" z1 = %g, %g, %g, %g\n",
         vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1  );

 
 v1.v4 = vq;
 printf(" q = %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3]  );
 printf(" q = %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q  );
 
 v2.v4 = vix;
 printf(" ix = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );
 printf(" ix = %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix  );
         
 v2.v4 = viy;
 printf(" iy = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iy = %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy  );

 v2.v4 = viz;
 printf(" iz = %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3]  );

 printf(" iz = %d, %d, %d, %d\n",
         vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz  );


 printf( "\n" );


 for( i = 0; i < 4; i ++ ) {
   CHECK(x0);   CHECK(x1);
   CHECK(y0);   CHECK(y1);
   CHECK(z0);   CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }

 printf( " ----------------------------------------------------------------------- \n" );

}

int main( int argc, char* argv[])
{
/* 
  test_load_store();

  test_splines();

  test_trim();
*/
  test_hp();
/*
  test_load_store_vp2d();
  
  test_load_store_vp3d();
*/
}
