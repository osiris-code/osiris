/*****************************************************************************************
  Tests for AVX specific code
*****************************************************************************************/

#include <math.h>
#include <stdio.h>

#include "vector-avx.h"
#include "splines-avx.h"
#include "os-spec-push-avx.h"


/*****************************************************************************************
test_load_store
Tests the load / store macros
*****************************************************************************************/

#define CHECKBUF_PS( vec, b, offset, stride ) { \
  int i; \
  for( i = 0; i < 8; i++ ) { \
     if ( vec.v[i] != b[ offset + i * stride ] ) { \
       printf( "(*error*) %g /= %g\n", vec.v[i], b[ offset + i * stride ] ); \
       exit(-1); \
     } \
  } \
}

void test_load_store()
{
 DECLARE_ALIGNED_32( float b[8*3] );
 fvec v1, v2, v3;
 __m256 s;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store single tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 8*3 ; i ++ ) {
   b[i] = i;
 }
 
 _MM256_LOAD8v3_PS(v1.v8, v2.v8, v3.v8, b)
 
 printf(" _MM256_LOAD8v3_PS test \n" );

 printf("                        v1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21]  ); 

 CHECKBUF_PS( v1, b, 0, 3 );

 printf("                        v2 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 4, 7,10,13,16,19,22] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22]  ); 

 CHECKBUF_PS( v2, b, 1, 3 );

 printf("                        v3 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7]  );
 printf("b[ 2, 5, 8,11,14,17,20,23] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23]  ); 

 CHECKBUF_PS( v3, b, 2, 3 );


 // Change the __m128d values and store them
 s  = _mm256_set1_ps( 10.0f );
 v1.v8 = _mm256_add_ps( v1.v8, s );
 v2.v8 = _mm256_add_ps( v2.v8, s );
 v3.v8 = _mm256_add_ps( v3.v8, s );

 _MM256_STORE8v3_PS(b, v1.v8, v2.v8, v3.v8);
 printf(" _MM256_STORE8v3_PS test \n" );

 printf("                        v1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21]  ); 

 CHECKBUF_PS( v1, b, 0, 3 );

 printf("                        v2 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 4, 7,10,13,16,19,22] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22]  ); 

 CHECKBUF_PS( v2, b, 1, 3 );

 printf("                        v3 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7]  );
 printf("b[ 2, 5, 8,11,14,17,20,23] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23]  ); 

 CHECKBUF_PS( v3, b, 2, 3 );
 
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 8*2 ; i ++ ) {
   b[i] = i;
 }

 _MM256_LOAD8v2_PS( v1.v8, v2.v8, b );
 printf(" _MM256_LOAD8v2_PS test \n" );

 printf("                        v1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14]  ); 

 CHECKBUF_PS( v1, b, 0, 2 );

 printf("                        v2 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15]  ); 

 CHECKBUF_PS( v2, b, 1, 2 );

 // Change the __m128d values and store them
 s = _mm256_set1_ps( 10.0 );
 v1.v8 = _mm256_add_ps( v1.v8, s );
 v2.v8 = _mm256_add_ps( v2.v8, s );

 _MM256_STORE8v2_PS(b, v1.v8, v2.v8);
 printf(" _MM256_STORE8v2_PS test \n" );
 printf("                        v1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14]  ); 

 CHECKBUF_PS( v1, b, 0, 2 );

 printf("                        v2 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15] = %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15]  ); 

 CHECKBUF_PS( v2, b, 1, 2 );

 printf( " ----------------------------------------------------------------------- \n" );


}


/*****************************************************************************************
test_load_store_int
Tests the load / store of integer values
*****************************************************************************************/

#define CHECKBUF_EPI32( vec, b, offset, stride ) { \
  int i; \
  for( i = 0; i < 8; i++ ) { \
     if ( vec.v[i] != b[ offset + i * stride ] ) { \
       printf( "(*error*) %d /= %d\n", vec.v[i], b[ offset + i * stride ] ); \
       exit(-1); \
     } \
  } \
}


void test_load_store_int()
{
 DECLARE_ALIGNED_32( int b[8*3] );
 ivec v1, v2, v3;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store integer tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 8*2 ; i ++ ) {
   b[i] = i+1;
 }

 _MM256_LOAD8v2_EPI32( v1.v8, v2.v8, b )

 printf(" _MM256_LOAD8v2_PS test \n" );
 printf("                        v1 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14]  ); 

 CHECKBUF_EPI32( v1, b, 0, 2 );

 printf("                        v2 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15]  ); 

 CHECKBUF_EPI32( v2, b, 1, 2 );
 
 v1.v8 = viadd( v1.v8, _mm256_set1_epi32(10) );
 v2.v8 = viadd( v2.v8, _mm256_set1_epi32(20) );
 
 _MM256_STORE8v2_EPI32( b, v1.v8, v2.v8 )
 printf(" _MM256_STORE8v2_PS test \n" );
 printf("                        v1 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14]  ); 

 CHECKBUF_EPI32( v1, b, 0, 2 );

 printf("                        v2 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15]  ); 

 CHECKBUF_EPI32( v2, b, 1, 2 );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 8*3 ; i ++ ) {
   b[i] = i+1;
 }

 _MM256_LOAD8v3_EPI32( v1.v8, v2.v8, v3.v8, b )

 printf(" _MM256_LOAD8v3_EPI32 test \n" );

 printf("                        v1 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21]  ); 

 CHECKBUF_EPI32( v1, b, 0, 3 );

 printf("                        v2 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 4, 7,10,13,16,19,22] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22]  ); 

 CHECKBUF_EPI32( v2, b, 1, 3 );

 printf("                        v3 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7]  );
 printf("b[ 2, 5, 8,11,14,17,20,23] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23]  ); 
         
 CHECKBUF_EPI32( v3, b, 2, 3 );
         
 v1.v8 = viadd( v1.v8, _mm256_set1_epi32(10) );
 v2.v8 = viadd( v2.v8, _mm256_set1_epi32(20) );
 v3.v8 = viadd( v3.v8, _mm256_set1_epi32(30) );
 
 _MM256_STORE8v3_EPI32( b, v1.v8, v2.v8, v3.v8 )

 printf(" _MM256_STORE8v3_PS test \n" );

 printf("                        v1 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21]  ); 

 CHECKBUF_EPI32( v1, b, 0, 3 );

 printf("                        v2 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf("b[ 1, 4, 7,10,13,16,19,22] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22]  ); 

 CHECKBUF_EPI32( v2, b, 1, 3 );

 printf("                        v3 = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7]  );
 printf("b[ 2, 5, 8,11,14,17,20,23] = %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23]  ); 

 CHECKBUF_EPI32( v3, b, 2, 3 );

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_splines
Tests the spline calculation
*****************************************************************************************/
void test_splines()
{
 __m256 v, s[5];
 int i, j;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Spline tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v = _mm256_setr_ps( -0.5, -0.4, -0.2, 0.1, 0.0, 0.1, 0.2, 0.4 );

 for( i = 1; i <= 4; i++ ) {
    fvec a,b ;
 
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

    printf(" vsplined_s%1d test \n", i );
    a.v8 = v;
    printf(" v   = %g, %g, %g, %g, %g, %g, %g, %g\n", 
            a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7] );
    
    b.v8 = _mm256_setzero_ps();
    
    for( j = 0; j < i+1; j++ ) {
      a.v8 = s[j];
      b.v8 = _mm256_add_ps( b.v8, a.v8 );
      printf( " s[%2d] = %g, %g, %g, %g, %g, %g, %g, %g\n", j-i+1, 
               a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7]  );
    }
    printf(" tot = %g, %g, %g, %g, %g, %g, %g, %g\n", 
            b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7] );
    printf( "\n" );

 }

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_trim
Tests the trimming calculations
*****************************************************************************************/

inline __m256 vmntrim( const __m256 vx )
{
   register __m256 va, vb;
   
   va = _mm256_cmp_ps( vx, _mm256_set1_ps( -0.5f ), _CMP_LT_OS );
   va = _mm256_and_ps( va, _mm256_set1_ps( +1.0f ) );
   
   vb = _mm256_cmp_ps( vx, _mm256_set1_ps( +0.5f ), _CMP_GE_OS );
   vb = _mm256_and_ps( vb, _mm256_set1_ps( +1.0f ) );
   
   return  _mm256_sub_ps( vb, va );
}

void test_trim()
{
 __m256 v1, v2;
 fvec a, b;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Trim tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v1 = _mm256_setr_ps( -0.6, -0.5, -0.4, -0.3, 0.3, 0.4, 0.5, 0.6 );
 v2 = vmntrim( v1 );
 
 a.v8 = v1; b.v8 = v2;
 
 printf( " v       = %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7] );
 printf( " vntrim  = %g, %g, %g, %g, %g, %g, %g, %g\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7] );

 printf( "\n" );

 printf( " ----------------------------------------------------------------------- \n" );

}




/*****************************************************************************************
test_hp
Tests the half-point calculations
*****************************************************************************************/

#define vget_hp_near(dx, ix, dxh, ixh, delta ) 	{   				                     \
   register __m256 cmp;                                                                  \
                                                                                         \
   cmp  = _mm256_cmp_ps( dx, _mm256_setzero_ps(), _CMP_LT_OS );                          \
   dxh =  _mm256_add_ps( dx, _mm256_blendv_ps( _mm256_set1_ps( -0.5f ),                  \
											   _mm256_set1_ps(  0.5f ), cmp ) );         \
   ixh = visub( ix, _mm256_castps_si256(                                                 \
					   _mm256_and_ps(                                                    \
						  cmp,                                                           \
						  _mm256_castsi256_ps( _mm256_set1_epi32( delta ) )              \
					   )                                                                 \
					) );                                                                 \
}


void test_hp()
{
 __m256 dx, dxh;
 __m256i ix, ixh;

 int delta;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Half-Point tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 dx   = _mm256_set_ps( -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4 );
 ix   = _mm256_set_epi32( 10, 20, 30, 40, 50, 60, 70, 80 );

 delta = 3;

 ixh   = _mm256_set1_epi32( -1 );

 vget_hp_near( dx, ix, dxh, ixh, delta );

 fvec a;
 ivec b;

 a.v8 = dx;
 printf( " dx      = %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7] );
 a.v8 = dxh;
 printf( " dxh     = %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7] );
 
 b.v8 = ix;
 printf( " ix      = %d, %d, %d, %d, %d, %d, %d, %d\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7] );
 
 printf( " delta   = %d \n", delta );
 
 b.v8 = ixh;
 printf( " ixh     = %d, %d, %d, %d, %d, %d, %d, %d\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7] );

 printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_load_store_vp
Tests the load/store virtual particle macros
*****************************************************************************************/



void test_load_store_vp2d()
{

 DECLARE_ALIGNED_32( t_vp2D vp[8] );

 __m256 vx0, vx1;
 __m256 vy0, vy1;
 __m256 vq, vvz;
 __m256i vix, viy;
 
 int i;
  
 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 2D VP macro tests tests \n" );
 printf( " size( t_vp2D ) = %zu \n", sizeof( vp[0] ) );
 printf( " size( t_vp2D ) / sizeof( float ) = %zu \n", sizeof( vp[0] ) / sizeof( float ) );
 printf( " ----------------------------------------------------------------------- \n" );
 
 for( i = 0; i < 8; i ++ ) {
   vp[i].x0 = -0.11 - 0.01*i; 
   vp[i].x1 =  0.11 + 0.01*i; 

   vp[i].y0 =  0.21 + 0.01*i;  
   vp[i].y1 = -0.21 - 0.01*i;  

   vp[i].q  = -0.01 - 0.01*i; 
   vp[i].vz =  0.91 + 0.01*i; 

   vp[i].ix = 10 + i;    
   vp[i].iy = 20 + i;
 }   


 LOAD8VP2D( ((float *) vp), vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

 printf( "LOAD8VP2D test\n" );

 fvec v1;
 ivec v2;

 v1.v8 = vx0;
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
         vp[4].x0, vp[5].x0, vp[6].x0, vp[7].x0  );


 v1.v8 = vx1;
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
         vp[4].x1, vp[5].x1, vp[6].x1, vp[7].x1  );

 v1.v8 = vy0;
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
         vp[4].y0, vp[5].y0, vp[6].y0, vp[7].y0  );
 
 v1.v8 = vy1;
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
         vp[4].y1, vp[5].y1, vp[6].y1, vp[7].y1  );
 
 v1.v8 = vq;
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q,
         vp[4].q, vp[5].q, vp[6].q, vp[7].q  );

v1.v8 = vvz;
 printf(" vvz = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" vvz = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz,
         vp[4].vz, vp[5].vz, vp[6].vz, vp[7].vz  );
 
 v2.v8 = vix;
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
         vp[4].ix, vp[5].ix, vp[6].ix, vp[7].ix  );
         
 v2.v8 = viy;
 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
         vp[4].iy, vp[5].iy, vp[6].iy, vp[7].iy  );

 printf( "\n" );

 #define CHECK(A) {\
   fvec a; \
   a.v8 = v##A; \
   if ( vp[i]. A !=  a.v[i] ) printf("(*error*) Bad %s !\n", #A ); }

 #define CHECKI(A)  {\
   ivec a; \
   a.v8 = v##A; \
   if ( vp[i]. A !=  a.v[i] ) printf("(*error*) Bad %s !\n", #A ); }


 for( i = 0; i < 8; i ++ ) {
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

 __m256 t = _mm256_set1_ps(0.001);
 vx0 = _mm256_add_ps( vx0, t );
 vx1 = _mm256_add_ps( vx1, t );
 vy0 = _mm256_add_ps( vy0, t );
 vy1 = _mm256_add_ps( vy1, t );
 vq  = _mm256_add_ps( vq, t );
 vvz = _mm256_add_ps( vvz, t );
 
 vix = viadd( vix, _mm256_set1_epi32( 100 ) );
 viy = viadd( viy, _mm256_set1_epi32( 200 ) );
 
 STORE8VP2D( ((float *) vp), vx0, vx1, vy0, vy1, vq, vvz, vix, viy );

 v1.v8 = vx0;
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
         vp[4].x0, vp[5].x0, vp[6].x0, vp[7].x0  );

 v1.v8 = vx1;
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
         vp[4].x1, vp[5].x1, vp[6].x1, vp[7].x1  );

 v1.v8 = vy0;
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
         vp[4].y0, vp[5].y0, vp[6].y0, vp[7].y0  );
 
 v1.v8 = vy1;
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
         vp[4].y1, vp[5].y1, vp[6].y1, vp[7].y1  );
 
 v1.v8 = vq;
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q,
         vp[4].q, vp[5].q, vp[6].q, vp[7].q  );

v1.v8 = vvz;
 printf(" vvz = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" vvz = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].vz, vp[1].vz, vp[2].vz, vp[3].vz,
         vp[4].vz, vp[5].vz, vp[6].vz, vp[7].vz  );
 
 v2.v8 = vix;
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
         vp[4].ix, vp[5].ix, vp[6].ix, vp[7].ix  );
         
 v2.v8 = viy;
 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
         vp[4].iy, vp[5].iy, vp[6].iy, vp[7].iy  );

 printf( "\n" );

 for( i = 0; i < 8; i ++ ) {
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

 t_vp3D vp[8];

 __m256 vx0, vx1;
 __m256 vy0, vy1;
 __m256 vz0, vz1;
 __m256 vq;
 __m256i vix, viy, viz;
 
 int i;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store 3D VP macro tests tests \n" );
 printf( " size( t_vp3D ) = %zu \n", sizeof( vp[0] ) );
 printf( " size( t_vp3D ) / sizeof( float ) = %zu \n", sizeof( vp[0] ) / sizeof( float ) );
 printf( " ----------------------------------------------------------------------- \n" );

  for( i = 0; i < 8; i ++ ) {
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
 

 LOAD8VP3D( ((float *)vp), vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );

 fvec v1;
 ivec v2;

 printf( "LOAD8VP3D test\n" );

v1.v8 = vx0;
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
         vp[4].x0, vp[5].x0, vp[6].x0, vp[7].x0  );

 v1.v8 = vx1;
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
         vp[4].x1, vp[5].x1, vp[6].x1, vp[7].x1  );

 v1.v8 = vy0;
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
         vp[4].y0, vp[5].y0, vp[6].y0, vp[7].y0  );
 
 v1.v8 = vy1;
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
         vp[4].y1, vp[5].y1, vp[6].y1, vp[7].y1  );

 v1.v8 = vz0;
 printf(" z0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" z0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0,
         vp[4].z0, vp[5].z0, vp[6].z0, vp[7].z0  );
 
 v1.v8 = vz1;
 printf(" z1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" z1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1,
         vp[4].z1, vp[5].z1, vp[6].z1, vp[7].z1  );

 
 v1.v8 = vq;
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q,
         vp[4].q, vp[5].q, vp[6].q, vp[7].q  );
 
 v2.v8 = vix;
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
         vp[4].ix, vp[5].ix, vp[6].ix, vp[7].ix  );
         
 v2.v8 = viy;
 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
         vp[4].iy, vp[5].iy, vp[6].iy, vp[7].iy  );

 v2.v8 = viz;
 printf(" iz = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iz = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz,
         vp[4].iz, vp[5].iz, vp[6].iz, vp[7].iz  );


 printf( "\n" );


 for( i = 0; i < 8; i ++ ) {
   CHECK(x0);   CHECK(x1);
   CHECK(y0);   CHECK(y1);
   CHECK(z0);   CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }

 
 printf( "STORE8VP3D test\n" );

 __m256 t = _mm256_set1_ps(0.001);
 vx0 = _mm256_add_ps( vx0, t );
 vx1 = _mm256_add_ps( vx1, t );
 vy0 = _mm256_add_ps( vy0, t );
 vy1 = _mm256_add_ps( vy1, t );
 vz0 = _mm256_add_ps( vz0, t );
 vz1 = _mm256_add_ps( vz1, t );
 vq  = _mm256_add_ps( vq, t );
 
 vix = viadd( vix, _mm256_set1_epi32( 100 ) );
 viy = viadd( viy, _mm256_set1_epi32( 200 ) );
 viz = viadd( viz, _mm256_set1_epi32( 300 ) );
 
 STORE8VP3D( ((float *)vp), vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );

 v1.v8 = vx0;
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x0, vp[1].x0, vp[2].x0, vp[3].x0,
         vp[4].x0, vp[5].x0, vp[6].x0, vp[7].x0  );

 v1.v8 = vx1;
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" x1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].x1, vp[1].x1, vp[2].x1, vp[3].x1,
         vp[4].x1, vp[5].x1, vp[6].x1, vp[7].x1  );

 v1.v8 = vy0;
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y0, vp[1].y0, vp[2].y0, vp[3].y0,
         vp[4].y0, vp[5].y0, vp[6].y0, vp[7].y0  );
 
 v1.v8 = vy1;
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" y1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].y1, vp[1].y1, vp[2].y1, vp[3].y1,
         vp[4].y1, vp[5].y1, vp[6].y1, vp[7].y1  );

 v1.v8 = vz0;
 printf(" z0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" z0 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].z0, vp[1].z0, vp[2].z0, vp[3].z0,
         vp[4].z0, vp[5].z0, vp[6].z0, vp[7].z0  );
 
 v1.v8 = vz1;
 printf(" z1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" z1 = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].z1, vp[1].z1, vp[2].z1, vp[3].z1,
         vp[4].z1, vp[5].z1, vp[6].z1, vp[7].z1  );

 
 v1.v8 = vq;
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7]  );
 printf(" q = %g, %g, %g, %g, %g, %g, %g, %g\n",
         vp[0].q, vp[1].q, vp[2].q, vp[3].q,
         vp[4].q, vp[5].q, vp[6].q, vp[7].q  );
 
 v2.v8 = vix;
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );
 printf(" ix = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].ix, vp[1].ix, vp[2].ix, vp[3].ix,
         vp[4].ix, vp[5].ix, vp[6].ix, vp[7].ix  );
         
 v2.v8 = viy;
 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iy = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iy, vp[1].iy, vp[2].iy, vp[3].iy,
         vp[4].iy, vp[5].iy, vp[6].iy, vp[7].iy  );

 v2.v8 = viz;
 printf(" iz = %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7]  );

 printf(" iz = %d, %d, %d, %d, %d, %d, %d, %d\n",
         vp[0].iz, vp[1].iz, vp[2].iz, vp[3].iz,
         vp[4].iz, vp[5].iz, vp[6].iz, vp[7].iz  );


 printf( "\n" );


 for( i = 0; i < 8; i ++ ) {
   CHECK(x0);   CHECK(x1);
   CHECK(y0);   CHECK(y1);
   CHECK(z0);   CHECK(z1);
   CHECK(q);
   CHECKI(ix);  CHECKI(iy);  CHECKI(iz);
 }


 printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_math
Tests additional math routines (i.e. integer math)
*****************************************************************************************/


void test_math()
{
 __m256i v1, v2, v3;
 ivec a, b, c;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Integer Math tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v1 = _mm256_set_epi32( 7, 6, 5, 4, 3, 2, 1, 0 );
 v2 = _mm256_set_epi32( 14, 12, 10, 8, 6, 4, 2, 0 );
 v3 = _mm256_set_epi32( 21, 18, 15, 12, 9, 6, 3, 0 );
  
 printf(  "vimul_s test:\n");
 printiv( "    v1 ", v1 );
 printiv( "v1 * 2 ", vimul_s( v1, 2 ) ); 
 printf( "\n" );

 printf(  "viadd test:\n");
 printiv( "     v1 ", v1 );
 printiv( "     v2 ", v2 );
 printiv( "v1 + v2 ", viadd( v1, v2 ) ); 
 printf( "\n" );

 printf(  "viadd3 test:\n");
 printiv( "          v1 ", v1 );
 printiv( "          v2 ", v2 );
 printiv( "          v3 ", v3 );
 printiv( "v1 + v2 + v3 ", viadd3( v1, v2, v3 ) ); 
 printf( "\n" );

 DECLARE_ALIGNED_32( int buf[8] );
 viadd3_store( buf, v1, v2, v3 );

 printf(  "viadd3 test:\n");
 printiv( "          v1 ", v1 );
 printiv( "          v2 ", v2 );
 printiv( "          v3 ", v3 );
 printf(  "v1 + v2 + v3  = %d, %d, %d, %d, %d, %d, %d, %d\n", 
          buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]); 
 printf( "\n" );

 printf( " ----------------------------------------------------------------------- \n" );

}


int main( int argc, char* argv[])
{
 
  test_load_store();

  test_load_store_int();

  test_splines();

  test_trim();

  test_hp();

  test_load_store_vp2d();
  
  test_load_store_vp3d();
  
  test_math();

}
