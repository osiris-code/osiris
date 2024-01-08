/*****************************************************************************************
  Tests for MIC specific code
*****************************************************************************************/

#include <math.h>
#include <stdio.h>

#include "vector-mic.h"
#include "splines-mic.h"
#include "os-spec-push-mic.h"


/*****************************************************************************************
test_load_store
Tests the load / store macros
*****************************************************************************************/

#define CHECKBUF_PS( vec, b, offset, stride ) { \
  int i; \
  for( i = 0; i < 16; i++ ) { \
     if ( vec.v[i] != b[ offset + i * stride ] ) { \
       printf( "(*error*) %g /= %g\n", vec.v[i], b[ offset + i * stride ] ); \
       exit(-1); \
     } \
  } \
}

void test_load_store()
{
 DECLARE_ALIGNED_64( float b[16*3] );
 fvec v1, v2, v3;
 __m512 s;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store single tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 16*3 ; i ++ ) {
   b[i] = i;
 }
 
 _MM512_LOAD16v3_PS(v1.v16, v2.v16, v3.v16, b)
 
 printf(" _MM512_LOAD16v3_PS test \n" );

 printf("                                                v1 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,42,45] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21],
         b[24], b[37], b[30], b[33],
         b[36], b[39], b[42], b[45]  ); 

 CHECKBUF_PS( v1, b, 0, 3 );

 printf("                                                v2 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 4, 7,10,13,16,19,22,25,28,31,34,37,40,43,46] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22],
         b[25], b[28], b[31], b[34],
         b[37], b[40], b[43], b[46]  ); 

 CHECKBUF_PS( v2, b, 1, 3 );

 printf("                                                v3 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7],
         v3.v[8], v3.v[9], v3.v[10], v3.v[11],
         v3.v[12], v3.v[13], v3.v[14], v3.v[15]  );
 printf("b[ 2, 5, 8,11,14,17,20,23,26,29,32,35,38,41,44,47] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23],
         b[26], b[29], b[32], b[35],
         b[38], b[41], b[44], b[47]  ); 

 CHECKBUF_PS( v3, b, 2, 3 );


 // Change the __m512 values and store them
 s  = _mm512_set1_ps( 10.0f );
 v1.v16 = _mm512_add_ps( v1.v16, s );
 v2.v16 = _mm512_add_ps( v2.v16, s );
 v3.v16 = _mm512_add_ps( v3.v16, s );

 _MM512_STORE16v3_PS(b, v1.v16, v2.v16, v3.v16);
 printf(" _MM512_STORE16v3_PS test \n" );

 printf("                                                v1 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,42,45] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21],
         b[24], b[37], b[30], b[33],
         b[36], b[39], b[42], b[45]  ); 

 CHECKBUF_PS( v1, b, 0, 3 );

 printf("                                                v2 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 4, 7,10,13,16,19,22,25,28,31,34,37,40,43,46] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22],
         b[25], b[28], b[31], b[34],
         b[37], b[40], b[43], b[46]  ); 

 CHECKBUF_PS( v2, b, 1, 3 );

 printf("                                                v3 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7],
         v3.v[8], v3.v[9], v3.v[10], v3.v[11],
         v3.v[12], v3.v[13], v3.v[14], v3.v[15]  );
 printf("b[ 2, 5, 8,11,14,17,20,23,26,29,32,35,38,41,44,47] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23],
         b[26], b[29], b[32], b[35],
         b[38], b[41], b[44], b[47]  ); 

 CHECKBUF_PS( v3, b, 2, 3 );
 
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 16*2 ; i ++ ) {
   b[i] = i;
 }

 _MM512_LOAD16v2_PS( v1.v16, v2.v16, b );
 printf(" _MM512_LOAD16v2_PS test \n" );

 printf("                                                v1 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14],
         b[16], b[18], b[20], b[22],
         b[24], b[26], b[28], b[30]  ); 

 CHECKBUF_PS( v1, b, 0, 2 );

 printf("                                                v2 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15],
         b[17], b[19], b[21], b[23],
         b[25], b[27], b[29], b[31]  ); 

 CHECKBUF_PS( v2, b, 1, 2 );

 // Change the __m128d values and store them
 s = _mm512_set1_ps( 10.0 );
 v1.v16 = _mm512_add_ps( v1.v16, s );
 v2.v16 = _mm512_add_ps( v2.v16, s );

 _MM512_STORE16v2_PS(b, v1.v16, v2.v16);

 printf("                                                v1 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14],
         b[16], b[18], b[20], b[22],
         b[24], b[26], b[28], b[30]  ); 

 CHECKBUF_PS( v1, b, 0, 2 );

 printf("                                                v2 = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15],
         b[17], b[19], b[21], b[23],
         b[25], b[27], b[29], b[31]  ); 

 CHECKBUF_PS( v2, b, 1, 2 );

 printf( " ----------------------------------------------------------------------- \n" );


}


/*****************************************************************************************
test_load_store_int
Tests the load / store of integer values
*****************************************************************************************/

#define CHECKBUF_EPI32( vec, b, offset, stride ) { \
  int i; \
  for( i = 0; i < 16; i++ ) { \
     if ( vec.v[i] != b[ offset + i * stride ] ) { \
       printf( "(*error*) %d /= %d\n", vec.v[i], b[ offset + i * stride ] ); \
       exit(-1); \
     } \
  } \
}


void test_load_store_int()
{
 DECLARE_ALIGNED_64( int b[16*3] );
 ivec v1, v2, v3;
 int i;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Load/store integer tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 16*2 ; i ++ ) {
   b[i] = i+1;
 }

 _MM512_LOAD16v2_EPI32( v1.v16, v2.v16, b )

 printf(" _MM512_LOAD16v2_EPI32 test \n" );
 printf("                                                v1 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11],
         v1.v[12], v1.v[13], v1.v[14], v1.v[15] );
         
 printf("b[ 0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14],
         b[16], b[18], b[20], b[22],
         b[24], b[26], b[28], b[30]  ); 

 CHECKBUF_EPI32( v1, b, 0, 2 );

 printf("                                                v2 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15] );
 printf("b[ 1, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15],
         b[17], b[19], b[21], b[23],
         b[25], b[27], b[29], b[31]  ); 

 CHECKBUF_EPI32( v2, b, 1, 2 );
 
 v1.v16 = _mm512_add_epi32( v1.v16, _mm512_set1_epi32(10) ); 
 v2.v16 = _mm512_add_epi32( v2.v16, _mm512_set1_epi32(20) ); 
 
 _MM512_STORE16v2_EPI32( b, v1.v16, v2.v16 )
 printf(" _MM512_STORE16v2_PS test \n" );
 printf("                                                v1 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3],
         v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11],
         v1.v[12], v1.v[13], v1.v[14], v1.v[15] );
         
 printf("b[ 0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[2], b[4], b[6],
         b[8], b[10], b[12], b[14],
         b[16], b[18], b[20], b[22],
         b[24], b[26], b[28], b[30]  ); 

 CHECKBUF_EPI32( v1, b, 0, 2 );

 printf("                                                v2 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15] );
 printf("b[ 1, 3, 5, 7, 9,11,13,15,17,19,21,23,25,27,29,31] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[3], b[5], b[7],
         b[9], b[11], b[13], b[15],
         b[17], b[19], b[21], b[23],
         b[25], b[27], b[29], b[31]  ); 

 printf( " ----------------------------------------------------------------------- \n" );

 // Initialize buffer values
 for( i = 0; i < 16*3 ; i ++ ) {
   b[i] = i+1;
 }

 _MM512_LOAD16v3_EPI32( v1.v16, v2.v16, v3.v16, b )

 printf(" _MM256_LOAD16v3_EPI32 test \n" );

 printf("                                                v1 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,42,45] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21],
         b[24], b[37], b[30], b[33],
         b[36], b[39], b[42], b[45]  ); 

 CHECKBUF_EPI32( v1, b, 0, 3 );

 printf("                                                v2 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 4, 7,10,13,16,19,22,25,28,31,34,37,40,43,46] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22],
         b[25], b[28], b[31], b[34],
         b[37], b[40], b[43], b[46]  ); 

 CHECKBUF_EPI32( v2, b, 1, 3 );

 printf("                                                v3 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7],
         v3.v[8], v3.v[9], v3.v[10], v3.v[11],
         v3.v[12], v3.v[13], v3.v[14], v3.v[15]  );
 printf("b[ 2, 5, 8,11,14,17,20,23,26,29,32,35,38,41,44,47] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23],
         b[26], b[29], b[32], b[35],
         b[38], b[41], b[44], b[47]  ); 

 CHECKBUF_EPI32( v3, b, 2, 3 );
         
 v1.v16 = _mm512_add_epi32( v1.v16, _mm512_set1_epi32(10) ); 
 v2.v16 = _mm512_add_epi32( v2.v16, _mm512_set1_epi32(20) ); 
 v3.v16 = _mm512_add_epi32( v3.v16, _mm512_set1_epi32(30) ); 
 
 _MM512_STORE16v3_EPI32( b, v1.v16, v2.v16, v3.v16 )

 printf(" _MM512_STORE16v3_EPI32 test \n" );

 printf("                                                v1 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v1.v[0], v1.v[1], v1.v[2], v1.v[3], v1.v[4], v1.v[5], v1.v[6], v1.v[7],
         v1.v[8], v1.v[9], v1.v[10], v1.v[11], v1.v[12], v1.v[13], v1.v[14], v1.v[15]  );
 printf("b[ 0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,42,45] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[0], b[3], b[6], b[9],
         b[12], b[15], b[18], b[21],
         b[24], b[37], b[30], b[33],
         b[36], b[39], b[42], b[45]  ); 

 CHECKBUF_EPI32( v1, b, 0, 3 );

 printf("                                                v2 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v2.v[0], v2.v[1], v2.v[2], v2.v[3],
         v2.v[4], v2.v[5], v2.v[6], v2.v[7],
         v2.v[8], v2.v[9], v2.v[10], v2.v[11],
         v2.v[12], v2.v[13], v2.v[14], v2.v[15]  );
 printf("b[ 1, 4, 7,10,13,16,19,22,25,28,31,34,37,40,43,46] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[1], b[4], b[7], b[10],
         b[13], b[16], b[19], b[22],
         b[25], b[28], b[31], b[34],
         b[37], b[40], b[43], b[46]  ); 

 CHECKBUF_EPI32( v2, b, 1, 3 );

 printf("                                                v3 = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         v3.v[0], v3.v[1], v3.v[2], v3.v[3],
         v3.v[4], v3.v[5], v3.v[6], v3.v[7],
         v3.v[8], v3.v[9], v3.v[10], v3.v[11],
         v3.v[12], v3.v[13], v3.v[14], v3.v[15]  );
 printf("b[ 2, 5, 8,11,14,17,20,23,26,29,32,35,38,41,44,47] = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
         b[2], b[5], b[8], b[11],
         b[14], b[17], b[20], b[23],
         b[26], b[29], b[32], b[35],
         b[38], b[41], b[44], b[47]  ); 

 CHECKBUF_EPI32( v3, b, 2, 3 );

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_splines
Tests the spline calculation
*****************************************************************************************/
void test_splines()
{
 __m512 v, s[5];
 int i, j;

 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Spline tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v = _mm512_setr_ps( -0.5, -0.4, -0.2, 0.1, 0.0, 0.1, 0.2, 0.4, 0., 0., 0., 0., 0., 0., 0., 0. );

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
    a.v16 = v;
    printf(" v   = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
            a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7],
            a.v[8], a.v[9], a.v[10], a.v[11], a.v[12], a.v[13], a.v[14], a.v[15] );
    
    b.v16 = _mm512_setzero_ps();
    
    for( j = 0; j < i+1; j++ ) {
      a.v16 = s[j];
      b.v16 = _mm512_add_ps( b.v16, a.v16 );
      printf( " s[%2d] = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", j-i+1, 
               a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7],
               a.v[8], a.v[9], a.v[10], a.v[11], a.v[12], a.v[13], a.v[14], a.v[15]  );
    }
    printf(" tot = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
            b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7],
            b.v[8], b.v[9], b.v[10], b.v[11], b.v[12], b.v[13], b.v[14], b.v[15] );
    printf( "\n" );

 }

 printf( " ----------------------------------------------------------------------- \n" );

}


/*****************************************************************************************
test_trim
Tests the trimming calculations
*****************************************************************************************/

#define VMNTRIM( vx, trim, cross ) { \
                                                                                         \
  trim = _mm512_setzero_ps();                                                            \
                                                                                         \
  __mmask16 mu = _mm512_cmp_ps_mask( vx, _mm512_set1_ps(  0.5f ), _MM_CMPINT_GE );       \
  __mmask16 ml = _mm512_cmp_ps_mask( vx, _mm512_set1_ps( -0.5f ), _MM_CMPINT_LT );       \
                                                                                         \
  trim = _mm512_mask_mov_ps( trim, mu, _mm512_set1_ps( +1.0f ) );                        \
  trim = _mm512_mask_mov_ps( trim, ml, _mm512_set1_ps( -1.0f ) );                        \
                                                                                         \
  cross = _mm512_kor( mu, ml );                                                          \
}

void test_trim()
{
 __m512 v1, v2;
 __mmask16 cross;
  fvec a, b;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Trim tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 v1 = _mm512_setr_ps( -0.6, -0.5, -0.4, -0.3, 0.3, 0.4, 0.5, 0.6, 
                       0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0 );
 
 VMNTRIM( v1, v2, cross );
 
 a.v16 = v1; b.v16 = v2;
 
 printf( " v       = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7],
         a.v[8], a.v[9], a.v[10], a.v[11], a.v[12], a.v[13], a.v[14], a.v[15] );
 printf( " vntrim  = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7],
         b.v[8], b.v[9], b.v[10], b.v[11], b.v[12], b.v[13], b.v[14], b.v[15] );

 printf( "\n" );

 printf( " ----------------------------------------------------------------------- \n" );

}




/*****************************************************************************************
test_hp
Tests the half-point calculations
*****************************************************************************************/
#define VGET_HP_NEAR(dx, ix, dxh, ixh, delta ) 	{   				                     \
 __mmask16 cmp = _mm512_cmp_ps_mask( dx, _mm512_setzero_ps(), _MM_CMPINT_LT );           \
                                                                                         \ 
 dxh = _mm512_add_ps( dx, _mm512_mask_mov_ps( _mm512_set1_ps( -0.5f ), cmp,              \
                                              _mm512_set1_ps( +0.5f ) ) );               \
 ixh = _mm512_mask_sub_epi32( ix, cmp, ix, _mm512_set1_epi32( delta ) );                 \
}




void test_hp()
{
 __m512 dx, dxh;
 __m512i ix, ixh;

 int delta;


 printf( " ----------------------------------------------------------------------- \n" );
 printf( " Half-Point tests \n" );
 printf( " ----------------------------------------------------------------------- \n" );

 dx   = _mm512_set_ps( -0.5, -0.4, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4,
                        0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0 );
 ix   = _mm512_set_epi32( 10, 20, 30, 40, 50, 60, 70, 80 ,
                           0,  0,  0,  0,  0,  0,  0,  0 );

 delta = 3;

 ixh   = _mm512_set1_epi32( -1 );

 VGET_HP_NEAR( dx, ix, dxh, ixh, delta );

 fvec a;
 ivec b;

 a.v16 = dx;
 printf( " dx      = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7],
         a.v[8], a.v[9], a.v[10], a.v[11], a.v[12], a.v[13], a.v[14], a.v[5] );
 a.v16 = dxh;
 printf( " dxh     = %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", 
         a.v[0], a.v[1], a.v[2], a.v[3], a.v[4], a.v[5], a.v[6], a.v[7],
         a.v[8], a.v[9], a.v[10], a.v[11], a.v[12], a.v[13], a.v[14], a.v[5] );
 
 b.v16 = ix;
 printf( " ix      = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7],
         b.v[8], b.v[9], b.v[10], b.v[11], b.v[12], b.v[13], b.v[14], b.v[15] );
 
 printf( " delta   = %d \n", delta );
 
 b.v16 = ixh;
 printf( " ixh     = %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
         b.v[0], b.v[1], b.v[2], b.v[3], b.v[4], b.v[5], b.v[6], b.v[7],
         b.v[8], b.v[9], b.v[10], b.v[11], b.v[12], b.v[13], b.v[14], b.v[15] );

 printf( " ----------------------------------------------------------------------- \n" );

}

/*****************************************************************************************
test_load_store_vp
Tests the load/store virtual particle macros
*****************************************************************************************/

/* These are not required in the MIC code because */




int main( int argc, char* argv[])
{
 
  test_load_store();

  test_load_store_int();

  test_splines();

  test_trim();

  test_hp();

}
