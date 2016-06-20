//Group 11 Members: Sagar Sharma (368058), Teodora Anitoaei (368014), Varun Gowtham  (368053)


#include <immintrin.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#define IDCT_SIZE         16
#define ITERATIONS        1000000
#define MAX_NEG_CROP      1024

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

static const short g_aiT16[16][16] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

static int64_t diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    int64_t d;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    d = temp.tv_sec*1000000000+temp.tv_nsec;
    return d;
}

static void compare_results(short *ref, short *res, const char *msg)
{
    int correct =1;

    printf("Comparing %s\n",msg);
    for(int j=0; j<IDCT_SIZE; j++)  {
        for(int i=0; i<IDCT_SIZE; i++){
            if(ref[j*IDCT_SIZE+i] != res[j*IDCT_SIZE+i]){
                correct=0;
                printf("failed at %d,%d\t ref=%d, res=%d\n ", i, j, ref[j*IDCT_SIZE+i],res[j*IDCT_SIZE+i]);
            }
        }
    }
    if (correct){
        printf("correct\n\n");
    }
}

// this function is for timing, do not change anything here
static void benchmark( void (*idct16)(short *, short *), short *input, short *output, const char *version )
{
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME,&start);

    for(int i=0;i<ITERATIONS;i++)
        idct16(input, output);

    clock_gettime(CLOCK_REALTIME,&end);
    double avg = (double) diff(start,end)/ITERATIONS;
    printf("%10s:\t %.3f ns\n", version, avg);
}

//scalar code for the inverse transform
static void partialButterflyInverse16(short *src, short *dst, int shift)
{
  int E[8],O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);

  for (int j=0; j<16; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (int k=0; k<8; k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ 16] + g_aiT16[ 3][k]*src[ 3*16] + g_aiT16[ 5][k]*src[ 5*16] + g_aiT16[ 7][k]*src[ 7*16] +
        g_aiT16[ 9][k]*src[ 9*16] + g_aiT16[11][k]*src[11*16] + g_aiT16[13][k]*src[13*16] + g_aiT16[15][k]*src[15*16];
    }
    //printf(" %d, %d, %d, %d, %d, %d, %d, %d ", O[0], O[1], O[2], O[3], O[4], O[5], O[6], O[7] );
    for (int k=0; k<4; k++)
    {
      EO[k] = g_aiT16[ 2][k]*src[ 2*16] + g_aiT16[ 6][k]*src[ 6*16] + g_aiT16[10][k]*src[10*16] + g_aiT16[14][k]*src[14*16];
    }
    EEO[0] = g_aiT16[4][0]*src[ 4*16 ] + g_aiT16[12][0]*src[ 12*16 ];
    EEE[0] = g_aiT16[0][0]*src[ 0    ] + g_aiT16[ 8][0]*src[  8*16 ];
    EEO[1] = g_aiT16[4][1]*src[ 4*16 ] + g_aiT16[12][1]*src[ 12*16 ];
    EEE[1] = g_aiT16[0][1]*src[ 0    ] + g_aiT16[ 8][1]*src[  8*16 ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    for (int k=0; k<2; k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }
    for (int k=0; k<4; k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }
    for (int k=0; k<8; k++)
    {
      dst[k]   = MAX( -32768, MIN( 32767, (E[k]   + O[k]   + add)>>shift ));
      dst[k+8] = MAX( -32768, MIN( 32767, (E[7-k] - O[7-k] + add)>>shift ));
    }
    src ++;
    dst += 16;
  }
}

static void simdPartialButterflyInverse16(short *src, short *dst, int shift)
{
  int E[8];
  int O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);
  __m128i *mm_O = (__m128i *) O;
  int temp_resO[8][4] __attribute__((aligned(16)));
  __m128i *temp_resmmO = (__m128i *) temp_resO;
  int temp_resEO[2][4] __attribute__((aligned(16)));
  __m128i *temp_resmmEO = (__m128i *) temp_resEO;
  int temp_resEEO[4] __attribute__((aligned(16)));
  __m128i *temp_resmmEEO = (__m128i *) temp_resEEO;
  int temp_resEE[4] __attribute__((aligned(16)));
  __m128i *temp_resmmEE = (__m128i *) temp_resEE;
  __m64 EEa;
  __m64 EEb;
  __m64 EEeven,EEodd;
  __m128i srcmm, srcmm1;
  __m128i tempEE;
  int final_res[8],final_res1[8];
  //__m128i temp_resmm, temp_resmm1;
  __m128i g_aiT16mm, g_aiT16mm1;
  __m128i mmsum,temp_mmsum, mmsum1;
  //__m128i *res_O = (__m128i *) O; 
  __m128i res_O, res_E, res_add, res_min, res_max,res_sum,res_sub,mmsub;
  __m128i *res_dst = (__m128i *) dst;
  int res_dst_sum[8] __attribute__((aligned(16)));
  int res_dst_sub[8] __attribute__((aligned(16)));
  __m128i *res_dst_summm = (__m128i *) res_dst_sum; 
  __m128i *res_dst_submm = (__m128i *) res_dst_sub;
  int max_val[4] __attribute__ ((aligned(32))) = {-32768,-32768,-32768,-32768};
  int min_val[4] __attribute__ ((aligned(32))) = {32767,32767,32767,32767};
  int add_val[4] __attribute__ ((aligned(32))) = {add,add,add,add};
  __m128i *a;
    a = (__m128i *) min_val;
    res_min = _mm_load_si128(&a[0]);  
    a = (__m128i *) max_val;       
    res_max = _mm_load_si128(&a[0]); 
    a = (__m128i *) add_val;
    res_add = _mm_load_si128(&a[0]); 
  ///////////////////////////////////////////
  int EO_vec[4][4] __attribute__((aligned(16)));
  __m128i *EO_vec_mm = (__m128i *) EO_vec;
  __m128i *EO_mm = (__m128i *) EO; 
  for (int j=0; j<16; j++)
  {  
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (int k=0; k<8; k++)
    { 
      srcmm = _mm_set_epi16 (src[16], src[3*16],src[5*16],src[7*16],src[9*16], src[11*16],src[13*16],src[15*16]);
      g_aiT16mm = _mm_set_epi16 (g_aiT16[1][k],g_aiT16[3][k],g_aiT16[5][k],g_aiT16[7][k],g_aiT16[9][k],g_aiT16[11][k],g_aiT16[13][k],g_aiT16[15][k]);
       _mm_store_si128(&temp_resmmO[k],_mm_madd_epi16 (srcmm, g_aiT16mm));
       O[k] = temp_resO[k][0]+ temp_resO[k][1]+ temp_resO[k][2]+ temp_resO[k][3];
    } 
     
    for (int k=0; k<2; k++)
    {
     srcmm = _mm_set_epi16 (src[2*16], src[6*16],src[10*16],src[14*16],src[2*16], src[6*16],src[10*16],src[14*16]);
     g_aiT16mm = _mm_set_epi16 (g_aiT16[2][2*k],g_aiT16[6][2*k],g_aiT16[10][2*k],g_aiT16[14][2*k],g_aiT16[2][1+(2*k)],g_aiT16[6][(2*k)+1],g_aiT16[10][(2*k)+1],g_aiT16[14][(2*k)+1]);
     _mm_store_si128(&temp_resmmEO[k],_mm_madd_epi16 (srcmm, g_aiT16mm));
     }
    
    srcmm = _mm_set_epi16 (src[4*16], src[12*16],src[4*16],src[12*16],src[0], src[8*16],src[0],src[8*16]);
    g_aiT16mm = _mm_set_epi16 (g_aiT16[4][0],g_aiT16[12][0],g_aiT16[4][1],g_aiT16[12][1],g_aiT16[0][0],g_aiT16[8][0],g_aiT16[0][1],g_aiT16[8][1]);
    _mm_store_si128(&temp_resmmEEO[0],_mm_madd_epi16 (srcmm, g_aiT16mm));
   
    
    E[0] = temp_resEEO[1] + temp_resEEO[3] + temp_resEO[0][2]+ temp_resEO[0][3];
    E[1] = temp_resEEO[0] + temp_resEEO[2] + temp_resEO[0][0]+ temp_resEO[0][1];
    E[2] = temp_resEEO[0] - temp_resEEO[2] + temp_resEO[1][2]+ temp_resEO[1][3];
    E[3] = temp_resEEO[1] - temp_resEEO[3] + temp_resEO[1][0]+ temp_resEO[1][1];
    E[4] = temp_resEEO[1] - temp_resEEO[3] - (temp_resEO[1][0]+ temp_resEO[1][1]);
    E[5] = temp_resEEO[0] - temp_resEEO[2] - (temp_resEO[1][2]+ temp_resEO[1][3]);
    E[6] = temp_resEEO[0] + temp_resEEO[2] - (temp_resEO[0][0]+ temp_resEO[0][1]);
    E[7] = temp_resEEO[1] + temp_resEEO[3] - (temp_resEO[0][2]+ temp_resEO[0][3]);

    __m128i *in_Emm = (__m128i *) E; 
    __m128i *in_Omm = (__m128i *) O;
    __m128i in_E, in_O, in_En, in_On;
    in_E = _mm_load_si128(&in_Emm[0]);  	
    in_O = _mm_load_si128(&in_Omm[0]);
    in_En = _mm_load_si128(&in_Emm[1]);  
    in_On = _mm_load_si128(&in_Omm[1]);
     
    mmsum = _mm_add_epi32 (in_O,in_E);
    mmsum = _mm_add_epi32 (mmsum,res_add);
    mmsum = _mm_srai_epi32 (mmsum, shift);
    mmsum = _mm_min_epi32 (mmsum , res_min);
    mmsum = _mm_max_epi32 (res_max , mmsum);
	
    mmsum1 = _mm_add_epi32(in_En, in_On);
    mmsum1 = _mm_add_epi32(mmsum1,res_add);
    mmsum1 =  _mm_srai_epi32 (mmsum1 , shift);
    mmsum1 = _mm_min_epi32 (mmsum1 , res_min);
    mmsum1 =	_mm_max_epi32 (res_max , mmsum1);
	
    mmsum =  _mm_packs_epi32(mmsum,mmsum1);
    _mm_store_si128(&res_dst[2*j],mmsum);
    in_En = _mm_shuffle_epi32 (in_En,27); 
    in_On = _mm_shuffle_epi32 (in_On,27);
    mmsum = _mm_sub_epi32(in_En,in_On); 
    mmsum = _mm_add_epi32(mmsum,res_add);
    mmsum =  _mm_srai_epi32 (mmsum , shift);

    mmsum = _mm_min_epi32 (mmsum , res_min);
    mmsum = _mm_max_epi32 (res_max , mmsum);

    in_E = _mm_shuffle_epi32 (in_E,27); 
    in_O = _mm_shuffle_epi32 (in_O,27);
    mmsum1 = _mm_sub_epi32(in_E,in_O);
    mmsum1 = _mm_add_epi32(mmsum1,res_add); 
    mmsum1 =  _mm_srai_epi32 (mmsum1 , shift);
    mmsum1 = _mm_min_epi32 (mmsum1 , res_min); 
    mmsum1 = _mm_max_epi32 (res_max , mmsum1);

    mmsum =  _mm_packs_epi32(mmsum,mmsum1);
    _mm_store_si128(&res_dst[j*2+1],mmsum);

    src ++;
  }
}
static void idct16_scalar(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  partialButterflyInverse16(pCoeff, tmp, 7);
  partialButterflyInverse16(tmp, pDst, 12);
}

/// CURRENTLY SAME CODE AS SCALAR !!
/// REPLACE HERE WITH SSE intrinsics
static void idct16_simd(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  simdPartialButterflyInverse16(pCoeff, tmp, 7);
  simdPartialButterflyInverse16(tmp, pDst, 12);
}

int main(int argc, char **argv)
{
    //allocate memory 16-byte aligned
    short *scalar_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *scalar_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    short *simd_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *simd_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    //initialize input
    printf("input array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
            short value = rand()%2 ? (rand()%32768) : -(rand()%32768) ;
            scalar_input[j*IDCT_SIZE+i] = value;
            simd_input  [j*IDCT_SIZE+i] = value;
	    printf("%d\t", value);
        }
        printf("\n");
    }

    idct16_scalar(scalar_input, scalar_output);
    idct16_simd  (simd_input  , simd_output);

    //check for correctness
    compare_results (scalar_output, simd_output, "scalar and simd");

    printf("output array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
	    printf("%d\t", scalar_output[j*IDCT_SIZE+i]);
        }
        printf("\n");
    }

    //Measure the performance of each kernel
    benchmark (idct16_scalar, scalar_input, scalar_output, "scalar");
    benchmark (idct16_simd, simd_input, simd_output, "simd");

    //cleanup
    free(scalar_input);    free(scalar_output);
    free(simd_input); free(simd_output);
}
