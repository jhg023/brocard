#include <cstdint>
#include <immintrin.h>

auto mulmod_vectorized( __m256i a, uint64_t b, __m256i m ) -> __m256i {
  __m256i d = _mm256_setzero_si256();
  __m256i b_vector = _mm256_set1_epi64x( b );
  __m256i mp2 = _mm256_srli_epi64( m, 1 );

  __m256i high_bit = _mm256_set1_epi64x( 0x8000000000000000ULL );

  for ( int i = 0; i < 64; ++i ) {
    __m256i d_mp2_gt = _mm256_cmpgt_epi64( d, mp2 );
    d = _mm256_slli_epi64( d, 1 );
    __m256i d_sub_m = _mm256_sub_epi64( d, m );
    d = _mm256_blendv_pd( d, d_sub_m, d_mp2_gt );

    __m256i a_high_bit = _mm256_and_si256( a, high_bit );
    __m256i d_add = _mm256_add_epi64( d, b_vector );
    d = _mm256_blendv_pd( d, d_add, a_high_bit );

    __m256i gt = _mm256_cmpgt_epi64( d, m );
    __m256i eq = _mm256_cmpeq_epi64( d, m );
    __m256i gt_eq = _mm256_or_si256( gt, eq );
    __m256i d_sub = _mm256_sub_epi64( d, m );
    d = _mm256_blendv_pd( d, d_sub, gt_eq );

    a = _mm256_slli_epi64( a, 1 );
  }
  
  return d;
}