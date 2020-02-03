#include <cstdint>
#include <cstdio>
#include <immintrin.h>

__m256i jacobi_vector( int64_t a, __m256i b ) {
  auto ones = _mm256_set1_epi64x(0xffffffffffffffff);
  auto a_v = _mm256_set1_epi64x(a);
  auto bit = _mm256_setzero_si256();
  b = _mm256_srli_epi64(b, 1);

  int64_t c_s = __builtin_ctzll( a );

  auto c = _mm256_set1_epi64x(c_s);
  auto b_shifted = _mm256_srli_epi64(b, 1);
  auto b_xor = _mm256_xor_si256(b, b_shifted);
  auto c_and = _mm256_and_si256(c, b_xor);
  bit = _mm256_xor_si256(bit, c_and);
  a_v = _mm256_srlv_epi64(a_v, c);
  a_v = _mm256_srli_epi64(a_v, 1);
  __m256i t;
  bool final1 = false;
  bool final2 = false;
  bool final3 = false;
  bool final4 = false;
  int bit1, bit2, bit3, bit4;
  int final_bit1, final_bit2, final_bit3, final_bit4;
  printf("%llu\n", a);
  do {
    t = _mm256_sub_epi64(a_v, b);
    auto bgta = _mm256_srli_epi64(t, 63);
    /*
    if( t == 0 ) {
      return 0;
    }
    */
    bit = _mm256_xor_si256(bit, _mm256_and_si256(bgta, _mm256_and_si256(a_v, b)));
    b = _mm256_add_epi64(b, _mm256_and_si256(bgta, t));
    a_v = _mm256_sub_epi64(_mm256_xor_si256(t, bgta), bgta);

    int64_t t_arr[4] __attribute__ ((aligned));
    _mm256_store_si256((__m256i*)t_arr, t);

    int64_t c1 = t_arr[0] == 0 ? 0 : __builtin_ctzll( t_arr[0] );
    int64_t c2 = t_arr[1] == 0 ? 0 : __builtin_ctzll( t_arr[1] );
    int64_t c3 = t_arr[2] == 0 ? 0 : __builtin_ctzll( t_arr[2] );
    int64_t c4 = t_arr[3] == 0 ? 0 : __builtin_ctzll( t_arr[3] );
    c = _mm256_set_epi64x(c4, c3, c2, c1);
    c = _mm256_add_epi64(c, _mm256_set1_epi64x(1ULL));
    b_shifted = _mm256_srli_epi64(b, 1);
    b_xor = _mm256_xor_si256(b, b_shifted);
    c_and = _mm256_and_si256(c, b_xor);
    bit = _mm256_xor_si256(bit, c_and);
    a_v = _mm256_srlv_epi64(a_v, c);

    int64_t b1, b2, b3, b4;
    b1 = _mm256_extract_epi64(b, 0);
    b2 = _mm256_extract_epi64(b, 1);
    b3 = _mm256_extract_epi64(b, 2);
    b4 = _mm256_extract_epi64(b, 3);
    bit1 = _mm256_extract_epi64(bit, 0);
    bit2 = _mm256_extract_epi64(bit, 1);
    bit3 = _mm256_extract_epi64(bit, 2);
    bit4 = _mm256_extract_epi64(bit, 3);
    if(b1 == 0 && !final1){
      final_bit1 = bit1;
      final1 = true;
    }
    if(b2 == 0 && !final2){
      final_bit2 = bit2;
      final2 = true;
    }
    if(b3 == 0 && !final3){
      final_bit3 = bit3;
      final3 = true;
    }
    if(b4 == 0 && !final4){
      final_bit4 = bit4;
      final4 = true;
    }
    printf("%d, %d, %d, %d\n", bit1, bit2, bit3, bit4);

  } while(!final1 || !final2 || !final3 || !final4);
  bit = _mm256_set_epi64x(bit4, bit3, bit2, bit1);

  return _mm256_sub_epi64(_mm256_set1_epi64x(1ULL), _mm256_slli_epi64(_mm256_and_si256(bit, _mm256_set1_epi64x(1ULL)), 1));
}

/*

int main() {
    auto a = _mm256_set1_epi64x(25ULL);
    auto b = _mm256_set_epi64x(29, 31, 37, 41);
    auto result = jacobi_vector(a, b);
    auto t1 = _mm256_extract_epi64(result, 0);
    auto t2 = _mm256_extract_epi64(result, 1);
    auto t3 = _mm256_extract_epi64(result, 2);
    auto t4 = _mm256_extract_epi64(result, 3);
    printf("%d, %d, %d, %d\n", t1, t2, t3, t4);
}
*/
