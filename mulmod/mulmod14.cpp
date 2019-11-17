#include <cstdint>

uint64_t mulmod_14( uint64_t a, uint64_t b, uint64_t m, uint64_t c, uint64_t t ) {
  __int128 a_128 = a;
  __int128 x = a_128 * b;

  __int128 orig_q = x >> t;
  uint64_t r = x - ( orig_q << t );
  uint64_t newq = ( orig_q * c ) >> t;
  r += ( orig_q * c ) - ( newq << t );
  uint64_t q = newq;
  newq = ( q * c ) >> t;
  r += ( q * c ) - ( newq << t );
  q = newq;
  while( r >= m ) {
    r -= m;
  }
  return r;
}
