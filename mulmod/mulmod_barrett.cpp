#include <cstdint>

auto mulmod_barrett( uint64_t a, uint64_t b, uint64_t n, int k, uint64_t r ) -> uint64_t {
  __int128 a_128 = a;
  __int128 b_128 = b;
  __int128 x = a_128 * b_128;
  uint64_t t = x - ( ( ( x * r ) >> k ) * n );
  return ( t < n ? t : t - n );
}
