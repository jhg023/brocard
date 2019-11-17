#include <cstdint>

uint64_t mulmod_naive( uint64_t i, uint64_t j, uint64_t k ) {
  uint64_t r = 0;
  while( j > 0 ) {
    if( j & 1 )
      r = ( r + i ) % k;
    i = ( i + i ) % k;
    j >>= 1;
  }
  return r;
}
