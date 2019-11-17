#include <cstdint>

uint64_t mulmod_mcqueen( uint64_t a, uint64_t b, uint64_t m ) {
  uint64_t res = 0;
  uint64_t temp_b;

  while( a != 0 ) {
    if( a & 1 ) {
      /* Add b to res, modulo m, without overflow */
      if( b >= m - res ) /* Equiv to if (res + b >= m), without overflow */
        res -= m;
      res += b;
    }
    a >>= 1;

    /* Double b, modulo m */
    temp_b = b;
    if( b >= m - b ) /* Equiv to if (2 * b >= m), without overflow */
      temp_b -= m;
    b += temp_b;
  }
  return res;
}

