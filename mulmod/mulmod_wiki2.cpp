#include <cstdint>

uint64_t mulmod_wiki2( uint64_t a, uint64_t b, uint64_t m ) {
  long double x;
  uint64_t c;
  int64_t r;
  if( a >= m )
    a %= m;
  if( b >= m )
    b %= m;
  x = a;
  c = x * b / m;
  r = ( int64_t )( a * b - c * m ) % ( int64_t ) m;
  return r < 0 ? r + m : r;
}
