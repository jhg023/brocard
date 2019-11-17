#include <cstdint>

auto mulmod_precomp( uint64_t a, uint64_t b, uint64_t n, double npre ) -> uint64_t {
  uint64_t quot;
  int64_t rem;

  quot = static_cast<uint64_t>( static_cast<double>( a ) * static_cast<double>( b ) * npre );
  rem = a * b - quot * n;
  if( rem < 0 ) {
    rem += n;
    if( rem < 0 ) {
      return rem + n;
    }
  } else if( rem >= n ) {
    return rem - n;
  }
  return rem;
}
