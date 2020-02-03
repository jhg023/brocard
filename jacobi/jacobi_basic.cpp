#include <cstdint>

auto jacobi_basic( int64_t a, int64_t n ) -> int {
  int t = 1;
  int64_t tmp;
  int64_t r;
  while( a != 0 ) {
    while( a % 2 == 0 ) {
      a /= 2;
      r = n % 8;
      if( r == 3 || r == 5 ) {
        t = -t;
      }
    }
    tmp = a;
    a = n;
    n = tmp;
    if( a % 4 == n % 4 && n % 4 == 3 ) {
      t = -t;
    }
    a %= n;
  }
  if( n == 1 ) {
    return t;
  }
  return 0;
}
