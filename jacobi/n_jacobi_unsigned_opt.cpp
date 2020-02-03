#include <cstdint>
#include <cstdio>

auto n_jacobi_unsigned_opt( int64_t b, int64_t a ) -> int {
  int64_t temp;
  int s;
  int exp = __builtin_ctzll( b );
  b >>= exp;

  if( ( ( exp * ( a * a - 1 ) ) & 8 ) != 0U ) {
    s = -1;
  } else {
    s = 1;
  }

  if( ( ( ( a - 1 ) * ( b - 1 ) ) & 4 ) != 0U ) {
    s = -s;
  }

  while( b != 1ULL ) {
    if( ( a >> 2 ) < b ) {
      temp = a - b;
      a = b;

      if( temp < b ) {
        b = temp;
      } else if( temp < ( b << 1 ) ) {
        b = temp - a;
      } else {
        b = temp - ( a << 1 );
      }
    } else {
      switch( b ) {
        case 3:
          temp = a % 3;
          break;
        case 5:
          temp = a % 5;
          break;
        default:
          temp = a % b;
          break;
      }
      a = b;
      b = temp;
    }

    if( __builtin_expect( static_cast<long>( b == 0ULL ), 0 ) != 0 ) {
      return 0;
    }

    exp = __builtin_ctzll( b );

    b >>= exp;

    if( ( ( exp * ( a * a - 1 ) ) & 8 ) != 0U ) {
      s = -s;
    }

    if( ( ( ( a - 1 ) * ( b - 1 ) ) & 4 ) != 0U ) {
      s = -s;
    }
  }

  return s;
}
