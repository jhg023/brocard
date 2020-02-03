#include <cstdint>
#include <cstdio>

int jacobi_unsigned_cmov( int64_t b, int64_t a ) {
  int64_t temp;
  int s;

  int32_t exp = __builtin_ctzll( b );
  b >>= exp;

  if( ( ( exp * ( a * a - 1 ) ) & 8 ) != 0U ) {
    s = -1;
  } else {
    s = 1;
  }

  if( ( ( ( a - 1 ) * ( b - 1 ) ) & 4 ) != 0U ) {
    s = -s;
  }

  while( b != 1 ) {
    int64_t a1, b1, a2, b2, temp1;
    a2 = b;
    b2 = a % b;
    temp1 = a - b;
    a1 = b;
    if( temp1 < b ) {
      b1 = temp1;
    } else {
      b1 = ( temp1 < ( b << 1 ) ) ? temp1 - a : temp1 - ( a << 1 );
    }

    a = ( ( a >> 2 ) < b ) ? a2 : a1;
    b = ( ( a >> 2 ) < b ) ? b2 : b1;

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

