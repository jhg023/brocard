#include <cstdint>
#define CNST_LIMB( C ) ( ( int64_t ) C##LL )
#define MP_LIMB_T_MAX ( ~( int64_t ) 0 )
#define GMP_LIMB_HIGHBIT ( MP_LIMB_T_MAX ^ ( MP_LIMB_T_MAX >> 1 ) )
#define LIMB_HIGHBIT_TO_MASK( n )                                                                                                                    \
  ( ( ( int64_t ) -1 >> 1 ) < 0 ? ( int64_t )( n ) >> ( 63 ) : ( n ) &GMP_LIMB_HIGHBIT ? MP_LIMB_T_MAX : CNST_LIMB( 0 ) )

int jacobi_gmp( int64_t a, int64_t b ) {
  int c, bit = 0;

  bit >>= 1;
  b >>= 1;
  c = __builtin_ctzll( a );
  bit ^= c & ( b ^ ( b >> 1 ) );
  a >>= c;
  a >>= 1;

  do {
    int64_t t = a - b;
    int64_t bgta = LIMB_HIGHBIT_TO_MASK( t );

    if( t == 0 ) {
      return 0;
    }
    bit ^= ( bgta & a & b );
    b += ( bgta & t );
    a = ( t ^ bgta ) - bgta;
    c = __builtin_ctzll( t );
    c++;
    bit ^= c & ( b ^ ( b >> 1 ) );
    a >>= c;
  } while( b > 0 );

  return 1 - 2 * ( bit & 1 );
}

