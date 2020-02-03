#include <cstdint>

#define FLINT_BITS 64

#define umul_ppmm( w1, w0, u, v )                                                                                      \
  __asm__( "mulq %3" : "=a"( w0 ), "=d"( w1 ) : "%0"( ( int64_t )( u ) ), "rm"( ( int64_t )( v ) ) )

#define add_ssaaaa( sh, sl, ah, al, bh, bl )                                                                           \
  do {                                                                                                                 \
    unsigned long int __x;                                                                                             \
    __x = ( al ) + ( bl );                                                                                             \
    ( sh ) = ( ah ) + ( bh ) + ( __x < ( al ) );                                                                       \
    ( sl ) = __x;                                                                                                      \
  } while( 0 )

#define r_shift( in, shift ) ( ( shift == FLINT_BITS ) ? 0LL : ( ( in ) >> ( shift ) ) )

int64_t mulmod_preinv( int64_t a, int64_t b, int64_t n, int64_t ninv ) {
  int64_t a_hi, a_lo;
  umul_ppmm( a_hi, a_lo, a, b );
  int64_t q0, q1, r;

  int norm = __builtin_clzll( n );
  n <<= norm;
  a_hi <<= norm;
  const int64_t u1 = a_hi + r_shift( a_lo, FLINT_BITS - norm );
  const int64_t u0 = ( a_lo << norm );

  umul_ppmm( q1, q0, ninv, u1 );
  add_ssaaaa( q1, q0, q1, q0, u1, u0 );

  r = ( u0 - ( q1 + 1 ) * n );

  if( r > q0 ) {
    r += n;
  }

  return ( r < n ) ? ( r >> norm ) : ( ( r - n ) >> norm );
}

int jacobi_unsigned_exponent( int64_t a, int64_t b, int64_t n, int64_t ninv ) {
  int64_t r = 1;
  while( b > 0 ) {
    if( ( b & 1 ) != 0 ) {
      r = mulmod_preinv( r, a, n, ninv );
    }

    b >>= 1;
    a = mulmod_preinv( a, a, n, ninv );
  }
  return r - n;
}
