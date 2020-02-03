static inline int jacobi_modified3( int64_t a, int64_t b, int64_t bit ) {
  if ( __builtin_expect( a == b, 0 ) ) {
    return 0;
  }

  b >>= 1;

  unsigned int c = __builtin_ctzll( a );
  bit &= c;

  a >>= c;
  a >>= 1;

  do {
    #pragma unroll( 1 )
    for ( unsigned int i = 0; i < 1; ++i ) {
      int64_t t = a - b;

      /* If b > a, invoke reciprocity */
      bit ^= ( a & ( a >= b ? 0 : b ) );

      /* b <-- min (a, b) */
      b = a < b ? a : b;

      /* a <-- |a - b| */
      a = ( ( t < 0 ) ? -t : t );

      c = __builtin_ctzll( a ) + 1;

      bit ^= c & ( b ^ ( b >> 1 ) );

      a >>= c;
    }
  } while( b > 0 );

  return bit & 1;
}
