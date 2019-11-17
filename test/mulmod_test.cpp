#include "gtest/gtest.h"
#include <flint/flint.h>
#include <flint/ulong_extras.h>

#define NUM_TESTS 53

// expected, a, b, m, c, t
uint64_t test_data[53][6] = {
    {1554, 1024, 1025, 2039, 9, 11},
    {1027, 2048, 2049, 4093, 3, 12},
    {6144, 4096, 4097, 8191, 1, 13},
    {4099, 8192, 8193, 16381, 3, 14},
    {8287, 16384, 16385, 32749, 19, 15},
    {16444, 32768, 32769, 65521, 15, 16},
    {98304, 65536, 65537, 131071, 1, 17},
    {196613, 131072, 131073, 262139, 5, 18},
    {393216, 262144, 262145, 524287, 1, 19},
    {262147, 524288, 524289, 1048573, 3, 20},
    {1572882, 1048576, 1048577, 2097143, 9, 21},
    {1048579, 2097152, 2097153, 4194301, 3, 22},
    {2097212, 4194304, 4194305, 8388593, 15, 23},
    {4194307, 8388608, 8388609, 16777213, 3, 24},
    {8388998, 16777216, 16777217, 33554393, 39, 25},
    {50331653, 33554432, 33554433, 67108859, 5, 26},
    {33554822, 67108864, 67108865, 134217689, 39, 27},
    {201327390, 134217728, 134217729, 268435399, 57, 28},
    {134217731, 268435456, 268435457, 536870909, 3, 29},
    {268435771, 536870912, 536870913, 1073741789, 35, 30},
    {1610612736, 1073741824, 1073741825, 2147483647, 1, 31},
    {3221225477, 2147483648, 2147483649, 4294967291, 5, 32},
    {6442450962, 4294967296, 4294967297, 8589934583, 9, 33},
    {12884902298, 8589934592, 8589934593, 17179869143, 41, 34},
    {8589934840, 17179869184, 17179869185, 34359738337, 31, 35},
    {51539607557, 34359738368, 34359738369, 68719476731, 5, 36},
    {103079215254, 68719476736, 68719476737, 137438953447, 25, 37},
    {206158430703, 137438953472, 137438953473, 274877906899, 45, 38},
    {137438953486, 274877906944, 274877906945, 549755813881, 7, 39},
    {274877908858, 549755813888, 549755813889, 1099511627689, 87, 40},
    {1649267441769, 1099511627776, 1099511627777, 2199023255531, 21, 41},
    {1099511627809, 2199023255552, 2199023255553, 4398046511093, 11, 42},
    {6597069767454, 4398046511104, 4398046511105, 8796093022151, 57, 43},
    {13194139533380, 8796093022208, 8796093022209, 17592186044399, 17, 44},
    {8796093022978, 17592186044416, 17592186044417, 35184372088777, 55, 45},
    {52776558133353, 35184372088832, 35184372088833, 70368744177643, 21, 46},
    {35184372092167, 70368744177664, 70368744177665, 140737488355213, 115, 47},
    {70368744178549, 140737488355328, 140737488355329, 281474976710597, 59, 48},
    {422212465067604, 281474976710656, 281474976710657, 562949953421231, 81, 49},
    {281474976710845, 562949953421312, 562949953421313, 1125899906842597, 27, 50},
    {1688849860268064, 1125899906842624, 1125899906842625, 2251799813685119, 129, 51},
    {1125899906843188, 2251799813685248, 2251799813685249, 4503599627370449, 47, 52},
    {2251799813688356, 4503599627370496, 4503599627370497, 9007199254740881, 111, 53},
    {13510798882111752, 9007199254740992, 9007199254740993, 18014398509481951, 33, 54},
    {9007199254741762, 18014398509481984, 18014398509481985, 36028797018963913, 55, 55},
    {54043195528445957, 36028797018963968, 36028797018963969, 72057594037927931, 5, 56},
    {108086391056891943, 72057594037927936, 72057594037927937, 144115188075855859, 13, 57},
    {72057594037928125, 144115188075855872, 144115188075855873, 288230376151711717, 27, 58},
    {144115188075856642, 288230376151711744, 288230376151711745, 576460752303423433, 55, 59},
    {864691128455137371, 576460752303423488, 576460752303423489, 1152921504606846883, 93, 60},
    {1729382256910270464, 1152921504606846976, 1152921504606846977, 2305843009213693951, 1, 61},
    {3458764513820541726, 2305843009213693952, 2305843009213693953, 4611686018427387847, 57, 62},
    {6917529027641082006, 4611686018427387904, 4611686018427387905, 9223372036854775783, 25, 63}};

// Computes 'a^b / n'
constexpr auto pow( uint64_t a, int b, uint64_t n ) -> uint64_t {
  __int128 res = 1;
  for( int i = 0; i < b; ++i ) {
    res *= a;
  }
  return res / n;
}

constexpr auto n_clog( uint64_t n ) -> uint64_t {
  uint64_t r = 0;
  while( n >>= 1 != 0u ) {
    ++r;
  }
  return r;
}

auto mulmod_barrett( uint64_t a, uint64_t b, uint64_t n, int k, uint64_t r ) -> uint64_t {
  __int128 a_128 = a;
  __int128 b_128 = b;
  __int128 n_128 = n;
  __int128 x = a_128 * b_128;
  __int128 t = x - ( ( ( x * r ) >> k ) * n_128 );
  return ( t < n ? t : t - n );
}

uint64_t mulmod_asm( uint64_t a, uint64_t b, uint64_t n ) {
  uint64_t d;
  uint64_t unused; // dummy output, unused, to tell GCC that RAX register is modified by this snippet
  asm( "mulq %3\n\t"
       "divq %4"
       : "=a"( unused ), "=&d"( d )
       : "a"( a ), "rm"( b ), "rm"( n )
       : "cc" );
  return d;
}

uint64_t mulmod_14( uint64_t a, uint64_t b, uint64_t m, uint64_t c, uint64_t t ) {
  __int128 a_128 = a;
  __int128 x = a_128 * b;

  __int128 orig_q = x >> t;
  uint64_t r = x - ( orig_q << t );
  uint64_t newq = ( orig_q * c ) >> t;
  r += ( orig_q * c ) - ( newq << t );
  uint64_t q = newq;
  newq = ( q * c ) >> t;
  r += ( q * c ) - ( newq << t );
  q = newq;
  while( r >= m ) {
    r -= m;
  }
  return r;
}

/*
#define FLINT_BITS 64
#define umul_ppmm( w1, w0, u, v )                                                                                      \
  __asm__( "mulq %3" : "=a"( w0 ), "=d"( w1 ) : "%0"( ( uint64_t )( u ) ), "rm"( ( uint64_t )( v ) ) )

#define count_leading_zeros( count, x )                                                                                \
  do {                                                                                                                 \
    uint64_t __cbtmp;                                                                                                  \
    __asm__( "bsrq %1,%0" : "=r"( __cbtmp ) : "rm"( ( uint64_t )( x ) ) );                                             \
    ( count ) = __cbtmp ^ ( uint64_t ) 63;                                                                             \
  } while( 0 )

#define add_ssaaaa( sh, sl, ah, al, bh, bl )                                                                           \
  __asm__( "addq %5,%q1\n\tadcq %3,%q0"                                                                                \
           : "=r"( sh ), "=&r"( sl )                                                                                   \
           : "0"( ( uint64_t )( ah ) ), "rme"( ( uint64_t )( bh ) ), "%1"( ( uint64_t )( al ) ),                       \
             "rme"( ( uint64_t )( bl ) ) )

#define r_shift( in, shift ) ( ( shift == FLINT_BITS ) ? 0LL : ( ( in ) >> ( shift ) ) )

double n_precompute_inverse( uint64_t n ) {
  return ( double ) 1 / ( double ) n;
}

uint64_t n_ll_mod_preinv( uint64_t a_hi, uint64_t a_lo, uint64_t n, uint64_t ninv ) {
  uint64_t q0, q1, r, norm;

  count_leading_zeros( norm, n );

  if( a_hi >= n ) {
    const uint64_t u1 = r_shift( a_hi, FLINT_BITS - norm );
    const uint64_t u0 = ( a_hi << norm );

    n <<= norm;

    umul_ppmm( q1, q0, ninv, u1 );
    add_ssaaaa( q1, q0, q1, q0, u1, u0 );

    a_hi = ( u0 - ( q1 + 1 ) * n );

    if( a_hi > q0 )
      a_hi += n;

    if( a_hi >= n )
      a_hi -= n;
  } else {
    n <<= norm;
    a_hi <<= norm;
  }

  {
    const uint64_t u1 = a_hi + r_shift( a_lo, FLINT_BITS - norm );
    const uint64_t u0 = ( a_lo << norm );

    umul_ppmm( q1, q0, ninv, u1 );
    add_ssaaaa( q1, q0, q1, q0, u1, u0 );

    r = ( u0 - ( q1 + 1 ) * n );

    if( r > q0 )
      r += n;

    return ( r < n ) ? ( r >> norm ) : ( ( r - n ) >> norm );
  }
}

uint64_t n_mulmod2_preinv( uint64_t a, uint64_t b, uint64_t n, uint64_t ninv ) {
  uint64_t p1, p2;

  umul_ppmm( p1, p2, a, b );
  return n_ll_mod_preinv( p1, p2, n, ninv );
}
*/
/*
uint64_t mulmod_14( uint64_t a, uint64_t b, uint64_t m, uint64_t c, uint64_t t ) {
  __int128 a_128 = a;
  __int128 x = a_128 * b;

  __int128 q = x >> t;
  uint64_t r = x - ( q << t );
  while( q > 0 ) {
    uint64_t newq = ( q * c ) >> t;
    r += ( q * c ) - ( newq << t );
    q = newq;
  }
  while( r >= m ) {
    r -= m;
  }
  return r;
}
*/

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

TEST( MulmodTest, mulmod_14 ) {
  for( uint64_t i = 0; i < NUM_TESTS; ++i ) {
    ASSERT_EQ( test_data[i][0],
               mulmod_14( test_data[i][1], test_data[i][2], test_data[i][3], test_data[i][4], test_data[i][5] ) );
  }
}

TEST( MulmodTest, mulmod_asm ) {
  for( uint64_t i = 0; i < NUM_TESTS; ++i ) {
    ASSERT_EQ( test_data[i][0], mulmod_asm( test_data[i][1], test_data[i][2], test_data[i][3] ) );
  }
}

TEST( MulmodTest, mulmod_barrett ) {
  for( uint64_t i = 0; i < NUM_TESTS; ++i ) {
    int k = n_clog( test_data[i][3] ) << 1;
    uint64_t r = pow( 4, k >> 1, test_data[i][3] );
    ASSERT_EQ( test_data[i][0], mulmod_barrett( test_data[i][1], test_data[i][2], test_data[i][3], k, r ) );
  }
}

TEST( MulmodTest, mulmod_precomp ) {
  for( uint64_t i = 0; i < NUM_TESTS; ++i ) {
    double npre = n_precompute_inverse( test_data[i][3] );
    ASSERT_EQ( test_data[i][0], mulmod_precomp( test_data[i][1], test_data[i][2], test_data[i][3], npre ) );
  }
}

TEST( MulmodTest, n_mulmod2_preinv ) {
  for( uint64_t i = 0; i < NUM_TESTS; ++i ) {
    auto pinv = n_preinvert_limb( test_data[i][3] );
    ASSERT_EQ( test_data[i][0], n_mulmod2_preinv( test_data[i][1], test_data[i][2], test_data[i][3], pinv ) );
  }
}

TEST( MulmodTest, mulmod_14_comprehensive ) {
  constexpr uint64_t m = 9223372036854775783;
  constexpr uint64_t c = 25;
  constexpr uint64_t t = 63;
  uint64_t result_asm = 1;
  uint64_t result_14 = 1;
  for( uint64_t i = 1; i < 1'000'000'000; ++i ) {
    result_14 = mulmod_14( result_14, i, m, c, t );
    result_asm = mulmod_asm( result_asm, i, m );
    ASSERT_EQ( result_asm, result_14 );
  }
}

auto main() -> int {
  testing::InitGoogleTest();
  RUN_ALL_TESTS();
  return 0;
}
