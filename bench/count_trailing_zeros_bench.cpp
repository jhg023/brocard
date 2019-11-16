#include <celero/Celero.h>

using celero::DoNotOptimizeAway;

#define START 1
#define END 1'000'000

#define NUM_RUNS 0
#define NUM_ITERATIONS 1

int ctz_naive( uint64_t x ) {
  if( x == 0 ) {
    return 64;
  }
  int count = 0;
  while( ( x & 1 ) == 0 ) {
    x = x >> 1;
    count++;
  }
  return count;
}

int ctz_lookup( uint32_t x ) {
  static const int lookup[] = { 32, 0,  1,  26, 2,  23, 27, 0,  3, 16, 24, 30, 28, 11, 0,  13, 4,  7, 17,
                                0,  25, 22, 31, 15, 29, 10, 12, 6, 0,  21, 14, 9,  5,  20, 8,  19, 18 };

  return lookup[( -x & x ) % 37];
}

#define ctz_asm( count, x ) __asm__( "bsfq %1,%q0" : "=r"( count ) : "rm"( ( uint64_t )( x ) ) );

int ctz_linear( uint32_t v ) {
  int c;
  if( v ) {
    v = ( v ^ ( v - 1 ) ) >> 1;
    for( c = 0; v; c++ ) {
      v >>= 1;
    }
  } else {
    c = 32;
  }
  return c;
}

int ctz_parallel( uint32_t v ) {
  unsigned int c = 32; // c will be the number of zero bits on the right
  v &= -signed( v );
  if( v )
    c--;
  if( v & 0x0000FFFF )
    c -= 16;
  if( v & 0x00FF00FF )
    c -= 8;
  if( v & 0x0F0F0F0F )
    c -= 4;
  if( v & 0x33333333 )
    c -= 2;
  if( v & 0x55555555 )
    c -= 1;
  return c;
}

int ctz_binary( uint32_t v ) {
  unsigned int c;
  if( v & 0x1 ) {
    c = 0;
  } else {
    c = 1;
    if( ( v & 0xffff ) == 0 ) {
      v >>= 16;
      c += 16;
    }
    if( ( v & 0xff ) == 0 ) {
      v >>= 8;
      c += 8;
    }
    if( ( v & 0xf ) == 0 ) {
      v >>= 4;
      c += 4;
    }
    if( ( v & 0x3 ) == 0 ) {
      v >>= 2;
      c += 2;
    }
    c -= v & 0x1;
  }
  return c;
}

int ctz_float( uint32_t v ) {
  int r;
  float f = ( float ) ( v & -v );
  r = ( *( uint32_t * ) &f >> 23 ) - 0x7f;
  return r;
}

int ctz_debruijn( uint32_t v ) {
  int r;
  static const int MultiplyDeBruijnBitPosition[32] = { 0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
                                                       31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9 };
  r = MultiplyDeBruijnBitPosition[( ( uint32_t )( ( v & -v ) * 0x077CB531U ) ) >> 27];
  return r;
}

BASELINE( CTZBench, ctz_naive, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_naive( i ) );
  }
}

BENCHMARK( CTZBench, __builtin_ctz, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = __builtin_ctz( i ) );
  }
}

BENCHMARK( CTZBench, __builtin_ctzll, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = __builtin_ctzll( i ) );
  }
}

BENCHMARK( CTZBench, ctz_asm, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  uint64_t tmp = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    ctz_asm( tmp, i );
    DoNotOptimizeAway( result = tmp );
  }
}

BENCHMARK( CTZBench, ctz_debruijn, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_debruijn( i ) );
  }
}

BENCHMARK( CTZBench, ctz_float, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_float( i ) );
  }
}

BENCHMARK( CTZBench, ctz_binary, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_binary( i ) );
  }
}

BENCHMARK( CTZBench, ctz_lookup, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_lookup( i ) );
  }
}

BENCHMARK( CTZBench, ctz_linear, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_linear( i ) );
  }
}

BENCHMARK( CTZBench, ctz_parallel, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = 0;
  for( uint64_t i = START; i <= END; ++i ) {
    DoNotOptimizeAway( result = ctz_parallel( i ) );
  }
}

CELERO_MAIN
