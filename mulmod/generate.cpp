#include <cstdint>
#include <flint/flint.h>
#include <flint/ulong_extras.h>

constexpr auto log2( uint64_t n ) -> uint64_t {
  if( n == 0 ) {
    return 0;
  }
  uint64_t r = 1;
  while( n >>= 1 != 0u ) {
    ++r;
  }
  return r;
}

auto generate_primes( uint64_t start ) -> uint64_t {
  uint64_t num_bits = log2( start );
  uint64_t result = 0;
  uint64_t prime = 0;
  n_primes_t iter;
  n_primes_init( iter );
  n_primes_jump_after( iter, ( 1ULL << num_bits ) - 1000 );

  while( prime < ( 1ULL << num_bits ) ) {
    result = prime;
    prime = n_primes_next( iter );
  }

  n_primes_clear( iter );

  return result;
}

int main() {
  for( uint64_t shift = 10; shift < 64; ++shift ) {
    uint64_t aStart = 1ULL << shift;
    uint64_t m = generate_primes( aStart );
    uint64_t c = ( 1ULL << ( log2( aStart ) ) ) - m;
    uint64_t t = log2( aStart );
    printf( "print((%llu * %llu) %% %llu, %llu, %llu, %llu, %llu, %llu)\n", aStart, aStart + 1, m, aStart, aStart + 1,
            m, c, t );
  }
  return 0;
}
