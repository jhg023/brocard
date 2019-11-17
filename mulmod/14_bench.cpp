#include <celero/Celero.h>
#include <cstdint>
#include <flint/flint.h>
#include <flint/ulong_extras.h>

using celero::DoNotOptimizeAway;

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

#define NUM_RUNS 15
#define NUM_ITERATIONS 10

class MulmodFixture : public celero::TestFixture {
public:
  MulmodFixture() = default;

  virtual std::vector<celero::TestFixture::ExperimentValue> getExperimentValues() const override {
    std::vector<celero::TestFixture::ExperimentValue> problemSpace;

    const int totalNumberOfTests = 50;

    for( int i = 2; i < totalNumberOfTests; i++ ) {
      problemSpace.push_back( {int64_t( pow( 2, i ) )} );
    }

    return problemSpace;
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

  virtual void setUp( const celero::TestFixture::ExperimentValue &experimentValue ) override {
    aStart = experimentValue.Value;
    aEnd = aStart + 1'000'000;
    mStart = generate_primes( aEnd );
    mEnd = mStart;
    c = ( 1ULL << ( log2( aEnd ) ) ) - mStart;
    t = log2( aEnd );
    //printf("aStart = %llu, aEnd = %llu, mStart = %llu, c = %llu, t = %llu\n", aStart, aEnd, mStart, c, t);
  }

  uint64_t aStart;
  uint64_t aEnd;
  uint64_t mStart;
  uint64_t mEnd;
  uint64_t c;
  uint64_t t;
};

BASELINE_F( MulmodBench, mulmod_14, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = mulmod_14( result, i, m, this->c, this->t ) );
    }
  }
}

CELERO_MAIN
