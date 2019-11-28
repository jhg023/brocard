#include <celero/Celero.h>
#include <cstdint>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <immintrin.h>
#include "../mulmod/mulmod_naive.cpp"
#include "../mulmod/mulmod_asm.cpp"
#include "../mulmod/mulmod_barrett.cpp"
#include "../mulmod/mulmod14.cpp"
#include "../mulmod/mulmod_precomp.cpp"
#include "../mulmod/mulmod_vectorized.cpp"

using celero::DoNotOptimizeAway;

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

    const int totalNumberOfTests = 40;

    for( int i = 20; i < totalNumberOfTests; i++ ) {
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
    mEnd = mStart + 3;
    k = n_clog( mStart ) << 1;
    r = pow( 4, k >> 1, mStart );
    pinv = n_preinvert_limb( mStart );
    npre = n_precompute_inverse( mStart );
    c = ( 1ULL << ( log2( aEnd ) ) ) - mStart;
    t = log2( aEnd );
  }

  uint64_t aStart;
  uint64_t aEnd;
  uint64_t mStart;
  uint64_t mEnd;
  int k;
  uint64_t r;
  uint64_t pinv;
  double npre;
  uint64_t c;
  uint64_t t;
};

BASELINE_F( MulmodBench, mulmod_precomp, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = mulmod_precomp( result, i, m, this->npre ) );
    }
  }
}

BENCHMARK_F( MulmodBench, mulmod_14, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = mulmod_14( result, i, m, this->c, this->t ) );
    }
  }
}

BENCHMARK_F( MulmodBench, n_mulmod2_preinv, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = n_mulmod2_preinv( result, i, m, this->pinv ) );
    }
  }
}

BENCHMARK_F( MulmodBench, mulmod_barrett, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = mulmod_barrett( result, i, m, this->k, this->r ) );
    }
  }
}

BENCHMARK_F( MulmodBench, mulmod_asm, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  uint64_t result = this->aStart - 1;
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; ++m ) {
      DoNotOptimizeAway( result = mulmod_asm( result, i, m ) );
    }
  }
}

BENCHMARK_F( MulmodBench, mulmod_vectorized, MulmodFixture, NUM_RUNS, NUM_ITERATIONS ) {
  __m256i result = _mm256_set1_epi64x( this->aStart - 1 );
  for( uint64_t i = this->aStart; i <= this->aEnd; ++i ) {
    for( uint64_t m = this->mStart; m <= this->mEnd; m += 4 ) {
      __m256i m_vector = _mm256_set_epi64x( m + 3, m + 2, m + 1, m );
      DoNotOptimizeAway( result = mulmod_vectorized( result, i, m_vector ) );
    }
  }
}

CELERO_MAIN
