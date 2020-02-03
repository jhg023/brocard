#include "../jacobi/jacobi_basic.cpp"
#include "../jacobi/jacobi_unsigned_exponent.cpp"
#include "../jacobi/n_jacobi_unsigned_opt.cpp"
#include "../jacobi/jacobi_unsigned_cmov.cpp"
#include "../jacobi/jacobi_gmp.cpp"
#include "../jacobi/jacobi1.cpp"
#include "../jacobi/jacobi2.cpp"
#include "../jacobi/jacobi3.cpp"
#include "../jacobi/jacobi4.cpp"
#include "../jacobi/jacobi5.cpp"
#include "../jacobi/jacobi6.cpp"
#include "../jacobi/jacobi7.cpp"
#include <celero/Celero.h>
#include <cstdint>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <vector>
#include "../jacobi/jacobi_vector.cpp"


using celero::DoNotOptimizeAway;
using std::vector;

#define NUM_RUNS 50
#define NUM_ITERATIONS 50


class JacobiFixture : public celero::TestFixture {
public:
  JacobiFixture() = default;

  auto getExperimentValues() const -> std::vector<celero::TestFixture::ExperimentValue> override {
    std::vector<celero::TestFixture::ExperimentValue> problemSpace;

    const int totalNumberOfTests = 16;

    for( int i = 12; i < totalNumberOfTests; i++ ) {
      problemSpace.emplace_back( int64_t( pow( 10, i ) ) );
    }

    return problemSpace;
  }

  auto generate_primes( uint64_t start ) -> uint64_t {
    n_primes_t iter;
    n_primes_init( iter );
    n_primes_jump_after( iter, start );

    uint64_t prime = n_primes_next( iter );
    n_primes_clear( iter );

    return prime;
  }

  void setUp( const celero::TestFixture::ExperimentValue &experimentValue ) override {
    aStart = experimentValue.Value - 1000;
    aEnd = experimentValue.Value;
    auto p1 = generate_primes(aEnd);
    auto p2 = generate_primes(p1);
    auto p3 = generate_primes(p2);
    auto p4 = generate_primes(p3);
    primes = {p1, p2, p3, p4};
    ninvs = {n_preinvert_limb( primes[0] ), n_preinvert_limb( primes[1]), n_preinvert_limb(primes[2]), n_preinvert_limb(primes[3])};
  }

  uint64_t aStart{};
  uint64_t aEnd{};
  vector<uint64_t> primes;
  vector<uint64_t> ninvs;
};

BASELINE_F( JacobiBench, jacobi1, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified1( a, p, 0 ) );
    }
  }
}

BENCHMARK_F( JacobiBench, jacobi2, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified2( a, p, 0 ) );
    }
  }
}

BENCHMARK_F( JacobiBench, jacobi3, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified3( a, p, 0 ) );
    }
  }
}
BENCHMARK_F( JacobiBench, jacobi4, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified4( a, p, 0 ) );
    }
  }
}
BENCHMARK_F( JacobiBench, jacobi5, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified5( a, p, 0 ) );
    }
  }
}
BENCHMARK_F( JacobiBench, jacobi6, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified6( a, p, 0 ) );
    }
  }
}
BENCHMARK_F( JacobiBench, jacobi7, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_modified7( a, p, 0 ) );
    }
  }
}
/*
BASELINE_F( JacobiBench, n_jacobi_unsigned, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( n_jacobi_unsigned( a, p ) );
    }
  }
}

BENCHMARK_F( JacobiBench, jacobi_vector, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    DoNotOptimizeAway( jacobi_vector( a, _mm256_set_epi64x(primes[0], primes[1], primes[2], primes[3]) ) );
  }
}

BENCHMARK_F( JacobiBench, jacobi_gmp, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_gmp( a, p ) );
    }
  }
}

BENCHMARK_F( JacobiBench, opt_n_jacobi_unsigned, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( n_jacobi_unsigned_opt( a, p ) );
    }
  }
}

BENCHMARK_F( JacobiBench, jacobi_unsigned_exponent, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( int i = 0; i < 4; ++i ) {
      auto p = this->primes[i];
      auto ninv = this->ninvs[i];
      DoNotOptimizeAway( jacobi_unsigned_exponent( a, ( p - 1 ) / 2, p, ninv ) );
    }
  }
}


BENCHMARK_F( JacobiBench, jacobi_unsigned_cmov, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_unsigned_cmov( a, p ) );
    }
  }
}

BENCHMARK_F( JacobiBench, jacobi_basic, JacobiFixture, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t a = this->aStart; a <= this->aEnd; ++a ) {
    for( auto &&p: primes ) {
      DoNotOptimizeAway( jacobi_basic( a, p ) );
    }
  }
}
*/
CELERO_MAIN
