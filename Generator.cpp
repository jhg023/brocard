#include <flint/ulong_extras.h>

// The number of threads to use when computing the initial factorial values.
// This is a memory hog, so a small amount of threads should be used here to
// avoid running out of memory.
constexpr int FACTORIAL_NUM_THREADS = 8;

// Milestone used for printing progress (normally 1 billion).
constexpr uint64_t MILESTONE = 100'000'000;

// If 'n - last_n[i] <= MULMOD_DIFFERENCE', then a more efficient method will be used
// to catch up 'last_n[i]' instead of repeatedly calling 'mulmod_preinv'.
constexpr int MULMOD_DIFFERENCE = 2'000'000;

// The number of primes to use when testing.
constexpr int NUM_PRIMES = 40;

// The amount of sub-ranges that the range should be partitioned into (normally 32).
constexpr int NUM_SUB_RANGES = 32;

// The start of the range (inclusive) to search for solutions in.
constexpr uint64_t RANGE_START = 1;

// The end of the range (inclusive) to search for solutions in.
constexpr uint64_t RANGE_END = 1'000'000'000;

// The name of the file to write potential solutions to.
constexpr char SOLUTION_FILE_NAME[] = "brocard_solutions.txt";

/**
 * Generates the inverses of each prime in 'primes'.
 */
void generate_pinvs( uint64_t *pinvs, const uint64_t *primes ) {
  for ( uint64_t i = 0; i < NUM_PRIMES; ++i ) {
    pinvs[i] = n_preinvert_limb( primes[i] );
  }
}

/**
 * Generates the first NUM_PRIMES primes after 'start + 1'.
 */
void generate_primes( uint64_t *primes, uint64_t start ) {
  n_primes_t iter;
  n_primes_init( iter );
  n_primes_jump_after( iter, start + 1 );

  for ( uint64_t i = 0; i < NUM_PRIMES; ++i ) {
      primes[i] = n_primes_next( iter );
  }

  n_primes_clear( iter );
}

void print_headers() {
    printf( "#include <flint/nmod_poly.h>\n" );
	printf( "#include <flint/ulong_extras.h>\n\n" );
}

void print_ll_mod_preinv() {
	printf( "static inline uint64_t ll_mod_preinv( uint64_t a_hi, uint64_t a_lo, uint64_t n, uint64_t ninv ) {\n" );
	printf( "    uint64_t q0, q1, r;\n\n" );
	printf( "    const int norm = __builtin_clzll( n );\n\n" );
	printf( "    n <<= norm;\n" );
	printf( "    a_hi <<= norm;\n\n" );
	printf( "    const uint64_t u1 = a_hi + ( a_lo >> ( %d - norm ) );\n", FLINT_BITS );
	printf( "    const uint64_t u0 = ( a_lo << norm );\n\n" );
	printf( "    umul_ppmm( q1, q0, ninv, u1 );\n" );
	printf( "    add_ssaaaa( q1, q0, q1, q0, u1, u0 );\n\n" );
	printf( "    r = ( u0 - ( q1 + 1 ) * n );\n\n" );
	printf( "    if ( r > q0 ) {\n" );
	printf( "        r += n;\n" );
	printf( "    }\n\n" );
	printf( "    if ( __builtin_expect( r < n, 1 ) ) {\n" );
	printf( "        return r >> norm;\n" );
	printf( "    }\n\n" );
	printf( "    return ( r - n ) >> norm;\n" );
	printf( "}\n\n" );
}

void print_mulmod_preinv() {
	printf( "static inline uint64_t mulmod_preinv( uint64_t a, uint64_t b, uint64_t n, uint64_t ninv ) {\n" );
	printf( "    uint64_t p1, p2;\n" );
	printf( "    umul_ppmm( p1, p2, a, b );\n" );
	printf( "    return ll_mod_preinv( p1, p2, n, ninv );\n" );
	printf( "}\n\n" );
}

void print_brocard( const int sub_range, const uint64_t sub_range_start, const uint64_t sub_range_end, 
	                const uint64_t *primes, const uint64_t *pinvs ) {
    printf( "void brocard_%d( uint64_t *factorials ) {\n", sub_range );
    printf( "    uint64_t last_n[%d] = { 0 };\n\n", NUM_PRIMES );
    printf( "    for ( uint64_t &i : last_n ) {\n" );
    printf( "        i = %lluULL;\n", sub_range_start - 1 );
    printf( "    }\n\n" );
    printf( "    int best_i = 25, i;\n" );
    printf( "    for ( uint64_t n = %lluULL; n <= %lluULL; ++n ) {\n", sub_range_start, sub_range_end );

    for ( int index = 0; index < NUM_PRIMES; ++index ) {
    	printf( "        if ( last_n[%d] == n - 1 ) {\n", index );
        printf( "            factorials[%d] = mulmod_preinv( factorials[%d], n, %lluULL, %lluULL );\n", index,
            index, primes[index], pinvs[index] );
        printf( "        } else if ( n - last_n[%d] <= %d ) {\n", index, MULMOD_DIFFERENCE );
        printf( "            for ( uint64_t j = last_n[%d] + 1; j <= n; ++j ) {\n", index );
        printf( "                factorials[%d] = mulmod_preinv( factorials[%d], j, %lluULL, %lluULL );\n",
            index, index, primes[index], pinvs[index] );
        printf( "            }\n" );
        printf( "        } else {\n" );
        printf( "            factorials[%d] = initialize_factorial( n, %lluULL, %lluULL );\n", index,
            primes[index], pinvs[index] );
        printf( "        }\n\n" );
        printf( "        last_n[%d] = n;\n\n", index );
        printf( "        if ( jacobi_unsigned( factorials[%d] + 1, %lluULL ) == -1 ) {\n", index, 
        	primes[index] );
        printf( "            if ( __builtin_expect( n %% %lluULL == 0, 0 ) != 0 ) {\n", MILESTONE );
        printf( "                printf( \"[Sub Range #%d] Progress: %%llu (%.2f%%%%)\\n\", n );\n",
            sub_range + 1, 100.0 * sub_range / NUM_SUB_RANGES );
        printf( "            } else if ( __builtin_expect( %d >= best_i, 0 ) != 0 ) {\n", index );
        printf( "                best_i = %d;\n", index );
        printf( "                printf( \"[Sub Range #%d] Progress: %%llu (%.2f%%%%), Tests Passed: %d\\n\", n );\n",
            sub_range + 1, 100.0 * sub_range / NUM_SUB_RANGES, index );
        printf( "            }\n\n" );
        printf( "            continue;\n" );
        printf( "        }\n\n" );
    }

    printf( "        printf( \"[Sub Range #%d] Potential Solution: %%llu, primes[0] = %llu, factorials[0] = %%llu\\n\", n, factorials[0] );\n",
        sub_range + 1, primes[0] );
    printf( "        FILE *fp = fopen( \"%s\", \"ae\" );\n", SOLUTION_FILE_NAME );
    printf( "        fprintf( fp, \"%%llu\\n\", n );\n" );
    printf( "        fclose( fp );\n" );
    printf( "    }\n" );
    printf(" }\n\n" );
}

void print_factorial_fast_mod2_preinv() {
	printf( "static inline uint64_t factorial_fast_mod2_preinv( uint64_t n, uint64_t p, uint64_t pinv ) {\n" );
	printf( "    slong i, m;\n" );
	printf( "    nmod_t mod;\n" );
	printf( "    mp_ptr t, u, v;\n" );
	printf( "    mp_limb_t r, s;\n\n" );
	printf( "    nmod_init( &mod, p );\n\n" );
	printf( "    m = n_sqrt( n );\n\n" );
	printf( "    t = _nmod_vec_init( m + 1 );\n" );
	printf( "    u = _nmod_vec_init( m + 1 );\n" );
	printf( "    v = _nmod_vec_init( m + 1 );\n\n" );
	printf( "    t[0] = UWORD( 0 );\n\n" );
	printf( "    for ( i = 1; i < m; i++ ) {\n" );
	printf( "        t[i] = n_submod( t[i - 1], UWORD( 1 ), p );\n" );
	printf( "    }\n\n" );
	printf( "    _nmod_poly_product_roots_nmod_vec( u, t, m, mod );\n\n" );
	printf( "    for ( i = 0; i < m; i++ ) {\n" );
	printf( "        t[i] = n_mod2_preinv( i * m + 1, p, pinv );\n" );
	printf( "    }\n\n" );
	printf( "    _nmod_poly_evaluate_nmod_vec_fast( v, u, m + 1, t, m, mod );\n\n" );
	printf( "    r = 1;\n\n" );
	printf( "    for ( i = 0; i < m; i++ ) {\n" );
	printf( "        r = mulmod_preinv( r, v[i], mod.n, mod.ninv );\n" );
	printf( "    }\n\n" );
	printf( "    for ( s = m * m + 1; s <= n; s++ ) {\n" );
	printf( "        r = mulmod_preinv( r, s, mod.n, mod.ninv );\n" );
	printf( "    }\n\n" );
	printf( "    _nmod_vec_clear( t );\n" );
	printf( "    _nmod_vec_clear( u );\n" );
	printf( "    _nmod_vec_clear( v );\n\n" );
	printf( "    return r;\n" );
	printf( "}\n\n" );
}

void print_initialize_factorial() {
	printf( "static inline uint64_t initialize_factorial( uint64_t n, uint64_t prime, uint64_t pinv ) {\n" );
	printf( "    if ( n < ( prime >> 1 ) ) {\n" );
	printf( "        return factorial_fast_mod2_preinv( n, prime, pinv );\n" );
	printf( "    }\n\n" );
	printf( "    uint64_t factorial = factorial_fast_mod2_preinv( prime - n - 1, prime, pinv );\n\n" );
	printf( "    factorial = n_invmod( factorial, prime );\n\n" );
	printf( "    if ( ( n & 1 ) == 0 ) {\n" );
	printf( "        factorial = -factorial + prime;\n" );
	printf( "    }\n\n" );
	printf( "    return factorial %% prime;\n" );
	printf( "}\n\n" );
}

void print_jacobi_unsigned() {
	printf( "static inline int jacobi_unsigned( uint64_t x, uint64_t y ) {\n" );
	printf( "    uint64_t b = x, a = y, temp;\n" );
	printf( "    int s;\n\n" );
	printf( "    int exp = __builtin_ctzll( b );\n" );
	printf( "    b >>= exp;\n\n" );
	printf( "    bool first = (( exp * ( a * a - 1 ) ) & 8) != 0;\n" );
	printf( "    bool second = (( ( a - 1 ) * ( b - 1 ) ) & 4) != 0;\n\n" );
	printf( "    if ( first != second ) {\n" );
	printf( "        s = -1;\n" );
	printf( "    } else {\n" );
	printf( "        s = 1;\n" );
	printf( "    }\n\n" );
	printf( "    while ( b != 1 ) {\n" );
	printf( "        if ( ( a >> 2 ) < b ) {\n" );
	printf( "            temp = a - b;\n" );
	printf( "            a = b;\n\n" );
	printf( "            if ( temp < b ) {\n" );
	printf( "                b = temp;\n" );
	printf( "            } else if ( temp < ( b << 1 ) ) {\n" );
	printf( "                b = temp - a;\n" );
	printf( "            } else {\n" );
	printf( "                b = temp - ( a << 1 );\n" );
	printf( "            }\n" );
	printf( "        } else {\n" );
	printf( "            temp = a %% b;\n" );
	printf( "            a = b;\n" );
	printf( "            b = temp;\n" );
	printf( "        }\n\n" );
	printf( "        if ( b == 0 ) {\n" );
	printf( "            return 0;\n" );
	printf( "        }\n\n" );
	printf( "        exp = __builtin_ctzll( b );\n" );
	printf( "        b >>= exp;\n\n" );
	printf( "        first = (( exp * ( a * a - 1 ) ) & 8) != 0;\n" );
	printf( "        second =  (( ( a - 1 ) * ( b - 1 ) ) & 4) != 0;\n\n" );
	printf( "        if ( first != second ) {\n" );
	printf( "            s = -s;\n" );
	printf( "        }\n" );
	printf( "    }\n\n" );
	printf( "    return s;\n" );
	printf( "}\n\n" );
}

void print_primes( uint64_t primes[][NUM_PRIMES] ) {
    printf( "    uint64_t primes[%d][%d] = {\n", NUM_SUB_RANGES, NUM_PRIMES );

	for ( int sub_range = 0; sub_range < NUM_SUB_RANGES; ++sub_range ) {
		printf( "        { " );

		for ( int index = 0; index < NUM_PRIMES; ++index ) {
			printf( "%lluULL", primes[sub_range][index] );

			if ( index < NUM_PRIMES - 1 ) {
				printf( ", " );
			}
		}

		printf( " }" );

		if (sub_range == NUM_SUB_RANGES - 1) {
			printf( "\n" );
		} else {
			printf( ",\n" );
		}
	}

	printf( "    };\n\n" );
}

void print_pinvs( const uint64_t pinvs[][NUM_PRIMES] ) {
    printf( "    uint64_t pinvs[%d][%d] = {\n", NUM_SUB_RANGES, NUM_PRIMES );

	for ( int sub_range = 0; sub_range < NUM_SUB_RANGES; ++sub_range ) {
		printf( "        { " );

		for ( int index = 0; index < NUM_PRIMES; ++index ) {
			printf( "%lluULL", pinvs[sub_range][index] );

			if ( index < NUM_PRIMES - 1 ) {
				printf( ", " );
			}
		}

		printf( " }" );

		if (sub_range == NUM_SUB_RANGES - 1) {
			printf( "\n" );
		} else {
			printf( ",\n" );
		}
	}

	printf( "    };\n\n" );
}

void print_main( const uint64_t partition_size, const uint64_t *sub_range_start, uint64_t primes[][NUM_PRIMES], 
	             const uint64_t pinvs[][NUM_PRIMES] ) {
	printf( "auto main() -> int {\n" );
	printf( "    remove( \"%s\" );\n\n", SOLUTION_FILE_NAME );
	printf( "    printf( \"Partition Size: %llu\\n\" );\n\n", partition_size );
	printf( "    uint64_t factorials[%d][%d];\n\n", NUM_SUB_RANGES, NUM_PRIMES );
	
	print_primes( primes );
	print_pinvs( pinvs );

	for ( int sub_range = 0; sub_range < NUM_SUB_RANGES; ++sub_range ) {
	    printf( "#pragma omp parallel for num_threads( %d ) schedule( dynamic )\n", FACTORIAL_NUM_THREADS );
	    printf( "    for ( int index = 0; index < %d; ++index ) {\n", NUM_PRIMES );
	    printf( "        factorials[%d][index] = initialize_factorial( %lluULL, primes[%d][index], pinvs[%d][index] );\n",
	        sub_range, sub_range_start[sub_range] - 1, sub_range, sub_range );
	    printf( "    }\n\n" );
	}

	printf( "#pragma omp parallel for schedule( dynamic )\n" );
	printf( "    for ( int sub_range = 0; sub_range < %d; ++sub_range ) {\n", NUM_SUB_RANGES );
	printf( "        switch ( sub_range ) {\n" );

	for ( int sub_range = 0; sub_range < NUM_SUB_RANGES; ++sub_range ) {
		printf( "            case %d:\n", sub_range );
		printf( "                brocard_%d( factorials[%d] );\n", sub_range, sub_range );
		printf( "                break;\n" );
	}

    printf( "        }\n" );
	printf( "    }\n" );
	printf( "}\n" );
}

int main() {
    if ( RANGE_START > RANGE_END ) {
        printf( "RANGE_START cannot be greater than RANGE_END!" );
        return 1;
    }

    const uint64_t partition_size = ( RANGE_END - RANGE_START ) / NUM_SUB_RANGES;

    print_headers();
    print_ll_mod_preinv();
    print_mulmod_preinv();
    print_factorial_fast_mod2_preinv();
    print_initialize_factorial();
    print_jacobi_unsigned();

    uint64_t sub_range_start[NUM_SUB_RANGES];
    uint64_t sub_range_end[NUM_SUB_RANGES];
    uint64_t primes[NUM_SUB_RANGES][NUM_PRIMES];
    uint64_t pinvs[NUM_SUB_RANGES][NUM_PRIMES];

    for ( int sub_range = 0; sub_range < NUM_SUB_RANGES; ++sub_range ) {
    	sub_range_start[sub_range] = RANGE_START + (partition_size * sub_range);
    	sub_range_end[sub_range] = sub_range == NUM_SUB_RANGES - 1 ? RANGE_END : sub_range_start[sub_range] +
    	    partition_size - 1;
        generate_primes( primes[sub_range], sub_range_end[sub_range] );
        generate_pinvs( pinvs[sub_range], primes[sub_range] );
        
        print_brocard( sub_range, sub_range_start[sub_range], sub_range_end[sub_range], primes[sub_range], 
        	pinvs[sub_range] );
    }

    print_main( partition_size, sub_range_start, primes, pinvs );
}










