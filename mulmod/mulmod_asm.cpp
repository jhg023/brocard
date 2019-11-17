#include <cstdint>

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
