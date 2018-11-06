#include "test_functions.h"

double timestamp()
{
  struct timeval tv;

  gettimeofday( &tv, ( struct timezone * ) 0 );
  return ( tv.tv_sec + (tv.tv_usec / 1000000.0) );
}
