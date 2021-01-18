#ifndef RANDOM_H_
#define RANDOM_H_

/**
 * This header contains wrapper methods to generate random numbers.
 * The implementation uses a mersenne twister.
 */

#include "mersennetwister.h"

/**
 * Sets the seed of the random number generator to <seed>
 */

static inline void setrandomseed(unsigned long seed)
{
	init_genrand(seed);
}

/**
 *  Returns a random double in [0,1)
 */
static inline double doublerand_1()
{
	return genrand_real2();
}


/**
 *  Returns a random double in [0,<maxVal>)
 */
static inline double doublerand(double maxVal)
{
	return doublerand_1() * maxVal;
}



/**
 * Returns a random integer value in [0,maxVal-1]
 */
static inline unsigned int intrand(unsigned int maxVal)
{
	return genrand_int32() % maxVal;
}

#endif /*RANDOM_H_*/
