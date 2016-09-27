#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <boost/random.hpp>
#include <ctime>
#include <cmath>
#include <vector>

//static boost::mt19937 generator(static_cast<unsigned> (std::time(0)));


// return a random number between 0-1
double random_uni01();

// return a random number between lowest(including) and highest(excluding)
unsigned int get_a_random_number(int lowest, int highest);

// return a random number between lowest(including) and highest(excluding) using boost
unsigned int boost_get_a_random_number(int lowest, int highest);

// return a random number with discrete prob.; pass the cum. distribution vector
unsigned int randomWithDiscreteProbability(const std::vector<double>& accum_prob_vec);

// return a random number with discrete prob.; pass the cum. distribution vector
unsigned int randomWithDiscreteProbability(const std::vector<int>& accum_prob_vec);
#endif
