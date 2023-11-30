/*!
 * \file rng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Class RNG provides wrappers for the GSL random number generate
 */

#ifndef __RNG_H__
#define __RNG_H__

#include <random>
#include <string>
#include <sstream>

namespace VMTutorial
{

	struct RNGState
	{
		std::string seed;
		std::string u_dist;
		std::string n_dist;
	};

	/*! Class handles random numbers in the system */
	class RNG
	{
	public:
		//! Constructor (initialize random number generator)
		RNG(unsigned int seed) : _generator(seed), _uniform_distribution(0.0, 1.0), _normal_distribution(0.0, 1.0) {}

		//! Destructor
		~RNG() {}

		//! Return random number between 0 and 1
		double drnd();

		//! Return random integer between 0 and N
		int lrnd(int);

		//! Return a Gaussian distributed number with a given standard deviation
		double gauss_rng();

		//! Return a Gaussian random number from a normal distribution with a given mean and standard deviation
		double normal(double, double);

		//! Set random number generator from a state
		void set(const RNGState &);

		// Get the random number generator state
		RNGState get_state();

	private:
		std::mt19937_64 _generator;									  //!< Mersenne Twister engine
		std::uniform_real_distribution<double> _uniform_distribution; // Uniform random numbers
		std::normal_distribution<double> _normal_distribution;		  // Gaussian distribution zero mean, unit variance
	};

}

#endif
