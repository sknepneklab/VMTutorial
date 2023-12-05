/*!
 * \file rng.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Wrappers for the GSL random number generate
 */

#include "rng.hpp"

namespace VMTutorial
{
	//! Get a random number between 0 and 1 drawn from an uniform distribution
	//! \return random number between 0 and 1
	double RNG::drnd()
	{
		return _uniform_distribution(_generator);
	}

	//! Return a random number from a Gaussian distribution with a given standard deviation
	//! \param sigma standard deviation
	double RNG::gauss_rng()
	{
		return _normal_distribution(_generator);
	}

	//! Return a random number from a Gaussian distribution with a given mean and standard deviation
	double RNG::normal(double avg, double std)
	{
		return avg + std * this->gauss_rng();
	}

	//! Get an integer random number between 0 and N drawn from an uniform distribution
	//! \param N upper bound for the interval
	//! \return integer random number between 0 and N
	int RNG::lrnd(int N)
	{
		return static_cast<int>(N * drnd());
	}

	// Set the random number generator from a state
	void RNG::set(const RNGState &state)
	{
		std::stringstream seed;
		std::stringstream u_dist;
		std::stringstream n_dist;

		seed.str(state.seed);
		u_dist.str(state.u_dist);
		n_dist.str(state.n_dist);

		seed >> _generator;
		u_dist >> _uniform_distribution;
		n_dist >> _normal_distribution;
	}

	RNGState RNG::get_state()
	{
		RNGState state;
		std::stringstream seedout;
		std::stringstream udistout;
		std::stringstream ndistout;

		seedout << _generator;
		udistout << _uniform_distribution;
		ndistout << _normal_distribution;

		state.seed = seedout.str();
		state.u_dist = udistout.str();
		state.n_dist = ndistout.str();

		return state;
	}

}
