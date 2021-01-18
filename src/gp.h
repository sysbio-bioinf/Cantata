/*
 * gp.h
 * The main Genetic Programming algorithm
 *
 *  Created on: Apr 16, 2010
 *      Author: muessel
 */

#ifndef GP_H_
#define GP_H_

#ifdef __GNUC__
// With gcc >= 4.2, the built-in tr1/unordered_map should work,
// so that streamlined versions of Boost can also be used.
// If gcc < 4.2 is installed, a full Boost installation is needed,
// and unordered_map is taken from Boost.
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 2))
	#include <unordered_map>
#else
	#include <boost/tr1/unordered_map.hpp>
#endif
#else
	#include <boost/tr1/unordered_map.hpp>
#endif

#include "individual.h"
#include "networkconstraint.h"
#include "constraintviolation.h"

#include <vector>
#include <list>
#include <set>

class Individual;

class NetworkConstraint;

class ConstraintViolation;

/**
 * The main wrapper class
 * for the Genetic Programming algorithm
 */
class GP
{
private:
	// The individuals of the current run
	std::vector<Individual *> individuals;

	// A collection of best individuals across multiple runs
	std::vector<Individual *> resultList;

	// The functions of the input network
	std::vector<BooleanTree *> originalFunctions;

	// A normalized version of the functions
	std::vector<BooleanTree *> normalizedFunctions;

	// The names of the genes in the input network
	std::vector<std::string> originalVars;

	// The rule sets used for rating the mutated networks
	std::vector<NetworkConstraint *> constraints;

	// A map mapping the gene names to their indices in the gene name list
	std::unordered_map<std::string, unsigned int> geneIndices;

	// The weights of the subscores in the topology fitness function
	std::vector<double> topologyWeights;

	// The maximum number of start states tested to verify the rule sets
	unsigned int maxStartStates;

	// The maximum number of transitions performed for each start state
	// to verify the rule sets
	unsigned int maxTransitions;

	// Switches heuristic and full check of rules
	bool heuristicCheck;

	// The epsilon for comparison of two fitness values
	double epsilon;

	/**
	 * Initializes the algorithm with <numIndividuals> individuals.
	 * <run> is the number of the current run which is stored for the individuals.
	 * <initialMutations> is the number of mutations to be applied to the first generation
	 */
	void initialize(unsigned int numIndividuals, unsigned int run, unsigned int initialMutations);

public:

	/**
	 * Perform a state transition with the original network for state <state>.
	 */
	void stateTransition(bool * state);

	/**
	 * Recalculates the fitness of individual <individual>
	 */
	void updateFitness(Individual &individual);

	/**
	 * Returns the vector of individuals in the current run
	 */
	inline std::vector<Individual *> getIndividuals() const
	{
		return individuals;
	}

	/**
	 * Returns the result list of the best individuals across
	 * multiple runs of the algorithm
	 */
	inline std::vector<Individual *> getResultList() const
	{
		return resultList;
	}

	/**
	 * Returns the number of variables in the original network
	 */
	inline unsigned int getNumberOfVars() const
	{
		return ((unsigned int)originalFunctions.size());
	}

	/**
	 * Returns the names of the variables in the original network
	 */
	inline std::string * getVarNames() const
	{
		return (std::string *)&originalVars[0];
	}

	/**
	 * Returns the number of objective functions used
	 */
	inline unsigned int getNumFitnessFunctions() const
	{
		return 3;
	}

	/**
	 * Starts the Genetic Programming algorithm. Parameters:
	 * <numIndividuals>: The number of individuals per generation
	 * <numOffspring>: The number of offspring in each generation
	 * <numIterations>: The number of generations
	 * <numStarts>: The number of independent starts of the algorithm
	 * <verbose>: Specifies whether additional status messages are printed.
	 */
	void start(unsigned int numIndividuals, unsigned int numOffspring,
			   unsigned int numIterations, unsigned int numStarts,
			   unsigned int initialMutations,
			   double injectionFraction, unsigned int negationFrac,
			   bool verbose);

	/**
	 * Initializes the GP algorithm.
	 * <originalVars>: The names of the variables in the input network
	 * <originalFunctions>: The transition functions of the original network
	 * <constraints>: A vector of rule sets used for validating and rating the networks
	 * <maxStartStates>: The maximum number of start states tested to verify the rule sets
	 * <topologyWeights>: The weights for the subscores in the topology fitness function, or NULL
	 * <calculatedTransitions>: The maximum number of transitions performed for each start state
	 * <epsilon> The epsilon for comparison of two fitness values
	 * <heuristicCheck>: Switches heuristic and full checking of rules
	 * to verify the rule sets
	 */
	GP(const std::vector<std::string> &originalVars, const std::vector<BooleanTree *> &originalFunctions,
	   const std::vector<NetworkConstraint *> &constraints,
	   const std::vector<double>  * topologyWeights,
	   unsigned int maxStartStates,
	   unsigned int calculatedTransitions,
	   double epsilon,
	   bool heuristicCheck);

	virtual ~GP();

};

#endif /* BOOLEANGP_H_ */
