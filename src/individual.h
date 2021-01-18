/*
 * Individual.h
 * Defines a wrapper class for an individual
 * of the GP approach
 *
 *  Created on: Apr 15, 2010
 *      Author: muessel
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#define TRACE_EVOLUTION

#include "networkconstraint.h"
#include "booleanformula.h"
#include "gp.h"
#include "random.h"
#include "helpers.h"
#include <vector>

class GP;
class SimpleCondition;

/**
 * Represents a single individual in the GP algorithm
 */
class Individual
{
private:
	friend class GP;
	friend class NetworkConstraint;
	friend void removeDuplicates(std::vector<Individual *> &population, bool freeDuplicates);

	// The names of the variables in the individual
	// (usually the same for all individuals, but can be used to
	// extend the algorithm to introduce unknown genes)
	std::vector<std::string> varNames;

	// The parent GP algorithm
	GP * parent;

	// Updates the cached normalized representation of gene <index>.
	void updateCachedNormalization(unsigned int index);

public:
	// The fitness vector of the individual
	double * fitness;

	// The run in which the individual was created
	unsigned int run;

	// The generation in which the individual was created
	unsigned int generation;

	// The Pareto front on which the individual is located
	unsigned int paretoFront;

	// The transition functions of the network represented by the individual
	std::vector<BooleanTree *> network;

	// A normalized cached version of <network>
	std::vector<BooleanTree *> normalizedNetwork;

#ifdef TRACE_EVOLUTION
	std::vector<Individual *> ancestors;
#endif

	/**
	 * Creates a new individual for the algorithm instance <parent> and the network <network> with
	 * variable names <varNames>.
	 * <run> denotes the number of the run in which the individual is created.
	 * <generation> denotes the number of the generation in which the individual is created.
	 * <startMutations> is the number of mutations applied to <network>.
	 */
	Individual(GP &parent, std::vector<BooleanTree *> &network, std::vector<std::string> &varNames,
			   unsigned int run, unsigned int generation,
			   unsigned int startMutations, bool replaceUnknown = true);

	/**
	 * Creates an exact copy of <clone>
	 */
#ifdef TRACE_EVOLUTION
	Individual(Individual &clone, bool copyAncestors = true);
#else
	Individual(Individual &clone);
#endif

	/**
	 * Returns an exact copy of the individual
	 */
	Individual * copy();

	/**
	 * Applies a mutation to the individual.
	 * If <allowNegation> is true, the negation operation
	 * is allowed as a mutation.
	 */
	void mutate(bool allowNegation = false);

	/**
	 * Simplifies the functions of an individual by removing
	 * duplicates, identifying tautologies, etc.
	 */
	void simplify();

	/**
	 * Performs a state transition of the network represented by the individual
	 * starting from state <initialState>.
	 */
	void stateTransition(bool * initialState, SimpleCondition * fixedGenes);

	/**
	 * Outputs the network represented by the individual.
	 * <separator> is the character that separates the targets from the functions.
	 * If <printFitness> is true, the fitness values are included in the output.
	 * If <printRun> is true, the run number is included in the output.
	 * Returns the string representation of the network.
	 */
	std::string toString(std::string separator = " = ",
						 bool printFitness = true,
						 bool printRun = true) const;
	/**
	 * Tests two individuals for equality.
	 */
	virtual bool operator==(const Individual &ind) const;

	virtual ~Individual();
};

/**
 * Returns a random Boolean operator with no children.
 * <parent> is an optional parent operator.
 */
inline BooleanOperator * randomOperator(BooleanOperator * parent = NULL)
{
	return new BooleanOperator((OperatorType)intrand(2),parent,false);//(bool)intrand(2));
}

/**
 * Returns a randomly chosen Literal.
 * <numberOfVariables> is the total number of variables in the network.
 * <parent> is an optional parent operator.
 */
inline Literal * randomLiteral(unsigned int numberOfVariables, BooleanOperator * parent = NULL)
{
	//return new Literal(intrand(numberOfVariables + 1) - 1,parent,(bool)intrand(2));
	return new Literal(intrand(numberOfVariables),parent,intrand(2) != 0);
}

#endif /* INDIVIDUAL_H_ */
