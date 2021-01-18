/*
 * networkconstraint.h
 * Provides classes to model the
 * rule sets for the specification
 * of expectations for the target networks
 *
 *  Created on: Nov 26, 2010
 *      Author: muessel
 */

#ifndef NETWORKCONSTRAINT_H_
#define NETWORKCONSTRAINT_H_

#ifdef __GNUC__
// With gcc >= 4.2, the built-in /unordered_map should work,
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
#include <vector>
#include "constraintviolation.h"
#include "individual.h"
#include "helpers.h"

class ConstraintViolation;

class Individual;

/**
 * Specifies a simple condition consisting
 * of a conjunction of expected variable values
 */
class SimpleCondition
{
private:

	friend class ConstraintViolation;
	friend void finishConstraint(NetworkConstraint * constraint);

	// the index of the condition in the parent list
	unsigned int index;

	// the alternative the condition belongs to
	unsigned int alternative;

	// a mapping of gene indices to their desired values
	std::unordered_map<unsigned int, bool > geneValues;
public:

	/**
	 * Create a SimpleCondition object which has index <index> in the
	 * parent list and belongs to alternative <alternative>
	 */
	SimpleCondition(unsigned int index, unsigned int alternative);

	/**
	 * Create a SimpleCondition by parsing <cond> and using the variable names supplied in <varNames>.
	 * The returned object has index <index> in the parent list and belongs to alternative <alternative>.
	 */
	SimpleCondition(std::string &cond, std::vector<std::string> &varNames, unsigned int index, unsigned int alternative);

	/**
	 * Determine the fraction of satisfied conditions in <state>
	 */
	double getFractionSatisfied(bool * state) const;

	/**
	 * Returns true if all conditions are satisfied by <state>,
	 * and false otherwise.
	 */
	bool isSatisfied(bool * state) const;

	/**
	 * Determines all conditions that are violated by <currentState> for matching <matching>,
	 * and adds them to the constraint violation object <res>.
	 * This requires the start state <startState>
	 * that yielded <currentState>.
	 */
	void getViolations(bool * currentState, bool * startState,
					   ConstraintViolation &res, unsigned int matching);

	/**
	 * Generates a state matching the conditions in the object,
	 * and write it to <result>. Here, <size> is the number of genes
	 * in <result>.
	 */
	void generateMatchingState(bool * result, unsigned int size);

	/**
	 * Generates the <number>-th state matching the conditions in the object,
	 * and write it to <result>. Here, <size> is the number of genes
	 * in <result>. This means filling the binary representation of <number>
	 * into the "free" (unspecified) genes.
	 */
	void generateStateNo(bool * result, unsigned int size, unsigned int number);

	/**
	 * Retrieves the index of the condition object in the parent list
	 */
	unsigned int getIndex()
	{
		return index;
	}

	/**
	 * Returns the alternative the object belongs to
	 */
	unsigned int getAlternative()
	{
		return alternative;
	}

	unsigned int getNumConditions() const;

	/**
	 * Copies the object
	 */
	SimpleCondition * copy();

	inline bool contains(unsigned int const&gene) const
	{
		return geneValues.count(gene) > 0;
	}

	inline const bool& operator[](unsigned int const&index)
	{
		return geneValues[index];
	}

};

/**
 * Specifies a full rule set for an attractor
 * or a chain with its preconditions
 */
class NetworkConstraint
{
private:
	friend void loadNetworkConstraints(std::string file, std::vector<std::string> &varNames, std::vector<NetworkConstraint *> &res);
	friend void finishConstraint(NetworkConstraint * constraint);
	friend class ConstraintViolation;

	// the specification of the initial state
	SimpleCondition * initialValues;

	// the specification of overexpressed/knocked-out genes
	SimpleCondition * fixedGenes;

	// a list of preconditions
	std::vector<SimpleCondition *> preconditions;

	// the state specification list
	std::vector<SimpleCondition *> stateSpecifications;

	// a vector specifying which entries of <stateSpecifications>
	// belong to which alternative
	std::vector<unsigned int> groupIndices;

	// the importance weight of the rule
	double importance;

	// is this an attractor rule or a chain rule?
	bool isAttractor;
public:

	/**
	 * Calculates the score of <ind> according to the constraint in a greedy manner.
	 * <maxStartStates> is the number of start states to use for the estimation.
	 * <calculatedTransitions is the maximum number of state transitions to calculate.
	 * Returns the score.
	 */
	double checkConstraintGreedy(Individual &ind, unsigned int maxStartStates, unsigned int calculatedTransitions) const;

	/**
	 * Calculates the score of <ind> according to the constraint using Dynamic Programming.
	 * <maxStartStates> is the number of start states to use for the estimation.
	 * <calculatedTransitions is the maximum number of state transitions to calculate.
	 * Returns the score.
	 */
	double checkConstraintDP(Individual &ind, unsigned int maxStartStates, unsigned int calculatedTransitions) const;

	/**
	 * Determines rule violations of <ind> according to the constraint in a greedy manner.
	 * <maxStartStates> is the number of start states to use for the estimation.
	 * <calculatedTransitions is the maximum number of state transitions to calculate.
	 * Returns a list of violations.
	 */
	ConstraintViolation * getConstraintViolationsGreedy(Individual &ind,
										   unsigned int maxStartStates,
										   unsigned int calculatedTransitions,
										   bool useRandomStates) const;

	/**
	 * Determines rule violations of <ind> according to the constraint using Dynamic Programming.
	 * <maxStartStates> is the number of start states to use for the estimation.
	 * <calculatedTransitions is the maximum number of state transitions to calculate.
	 * Returns a list of violations.
	 */
	ConstraintViolation * getConstraintViolationsDP(Individual &ind,
										   unsigned int maxStartStates,
										   unsigned int calculatedTransitions,
										   bool useRandomStates) const;

	virtual ~NetworkConstraint();

};

#endif /* NETWORKCONSTRAINT_H_ */
