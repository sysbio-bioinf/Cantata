/*
 * constraintviolation.h
 * Provides classes to monitor
 * violations of rule sets for
 * manual analysis
 *
 *  Created on: Nov 25, 2010
 *      Author: muessel
 */

#ifndef CONSTRAINTVIOLATION_H_
#define CONSTRAINTVIOLATION_H_

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
#include <map>

#include "networkconstraint.h"
#include "helpers.h"
#include <vector>

class SimpleCondition;

/**
 * Stores the violation of a single variable restriction
 * in a SimpleCondition rule
 */
class ConditionViolation
{
private:
	// The condition object that is violated
	SimpleCondition * condition;

	// The index of the specification in the parent specification list
	unsigned int specificationEntry;

	// The index of the restriction in <condition> that is violated
	unsigned int conditionIndex;
public:


	/**
	 * Creates a new ConditionViolation object with parameters:
	 * <condition>: The condition object that is violated
	 * <specificationEntry>: The index of the specification in the parent specification list
	 * <conditionIndex>: The index of the restriction in <condition> that is violated
	 */
	ConditionViolation(SimpleCondition &condition, unsigned int specificationEntry, unsigned int conditionIndex)
	{
		this->condition = &condition;
		this->conditionIndex = conditionIndex;
		this->specificationEntry = specificationEntry;
	}

	/**
	 * Retrieves the alternative that is violated
	 */
	unsigned int getAlternative() const;

	/**
	 * Retrieves the index of the restriction that is violated
	 */
	unsigned int getConditionIndex() const;

	/**
	 * Retrieves the index of the specification in the parent specification list
	 */
	unsigned int getSpecificationEntry() const;

	/**
	 * Retrieves the condition object that contains the violated condition
	 */
	SimpleCondition * getCondition() const;
};

/**
 * Hash calculation for condition violations
 */
struct HashConditionViolations
{

	size_t operator() (const ConditionViolation &val) const
	{
		return val.getConditionIndex() + val.getAlternative() + val.getSpecificationEntry();
	}
};

/**
 * Test for equality of two condition violations for hashing
 */
struct CompareConditionViolations
{

	bool operator() (const ConditionViolation &cond1, const ConditionViolation &cond2) const
	{
		return (cond1.getConditionIndex() == cond2.getConditionIndex()
				&& cond1.getAlternative() == cond2.getAlternative()
				&& cond1.getSpecificationEntry() == cond2.getSpecificationEntry());
	}
};

class CompareStatesLess;
typedef std::map<bool *, int, CompareStatesLess> StateSet;
typedef std::unordered_map<ConditionViolation, StateSet *, HashConditionViolations, CompareConditionViolations> ViolationMap;

class NetworkConstraint;

/**
 * A class that stores all violations of a set of rules
 */
class ConstraintViolation
{
private:

	// The number of genes in the network
	unsigned int numGenes;

	// The constraint for which violations are listed
	NetworkConstraint * constraint;

	// A list of attractors identified in the network
	std::vector<std::vector<bool *> > attractors;

	// A list of matchings of attractor states and rules
	std::vector<unsigned int *> attractorMatchings;

	// A list storing which rule alternative was chosen for the matching
	std::vector<unsigned int> alternatives;

	// A list of maps holding the violations for the attractor matchings
	// (or in case of a chain matching, a vector with a single entry)
	// and a vector of initial states that cause the violation
	std::vector<ViolationMap> violations;

	CompareStatesLess * comp;

public:

	/**
	 * Creates a new ConstraintViolation object
	 * and initializes the number of genes with <numGenes>.
	 * If <singleViolationList> is true, only one violation list is maintained
	 * (not a separate one for each matching).
	 */
	ConstraintViolation(NetworkConstraint * constraint, unsigned int numGenes,
						bool singleViolationList);

	virtual ~ConstraintViolation();

	/**
	 * Add a new violation to the list.
	 * <matching> is the index of the matching in the list
	 * of attractor matchings, or 0 for a chain matching
	 * <state> is the state that causes the violation.
	 * <violation> is a description of the violation.
	 */
	void addViolation(unsigned int matching, bool * state, ConditionViolation &violation);

	/**
	 * Adds an assignment of attractor states to rule entries
	 * to the internal list.
	 * <attractor> is a vector of states describing an attractor
	 * <assignment> is a vector describing which states of the attractor
	 * are matched with which rules.
	 * <alternative> is the rule set chosen for the matching.
	 */
	void addAttractorMatching(std::vector<bool *> attractor, unsigned int * assignment,
										unsigned int alternative);

	/**
	 * Returns a string representation of the violations based on the
	 * variable names in <varNames>.
	 * <collapse> specifies the maximum number of initial states that are printed
	 * for a violation. All further states are skipped.
	 */
	std::string toString(std::string * varNames, unsigned int collapse);

	/**
	 * Returns the number of matchings in the violation list.
	 */
	inline unsigned int getNumMatchings()
	{
		return (unsigned int)violations.size();
	}
};

#endif /* CONSTRAINTVIOLATION_H_ */
