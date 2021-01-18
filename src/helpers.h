/*
 * helpers.h
 * Helper functions for loading files,
 * non-dominated sorting, string processing, etc.
 *
 *  Created on: Nov 26, 2010
 *      Author: muessel
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include <vector>
#include <list>
#include <cmath>
#include <set>
#include "booleanformula.h"
#include "networkconstraint.h"
#include "individual.h"

/**
 * Allocates a 2-dimensional <dim1>x<dim2> array of type <type>
 * in variable <name>
 */
#define ALLOC_2D_ARRAY(name, type, dim1, dim2)	name = new type * [dim1]; \
												for (unsigned int i = 0; i < dim1; ++i) \
													name[i] = new type[dim2];

/**
 * Allocates a 3-dimensional <dim1>x<dim2>x<dim3> array of type <type>
 * in variable <name>
 */
#define ALLOC_3D_ARRAY(name, type, dim1, dim2, dim3)	name = new type ** [dim1]; \
												for (unsigned int i = 0; i < dim1; ++i) \
												{\
													name[i] = new type * [dim2];\
													for (unsigned int j = 0; j < dim2; ++j)\
														name[i][j] = new type[dim3];\
												}

/**
 * Deletes a 3-dimensional array in <name> which has <dim1> elements
 * in the first dimension and <dim2> elements in the second dimension
 */
#define DELETE_3D_ARRAY(name, dim1, dim2) for (unsigned int i = 0; i < dim1; ++i) \
		{\
			for (unsigned int j = 0; j < dim2; ++j)\
				delete [] name[i][j];\
		delete [] name[i]; \
		}\
	delete [] name;

/**
 * Deletes a 2-dimensional array in <name> which has <dim1> elements
 * in the first dimension
 */
#define DELETE_2D_ARRAY(name, dim1) for (unsigned int i = 0; i < dim1; ++i) \
										delete [] name[i]; \
									delete [] name;

class NetworkConstraint;

class Individual;

/**
 * Rounds <value> to <digits> digits after the decimal point.
 */
inline double roundTo(double value, unsigned int digits)
{
	return floor(((value*pow(10.0,(int)digits)) / pow(10.0,(int)digits)) + 0.5);
}

/**
 * Replaces the format sign "%d" in <fmt> by <num> and return
 * the corresponding string.
 */
std::string insertNumber(std::string fmt, unsigned int num);

/**
 * Removes leading and trailing whitespaces in <str>
 */
void trim(std::string& str);

/**
 * Loads the rule sets to rate the target networks.
 * <varNames> is the set of known variable names.
 * <res> is a vector of rule sets that receives the results.
 */
void loadNetworkConstraints(std::string file, std::vector<std::string> &varNames, std::vector<NetworkConstraint *> &res);

/**
 * Loads additional genetic dependencies -- currently not used
 */
void loadGeneDependencies(std::string file, std::vector<std::string> &varNames, double ** pos, double ** neg);

/**
 * Loads a network from <file>.
 * <varNames> is a vector that receives the names of the variables in the network.
 * <res> is a vector that receives the transition functions of the network.
 * If <simplify> is true, the rules of the network are simplified.
 */
void loadNetwork(std::string file, std::vector<std::string> &varNames, std::vector<BooleanTree *> &res, bool simplify);

/**
 * Removes duplicate networks from a list of individuals.
 * If <freeDuplicates> is true, the eliminated duplicate individuals are freed.
 */
void removeDuplicates(std::vector<Individual *> &population, bool freeDuplicates = true);

/**
 * Extracts non-dominated individuals from <nonDom>.
 * <dom> receives a list of dominated individuals, whereas the non-dominated individuals remain in <nonDom>.
 * <numFitnessFunctions> is the number of fitness functions in the optimization.
 * If <weak> is true, weak domination is used. Otherwise, "normal" domination is used.
 */
template <typename Ind > inline void dominationSort(std::list<Ind *> &nonDom, std::list<Ind *> &dom, unsigned int numFitnessFunctions, bool weak)
{
	for (typename std::list<Ind *>::iterator it1 = nonDom.begin(); it1 != nonDom.end(); )
	{
		bool dominated = false;
		for (typename std::list<Ind *>::iterator it2 = nonDom.begin(); it2 != nonDom.end(); ++it2)
			if (it2 != it1)
			{
				bool weaklyDominated = true;
				for (unsigned int i = 0; i < numFitnessFunctions; ++i)
				{
					if ((*it1)->fitness[i] < (*it2)->fitness[i])
					{
						weaklyDominated=false;
						break;
					}
				}
				if (weak && weaklyDominated)
				{
					dominated = true;
					continue;
				}
				if (weaklyDominated)
				{
					for (unsigned int i = 0; i < numFitnessFunctions; ++i)
					{
						if ((*it1)->fitness[i] > (*it2)->fitness[i])
						{
							dominated = true;
							break;
						}
					}
					if (dominated)
					{
						break;
					}
				}
			}
			if (dominated)
			{
				dom.push_back(*it1);
				it1 = nonDom.erase(it1);
			}
			else
				++it1;
		}
}

/**
 * Calculates a hash value for a state by encoding it
 * as a decimal number
 */
class HashStates
{
private:
	unsigned int numGenes;
public:

	/**
	 * Creates a hash object for states with <numGenes> genes.
	 */
	HashStates(unsigned int numGenes)
	{
		this->numGenes = numGenes;
	}

	/**
	 * Calculates the hash value for state <val>.
	 */
	size_t operator() (const bool * val) const
	{
		unsigned long h = 0;
		for (unsigned int i = 0; i < numGenes; ++i)
			h |= val[i] << i;
		return h;
	}
};

/**
 * Checks the equality of two states
 */
class CompareStatesEqual
{
private:
	unsigned int numGenes;
public:

	/**
	 * Creates a comparator for states with <numGenes> genes.
	 */
	CompareStatesEqual(unsigned int numGenes)
	{
		this->numGenes = numGenes;
	}

	/**
	 * Checks for equality of <state1> and <state2>
	 */
	bool operator() (const bool * state1, const bool * state2) const
	{
		for (unsigned int i = 0; i < numGenes; ++i)
		{
			if (state1[i] != state2[i])
				return false;
		}
		return true;
	}
};

/**
 * Obtains an ordering of states
 * by sorting them according to their decimal representation
 */
class CompareStatesLess
{
private:
	unsigned int numGenes;
public:

	/**
	 * Creates a comparator for states with <numGenes> genes.
	 */
	CompareStatesLess(unsigned int numGenes)
	{
		this->numGenes = numGenes;
	}

	/**
	 * Checks whether <state1> precedes <state2>.
	 */
	const bool operator() (const bool * state1, const bool * state2) const
	{
		for (unsigned int i = 0; i < numGenes; ++i)
		{
			if (!state1[i] && state2[i])
				return true;
			else
			if (state1[i] != state2[i])
				return false;
		}
		return false;
	}

};

extern std::set<std::string> warnings;

extern std::vector<double> * decodeVector(std::string vec);

#endif /* HELPERS_H_ */
