/*
 * gp.cpp
 * The main Genetic Programming algorithm
 *
 *  Created on: Apr 16, 2010
 *      Author: muessel
 * 	Edited 2020 by schwab
 */

#include "gp.h"
#include "random.h"
#include "formulaparser.h"
#include "networkconstraint.h"
#include "constraintviolation.h"
#include "helpers.h"
#include <map>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define _HIERARCH
//#define _NONDOMSORT
//#define _TOURNAMENT
//#define _ROUND
#define _EPSILON_COMPARISON
//#define _STRUCTURE_FIRST

using namespace std;


#ifdef _EPSILON_COMPARISON
class SortByAccuracy
{
	private:
		double epsilon1;
		double epsilon2;
		double epsilon3;
	public:

	SortByAccuracy(double epsilon1, double epsilon2, double epsilon3)
	{
		this->epsilon1 = epsilon1;
		this->epsilon2 = epsilon2;
		this->epsilon3 = epsilon3;
	}

	bool operator() (const Individual * i1, const Individual * i2) const
	{
		if (abs(i1->fitness[0] - i2->fitness[0]) <= epsilon1)
		{
//			if (abs(i1->fitness[1] - i2->fitness[1]) < epsilon2 && i1->fitness[1] != i2->fitness[1])
//				cout << i1->fitness[1] << " " << i2->fitness[1] << " " << epsilon2 << endl;
			if (abs(i1->fitness[1] - i2->fitness[1]) <= epsilon2)
				return (i1->fitness[2] < i2->fitness[2] - epsilon3);
			else
				return (i1->fitness[1] < i2->fitness[1]);
		}
		return (i1->fitness[0] < i2->fitness[0]);
	}
};
#else
#ifdef _STRUCTURE_FIRST
/**
 * Comparison operator sorting two individuals
 * by their objective values, starting with
 * objective 1
 */
struct SortByAccuracy
{
	bool operator() (const Individual * i1, const Individual * i2) const
	{
		if (i1->fitness[0] == i2->fitness[0])
		{
			if (i1->fitness[1] == i2->fitness[1])
				return (i1->fitness[2] < i2->fitness[2]);
			else
				return (i1->fitness[1] < i2->fitness[1]);
		}
		return (i1->fitness[0] < i2->fitness[0]);
	}
};
#else
/**
 * Comparison operator sorting two individuals
 * by their objective values, starting with
 * objective 1
 */
struct SortByAccuracy
{
	bool operator() (const Individual * i1, const Individual * i2) const
	{
		if (i1->fitness[0] == i2->fitness[0])
		{
			if (i1->fitness[2] == i2->fitness[2])
				return (i1->fitness[1] < i2->fitness[1]);
			else
				return (i1->fitness[2] < i2->fitness[2]);
		}
		return (i1->fitness[0] < i2->fitness[0]);
	}
};
#endif
#endif

#ifdef _TOURNAMENT
struct Tournament
{
	bool operator()(const Individual * i1, const Individual * i2) const
	{
		if (i1->paretoFront == i2->paretoFront)
		{
			if (i1->fitness[0] == i2->fitness[0])
			{
				if (i1->fitness[1] == i2->fitness[1])
					return (i1->fitness[2] < i2->fitness[2]);
				else
					return (i1->fitness[1] < i2->fitness[1]);
			}
			return (i1->fitness[0] < i2->fitness[0]);
		}
		return (i1->paretoFront < i2->paretoFront);
	}
};
#endif

GP::GP(const vector<string> &originalVars, const vector<BooleanTree *> &originalFunctions,
	   const vector<NetworkConstraint *> &constraints,
	   const vector<double> * topologyWeights,
	   unsigned int maxStartStates,
	   unsigned int calculatedTransitions,
	   double epsilon,
	   bool heuristicCheck)
{
	this->maxStartStates = maxStartStates;
	this->maxTransitions = calculatedTransitions;
	this->epsilon = epsilon;
	this->originalVars.assign(originalVars.begin(), originalVars.end());
	this->heuristicCheck = heuristicCheck;

	for (unsigned int i = 0; i < originalVars.size(); ++i)
		geneIndices[originalVars[i]] = i;

	this->constraints.assign(constraints.begin(), constraints.end());

	if (topologyWeights != NULL)
		this->topologyWeights.assign(topologyWeights->begin(), topologyWeights->end());
//	for (unsigned int i = 0; i < maxAdditionalVars; ++i)
//	{
//		this->originalVars.push_back("UNK_"+i);
//	}

	this->originalFunctions.assign(originalFunctions.begin(), originalFunctions.end());

	for (unsigned int i = 0; i < this->originalFunctions.size(); ++i)
		normalizedFunctions.push_back(this->originalFunctions[i]->getNormalizedCopy());
}

void GP::initialize(unsigned int numIndividuals, unsigned int run, unsigned int initialMutations)
{
	// remove old individuals
	for (vector<Individual *>::iterator it = individuals.begin(); it != individuals.end(); ++it)
	{
		delete *it;
	}
	individuals.clear();

	// create new individuals
	for (unsigned int i = 0; i < numIndividuals; ++i)
	{
		Individual * fi = new Individual(*this,this->originalFunctions,this->originalVars,run,0,initialMutations,true);
		individuals.push_back(fi);
		// calculate initial fitness
		updateFitness(*fi);
	}
}

GP::~GP()
{

	for (vector<Individual *>::iterator it = individuals.begin(); it != individuals.end(); ++it)
	{
		delete *it;
	}
	for (vector<Individual *>::iterator it = resultList.begin(); it != resultList.end(); ++it)
	{
		delete *it;
	}
	for (vector<NetworkConstraint *>::iterator it = constraints.begin(); it != constraints.end(); ++it)
	{
		delete *it;
	}
	for (vector<BooleanTree *>::iterator it = originalFunctions.begin(); it != originalFunctions.end(); ++it)
	{
		delete *it;
	}
	for (vector<BooleanTree *>::iterator it = normalizedFunctions.begin(); it != normalizedFunctions.end(); ++it)
	{
		delete *it;
	}
}

void GP::start(unsigned int numIndividuals, unsigned int numOffspring,
			   unsigned int numIterations, unsigned int numStarts,
			   unsigned int initialMutations,
			   double injectionFraction, unsigned int negationFrac,
			   bool verbose)
{

	unsigned int numNew = (unsigned int)(numOffspring * injectionFraction);
	unsigned int numMutated = numOffspring - numNew;

	if (verbose)
		cout << "Starting optimization..." << endl;

	#ifdef _NONDOMSORT
	cout << "Non-dominated sorting" << endl;
	#endif

	#ifdef _HIERARCH
	cout << "Hierarchical optimization" << endl;
	#endif

	#ifdef _TOURNAMENT
	cout << "Tournament selection" << endl;
	#endif

	for (unsigned int s = 0; s < numStarts; ++s)
	{
		initialize(numIndividuals, s, initialMutations);
		for (unsigned int i = 0; i < numIterations; ++i)
		{
			if (verbose && i % 50 == 49)
				cout << "Run " << s + 1 << " Generation " << i + 1 << endl;

			list<Individual *> domList;
			#ifndef _HIERARCH
			list<Individual *> nonDomList;
			#else
			vector<Individual *> candidates;
			#endif

			candidates.assign(individuals.begin(),individuals.end());

			// create mutated offspring of individuals
			for (unsigned int j = 0; j < numMutated; ++j)
			{
					Individual * newInd = new Individual(*individuals[intrand((unsigned int)individuals.size())]);
											//individuals[intrand((unsigned int)individuals.size())]->copy();
					newInd->generation = i+1;
					newInd->mutate(j % negationFrac == 0);
					updateFitness(*newInd);
					candidates.push_back(newInd);
			}

			// create mutated copies or original network
			for (unsigned int j = 0; j < numNew; ++j)
			{
				Individual * newInd = new Individual(*this,this->originalFunctions,this->originalVars,s,i+1,initialMutations,true);
				candidates.push_back(newInd);
				for (unsigned int k = 0; k < 2; ++k)
				{
					newInd->mutate(j % negationFrac == 0);
				}
				updateFitness(*newInd);
			}

			individuals.clear();

			#ifdef _TOURNAMENT
			nonDomList.assign(candidates.begin(),candidates.end());
			unsigned int front = 0;
			// non-dominated sorting
			do
			{
				dominationSort(nonDomList,domList,getNumFitnessFunctions(),false);
				for (list<Individual *>::iterator it = nonDomList.begin(); it != nonDomList.end(); ++it)
				{
					//individuals.push_back(*it);
					(*it)->paretoFront = front;
				}
				nonDomList.assign(domList.begin(),domList.end());
				domList.clear();
				++front;
			}
			while (nonDomList.size() != 0);

			for (unsigned int j = 0; j < numIndividuals; ++j)
			{
				vector<Individual *> tourn;

				for (unsigned int k = 0; k < 5; ++k)
					tourn.push_back(candidates[intrand(candidates.size())]);
				sort(tourn.begin(),tourn.end(),Tournament());
				individuals.push_back(new Individual(*tourn[0]));
			}

			// free remaining bad candidates
			for (vector<Individual *>::iterator it = candidates.begin(); it != candidates.end(); ++it)
			{
				delete *it;
			}
			#endif

			#ifdef _NONDOMSORT
			// non-dominated sorting
			nonDomList.assign(candidates.begin(),candidates.end());
			do
			{
				dominationSort(nonDomList,domList,getNumFitnessFunctions(),false);
				for (list<Individual *>::iterator it = nonDomList.begin(); it != nonDomList.end(); ++it)
				{
					individuals.push_back(*it);
				}
				nonDomList.assign(domList.begin(),domList.end());
				domList.clear();
			}
			while (individuals.size() < numIndividuals);

			// free remaining bad candidates
			for (list<Individual *>::iterator it = nonDomList.begin(); it != nonDomList.end(); ++it)
			{
				delete *it;
			}

			// free last generation
			for (unsigned int i = numIndividuals; i < individuals.size(); ++i)
			{
				delete individuals[i];
			}

			// too many candidates were selected by non-dominated sorting => free the rest
			if (numIndividuals < individuals.size())
				individuals.erase(individuals.begin() + numIndividuals, individuals.end());
			#endif

			#ifdef _HIERARCH
			#ifdef _EPSILON_COMPARISON
//			sort(candidates.begin(),candidates.end(),
//				 SortByAccuracy(epsilon,epsilon,epsilon));
			nth_element(candidates.begin(),candidates.begin() + numIndividuals, candidates.end(),
						SortByAccuracy(epsilon,epsilon,epsilon));
			#else
//			sort(candidates.begin(),candidates.end(),SortByAccuracy());
			nth_element(candidates.begin(),candidates.begin() + numIndividuals, candidates.end(),
						SortByAccuracy());
			#endif
			unsigned int j = 0;
			for (vector<Individual *>::iterator it = candidates.begin(); it != candidates.end(); ++it, ++j)
			{
				if (j < numIndividuals)
					individuals.push_back(*it);
				else
					delete *it;
			}
			#endif
		}

		// add results of the current run to the result list
		for (vector<Individual *>::iterator it = individuals.begin(); it != individuals.end(); ++it)
		{
			resultList.push_back(new Individual(**it));
		}
	}

	// simplify formulae
	for (vector<Individual *>::iterator it = resultList.begin(); it != resultList.end(); ++it)
	{
		(*it)->simplify();
	}

	// remove duplicate networks from result list
	removeDuplicates(resultList);

	// sort networks by the first two fitness functions
	sort(resultList.begin(), resultList.end(), SortByAccuracy(0,0,0));
}

void GP::stateTransition(bool * state)
{
	bool * oldState = new bool[originalVars.size()];
	memcpy(oldState,state,sizeof(bool) * originalVars.size());

	// calculate next state using the original functions
	for (unsigned int i = 0; i < originalVars.size(); ++i)
	{
		state[i] = originalFunctions[i]->evaluate(oldState);
	}
	delete [] oldState;
}

//inline double netDistance(const std::vector<BooleanTree *> & network1, const std::vector<BooleanTree *> & network2)
//{
//	LiteralMap * literals1;
//	LiteralMap * literals2;
//
//	Literal unknownSymb(-2,NULL,false);
//
//	double res = 0.0;
//	double divisor = 0.0;
//	unsigned int netSize = min((unsigned int)network1.size(), (unsigned int)network2.size());
//
//	for (unsigned int i = 0; i < netSize; ++i)
//	{
//		// extract literals of the current function
//		literals1 = network1[i]->getLiterals();
//		literals2 = network2[i]->getLiterals();
//
////		cout << literals1->size() << " " << literals2->size() << " ";
////		for (LiteralMap::iterator lit = literals1->begin(); lit != literals1->end(); ++lit)
////		{
////			cout << lit->first->literalIndex << endl;
////		}
//
//		// if a function is UNKNOWN, it is excluded from the distance calculation
//		if (literals1->count(&unknownSymb) > 0 || literals2->count(&unknownSymb) > 0)
//		{
//			continue;
//		}
//
//		unsigned int intersectCount = 0;
//		double numOrigLits = 0;
//
//		// calculate intersection size
//		for (LiteralMap::iterator lit = literals2->begin(); lit != literals2->end(); ++lit)
//		{
//			if (literals1->count(lit->first) != 0)
//				//++intersectCount;
////				intersectCount += (lit->second + (*literals2)[lit->first]);
//				intersectCount += min(lit->second,(*literals1)[lit->first]);
//			numOrigLits += lit->second;
//		}
//		numOrigLits *= numOrigLits;
//		divisor += numOrigLits;
//
////		unsigned int unionCount = literals2.size();
////		for (LiteralMap::iterator lit = literals1.begin(); lit != literals1.end(); ++lit)
////		{
////			if (literals2.count(lit->first) == 0)
////				++unionCount;
////		}
//
//		// calculate union size
//		unsigned int unionCount = 0;
//		for (LiteralMap::iterator lit = literals1->begin(); lit != literals1->end(); ++lit)
//		{
//			unionCount += lit->second;
//		}
//		for (LiteralMap::iterator lit = literals2->begin(); lit != literals2->end(); ++lit)
//		{
//			unionCount += lit->second;
//		}
//
//		// calculate Jaccard index with squared weights
//		res += intersectCount/((double)unionCount) * numOrigLits;
//	}
//	//res /= (double)(max(network1.size(), network2.size()));
//
//	//normalization
//	if (divisor != 0)
//		res /= divisor;
//	return res;
//}

inline double netDistance(const std::vector<BooleanTree *> & network1, const std::vector<BooleanTree *> & network2)
{
	LiteralMap * literals1;
	LiteralMap * literals2;

	Literal unknownSymb(-2,NULL,false);

	unsigned int intersectSum = 0, n1Sum = 0, n2Sum = 0;

	unsigned int netSize = min((unsigned int)network1.size(), (unsigned int)network2.size());

	for (unsigned int i = 0; i < netSize; ++i)
	{
		// extract literals of the current function
		literals1 = network1[i]->getLiterals();
		literals2 = network2[i]->getLiterals();


//		cout << literals1->size() << " " << literals2->size() << " ";
//		for (LiteralMap::iterator lit = literals1->begin(); lit != literals1->end(); ++lit)
//		{
//			cout << lit->first->literalIndex << endl;
//		}

		// if a function is UNKNOWN, it is excluded from the distance calculation
		if (literals1->count(&unknownSymb) > 0 || literals2->count(&unknownSymb) > 0)
		{
			continue;
		}

		// calculate intersection size

		for (LiteralMap::iterator lit = literals1->begin(); lit != literals1->end(); ++lit)
		{
			n1Sum += lit->second;
		}

		for (LiteralMap::iterator lit = literals2->begin(); lit != literals2->end(); ++lit)
		{
			n2Sum += lit->second;
			if (literals1->count(lit->first) != 0)
				intersectSum += min(lit->second,(*literals1)[lit->first]);
		}

	}
	double res = intersectSum;

	//normalization
	if (n1Sum != 0 && n2Sum != 0)
		res /= sqrt((double)(n1Sum * n2Sum));

	return res;
}

double calculateFitnessSharing(const Individual &individual, const GP &gp)
{
	double maxDifferent = 0.1;//15;
	//LiteralSet sets[gp.getIndividuals().size()];

//	for (vector<FormulaIndividual *>::iterator it = gp.getIndividuals().begin(); it != gp.getIndividuals().end(); ++it, ++i)
//	{
//		(*it)->getTree()->getLiterals(sets[i]);
//	}
//
//	i = 0;
	LiteralMap * individualSet;
	double fitnessSharingValue = 0.000001;

	for (unsigned int i = 0; i < individual.network.size(); ++i)
	{
		individualSet = individual.network[i]->getLiterals();

		//(*it)->fitnessSharingValue = 0;
		for (unsigned int j = 0; j < gp.getIndividuals().size(); ++j)
		{
			Individual * ind2 = gp.getIndividuals()[j];
			if (ind2 != &individual && i < ind2->network.size())
			{
				LiteralMap * map;
				map = ind2->network[i]->getLiterals();

				unsigned int intersectCount = 0;
				for (LiteralMap::iterator lit = individualSet->begin(); lit != individualSet->end(); ++lit)
				{
					if (map->count(lit->first) != 0)
//						++intersectCount;
//						intersectCount += (lit->second + (*map)[lit->first]);
						intersectCount += min(lit->second,(*map)[lit->first]);
				}

//				unsigned int unionCount = map->size();
//				for (LiteralMap::iterator lit = individualSet->begin(); lit != individualSet->end(); ++lit)
//				{
//					if (map->count(lit->first) == 0)
//						++unionCount;
//				}

				unsigned int unionCount = 0;
				for (LiteralMap::iterator lit = individualSet->begin(); lit != individualSet->end(); ++lit)
				{
					unionCount += lit->second;
				}
				for (LiteralMap::iterator lit = map->begin(); lit != map->end(); ++lit)
				{
					unionCount += lit->second;
				}

				double singleFSVal = intersectCount/((double)unionCount);
				//double singleFSVal = (map.size()+individualSet.size()-2*equalCount) / (double)(map.size()+individualSet.size());

				if (singleFSVal <= maxDifferent)
				{
					singleFSVal /= maxDifferent;
					singleFSVal = 1 - pow(singleFSVal,2);
					fitnessSharingValue += singleFSVal;
				}
			}
		}
	}
	return fitnessSharingValue;
}

void GP::updateFitness(Individual &individual)
{
	individual.fitness[0] = 1;

	// verify the rules
	for (unsigned int i = 0; i < constraints.size(); ++i)
	{
		if (heuristicCheck)
		{
			individual.fitness[0] -= constraints[i]->checkConstraintGreedy(individual,
																		   maxStartStates,
																		   maxTransitions);
		}
		else
		{
			individual.fitness[0] -= constraints[i]->checkConstraintDP(individual,
																	   maxStartStates,
																	   maxTransitions);
		}
	}
	//if (individual.fitness[0] < 1e-10)
	//	individual.fitness[0] = 0;
	#ifdef _ROUND
	individual.fitness[0] = roundTo(individual.fitness[0],3);
	#else
	if (individual.fitness[0] < 1e-10)
		individual.fitness[0] = 0;
	#endif

	// calculate the size of the network
	individual.fitness[1] = 0.0;
	for (unsigned int i = 0; i < individual.network.size(); ++i)
	{
		individual.fitness[1] += individual.network[i]->getNumberOfElements();
	}

	// calculate the number of new constant genes
	unsigned int newConst = 0;
	for (unsigned int i = 0; i < individual.varNames.size(); ++i)
	{
		if (typeid(*individual.network[i]->getRoot()) == typeid(Literal))
		{
			if (((Literal *)individual.network[i]->getRoot())->literalIndex == -1 &&
				((typeid(*originalFunctions[i]->getRoot()) != typeid(Literal) ||
				((Literal *)originalFunctions[i]->getRoot())->literalIndex >= 0)))
				++newConst;
		}
	}

	// calculate weighted sum of new constant genes,
	// network size and distance to the original network

	#ifdef _ROUND
	individual.fitness[1] =  roundTo(
							1.0
							- topologyWeights[2]*(1 - newConst/((double)individual.network.size()))
							- topologyWeights[0]/(1.0 + individual.fitness[1])
							- topologyWeights[1]*netDistance(individual.normalizedNetwork, this->normalizedFunctions), 3);
	#else
	individual.fitness[1] = 1.0
							- topologyWeights[2]*(1 - newConst/((double)individual.network.size()))
							- topologyWeights[0]/(1.0 + individual.fitness[1])
							- topologyWeights[1]*netDistance(individual.normalizedNetwork, this->normalizedFunctions);
	#endif

	// calculate robustness
	bool * state1 = new bool[individual.varNames.size()];
	bool * state2 =  new bool[individual.varNames.size()];

	unsigned int dist = 0;
	unsigned int numTests = maxStartStates;//50;
	for (unsigned int i = 0; i < numTests; ++i)
	{
		// create two start states with Hamming distance 1
		for (unsigned int j = 0; j < individual.varNames.size(); ++j)
		{
			state1[j] = (intrand(2) != 0);
			state2[j] = state1[j];
		}

		unsigned int flipIdx = intrand((unsigned int)individual.varNames.size());
		state2[flipIdx] = !state2[flipIdx];

		// calculate state transitions
		individual.stateTransition(state1,NULL);
		individual.stateTransition(state2,NULL);

		// measure Hamming distances of successor states
		for (unsigned int j = 0; j < individual.varNames.size(); ++j)
		{
			if (state1[j] != state2[j])
				++dist;
		}
	}
	delete [] state1;
	delete [] state2;

	// average over the normalized Hamming distances
	#ifdef _ROUND
	individual.fitness[2] = roundTo(dist/((double)individual.varNames.size() * numTests),3);
	#else
	individual.fitness[2] = dist/((double)individual.varNames.size() * numTests);
	#endif

//	double fitnessSharing = calculateFitnessSharing(individual, *this);
//	individual.fitness[2] = (1 + individual.fitness[0]) * fitnessSharing;

	//individual.fitness[1] *= (1+individual.fitness[0]);// * fitnessSharing;

//	double satisfiedDependencies = 0.0;
//	for (unsigned int i = 0; i < originalVars.size(); ++i)
//	{
//		double positiveWeight = 0.0, negativeWeight = 0.0;
//		LiteralMap map;
//		individual.network[i]->getLiterals(map);
//		for (LiteralMap::iterator it = map.begin(); it != map.end(); ++it)
//		{
//			if (it->first->literalIndex >= 0 && it->first->literalIndex < (int)originalVars.size())
//			{
////				cout << it->first->literalIndex << " " << i << endl;
////				cout << this->posGeneDependencies[it->first->literalIndex][i] << " " << this->negGeneDependencies[it->first->literalIndex][i] << endl;
//				if (it->first->literalIndex < (int)this->originalVars.size())
//				{
//					positiveWeight += this->posGeneDependencies[it->first->literalIndex][i];
//					negativeWeight += this->negGeneDependencies[it->first->literalIndex][i];
//				}
//			}
//		}
////		cout << posSums[i] << " " << negSums[i] << endl;
//		if (posSums[i] != 0.0)
//			positiveWeight /= posSums[i];
////			weightedJaccardPos /= originalVars.size();
////		else
////			positiveWeight = 1.0;
//		if (negSums[i] != 0.0)
//			negativeWeight /= negSums[i];
////			weightedJaccardNeg /= originalVars.size();
////		else
////			negativeWeight = 1.0;
//		satisfiedDependencies += 1 + positiveWeight - negativeWeight;
//	}
//	cout << satisfiedDependencies << endl;
//	individual.fitness[1] /= (1.0+satisfiedDependencies/originalVars.size());
//	individual.fitness[1] = (1.0 + netDistance(individual.network, this->originalFunctions)) - 1.5/(1+individual.fitness[1]);

}
