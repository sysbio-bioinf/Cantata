/*
 * networkconstraint.cpp
 * Provides classes to model the
 * rule sets for the specification
 * of expectations for the target networks
 *  Created on: Nov 26, 2010
 *      Author: muessel
 * Edited 2020 by schwab
 */

#include "networkconstraint.h"
#include "helpers.h"
#include <cstring>
#include <algorithm>
#include <cmath>

using namespace std;

SimpleCondition::SimpleCondition(unsigned int index, unsigned int alternative)
{
	this->index = index;
	this->alternative = alternative;
}

SimpleCondition::SimpleCondition(std::string &cond, std::vector<std::string> &varNames, unsigned int index, unsigned int alternative)
{
	this->index = index;
	this->alternative = alternative;

	// read the variable names and their negation status
	size_t lastSep = 0, nextSep = 0;
	while (true)
	{
		if (nextSep == cond.size())
			break;

		// find next separator (whitespace)
		nextSep = cond.find(' ',lastSep);
		if (nextSep == string::npos)
			nextSep = cond.size();

		std::string val = cond.substr(lastSep,nextSep-lastSep);
		trim(val);

		// is the literal negated?
		bool neg = false;
		if (val[0] == '!')
		{
			neg = true;
			val = val.substr(1);
		}

		// check for gene in variable list
		std::vector<std::string>::iterator varIndex = find(varNames.begin(),varNames.end(),val);
		if (varIndex == varNames.end())
			throw string("Unknown symbol " + val + "!");
		geneValues[(unsigned int)(varIndex - varNames.begin())] = !neg;
		lastSep = nextSep+1;
	}
}

double SimpleCondition::getFractionSatisfied(bool * state) const
{
	double res = 0.0;
	for (std::unordered_map<unsigned int, bool >::const_iterator it = geneValues.begin();
		 it != geneValues.end(); ++it)
	{
		if (state[it->first] == it->second)
			++res;
	}
	return (res/geneValues.size());
}

bool SimpleCondition::isSatisfied(bool * state) const
{
	for (std::unordered_map<unsigned int, bool >::const_iterator it = geneValues.begin();
		 it != geneValues.end(); ++it)
	{
		if (state[it->first] != it->second)
			return false;
	}
	return true;
}

void SimpleCondition::getViolations(bool * currentState, bool * startState,
									ConstraintViolation &res, unsigned int matching)
{
	for (std::unordered_map<unsigned int, bool >::const_iterator it = geneValues.begin();
			 it != geneValues.end(); ++it)
	{
		if (currentState[it->first] != it->second)
		{
			ConditionViolation cond(*this,index,it->first);
			res.addViolation(matching, startState, cond);
		}
	}
}

void SimpleCondition::generateMatchingState(bool * result, unsigned int size)
{
	for (unsigned int j = 0; j < size; ++j)
	{
		if (geneValues.count(j) != 0)
		// this gene is fixed => use the corresponding value
		{
			result[j] = geneValues[j];
		}
		else
		// this gene is not specified => generate it
		{
			result[j] = (intrand(2) != 0);
		}
	}
}

void SimpleCondition::generateStateNo(bool * result, unsigned int size, unsigned int number)
{
	unsigned int bit = 0;
	for (unsigned int j = 0; j < size; ++j)
	{
		if (geneValues.count(j) != 0)
		// this gene is fixed => use the corresponding value
		{
			result[j] = geneValues[j];
		}
		else
		// this gene is not specified => use the value for the next bit in <number>
		{
			result[j] = ((number & (1 << bit)) > 0);
			++bit;
		}
	}
}

unsigned int SimpleCondition::getNumConditions() const
{
	return (unsigned int)geneValues.size();
}

SimpleCondition * SimpleCondition::copy()
{
	SimpleCondition * res = new SimpleCondition(index, alternative);
	for (std::unordered_map<unsigned int, bool >::iterator it = geneValues.begin();
		 it != geneValues.end(); ++it)
	{
		res->geneValues[it->first] = it->second;
	}
	return res;
}

double NetworkConstraint::checkConstraintGreedy(Individual &ind, unsigned int maxStartStates, unsigned int calculatedTransitions) const
{
	unsigned int numGroups = (unsigned int)groupIndices.size() - 1;
	double res = 0.0;
	bool * state = new bool[ind.varNames.size()];

	unsigned int possibleStates = ((unsigned int)ind.varNames.size() - initialValues->getNumConditions());
	if (possibleStates > 31)
		possibleStates = ~0;
	else
		possibleStates = (1 << possibleStates);

	bool random = true;
	if (maxStartStates >= possibleStates)
	{
		random = false;
		maxStartStates = possibleStates;
	}

	unsigned int checkedStates = maxStartStates;

	HashStates h1((unsigned int)ind.varNames.size());
	CompareStatesEqual s1((unsigned int)ind.varNames.size());
	unordered_map<bool*, double, HashStates, CompareStatesEqual> savedScores(8,h1,s1,std::allocator<std::pair<bool*, double> >());

	for (unsigned int i = 0; i < maxStartStates; )
	{
//		cout << i << endl;
		if (random)
			initialValues->generateMatchingState(state, (unsigned int)ind.varNames.size());
		else
			initialValues->generateStateNo(state, (unsigned int)ind.varNames.size(), i);

		unsigned int preconditionsVerified = 0;
		if (preconditions.size() > 0)
		{
			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
				if (preconditions[preconditionsVerified]->isSatisfied(state))
				{
					++preconditionsVerified;
					if (preconditionsVerified == preconditions.size())
						break;
				}
				ind.stateTransition(state, fixedGenes);
			}
			if (preconditionsVerified != preconditions.size())
			{
				if (!random)
				{
					++i;
					--checkedStates;
				}
				continue;
			}
		}

		if (isAttractor)
		{
			HashStates h2((unsigned int)ind.varNames.size());
			CompareStatesEqual s2((unsigned int)ind.varNames.size());
			unordered_map<bool*, bool*, HashStates, CompareStatesEqual> states(8,h2,s2,std::allocator<std::pair<bool*, bool*> >());

			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
//				cout << "k=" << k  << endl;
				if (states.count(state) == 0)
				{
					bool * startState = new bool[ind.varNames.size()];
					memcpy(startState, state, sizeof(bool) * ind.varNames.size());

					bool * nextState = new bool[ind.varNames.size()];
					ind.stateTransition(state, fixedGenes);
					memcpy(nextState, state, sizeof(bool) * ind.varNames.size());

					states[startState] = nextState;
				}
				else
				{

//					cout << "Found attractor after " << k << " transitions."  << endl;

					if (savedScores.count(state) > 0)
					{
						res += savedScores[state];
						break;
					}

					vector<bool * > attractor;
					bool * startState = new bool[ind.varNames.size()];
					memcpy(startState, state, sizeof(bool) * ind.varNames.size());
					bool * curState = state;
					do
					{
						attractor.push_back(curState);
						curState = states[curState];
					}
					while (memcmp(curState,startState,sizeof(bool) * ind.varNames.size()) != 0);

					delete [] startState;
//					cout << "Attractor size: " << attractor.size() << endl;
//					for (unsigned int l = 0; l < attractor.size(); ++l)
//					{
//						for (unsigned int p = 0; p < ind.varNames.size(); ++p)
//						{
//							cout << ind.varNames[p] << "=" << attractor[l][p] << " ";
//						}
//						cout << endl;
//					}

					double bestGroupScore = 0.0;
					for (unsigned int r = 0; r < numGroups; ++r)
					{
						unsigned int specSize = groupIndices[r+1] - groupIndices[r];

//						cout << "r: " << r << " Rule size: " << ruleSize << " " << groupIndices[r+1] << " " << groupIndices[r] << endl;
						double maxScore = 0.0;
//						unsigned int conditionsVerified = 0;
						for (unsigned int l = 0; l < attractor.size(); ++l)
						{
	//						cout << "l=" << l  << endl;
							unsigned int stateNo = 0;
							double score = 0.0;
							int currentCondition = -1;
							unsigned int maxTrans = (unsigned int)attractor.size();//max((size_t)specSize,attractor.size());
							for (stateNo = 0; stateNo < maxTrans; ++stateNo)

							{
								double subscore1 = 0, subscore2;

//								cout << "CC " << groupIndices[r] + (currentCondition % ruleSize) << " "
//									<< groupIndices[r] + ((currentCondition+ 1) % ruleSize) << endl;

								if (currentCondition != -1)
									subscore1 = stateSpecifications[groupIndices[r] + (currentCondition % specSize)]->getFractionSatisfied(attractor[(stateNo+l) % attractor.size()]);

								subscore2 = stateSpecifications[groupIndices[r] + ((currentCondition + 1) % specSize)]->getFractionSatisfied(attractor[(stateNo+l) % attractor.size()]);
	//							cout << "S " << subscore2 << " " << subscore1 << endl;

								if (subscore2 >= subscore1)
								{
	//								cout << (currentCondition+1) % rules.size() << " " << (stateNo+l) % attractor.size() << endl;
									score += subscore2;
									currentCondition = currentCondition + 1;
								}
								else
								{
	//								cout << currentCondition % rules.size() << " " << (stateNo+l) % attractor.size() << endl;
									score += subscore1;
								}
	//							cout << "Score " << score << endl;
	//							++stateNo;
	//							cout << stateNo << endl;
							}
							score /= maxTrans;

							if ((int)specSize > currentCondition + 1)
								score *= (currentCondition + 1)/(double)specSize;
							if (score > maxScore)
							{
								maxScore = score;
//								conditionsVerified = currentCondition + 1;
							}
						}

//						cout << "MaxScore: " << maxScore << endl;

						if (maxScore > bestGroupScore)
							bestGroupScore = maxScore;
					}
					res += bestGroupScore;

					for (vector<bool *>::iterator it = attractor.begin(); it != attractor.end(); ++it)
					{
						bool * state = new bool[ind.varNames.size()];
						memcpy(state,*it,sizeof(bool)*ind.varNames.size());
						savedScores[state] = bestGroupScore;
					}
					break;
				}
			}
			for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
				it != states.end(); ++it)
			{
				delete [] it->first;
				delete [] it->second;
			}
		}
		else
		{
			bool * startState = new bool[ind.varNames.size()];
			memcpy(startState, state, sizeof(bool) * ind.varNames.size());
			double bestGroupScore = 0.0;
			bool exitLoop = false;

			for (unsigned int r = 0; r < numGroups; ++r)
			{
				if (exitLoop)
					break;
				memcpy(state, startState, sizeof(bool) * ind.varNames.size());

				double maxScore = 0.0;
				unsigned int conditionVerified = 0;
				unsigned int specSize = groupIndices[r+1] - groupIndices[r];

				for (unsigned int k = 0; k < calculatedTransitions; ++k)
				{

					double score = stateSpecifications[groupIndices[r] + conditionVerified]->getFractionSatisfied(state);
					if (score > maxScore)
					{
						maxScore = score;
					}
					if (score == 1.0)
					{
						++conditionVerified;
						maxScore = 0.0;
					}
					if (conditionVerified==specSize)
					{
						exitLoop = true;
						break;
					}
					ind.stateTransition(state, fixedGenes);
				}
				double specScore = (conditionVerified + maxScore)/((double)specSize);

				if (specScore > bestGroupScore)
					bestGroupScore = specScore;

			}
			delete [] startState;
			res += bestGroupScore;
		}
		++i;
	}

	delete [] state;
	for (unordered_map<bool*, double, HashStates, CompareStatesEqual>::iterator it = savedScores.begin();
		it != savedScores.end(); ++it)
	{
		delete [] it->first;
	}

	return res / checkedStates * importance;
}

double NetworkConstraint::checkConstraintDP(Individual &ind, unsigned int maxStartStates, unsigned int calculatedTransitions) const
{
	unsigned int numGroups = (unsigned int)groupIndices.size() - 1;
	double res = 0.0;
	bool * state = new bool[ind.varNames.size()];

	unsigned int possibleStates = ((unsigned int)ind.varNames.size() - initialValues->getNumConditions());
	if (possibleStates > 31)
		possibleStates = ~0;
	else
		possibleStates = (1 << possibleStates);

	bool random = true;
	if (maxStartStates >= possibleStates)
	{
		random = false;
		maxStartStates = possibleStates;
	}

	// maintain a hash table of attractor states for which a score has already been calculated
	HashStates h1((unsigned int)ind.varNames.size());
	CompareStatesEqual s1((unsigned int)ind.varNames.size());
	unordered_map<bool*, double, HashStates, CompareStatesEqual> savedScores(8,h1,s1,std::allocator<std::pair<bool*, double> >());

	unsigned int maxSpecSize = 0;
	for (unsigned int r = 0; r < numGroups; ++r)
	{
		unsigned int specSize = groupIndices[r+1] - groupIndices[r];
		if (specSize > maxSpecSize)
			maxSpecSize = specSize;
	}

	unsigned int checkedStates = maxStartStates;

	for (unsigned int i = 0; i < maxStartStates; )
	{
//		cout << i << endl;

		// generate a start state
		if (random)
			initialValues->generateMatchingState(state, (unsigned int)ind.varNames.size());
		else
			initialValues->generateStateNo(state, (unsigned int)ind.varNames.size(), i);

		// verify preconditions by repeated state calculations
		unsigned int preconditionsVerified = 0;
		if (preconditions.size() > 0)
		{
			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
				if (preconditions[preconditionsVerified]->isSatisfied(state))
				{
					++preconditionsVerified;
					if (preconditionsVerified == preconditions.size())
						break;
				}
				ind.stateTransition(state, fixedGenes);
			}
			if (preconditionsVerified != preconditions.size())
			{
				if (!random)
				{
					++i;
					--checkedStates;
				}
				continue;
			}
		}

		// save sequence of states or attractor in a hash table
		HashStates h2((unsigned int)ind.varNames.size());
		CompareStatesEqual s2((unsigned int)ind.varNames.size());
		unordered_map<bool*, bool*, HashStates, CompareStatesEqual> states(8,h2,s2,std::allocator<std::pair<bool*, bool*> >());
		vector<bool * > stateList;

		unsigned int k;
		for (k = 0; k < calculatedTransitions; ++k)
		{
//				cout << "k=" << k  << endl;
			if (states.count(state) == 0)
			// this state has not been reached previously
			{
				bool * startState = new bool[ind.varNames.size()];
				memcpy(startState, state, sizeof(bool) * ind.varNames.size());

				bool * nextState = new bool[ind.varNames.size()];
				ind.stateTransition(state, fixedGenes);
				memcpy(nextState, state, sizeof(bool) * ind.varNames.size());

				states[startState] = nextState;

				if (!isAttractor && states.count(nextState) == 0)
				// if this is a chain, add the state to the state list
				{
					stateList.push_back(nextState);
				}
			}
			else
			if (!isAttractor && stateList.size() < maxSpecSize)
			// if the state list is too short, add the attractor states as well
			{
				bool * nextState = states[state];
				stateList.push_back(nextState);
			}
			else
				break;
		}

		if (k == calculatedTransitions)
		// The number of transitions was not sufficient to reach the attractor
		{
			warnings.insert("The maximum number of transitions is not sufficient to reach the attractor! It is recommended to increase \"-mt\"!");
			for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
					it != states.end(); ++it)
			{
				delete [] it->first;
				delete [] it->second;
			}
			continue;
		}

		if (isAttractor)
		// this is an attractor rule
		{

			if (savedScores.count(state) > 0)
			// the score for this attractor has already been calculated
			// => look it up in the hash table
			{
				res += savedScores[state];
				++i;
				for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
						it != states.end(); ++it)
				{
					delete [] it->first;
					delete [] it->second;
				}
				continue;
			}

			// determine the attractor
			bool * startState = new bool[ind.varNames.size()];
			memcpy(startState, state, sizeof(bool) * ind.varNames.size());
			bool * curState = state;
			do
			{
				stateList.push_back(curState);
				curState = states[curState];
			}
			while (memcmp(curState,startState,sizeof(bool) * ind.varNames.size()) != 0);
			delete [] startState;

//			cout << "Attractor size: " << stateList.size() << endl;
//			for (unsigned int l = 0; l < stateList.size(); ++l)
//			{
//				for (unsigned int p = 0; p < ind.varNames.size(); ++p)
//				{
//					cout << ind.varNames[p] << "=" << stateList[l][p] << " ";
//				}
//				cout << endl;
//			}

			// calculate the scores for all alternatives and choose the best one
			double bestGroupScore = 0.0;

			for (unsigned int r = 0; r < numGroups; ++r)
			{
				unsigned int specSize = groupIndices[r+1] - groupIndices[r];

				for (unsigned int l = 0; l < stateList.size(); ++l)
				// iterate over attractor state offsets
				{
					double *** ffTable;
					ALLOC_3D_ARRAY(ffTable, double, specSize, stateList.size() + 1, specSize + 1);
					for (unsigned int u = 0; u < specSize; ++u)
					// iterate over the maximum number of uncovered specifications
					{

						// initialize Dynamic Programming matrix

						ffTable[u][0][0] = 0;

						for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
						{
							ffTable[u][stateIdx][0] = -1e10;
						}

						for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
						{
							if (specIdx <= u)
								ffTable[u][0][specIdx] = 0;
							else
								ffTable[u][0][specIdx] = -1e10;
						}

						// main DP loops
						for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
						{
							for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
							{
								double fractionFulfilled = stateSpecifications[groupIndices[r] + specIdx - 1]->getFractionSatisfied(stateList[(stateIdx - 1 + l) % stateList.size()]);

								double scores[3];
								unsigned int maxIdx = 0;

//									cout << "FF " << ffTable[stateIdx][specIdx - 1] << " " << coveredSpecTable[stateIdx][specIdx-1] << endl;
								scores[0] = (fractionFulfilled + ffTable[u][stateIdx - 1][specIdx - 1]);
								scores[1] = (fractionFulfilled + ffTable[u][stateIdx - 1][specIdx]);
                                
                                //improve prediction, go one field to the left. 
                                if (u > 0)
									scores[2] = ffTable[u-1][stateIdx][specIdx - 1];
								else
									scores[2] = -1e10;

//									cout << stateIdx << " " << specIdx << ": " << scores[0] << " " << scores[1] << " " << scores[2] << endl;

								if (scores[0] >= scores[1])
								{
									if (scores[0] >= scores[2])
										maxIdx = 0;
									else
										maxIdx = 2;
								}
								else
								{
									if (scores[1] >= scores[2])
										maxIdx = 1;
									else
										maxIdx = 2;
								}

								ffTable[u][stateIdx][specIdx] = scores[maxIdx];
							}
						}

						// calculate score from DP tables
						double score = (ffTable[u][stateList.size()][specSize] * (specSize - u))/(stateList.size() * specSize);
//						cout << "u=" << u << " " << score << endl;
						if (score > bestGroupScore)
							bestGroupScore = score;
					}
					DELETE_3D_ARRAY(ffTable, specSize, stateList.size() + 1);
				}

//				cout << "Score: " << bestGroupScore << endl;
				res += bestGroupScore;

				// memorize score in hash table
				for (vector<bool *>::iterator it = stateList.begin(); it != stateList.end(); ++it)
				{
					bool * state = new bool[ind.varNames.size()];
					memcpy(state,*it,sizeof(bool)*ind.varNames.size());
					savedScores[state] = bestGroupScore;
				}
			}
		}
		else
		// this is a chain rule
		{
			double bestGroupScore = 0.0;

			// iterate over alternatives and choose score for best alternative
			for (unsigned int r = 0; r < numGroups; ++r)
			{
				unsigned int specSize = groupIndices[r+1] - groupIndices[r];

				// initialize Dynamic Programming matrix
				double ** scoreTable;
				ALLOC_2D_ARRAY(scoreTable, double, stateList.size() + 1, specSize + 1);

				for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
				{
					scoreTable[0][specIdx] = -1e10;
				}

				for (unsigned int stateIdx = 0; stateIdx <= stateList.size(); ++stateIdx)
				{
					scoreTable[stateIdx][0] = 0;
				}

				// main DP loops
				for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
				{
					for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
					{
						double fractionFulfilled = stateSpecifications[groupIndices[r] + specIdx - 1]->getFractionSatisfied(stateList[stateIdx - 1]);

						double scores[2];

						scores[0] = fractionFulfilled + scoreTable[stateIdx - 1][specIdx - 1];
						scores[1] = scoreTable[stateIdx - 1][specIdx];

						scoreTable[stateIdx][specIdx] = max(scores[0],scores[1]);
					}
				}

				double score = scoreTable[stateList.size()][specSize]/specSize;

				if (score > bestGroupScore)
				// this is the currently best alternative
					bestGroupScore = score;

				DELETE_2D_ARRAY(scoreTable, stateList.size() + 1);
			}
			res += bestGroupScore;
		}

		for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
				it != states.end(); ++it)
		{
			delete [] it->first;
			delete [] it->second;
		}

		++i;
	}

	delete [] state;

	for (unordered_map<bool*, double, HashStates, CompareStatesEqual>::iterator it = savedScores.begin();
					it != savedScores.end(); ++it)
	{
		delete [] it->first;
	}
//	cout << "Fitness: " << res / checkedStates << endl;
	return res / checkedStates * importance;
}


ConstraintViolation * NetworkConstraint::getConstraintViolationsGreedy(Individual &ind,
														  unsigned int maxStartStates,
														  unsigned int calculatedTransitions,
														  bool useRandomStates) const
{
	unsigned int numGroups = (unsigned int)groupIndices.size() - 1;
	ConstraintViolation * res = new ConstraintViolation((NetworkConstraint*)this, (unsigned int)ind.varNames.size(), true);

	bool * state =  new bool[ind.varNames.size()];
	bool * initialState =  new bool[ind.varNames.size()];

	unsigned int possibleStates = ((unsigned int)ind.varNames.size() - initialValues->getNumConditions());
	if (possibleStates > 31)
		possibleStates = ~0;
	else
		possibleStates = (1 << possibleStates);

	if (maxStartStates >= possibleStates)
	{
		useRandomStates = false;
		maxStartStates = possibleStates;
	}

	unsigned int checkedStates = maxStartStates;

	for (unsigned int i = 0; i < maxStartStates; )
	{
		//cout << i << endl;
		if (useRandomStates)
			initialValues->generateMatchingState(initialState, (unsigned int)ind.varNames.size());
		else
			initialValues->generateStateNo(initialState, (unsigned int)ind.varNames.size(), i);

		memcpy(state, initialState, sizeof(bool) * ind.varNames.size());

		unsigned int preconditionsVerified = 0;
		if (preconditions.size() > 0)
		{
			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
				if (preconditions[preconditionsVerified]->isSatisfied(state))
				{
					++preconditionsVerified;
					if (preconditionsVerified == preconditions.size())
						break;
				}
				ind.stateTransition(state, fixedGenes);
			}
			if (preconditionsVerified != preconditions.size())
			{
				if (!useRandomStates)
				{
					++i;
					--checkedStates;
				}
				continue;
			}
		}

		if (isAttractor)
		{
			HashStates h((unsigned int)ind.varNames.size());
			CompareStatesEqual s((unsigned int)ind.varNames.size());
			unordered_map<bool*, bool*, HashStates, CompareStatesEqual> states(8,h,s,std::allocator<std::pair<bool*, bool*> >());

			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
//				cout << "k=" << k  << endl;
				if (states.count(state) == 0)
				{
					bool * startState = new bool[ind.varNames.size()];
					memcpy(startState, state, sizeof(bool) * ind.varNames.size());

					bool * nextState = new bool[ind.varNames.size()];
					ind.stateTransition(state, fixedGenes);
					memcpy(nextState, state, sizeof(bool) * ind.varNames.size());

					states[startState] = nextState;
				}
				else
				{
//					cout << "Found attractor after " << k << " transitions."  << endl;
					vector<bool * > attractor;
					bool * startState = new bool[ind.varNames.size()];
					memcpy(startState, state, sizeof(bool) * ind.varNames.size());
					bool * curState = state;
					do
					{
						attractor.push_back(curState);
						curState = states[curState];
					}
					while (memcmp(curState,startState,sizeof(bool) * ind.varNames.size()) != 0);

					delete [] startState;

					cout << "Attractor size: " << attractor.size() << endl;
					for (unsigned int l = 0; l < attractor.size(); ++l)
					{
						for (unsigned int p = 0; p < ind.varNames.size(); ++p)
						{
							cout << ind.varNames[p] << "=" << attractor[l][p] << " ";
						}
						cout << endl;
					}
					double bestGroupScore = 0.0;
					unsigned int bestGroup = 0;
					unsigned int * bestGroupAssignment = new unsigned int[attractor.size()];

					for (unsigned int r = 0; r <  numGroups; ++r)
					{
						double maxScore = 0.0;
						unsigned int * stateAssignment = new unsigned int[attractor.size()];
//						unsigned int conditionsVerified = 0;
						unsigned int specSize = groupIndices[r+1] - groupIndices[r];
//						cout << "specSize " << specSize << endl;
						for (unsigned int l = 0; l < specSize; ++l)
						{

							unsigned int * tempAssignment = new unsigned int[attractor.size()];
							int currentCondition = -1;
							unsigned int stateNo = 0;
							double score = 0.0;
							unsigned int maxTrans = (unsigned int)attractor.size();//max((size_t)specSize,attractor.size());
//							cout << "Start" << endl;
							for (stateNo = 0; stateNo < maxTrans; ++stateNo)

							{
								double subscore1 = 0, subscore2;
//								cout << "CC " << currentCondition+1 << endl;
	//							cout << "Curr: " << currentCondition << "/" << rules.size() << endl;
//								cout << "StateNo: " << stateNo << "/" << attractor.size() << endl;
								if (currentCondition != -1)
									subscore1 = stateSpecifications[groupIndices[r] +  (currentCondition % specSize)]->getFractionSatisfied(attractor[(stateNo+l) % attractor.size()]);

								subscore2 = stateSpecifications[groupIndices[r] + ((currentCondition + 1) % specSize)]->getFractionSatisfied(attractor[(stateNo+l) % attractor.size()]);
//								cout << "S " << subscore2 << " " << subscore1 << endl;

								if (subscore2 >= subscore1)
								{
	//								cout << (currentCondition+1) % rules.size() << " " << (stateNo+l) % attractor.size() << endl;
									score += subscore2;
									tempAssignment[(stateNo+l) % attractor.size()] = ((currentCondition + 1) % specSize);
									currentCondition = currentCondition + 1;
								}
								else
								{
	//								cout << currentCondition % rules.size() << " " << (stateNo+l) % attractor.size() << endl;
									score += subscore1;
									tempAssignment[(stateNo+l) % attractor.size()] = currentCondition % specSize;
								}
	//							cout << "Score " << score << endl;
	//							++stateNo;
	//							cout << stateNo << endl;
							}
							score /= maxTrans;

							if ((int)specSize > currentCondition + 1)
								score *= (currentCondition + 1)/(double)specSize;

							if (score > maxScore)
							{
								maxScore = score;
//								conditionsVerified = currentCondition + 1;
								memcpy(stateAssignment, tempAssignment, sizeof(unsigned int) * attractor.size());
//								for (unsigned int i = 0; i < attractor.size(); ++i)
//									stateAssignment[i].assign(tempAssignment[i].begin(), tempAssignment[i].end());
							}
							delete [] tempAssignment;
						}
//						if (ruleSize > conditionsVerified)
//						{
//							maxScore *= conditionsVerified/(double)ruleSize;
//						};

						if (bestGroupScore < maxScore)
						{
							bestGroupScore = maxScore;
							bestGroup = r;
							memcpy(bestGroupAssignment, stateAssignment, sizeof(unsigned int) * attractor.size());
//							for (unsigned int i = 0; i < attractor.size(); ++i)
//								bestGroupAssignment[i].assign(stateAssignment[i].begin(), stateAssignment[i].end());
						}
						delete [] stateAssignment;
					}
//					cout << "Score " << bestGroupScore << endl;
//					cout << "Matching:\n";
					for (unsigned int v = 0; v < attractor.size(); ++v)
					{
//						cout << v << " => " << bestGroupAssignment[v] << endl;
						stateSpecifications[groupIndices[bestGroup] + bestGroupAssignment[v]]->getViolations(attractor[v], initialState, *res, 0);
					}
					delete [] bestGroupAssignment;
					break;
				}
			}
			for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
				it != states.end(); ++it)
			{
				delete [] it->first;
				delete [] it->second;
			}
		}
		else
		{
			bool * startState = new bool[ind.varNames.size()];
			memcpy(startState, state, sizeof(bool) * ind.varNames.size());
			double bestGroupScore = 0.0;
			unsigned int bestGroup = 0;
			unsigned int bestConditionVerified = 0;
			bool * bestGroupErrorState = new bool[ind.varNames.size()];
			bool exitLoop = false;

			for (unsigned int r = 0; r < numGroups; ++r)
			{
				if (exitLoop)
					break;

				memcpy(state, startState, sizeof(bool) * ind.varNames.size());

				double maxScore = 0.0;
				bool * maxState = new bool[ind.varNames.size()];
				unsigned int conditionVerified = 0;
				unsigned int ruleSize = groupIndices[r+1] - groupIndices[r];

				for (unsigned int k = 0; k < calculatedTransitions; ++k)
				{
					double score = stateSpecifications[groupIndices[r] + conditionVerified]->getFractionSatisfied(state);
					if (score > maxScore)
					{
						maxScore = score;
						memcpy(maxState, state, sizeof(bool) * ind.varNames.size());
					}
					if (score == 1.0)
					{
						++conditionVerified;
						maxScore = 0.0;
					}
					if (conditionVerified==ruleSize)
					{
						exitLoop = true;
						break;
					}
					ind.stateTransition(state, fixedGenes);
				}

				double ruleScore = (conditionVerified + maxScore)/((double)ruleSize);

				if (ruleScore > bestGroupScore)
				{
					bestGroupScore = ruleScore;
					bestGroup = r;
					bestConditionVerified = conditionVerified;
					memcpy(bestGroupErrorState, maxState, sizeof(bool) * ind.varNames.size());
				}
				delete [] maxState;
			}
			if (bestConditionVerified != (groupIndices[bestGroup+1] - groupIndices[bestGroup]))
			{
				stateSpecifications[groupIndices[bestGroup] + bestConditionVerified]->getViolations(bestGroupErrorState, initialState, *res, 0);
			}
			delete [] startState;
			delete [] bestGroupErrorState;
		}
		++i;
	}
	delete [] state;
	delete [] initialState;
	return res;
}

ConstraintViolation * NetworkConstraint::getConstraintViolationsDP(Individual &ind,
														  unsigned int maxStartStates,
														  unsigned int calculatedTransitions,
														  bool useRandomStates) const
{

	unsigned int numGroups = (unsigned int)groupIndices.size() - 1;
	ConstraintViolation * res = new ConstraintViolation((NetworkConstraint*)this,
														(unsigned int)ind.varNames.size(),
														!isAttractor);

	bool * state = new bool[ind.varNames.size()];
	bool * initialState = new bool[ind.varNames.size()];

	unsigned int possibleStates = ((unsigned int)ind.varNames.size() - initialValues->getNumConditions());
	if (possibleStates > 31)
		possibleStates = ~0;
	else
		possibleStates = (1 << possibleStates);

	if (maxStartStates >= possibleStates)
	{
		useRandomStates = false;
		maxStartStates = possibleStates;
	}

	HashStates h1((unsigned int)ind.varNames.size());
	CompareStatesEqual s1((unsigned int)ind.varNames.size());
	unordered_map<bool*, unsigned int, HashStates, CompareStatesEqual> savedAssignments(8,h1,s1,std::allocator<std::pair<bool*, unsigned int> >());

	unsigned int maxSpecSize = 0;
	for (unsigned int r = 0; r < numGroups; ++r)
	{
		unsigned int specSize = groupIndices[r+1] - groupIndices[r];
		if (specSize > maxSpecSize)
			maxSpecSize = specSize;
	}

	unsigned int checkedStates = maxStartStates;

	for (unsigned int i = 0; i < maxStartStates; )
	{
//		cout << i << endl;

		// generate a start state
		if (useRandomStates)
			initialValues->generateMatchingState(state, (unsigned int)ind.varNames.size());
		else
			initialValues->generateStateNo(state, (unsigned int)ind.varNames.size(), i);

		memcpy(initialState, state, sizeof(bool) * ind.varNames.size());

		// verify preconditions by repeated state calculations
		unsigned int preconditionsVerified = 0;
		if (preconditions.size() > 0)
		{
			for (unsigned int k = 0; k < calculatedTransitions; ++k)
			{
				if (preconditions[preconditionsVerified]->isSatisfied(state))
				{
					++preconditionsVerified;
					if (preconditionsVerified == preconditions.size())
						break;
				}
				ind.stateTransition(state, fixedGenes);
			}
			if (preconditionsVerified != preconditions.size())
			{
				if (!useRandomStates)
				{
					++i;
					--checkedStates;
				}
				continue;
			}
		}

		// save sequence of states or attractor in a hash table
		HashStates h2((unsigned int)ind.varNames.size());
		CompareStatesEqual s2((unsigned int)ind.varNames.size());
		unordered_map<bool*, bool*, HashStates, CompareStatesEqual> states(8,h2,s2,std::allocator<std::pair<bool*, bool*> >());
		vector<bool * > stateList;

		unsigned int k;
		for (k = 0; k < calculatedTransitions; ++k)
		{
//				cout << "k=" << k  << endl;
			if (states.count(state) == 0)
			// this state has not been reached previously
			{
				bool * startState = new bool[ind.varNames.size()];
				memcpy(startState, state, sizeof(bool) * ind.varNames.size());

				bool * nextState = new bool[ind.varNames.size()];
				ind.stateTransition(state, fixedGenes);
				memcpy(nextState, state, sizeof(bool) * ind.varNames.size());

				states[startState] = nextState;

				if (!isAttractor && states.count(nextState) == 0)
				// if this is a chain, add the state to the state list
				{
					stateList.push_back(nextState);
				}
			}
			else
			if (!isAttractor && stateList.size() < maxSpecSize)
			// if the state list is too short, add the attractor states as well
			{
				bool * nextState = states[state];
				stateList.push_back(nextState);
			}
			else
				break;
		}

		if (k == calculatedTransitions)
		// The number of transitions was not sufficient to reach the attractor
		{
			warnings.insert("The maximum number of transitions is not sufficient to reach the attractor! It is recommended to increase \"-mt\"!");
			for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
					it != states.end(); ++it)
			{
				delete [] it->first;
				delete [] it->second;
			}
			continue;
		}

		if (isAttractor)
		// this is an attractor rule
		{

			unsigned int matchIdx;
			if (savedAssignments.count(state) == 0)
				matchIdx = res->getNumMatchings();
			else
				matchIdx = savedAssignments[state];

			// determine the attractor
			bool * startState = new bool[ind.varNames.size()];
			memcpy(startState, state, sizeof(bool) * ind.varNames.size());
			bool * curState = state;
			do
			{
				stateList.push_back(curState);
				curState = states[curState];
			}
			while (memcmp(curState,startState,sizeof(bool) * ind.varNames.size()) != 0);

			delete [] startState;

//			cout << "Attractor size: " << stateList.size() << endl;
//			for (unsigned int l = 0; l < stateList.size(); ++l)
//			{
//				for (unsigned int p = 0; p < ind.varNames.size(); ++p)
//				{
//					cout << ind.varNames[p] << "=" << stateList[l][p] << " ";
//				}
//				cout << endl;
//			}

			// calculate the scores for all alternatives and choose the best one
			double bestGroupScore = 0.0;
			unsigned int bestGroup = 0;
			// remember the state assignment in a list
			unsigned int * bestGroupAssignment = new unsigned int[stateList.size()];

			for (unsigned int r = 0; r < numGroups; ++r)
			{
				unsigned int specSize = groupIndices[r+1] - groupIndices[r];

				for (unsigned int l = 0; l < stateList.size(); ++l)
				// iterate over attractor state offsets
				{
					double *** ffTable;
					unsigned char *** choiceTable;
					ALLOC_3D_ARRAY(ffTable, double, specSize, stateList.size() + 1, specSize + 1);
					ALLOC_3D_ARRAY(choiceTable, unsigned char, specSize, stateList.size() + 1, specSize + 1);

					for (unsigned int u = 0; u < specSize; ++u)
					// iterate over the maximum number of uncovered specifications
					{
						// initialize Dynamic Programming matrix

						ffTable[u][0][0] = 0;
						choiceTable[u][0][0] = 4;

						for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
						{
							ffTable[u][stateIdx][0] = -1e10;
						}

						for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
						{
							if (specIdx <= u)
								ffTable[u][0][specIdx] = 0;
							else
								ffTable[u][0][specIdx] = -1e10;

							choiceTable[u][0][specIdx] = 2;
						}

						// main DP loops
						for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
						{
							for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
							{
								double fractionFulfilled = stateSpecifications[groupIndices[r] + specIdx - 1]->getFractionSatisfied(stateList[(stateIdx - 1 + l) % stateList.size()]);

								double scores[3];
								unsigned int maxIdx = 0;

//									cout << "FF " << ffTable[stateIdx][specIdx - 1] << " " << coveredSpecTable[stateIdx][specIdx-1] << endl;
								scores[0] = (fractionFulfilled + ffTable[u][stateIdx - 1][specIdx - 1]);
								scores[1] = (fractionFulfilled + ffTable[u][stateIdx - 1][specIdx]);

								if (u > 0)
									scores[2] = ffTable[u-1][stateIdx][specIdx - 1];
								else
									scores[2] = -1e10;

//									cout << stateIdx << " " << specIdx << ": " << scores[0] << " " << scores[1] << " " << scores[2] << endl;

								if (scores[0] >= scores[1])
								{
									if (scores[0] >= scores[2])
										maxIdx = 0;
									else
										maxIdx = 2;
								}
								else
								{
									if (scores[1] >= scores[2])
										maxIdx = 1;
									else
										maxIdx = 2;
								}

								// also remember which choice was taken (i.e. if a specification entry was assigned
								// to the current state, and which
								choiceTable[u][stateIdx][specIdx] = maxIdx;
								ffTable[u][stateIdx][specIdx] = scores[maxIdx];
								//cout << "Res: " << ffTable[stateIdx][specIdx] << " " << coveredSpecTable[stateIdx][specIdx] << endl;
							}
						}

						// calculate score from DP tables
						double score = (ffTable[u][stateList.size()][specSize] * (specSize - u))/(stateList.size() * specSize);

//							cout << specSize << " " << attractor.size() << " " << score << " " << coveredSpecTable[attractor.size()][specSize] << " " <<  ffTable[attractor.size()][specSize] << endl;

						if (score > bestGroupScore)
						{
							bestGroup = r;
							bestGroupScore = score;
//									cout << "Best score: " << bestGroupScore << endl;

							int stateIdx = (int)stateList.size();
							int specIdx = specSize;
							int uIdx = u;

							// backtrace the choice table to determine
							// the assignment of states to specification entries
							while(true)
							{
								if (choiceTable[uIdx][stateIdx][specIdx] == 4)
								    break;
								if (choiceTable[uIdx][stateIdx][specIdx] == 0)
								{
									bestGroupAssignment[(stateIdx-1 + l) % stateList.size()] = (specIdx-1);
									--stateIdx;
									--specIdx;
								}
								else
								if (choiceTable[uIdx][stateIdx][specIdx] == 1)
								{
									bestGroupAssignment[(stateIdx-1 + l) % stateList.size()] = (specIdx-1);
									--stateIdx;
								}
								else
								{
								  --specIdx;
								  --uIdx;
								}

							 }
						}

					}
					DELETE_3D_ARRAY(ffTable, specSize, stateList.size() + 1);
					DELETE_3D_ARRAY(choiceTable, specSize, stateList.size() + 1);
				}
//				cout << "Score " << bestGroupScore << endl;
//				cout << "Matching:\n";

			}

			if (savedAssignments.count(stateList[0]) == 0)
			{
				res->addAttractorMatching(stateList, bestGroupAssignment, bestGroup);
				// memorize score in hash table
				for (vector<bool *>::iterator it = stateList.begin(); it != stateList.end(); ++it)
				{
					bool * state = new bool[ind.varNames.size()];
					memcpy(state,*it,sizeof(bool)*ind.varNames.size());
					savedAssignments[state] = matchIdx;
				}
			}

			//calculate mismatches of best matching
			for (unsigned int v = 0; v < stateList.size(); ++v)
			{
//					for (unsigned int p = 0; p < ind.varNames.size(); ++p)
//					{
//						cout << stateList[v][p] << " ";
//					}
//					cout << "=> ";
//					for (std::tr1::unordered_map<unsigned int, bool >::iterator it = stateSpecifications[groupIndices[r] + bestGroupAssignment[v]]->geneValues.begin();
//							it != stateSpecifications[groupIndices[r] + bestGroupAssignment[v]]->geneValues.end(); ++it)
//					{
//						cout << it->second << " ";
//						//cout << stateList[bestGroupAssignment[v]][p] << " ";
//					}
//					cout << endl;
				stateSpecifications[groupIndices[bestGroup] + bestGroupAssignment[v]]->getViolations(stateList[v],
																									 initialState, *res, matchIdx);
			}

			delete [] bestGroupAssignment;
		}
		else
		// this is a chain rule
		{
//			cout << "Chain size: " << stateList.size() << endl;
//			for (unsigned int l = 0; l < stateList.size(); ++l)
//			{
//				for (unsigned int p = 0; p < ind.varNames.size(); ++p)
//				{
//					cout << ind.varNames[p] << "=" << stateList[l][p] << " ";
//				}
//				cout << endl;
//			}

			// remember the score and the assignment of the states to specification entries
			// for the best alternative
			double bestGroupScore = 0.0;
			unsigned int bestGroup = 0;
			int * bestGroupAssignment = new int[stateList.size()];
			for (unsigned int v = 0; v < stateList.size(); ++v)
				bestGroupAssignment[v] = -1;

			// iterate over alternatives and choose score for best alternative
			for (unsigned int r = 0; r < numGroups; ++r)
			{
				unsigned int specSize = groupIndices[r+1] - groupIndices[r];
//				cout << "SpecSize " << specSize << endl;

				// initialize Dynamic Programming matrix
				double ** scoreTable;
				ALLOC_2D_ARRAY(scoreTable, double, stateList.size() + 1, specSize + 1);
				unsigned int ** choiceTable;
				ALLOC_2D_ARRAY(choiceTable, unsigned int, stateList.size() + 1, specSize + 1);

				for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
				{
					scoreTable[0][specIdx] = -1e10;
					choiceTable[0][specIdx] = 2;
				}

				for (unsigned int stateIdx = 0; stateIdx <= stateList.size(); ++stateIdx)
				{
					scoreTable[stateIdx][0] = 0;
					choiceTable[stateIdx][0] = 1;
				}

				// main DP loops
				for (unsigned int stateIdx = 1; stateIdx <= stateList.size(); ++stateIdx)
				{
					for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
					{
						double fractionFulfilled = stateSpecifications[groupIndices[r] + specIdx - 1]->getFractionSatisfied(stateList[stateIdx - 1]);

						double scores[2];

						scores[0] = fractionFulfilled + scoreTable[stateIdx - 1][specIdx - 1];
						scores[1] = scoreTable[stateIdx - 1][specIdx];

						if (scores[0] >= scores[1])
						{
							scoreTable[stateIdx][specIdx] = scores[0];
							choiceTable[stateIdx][specIdx] = 0;
						}
						else
						{
							scoreTable[stateIdx][specIdx] = scores[1];
							choiceTable[stateIdx][specIdx] = 1;
						}
					}
				}

				double score = scoreTable[stateList.size()][specSize]/specSize;
//				cout << "ChainScore " << score << endl;

				if (score > bestGroupScore)
				// this is the currently best alternative
				{
					bestGroup = r;
					bestGroupScore = score;
					int stateIdx = (int)stateList.size();
					int specIdx = specSize;

					// backtrace the choice table to determine
					// the assignment of states to specification entries
					while(stateIdx > 0)
					{
						if (choiceTable[stateIdx][specIdx] == 0)
						{
							bestGroupAssignment[stateIdx - 1] = (specIdx - 1);
							--stateIdx;
							--specIdx;
						}
						else
						if (choiceTable[stateIdx][specIdx] == 0)
						{
							--stateIdx;
						}
						else
							throw string("Unknown choice: " + std::to_string(choiceTable[stateIdx][specIdx]));
					 }
				}
				DELETE_2D_ARRAY(scoreTable, stateList.size() + 1);
				DELETE_2D_ARRAY(choiceTable, stateList.size() + 1);
			}
//			cout << "Score " << bestGroupScore << endl;
//			cout << "Matching:\n";

			// determine the mismatches of the optimal assignment
			for (unsigned int v = 0; v < stateList.size(); ++v)
			{
				if (bestGroupAssignment[v] != -1)
				{
//					cout << v << " => " << bestGroupAssignment[v] << endl;
					stateSpecifications[groupIndices[bestGroup] + bestGroupAssignment[v]]->getViolations(stateList[v], initialState, *res, 0);
				}
			}

			delete [] bestGroupAssignment;
		}

		for (unordered_map<bool*, bool*, HashStates, CompareStatesEqual>::iterator it = states.begin();
				it != states.end(); ++it)
		{
			delete [] it->first;
			delete [] it->second;
		}

		++i;
	}
	delete [] state;
	delete [] initialState;

	for (unordered_map<bool*, unsigned int, HashStates, CompareStatesEqual>::iterator it = savedAssignments.begin();
						it != savedAssignments.end(); ++it)
	{
		delete [] it->first;
	}

	return res;


//	unsigned int numGroups = groupIndices.size() - 1;
//	ConstraintViolation * res = new ConstraintViolation(ind.varNames.size());
//
//	bool state[ind.varNames.size()];
//	bool initialState[ind.varNames.size()];
//
//	unsigned int possibleStates = (ind.varNames.size() - initialValues->getNumConditions());
//	if (possibleStates > 31)
//		possibleStates = ~0;
//	else
//		possibleStates = (1 << possibleStates);
//
//	if (maxStartStates > possibleStates)
//	{
//		useRandomStates = false;
//		maxStartStates = possibleStates;
//	}
//
//	unsigned int checkedStates = maxStartStates;
//
//	for (unsigned int i = 0; i < maxStartStates; )
//	{
//		//cout << i << endl;
//		if (useRandomStates)
//			initialValues->generateMatchingState(initialState, ind.varNames.size());
//		else
//			initialValues->generateStateNo(initialState, ind.varNames.size(), i);
//
//		memcpy(state, initialState, sizeof(bool) * ind.varNames.size());
//
//		unsigned int preconditionsVerified = 0;
//		if (preconditions.size() > 0)
//		{
//			for (unsigned int k = 0; k < calculatedTransitions; ++k)
//			{
//				if (preconditions[preconditionsVerified]->isSatisfied(state))
//				{
//					++preconditionsVerified;
//					if (preconditionsVerified == preconditions.size())
//						break;
//				}
//				ind.stateTransition(state);
//			}
//			if (preconditionsVerified != preconditions.size())
//			{
//				if (!useRandomStates)
//				{
//					++i;
//					--checkedStates;
//				}
//				continue;
//			}
//		}
//
//		if (isAttractor)
//		{
//			HashStates h(ind.varNames.size());
//			CompareStates s(ind.varNames.size());
//			::unordered_map<bool*, bool*, HashStates, CompareStates> states(8,h,s,std::allocator<std::pair<bool*, bool*> >());
//
//			for (unsigned int k = 0; k < calculatedTransitions; ++k)
//			{
////				cout << "k=" << k  << endl;
//				if (states.count(state) == 0)
//				{
//					bool * startState = new bool[ind.varNames.size()];
//					memcpy(startState, state, sizeof(bool) * ind.varNames.size());
//
//					bool * nextState = new bool[ind.varNames.size()];
//					ind.stateTransition(state);
//					memcpy(nextState, state, sizeof(bool) * ind.varNames.size());
//
//					states[startState] = nextState;
//				}
//				else
//				{
////					cout << "Found attractor after " << k << " transitions."  << endl;
//					vector<bool * > attractor;
//					bool startState[ind.varNames.size()];
//					memcpy(startState, state, sizeof(bool) * ind.varNames.size());
//					bool * curState = state;
//					do
//					{
//						attractor.push_back(curState);
//						curState = states[curState];
//					}
//					while (memcmp(curState,startState,sizeof(bool) * ind.varNames.size()) != 0);
//
////					cout << "Attractor size: " << attractor.size() << endl;
////					for (unsigned int l = 0; l < attractor.size(); ++l)
////					{
////						for (unsigned int p = 0; p < ind.varNames.size(); ++p)
////						{
////							cout << ind.varNames[p] << "=" << attractor[l][p] << " ";
////						}
////						cout << endl;
////					}
//
//					double bestGroupScore = -1;
//					unsigned int bestGroup = 0;
//					unsigned int bestGroupAssignment[attractor.size()];
//
//					for (unsigned int r = 0; r < numGroups; ++r)
//					{
//						unsigned int specSize = groupIndices[r+1] - groupIndices[r];
//
//						int maxCovered = min((unsigned int)attractor.size(), specSize);
//
//						for (unsigned int l = 0; l < attractor.size(); ++l)
//						{
//							for (int c = 1; c <= maxCovered; ++c)
//							{
//	//							cout << "l = " << l << endl;
//								double ffTable[attractor.size() + 1][specSize + 1];
//								int coveredSpecTable[attractor.size() + 1][specSize + 1];
//								unsigned int choiceTable[attractor.size() + 1][specSize + 1];
//
//								ffTable[0][0] = 0;
//								coveredSpecTable[0][0] = specSize;
//
//								for (unsigned int stateIdx = 1; stateIdx <= attractor.size(); ++stateIdx)
//								{
//									ffTable[stateIdx][0] = -1e10;
//									coveredSpecTable[stateIdx][0] = specSize;
//								}
//
//								for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
//								{
//									ffTable[0][specIdx] = -1e10;
//									coveredSpecTable[0][specIdx] = specSize - specIdx;
//								}
//
//								for (unsigned int stateIdx = 1; stateIdx <= attractor.size(); ++stateIdx)
//								{
//									for (unsigned int specIdx = 1; specIdx <= specSize; ++specIdx)
//									{
//										double fractionFulfilled = stateSpecifications[groupIndices[r] + specIdx - 1]->getFractionSatisfied(attractor[(stateIdx - 1 + l) % attractor.size()]);
//
//										double scores[3];
//										unsigned int maxIdx = 0;
//
//	//									cout << "FF " << ffTable[stateIdx][specIdx - 1] << " " << coveredSpecTable[stateIdx][specIdx-1] << endl;
//										scores[0] = (fractionFulfilled + ffTable[stateIdx - 1][specIdx]) * c;
//										scores[1] = (fractionFulfilled + ffTable[stateIdx - 1][specIdx - 1]) * c;
//										scores[2] = (ffTable[stateIdx][specIdx - 1]) * c;
//
//	//									cout << stateIdx << " " << specIdx << ": " << scores[0] << " " << scores[1] << " " << scores[2] << endl;
//
//										if (scores[0] > scores[1])
//										{
//											if (scores[0] > scores[2] || coveredSpecTable[stateIdx][specIdx-1] == c)
//												maxIdx = 0;
//											else
//												maxIdx = 2;
//										}
//										else
//										{
//											if (scores[1] > scores[2] || coveredSpecTable[stateIdx][specIdx-1] == c)
//												maxIdx = 1;
//											else
//												maxIdx = 2;
//										}
//
//										choiceTable[stateIdx][specIdx] = maxIdx;
//										switch(maxIdx)
//										{
//											case 0:
//												ffTable[stateIdx][specIdx] = fractionFulfilled + ffTable[stateIdx - 1][specIdx];
//												coveredSpecTable[stateIdx][specIdx] = coveredSpecTable[stateIdx-1][specIdx];
//												break;
//											case 1:
//												ffTable[stateIdx][specIdx] = fractionFulfilled + ffTable[stateIdx - 1][specIdx - 1];
//												coveredSpecTable[stateIdx][specIdx] = coveredSpecTable[stateIdx-1][specIdx-1];
//												break;
//											case 2:
//												ffTable[stateIdx][specIdx] = ffTable[stateIdx][specIdx - 1];
//												coveredSpecTable[stateIdx][specIdx] = max(0,coveredSpecTable[stateIdx][specIdx-1] - 1);
//												break;
//										}
//										//cout << "Res: " << ffTable[stateIdx][specIdx] << " " << coveredSpecTable[stateIdx][specIdx] << endl;
//									}
//								}
//								double score = (ffTable[attractor.size()][specSize] * coveredSpecTable[attractor.size()][specSize])/(attractor.size() * specSize);
//
//	//							cout << specSize << " " << attractor.size() << " " << score << " " << coveredSpecTable[attractor.size()][specSize] << " " <<  ffTable[attractor.size()][specSize] << endl;
//
//								if (score > bestGroupScore)
//								{
//									bestGroupScore = score;
////									cout << "Best score: " << bestGroupScore << endl;
//
//									int stateIdx = attractor.size();
//									int specIdx = specSize;
//									while(stateIdx > 0)
//									{
//										if (choiceTable[stateIdx][specIdx] == 0)
//										{
//											bestGroupAssignment[(stateIdx-1 + l) % attractor.size()] = (specIdx-1);
//											--stateIdx;
//										}
//										else
//										if (choiceTable[stateIdx][specIdx] == 1)
//										{
//											bestGroupAssignment[(stateIdx-1 + l) % attractor.size()] = (specIdx-1);
//											--stateIdx;
//											--specIdx;
//										}
//										else
//										{
//										  --specIdx;
//										}
//
//									 }
//								}
//								if (coveredSpecTable[attractor.size()][specSize] > c)
//								{
//									c = coveredSpecTable[attractor.size()][specSize];
//								}
//							}
//						}
//					}
////					cout << "Score " << bestGroupScore << endl;
////					cout << "Matching:\n";
//					for (unsigned int v = 0; v < attractor.size(); ++v)
//					{
////						cout << v << " => " << bestGroupAssignment[v] << endl;
//						stateSpecifications[groupIndices[bestGroup] + bestGroupAssignment[v]]->getViolations(attractor[v], initialState, *res);
//					}
//					break;
//				}
//			}
//			for (tr1::unordered_map<bool*, bool*, HashStates, CompareStates>::iterator it = states.begin();
//				it != states.end(); ++it)
//			{
//				delete [] it->first;
//				delete [] it->second;
//			}
//		}
//		else
//		{
//			bool startState[ind.varNames.size()];
//			memcpy(startState, state, sizeof(bool) * ind.varNames.size());
//			double bestGroupScore = 0.0;
//			unsigned int bestGroup = 0;
//			unsigned int bestConditionVerified = 0;
//			bool * bestGroupErrorState = new bool[ind.varNames.size()];
//			bool exitLoop = false;
//
//			for (unsigned int r = 0; r < numGroups; ++r)
//			{
//				if (exitLoop)
//					break;
//
//				memcpy(state, startState, sizeof(bool) * ind.varNames.size());
//
//				double maxScore = 0.0;
//				bool maxState[ind.varNames.size()];
//				unsigned int conditionVerified = 0;
//				unsigned int ruleSize = groupIndices[r+1] - groupIndices[r];
//
//				for (unsigned int k = 0; k < calculatedTransitions; ++k)
//				{
//					double score = stateSpecifications[groupIndices[r] + conditionVerified]->getFractionSatisfied(state);
//					if (score > maxScore)
//					{
//						maxScore = score;
//						memcpy(maxState, state, sizeof(bool) * ind.varNames.size());
//					}
//					if (score == 1.0)
//					{
//						++conditionVerified;
//						maxScore = 0.0;
//					}
//					if (conditionVerified==ruleSize)
//					{
//						exitLoop = true;
//						break;
//					}
//					ind.stateTransition(state);
//				}
//
//				double ruleScore = (conditionVerified + maxScore)/((double)ruleSize);
//
//				if (ruleScore > bestGroupScore)
//				{
//					bestGroupScore = ruleScore;
//					bestGroup = r;
//					bestConditionVerified = conditionVerified;
//					memcpy(bestGroupErrorState, maxState, sizeof(bool) * ind.varNames.size());
//				}
//			}
//			if (bestConditionVerified != (groupIndices[bestGroup+1] - groupIndices[bestGroup]))
//			{
//				stateSpecifications[groupIndices[bestGroup] + bestConditionVerified]->getViolations(bestGroupErrorState, initialState, *res);
//			}
//			delete bestGroupErrorState;
//		}
//		++i;
//	}
//	return res;
}


NetworkConstraint::~NetworkConstraint()
{
	for (vector<SimpleCondition *>::iterator it = preconditions.begin(); it != preconditions.end(); ++it)
	{
		delete *it;
	}

	for (vector<SimpleCondition*>::iterator it1 = stateSpecifications.begin(); it1 != stateSpecifications.end(); ++it1)
	{
		delete *it1;
	}

	delete initialValues;
	delete fixedGenes;
}
