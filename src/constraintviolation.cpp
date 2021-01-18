/*
 * constraintviolation.cpp
 * Provides classes to monitor
 * violations of rule sets for
 * manual analysis
 *
 *  Created on: Nov 25, 2010
 *      Author: muessel
 * 	Edited 2020 by schwab
 */
#include "helpers.h"
#include "constraintviolation.h"
#include <cstring>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <iostream>

using namespace std;

unsigned int ConditionViolation::getAlternative() const
{
	return condition->getAlternative();
}

unsigned int ConditionViolation::getConditionIndex() const
{
	return conditionIndex;
}

unsigned int ConditionViolation::getSpecificationEntry() const
{
	return specificationEntry;
}

SimpleCondition * ConditionViolation::getCondition() const
{
	return condition;
}


ConstraintViolation::ConstraintViolation(NetworkConstraint * constraint, unsigned int numGenes,
										 bool singleViolationList)
{
	this->numGenes = numGenes;
	this->constraint = constraint;
	this->comp = new CompareStatesLess(numGenes);
	if (singleViolationList)
		this->violations.push_back(ViolationMap());
}

void ConstraintViolation::addViolation(unsigned int matching, bool * state, ConditionViolation &violation)
{
	ConditionViolation v = violation;
	if (violations[matching].count(v) == 0)
	{
		violations[matching][v] = new StateSet(*comp);
	}
	if (violations[matching][v]->count(state) == 0)
	{
		bool * nstate = new bool[numGenes];
		memcpy(nstate, state, sizeof(bool) * numGenes);
		(*violations[matching][v])[nstate] = 1;
	}
	else
		++(*violations[matching][v])[state];
}

void ConstraintViolation::addAttractorMatching(vector<bool *> attractor, unsigned int * assignment,
												unsigned int alternative)
{
	vector<bool *> attr;
	for (vector<bool *>::iterator it = attractor.begin(); it != attractor.end(); ++it)
	{
		bool * state = new bool[numGenes];
		memcpy(state, *it, sizeof(bool) * numGenes);
		attr.push_back(state);
	}
	unsigned int * ass = new unsigned int[attractor.size()];
	memcpy(ass, assignment, sizeof(unsigned int) * attractor.size());

	attractors.push_back(attr);
	attractorMatchings.push_back(ass);
	alternatives.push_back(alternative);
	violations.push_back(ViolationMap());
}

ConstraintViolation::~ConstraintViolation()
{
	delete comp;
	for (vector<ViolationMap>::iterator it1 = violations.begin(); it1 != violations.end(); ++it1)
	{
		for (ViolationMap::iterator it2 = it1->begin();
			 it2 != it1->end(); ++it2)
		{
			for (StateSet::iterator it3 = it2->second->begin(); it3 != it2->second->end(); ++it3)
				delete [] it3->first;
			delete it2->second;
		}
	}

	for (vector<vector<bool*> >::iterator it1 = attractors.begin(); it1 != attractors.end(); ++it1)
	{
		for (vector<bool *>::iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
			delete [] *it2;
	}

	for (vector<unsigned int*>::iterator it = attractorMatchings.begin(); it != attractorMatchings.end(); ++it)
		delete [] *it;
}

string ConstraintViolation::toString(std::string * varNames, unsigned int collapse)
{
	if (violations.size() == 0  || violations[0].size() == 0)
		return "(no violations)";

	ostringstream res;

	if (attractorMatchings.size() > 0)
	// the attractor matchings have been logged => print detailed matching information
	{
		for (unsigned int i = 0; i < attractorMatchings.size(); ++i)
		{
			res << "\nAttractor matching " << (i + 1) << " (using alternative " << (alternatives[i] + 1) << "):\n";
			res << setw(numGenes * 2) << left << "Attractor" << "\t=>\t" << setw(numGenes * 2) << left << "Specifications" << endl;

			unsigned int constraintSize = constraint->groupIndices[alternatives[i] + 1] - constraint->groupIndices[alternatives[i]];

			for (unsigned int j = 0; j < attractors[i].size(); ++j)
			{
				// print attractor state
				for (unsigned int k = 0; k < numGenes; ++k)
				{
					res << attractors[i][j][k] << " ";
				}
				res << "\t=>\t";

				SimpleCondition * spec = constraint->stateSpecifications[constraint->groupIndices[alternatives[i]] + attractorMatchings[i][j]];

				// print the coresponding state specification
				for (unsigned int k = 0; k < numGenes; ++k)
				{

					if (spec->geneValues.count(k) > 0)
						res << spec->geneValues[k];
					else
						res << "*";
					res << " ";
				}
				res << " (State spec. " << (attractorMatchings[i][j] + 1) << ")\n";

				int expectedNext = (attractorMatchings[i][j] + 1) % constraintSize;
				int trueNext = attractorMatchings[i][(j + 1) % attractors[i].size()];

				if (trueNext < expectedNext)
					trueNext += constraintSize;

				if (expectedNext != trueNext)
				// Some specification entries have been left out
				// => print them
				{
					for (int l = expectedNext; l < trueNext; ++l)
					{
						res << setw(numGenes * 2) << "" << "\t=>\t";

						unsigned int specIdx = l % constraintSize;
						SimpleCondition * spec = constraint->stateSpecifications[constraint->groupIndices[alternatives[i]] + specIdx];

						for (unsigned int k = 0; k < numGenes; ++k)
						{

							if (spec->geneValues.count(k) > 0)
								res << spec->geneValues[k];
							else
								res << "*";
							res << " ";
						}
						res << " (State spec. " << (specIdx + 1) << ")\n";
					}
				}
			}

			// print the violations
			res << "\nViolations caused by this matching:\n";
			for (ViolationMap::iterator it = violations[i].begin();
						 it != violations[i].end(); ++it)
			{
				res << "State specification " << (it->first.getSpecificationEntry()+1) << ": ";
				res << varNames[it->first.getConditionIndex()] << " != " << it->first.getCondition()->geneValues[it->first.getConditionIndex()] << "\n";
			}

			res << "\nStart states yielding this matching:\n";
			//unsigned int i = 0;

			// As all violations are caused by the same start states here, print them
			// only once using the first violation as a representative
			ViolationMap::iterator representative = violations[i].begin();
			for (StateSet::iterator it2 = representative->second->begin(); it2 != representative->second->end() && i < collapse; ++it2, ++i)
			{
				res << "\t";
				for (unsigned int j = 0; j < numGenes; ++j)
					res << it2->first[j] << " ";
				res << "\n";
			}
			if (representative->second->size() > collapse)
			{
				res << "\t(... further " << representative->second->size() - collapse << " states ...)\n";
			}
		}
	}
	else
	{
		res << "\nViolations of the chain rule:\n";
		for (ViolationMap::iterator it = violations[0].begin();
				 it != violations[0].end(); ++it)
		{
			res << "\nAlternative " <<  (it->first.getAlternative() + 1) << ", State specification " << (it->first.getSpecificationEntry()+1) << ": ";
			res << varNames[it->first.getConditionIndex()] << " != " << it->first.getCondition()->geneValues[it->first.getConditionIndex()] << " for start states: \n";

			unsigned int i = 0;
			for (StateSet::iterator it2 = it->second->begin(); it2 != it->second->end() && i < collapse; ++it2, ++i)
			{
				res << "\t";
				for (unsigned int j = 0; j < numGenes; ++j)
					res << it2->first[j] << " ";
				res << "\n";
			}
			if (it->second->size() > collapse)
			{
				res << "\t(... further " << it->second->size() - collapse << " states ...)\n";
			}
		}
	}
	return res.str();
}
