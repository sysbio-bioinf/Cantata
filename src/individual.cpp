/*
 * Individual.cpp
 * Defines a wrapper class for an individual
 * of the GP approach
 *
 *  Created on: Apr 15, 2010
 *      Author: muessel
 * Edited 2020 by schwab
 */

#include "individual.h"
#include "random.h"
#include <sstream>
#include <cstring>
#include <iostream>

using namespace std;

Individual::Individual(GP &parent, std::vector<BooleanTree *> &network, std::vector<std::string> &varNames,
					   unsigned int run, unsigned int generation,
					   unsigned int startMutations, bool replaceUnknown)
{
	this->parent = &parent;
	this->fitness = new double[parent.getNumFitnessFunctions()];
	this->run = run;
	this->generation = generation;

	unsigned int newGenes = 0;//intrand(6);

//	if (newGenes < 4)
//		newGenes = 0;
//	else
//		newGenes -= 3;

	this->network.reserve(network.size() + newGenes);
	this->varNames.assign(varNames.begin(),varNames.end());
	this->network.reserve(network.size() + newGenes);
	this->normalizedNetwork.reserve(network.size() + newGenes);

	for (vector<BooleanTree *>::iterator it = network.begin(); it != network.end(); ++it)
	{
		this->network.push_back((*it)->copy());
		this->normalizedNetwork.push_back(NULL);
	}

	if (replaceUnknown)
	{
		for (vector<BooleanTree *>::iterator it = this->network.begin(); it != this->network.end(); ++it)
		{
			if (typeid(*(*it)->getRoot()) == typeid(Literal))
			{
				Literal * lit = (Literal *)(*it)->getRoot();

				if (lit->literalIndex == -2)
				// replace "UNKNOWN" symbols by random literals
				{
					lit->literalIndex = intrand((unsigned int)varNames.size());
					lit->negated = (intrand(2) != 0);
				}

			}
		}
	}

	for (unsigned int i = 0; i < startMutations; ++i)
	// apply a mutation to the new individual
	{
//		for (unsigned int i = 0; i < newGenes; ++i)
//		{
//			ostringstream name;
//			name << "Gene " << this->varNames.size();
//			this->varNames.push_back(name.str());
//		}
//		for (unsigned int i = 0; i < newGenes; ++i)
//		{
//			BooleanTree * tree = new BooleanTree(randomOperator());
//			((BooleanOperator *)tree->getRoot())->addNewOperand(randomLiteral(this->varNames.size()));
//			this->network.push_back(tree);
//		}
//		for (unsigned int i = 0; i < newGenes + 1; ++i)

		this->mutate();
	}

	// calculate cached normalized versions
	for (unsigned int i = 0; i < normalizedNetwork.size(); ++i)
	{
		if (normalizedNetwork[i] == NULL)
			updateCachedNormalization(i);
	}
}

#ifdef TRACE_EVOLUTION
Individual::Individual(Individual &clone, bool copyAncestors)
#else
Individual::Individual(Individual &clone)
#endif
{
	this->parent = clone.parent;

	this->fitness = new double[parent->getNumFitnessFunctions()];
	memcpy(this->fitness,clone.fitness,parent->getNumFitnessFunctions()*sizeof(double));

	this->run = clone.run;
	this->generation = clone.generation;

	this->network.reserve(clone.network.size());
	this->normalizedNetwork.reserve(clone.network.size());
	this->varNames.assign(clone.varNames.begin(),clone.varNames.end());

	for (vector<BooleanTree *>::iterator it = clone.network.begin(); it != clone.network.end(); ++it)
	{
		this->network.push_back((*it)->copy());
	}

	for (vector<BooleanTree *>::iterator it = clone.normalizedNetwork.begin(); it != clone.normalizedNetwork.end(); ++it)
	{
		this->normalizedNetwork.push_back((*it)->copy());
	}

#ifdef TRACE_EVOLUTION
	if (copyAncestors)
	{
		for (vector<Individual *>::iterator it = clone.ancestors.begin(); it != clone.ancestors.end(); ++it)
		{
			this->ancestors.push_back(new Individual(**it, false));
		}
	}
#endif
}


Individual::~Individual()
{
	for (vector<BooleanTree *>::iterator it = network.begin(); it != network.end(); ++it)
	{
		delete *it;
	}

	for (vector<BooleanTree *>::iterator it = normalizedNetwork.begin(); it != normalizedNetwork.end(); ++it)
	{
		delete *it;
	}

#ifdef TRACE_EVOLUTION
	for (vector<Individual *>::iterator it = ancestors.begin(); it != ancestors.end(); ++it)
	{
		delete *it;
	}
#endif

	delete [] fitness;
}

Individual * Individual::copy()
{
	return new Individual(*this->parent,this->network,this->varNames,
						  this->run,this->generation,false,false);
}

void Individual::mutate(bool allowNegation)
{
#ifdef TRACE_EVOLUTION
	if (ancestors.size() == 0 || !(*ancestors[ancestors.size() - 1] == *this))
		ancestors.push_back(this->copy());//new Individual(*this, false));
#endif

//	return;

	// calculate the CDF of the uncertainties of the functions
	double * uncertainty = new double[network.size()];
	double uncertaintySum = 0.0;

	for (unsigned int i = 0; i < network.size(); ++i)
	{
		uncertainty[i] = network[i]->getUncertainty();
		uncertaintySum += uncertainty[i];
	}
//	cout << "US " << uncertaintySum << endl;

	unsigned int chosenTree;
	double cdf = 0.0;
	double prob = doublerand(uncertaintySum);

	// determine the function to modify according to the drawn random
	// value and the CDF
	for (chosenTree = 0; chosenTree < network.size(); ++chosenTree)
	{
		cdf += uncertainty[chosenTree];
		if (cdf > prob)
			break;
	}
	delete [] uncertainty;

//	cout << chosenTree << endl;

	BooleanTree * tree = network[chosenTree];

	unsigned int chosenIndex = intrand(tree->getNumberOfElements());
	TreePosition chosenPosition = tree->getElementByIndex(chosenIndex);
//	cdf = 0.0;
//	prob = doublerand(uncertainty[chosenTree]);
//	TreePosition chosenPosition = tree->getElementByIndex(0);
//
//	while (chosenPosition.element != NULL)
//	{
//		cdf += chosenPosition.element->uncertainty;
//
//		if (cdf > prob)
//			break;
//
//		chosenPosition = tree->getNextElement(chosenPosition);//Operator(chosenPosition);
//	}

	//cout << "\nChosen: " << chosenPosition.element->toString(parent->getVarNames()) << endl;
	if (typeid(*chosenPosition.element) == typeid(BooleanOperator))
	// an operator is changed
	{
		BooleanOperator * op = (BooleanOperator *)chosenPosition.element;
		BooleanOperator * subFormula;
		unsigned int mutationType, childIdx;

		if (allowNegation)
			mutationType = intrand(5);
		else
			mutationType = intrand(4);
		switch(mutationType)
		{
			case 0:
				//cout << "Change operator" << endl;
				// change type of operator
				op->operatorType = (OperatorType)intrand(2);
				break;
			case 1:
				// delete an element of the operator
				if (op->operands.size() == 2 && op != tree->getRoot())
				// only one element left => remove operator
				{
					BooleanFormula * replacement = op->operands[intrand((unsigned int)op->operands.size())]->copy(NULL);
					replacement->negated = (op->negated && !replacement->negated)
											|| (!op->negated && replacement->negated);
					tree->replaceElement(chosenPosition, replacement);
				}
				else
				if (op->operands.size() > 2 || (op == tree->getRoot() && op->operands.size() == 2))
				{
					op->deleteOperand(intrand((unsigned int)op->operands.size()));
				}
				break;
			case 2:
//				for (unsigned int i = 0; i < 20; ++i)
//				{
//					litIdx = intrand(this->varNames.size());
//					bool found = false;
//					for (std::vector<BooleanFormula*>::iterator it = op->operands.begin();
//						it != op->operands.end(); ++it)
//					{
//						if (typeid(**it) == typeid(Literal))
//						{
//							Literal * lit = (Literal *)*it;
//							if (lit->literalIndex == litIdx)
//							{
//								found = true;
//								break;
//							}
//						}
//					}
//					if (!found)
//					{
//						new Literal(litIdx,op,(bool)intrand(2));
//						break;
//					}
//				}
				// insert a new literal
				randomLiteral((unsigned int)this->varNames.size(),op);
				break;
			case 3:
				// insert new subformula
				childIdx = intrand((unsigned int)op->operands.size());

//				if (typeid(*op->operands[childIdx]) == typeid(Literal) && intrand(2) == 1)
//				{
//					((Literal *)op->operands[childIdx])->literalIndex = intrand(this->varNames.size());
//				}
//				else
				{
					subFormula = randomOperator(NULL);
					subFormula->addNewOperand(op->operands[childIdx]);
					op->operands[childIdx]->setParent(subFormula);
					randomLiteral((unsigned int)this->varNames.size(),subFormula);

					tree->replaceElement(TreePosition(op,childIdx,op->operands[childIdx]),subFormula,false);
				}
				break;
			case 4:
//				if (intrand(2) == 1)
//				{
//					childIdx = intrand(op->operands.size());
//					op->operands[childIdx]->negated = !op->operands[childIdx]->negated;
//				}
//				else
				// negate an operator
				op->negated = !op->negated;
				break;
		}
	}
	else
	// a literal is changed
	{
		unsigned int mutationType;

		if (allowNegation)
			mutationType = intrand(4);
		else
			mutationType = intrand(3);
		Literal * lit = (Literal *)chosenPosition.element;
		BooleanOperator * op;
		switch(mutationType)
		{
			case 0:
				// exchange literals
				//cout << "Change literal" << endl;
				// choose different literal
				//lit->literalIndex = intrand(parent->getNumberOfVars() + 1) - 1;
				lit->literalIndex = intrand((unsigned int)this->varNames.size());
				break;
			case 1:
				// delete the literal
				//break;
				//cout << "Delete literal 2" << endl;
				// delete literal
				if (chosenPosition.parent != NULL)
				{
					if (chosenPosition.parent->operands.size() == 2)
					{
						BooleanFormula * replacement = lit->copy(NULL);
						replacement->negated = (chosenPosition.parent->negated && !lit->negated)
												|| (!chosenPosition.parent->negated && lit->negated);

						BooleanOperator * parParent = chosenPosition.parent->parent;
						unsigned int parentIndex = 0;
						if (parParent != NULL)
						{
							for (;parParent->operands[parentIndex] != chosenPosition.parent; ++parentIndex);
						}
						//cout << parentIndex << endl;
						TreePosition parPos(parParent,parentIndex,chosenPosition.parent);

						tree->replaceElement(parPos, replacement);
					}
					else
					if (chosenPosition.parent->operands.size() > 2)
					{
						//cout << chosenPosition.index << endl;
						chosenPosition.parent->deleteOperand(chosenPosition.index);
						//cout << chosenPosition.parent->toString(parent->getVarNames()) << endl;
					}
				}
				break;
			case 2:
				// create a new operator and add this literal as well as a new one
				//cout << "New operator" << endl;
				// include literal in sub-formula
				op = randomOperator(NULL);
				//lit->copy(op);
				op->addNewOperand(lit);
				randomLiteral((unsigned int)this->varNames.size(),op);
				tree->replaceElement(chosenPosition,op,false);
				break;
			case 3:
				// negate the literal
				//cout << "Negate literal" << endl;
				// negate literal
				lit->negated = !lit->negated;
				break;
		}
	}
	// simplify the new expression, and mark it as changed
	tree->simplify();
	tree->touch();

	// update cached normalized copy
	updateCachedNormalization(chosenTree);
}

void Individual::updateCachedNormalization(unsigned int index)
{
	if (normalizedNetwork[index] != NULL)
		delete normalizedNetwork[index];
	normalizedNetwork[index] = network[index]->getNormalizedCopy();
}

void Individual::stateTransition(bool * state, SimpleCondition * fixedGenes)
{
	bool * oldState = new bool[varNames.size()];
	memcpy(oldState,state,sizeof(bool) * varNames.size());
	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		if (fixedGenes != NULL && fixedGenes->contains(i))
			state[i] = (*fixedGenes)[i];
		else
			state[i] = network[i]->evaluate(oldState);
	}
	delete [] oldState;
}

void Individual::simplify()
{
	// simplify all functions of the network
	for (unsigned int i = 0; i < network.size(); ++i)
	{
		network[i]->simplify();
		updateCachedNormalization(i);
	}
	parent->updateFitness(*this);
}

std::string Individual::toString(std::string separator, bool printFitness, bool printInfo) const
{
	ostringstream str;
	for (unsigned int i = 0; i < network.size(); ++i)
	{
		str << varNames[i] << separator << network[i]->toString(&varNames[0]) << "\n";
	}
	if (printFitness)
	{
		str << "Fitness: ";
		for (unsigned int i = 0; i < parent->getNumFitnessFunctions(); ++i)
		{
			if (i != 0)
				str << " ";
			str << fitness[i];
		}
		if (printInfo)
			str << " ";
	}
	if (printInfo)
		str << "Run: " << (run+1) << " Generation: " << generation;
	return str.str();
}

bool Individual::operator==(const Individual &ind) const
{
	if (ind.network.size() != this->network.size())
		return false;

	// compare the functions of the network
	for (unsigned int i = 0; i < network.size(); ++i)
	{
		if (!(*this->network[i] == *ind.network[i]))
			return false;
	}
	return true;
}
