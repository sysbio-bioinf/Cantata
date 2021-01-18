/*
 * randomnet.h
 * Utilities to randomly create and mutate networks
 *
 *  Created on: Oct 29, 2010
 *      Author: muessel
 */
#include "individual.h"
#include "random.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>

/**
 * Applies a maximum of <numMutations> mutations to the expression <formula>.
 * <numVariables> is the total number of variables in the network.
 * If <deletionsAllowed> is true, sub-expressions may be deleted from the expression.
 * If <completeRemovalAllowed>, a formula can be UNKNOWN after deletion.
 * If <alterationsAllowed> is true, sub-expressions may be inserted or modified.
 * Returns the actual number of mutations applied to the expression.
 */
unsigned int mutateFormula(BooleanTree * formula, unsigned int numVariables, unsigned int numMutations,
				  double insertionProb, bool completeRemovalAllowed)
{
	unsigned int res = 0;
	for (unsigned int i = 0; i < numMutations; ++i)
	{
		if (typeid(*formula->getRoot()) == typeid(BooleanOperator))
			((BooleanOperator *)formula->getRoot())->recalculateNumberOfElements();
		unsigned int idx = intrand(formula->getNumberOfElements());
		TreePosition el = formula->getElementByIndex(idx);

		double r = doublerand_1();

		if (typeid(*el.element) == typeid(BooleanOperator))
		{
			BooleanOperator * op = (BooleanOperator * )el.element;

			if (r < insertionProb)
			{
				if (intrand(2) == 0)
				{
					unsigned int elIdx = intrand(op->operands.size());
					TreePosition pos(op, elIdx, op->operands[elIdx]);
					BooleanOperator * subOp = randomOperator();
					subOp->addNewOperand(op->operands[elIdx]->copy());
					subOp->addNewOperand(randomLiteral(numVariables));
					formula->replaceElement(pos,subOp,true);
					++res;
				}
				else
				{
					op->addNewOperand(randomLiteral(numVariables));
					++res;
				}
			}
			else
			{
				if (intrand(2) == 0)
				{
					if (op->operands.size() > 1)
					{
						op->deleteOperand(intrand((unsigned int)op->operands.size()));
						++res;
					}
				}
				else
				{
					if (el.parent != NULL && (el.parent->operands.size() > 1 || completeRemovalAllowed))
					{
						el.parent->deleteOperand(el.index);
						++res;
					}
				}
			}
		}
		else
		{
			if (r < insertionProb)
			{
				Literal * l = (Literal *)el.element->copy();
				BooleanOperator * op = randomOperator();
				op->addNewOperand(l);
				op->addNewOperand(randomLiteral(numVariables));
				formula->replaceElement(el,op,true);
				++res;
			}
			else
			{
				if (el.parent != NULL && (el.parent->operands.size() > 1 || completeRemovalAllowed))
				{
					el.parent->deleteOperand(el.index);
					++res;
				}
				else
				if (el.parent == NULL && completeRemovalAllowed)
				{
					formula->replaceElement(el,new BooleanOperator(OP_AND),true);
					++res;
				}
			}
		}

	}
	formula->simplify();
	return res;
}

///**
// * Applies a maximum of <numMutations> mutations to the expression <formula>.
// * <numVariables> is the total number of variables in the network.
// * If <deletionsAllowed> is true, sub-expressions may be deleted from the expression.
// * If <completeRemovalAllowed>, a formula can be UNKNOWN after deletion.
// * If <alterationsAllowed> is true, sub-expressions may be inserted or modified.
// * Returns the actual number of mutations applied to the expression.
// */
//unsigned int mutateFormula(BooleanTree * formula, unsigned int numVariables, unsigned int numMutations,
//				  bool deletionsAllowed, bool completeRemovalAllowed, bool alterationsAllowed)
//{
//	unsigned int res = 0;
//	for (unsigned int i = 0; i < numMutations; ++i)
//	{
//		if (typeid(*formula->getRoot()) == typeid(BooleanOperator))
//			((BooleanOperator *)formula->getRoot())->recalculateNumberOfElements();
//		unsigned int idx = intrand(formula->getNumberOfElements());
//		TreePosition el = formula->getElementByIndex(idx);
//		unsigned int r;
//
//		if (typeid(*el.element) == typeid(BooleanOperator))
//		{
//			BooleanOperator * op;
//
//			if (deletionsAllowed)
//				r = intrand(8) + 1;//intrand(9);
//			if (!alterationsAllowed)
//				r = 5 + intrand(4);
//			else
//				r = intrand(5);
//
//			switch(r)
//			{
//				case 0:
//					el.element->negated = !el.element->negated;
//					++res;
//					break;
//				case 1:
//				case 2:
//					op = randomOperator();
//					op->addNewOperand(randomLiteral(numVariables));
//					op->addNewOperand(randomLiteral(numVariables));
//					((BooleanOperator *)el.element)->addNewOperand(op);
//					++res;
//					break;
//				case 3:
//				case 4:
//					((BooleanOperator *)el.element)->addNewOperand(randomLiteral(numVariables));
//					++res;
//					break;
//				case 5:
//				case 6:
//					if (((BooleanOperator *)el.element)->operands.size() > 1)
//					{
//						((BooleanOperator *)el.element)->deleteOperand(intrand((unsigned int)((BooleanOperator *)el.element)->operands.size()));
//						++res;
//					}
//					break;
//				case 7:
//				case 8:
//					if (el.parent != NULL && (el.parent->operands.size() > 1 || completeRemovalAllowed))
//					{
//						el.parent->deleteOperand(el.index);
//						++res;
//					}
//					break;
//			}
//		}
//		else
//		{
//			r = intrand(2);
//			if (!alterationsAllowed)
//				r = 1;
//
//			Literal * l;
//			BooleanOperator * op;
//
//			switch(r)
//			{
//				case 0:
//					l = (Literal *)el.element->copy();
//					op = randomOperator();
//					op->addNewOperand(l);
//					op->addNewOperand(randomLiteral(numVariables));
//					formula->replaceElement(el,op,true);
//					++res;
//					break;
//				case 1:
//					if (deletionsAllowed && el.parent != NULL && (el.parent->operands.size() > 1 || completeRemovalAllowed))
//					{
//						el.parent->deleteOperand(el.index);
//						++res;
//					}
//					else
//					if (deletionsAllowed && el.parent == NULL && completeRemovalAllowed)
//					{
//						formula->replaceElement(el,new BooleanOperator(OP_AND),true);
//						++res;
//					}
//			}
//			//el.element->negated = !el.element->negated;
//		}
//
//	}
//	formula->simplify();
//	return res;
//}

/**
 * Creates a random formula with <numVariables> potential input variables
 * by applying <numMutations> random mutations.
 */
BooleanTree * randomFormula(unsigned int numVariables, unsigned int numMutations)
{
	BooleanTree * res = new BooleanTree(randomOperator());
	BooleanOperator * root = ((BooleanOperator *)res->getRoot());
	root->addNewOperand(randomLiteral(numVariables));
	root->addNewOperand(randomLiteral(numVariables));


	numMutations = intrand(numMutations);
	mutateFormula(res,numVariables,numMutations,1.0,false);

	return res;
}

/**
 * Writes a randomly created network with <numGenes> genes to <outputFile>.
 * The network is created by applying <numOperations> mutations.
 */
void writeRandomNetworks(std::string outputFile,
						 unsigned int numGenes,
						 unsigned int numOperations)
{
	std::string * geneNames = new std::string[numGenes];

	for (unsigned int i = 0; i < numGenes; ++i)
	{
		std::ostringstream oss;
		oss << "gene" << i;
		geneNames[i] = oss.str();
	}

	std::fstream output(outputFile.c_str(),std::ios::out);

	output << "targets, factors\n";

	for (unsigned int i = 0; i < numGenes; ++i)
	{
		BooleanTree * formula = randomFormula(numGenes, numOperations);
		output << "gene" << i << ", " << formula->toString(geneNames) << std::endl;
		delete formula;
	}
	output.close();

	delete [] geneNames;
}

void internal_insertDependency(std::vector<BooleanTree *> network, int targetGene, int dependency)
{
	BooleanTree * formula = network[targetGene];
	if (typeid(*formula->getRoot()) == typeid(BooleanOperator))
				((BooleanOperator *)formula->getRoot())->recalculateNumberOfElements();
	unsigned int insertionIdx = intrand(formula->getNumberOfElements());
	TreePosition el = formula->getElementByIndex(insertionIdx);

	if (typeid(*el.element) == typeid(BooleanOperator))
	{
		BooleanOperator * op = (BooleanOperator * )el.element;


		if (intrand(2) == 0)
		{
			unsigned int elIdx = intrand(op->operands.size());
			TreePosition pos(op, elIdx, op->operands[elIdx]);
			BooleanOperator * subOp = randomOperator();
			op->operands[elIdx]->copy(subOp);
			new Literal(dependency,subOp,(bool)intrand(2));
			formula->replaceElement(pos,subOp,true);
		}
		else
		{
			new Literal(dependency,op,(bool)intrand(2));
		}
	}
	else
	{
		BooleanOperator * op = randomOperator();
		el.element->copy(op);
		new Literal(dependency,op,(bool)intrand(2));
		formula->replaceElement(el,op,true);
	}
}

void internal_removeDependency(std::vector<BooleanTree *> network, int targetGene, int dependency)
{
	std::list<TreePosition> queue;
	queue.push_back(TreePosition(NULL,0,network[targetGene]->getRoot()));
	while(queue.size() != 0)
	{
		TreePosition pos = *queue.begin();
		if (typeid(*pos.element) == typeid(BooleanOperator))
		{
			BooleanOperator * op = ((BooleanOperator *)pos.element);
			unsigned int i = 0;
			for (std::vector<BooleanFormula*>::iterator it = op->operands.begin();
								 it != op->operands.end(); ++it, ++i)
			{
				queue.push_back(TreePosition(op,i,*it));
			}
		}
		else
		{
			Literal * lit = ((Literal *)pos.element);
			if (lit->literalIndex >= 0)
			{
				if (lit->literalIndex == dependency)
				{
					if (pos.parent != NULL)
					{
						pos.parent->deleteOperand(pos.index);

						BooleanOperator * parent = pos.parent;

						while (parent->operands.size() == 0 && parent->parent != NULL)
						{
							unsigned int i = 0;
							for (std::vector<BooleanFormula*>::iterator it = parent->parent->operands.begin();
															 it != parent->parent->operands.end(); ++it, ++i)
							{
								if (*it == parent)
									break;
							}
							parent = parent->parent;
							parent->deleteOperand(i);
						}
					}
				}
			}
		}
		queue.pop_front();
	}
}

/**
 * Loads an existing network from <originalFile>, applies <numMutations> deletions,
 * and saves the truncated network to <mutatedFile>.
 * If <completeRemovalAllowed>, a formula can be UNKNOWN after deletion
 * If <fullDependencyChanges>, deletions remove all occurences of a dependency, and insertions
 * only insert new dependencies.
 */
void mutateNetwork(std::string originalFile, std::string mutatedFile, unsigned int numMutations,
				  double insertionProb, bool completeRemovalAllowed, bool fullDependencyChanges)
{
	std::vector<BooleanTree *> network;
	std::vector<std::string> varNames;

	loadNetwork(originalFile,varNames,network,true);

	if (fullDependencyChanges)
	{
		for (unsigned int z = 0; z < numMutations; ++z)
		{
			if (doublerand_1() < insertionProb)
			{
				int targetGene = (int)intrand((unsigned int)network.size());
				int dependency = -1;
				LiteralMap * literals = network[targetGene]->getLiterals();

				bool cont = true;
				while (cont)
				{
					dependency = (int)intrand(network.size());
					cont = false;
					for (LiteralMap::const_iterator it = literals->begin(); it != literals->end(); ++it)
					{
						if (it->first->literalIndex == dependency)
						{
							cont = true;
							break;
						}
					}
				}
				internal_insertDependency(network, targetGene, dependency);
				network[targetGene]->touch();
			}
			else
			{
				bool onlyLiteral = true;
				int targetGene = -1, dependency = -1;
				while (onlyLiteral && !completeRemovalAllowed)
				{
					targetGene = (int)intrand((unsigned int)network.size());
					dependency = -1;
					LiteralMap * literals = network[targetGene]->getLiterals();

					unsigned int idx = intrand(literals->size()), i = 0;
					for (LiteralMap::const_iterator it = literals->begin(); it != literals->end(); ++it, ++i)
					{
						if (i == idx)
						{
							dependency = it->first->literalIndex;
							break;
						}
					}
					for (LiteralMap::const_iterator it = literals->begin(); it != literals->end(); ++it)
					{
						if (it->first->literalIndex != dependency)
						{
							onlyLiteral = false;
							break;
						}
					}
				}
				internal_removeDependency(network, targetGene, dependency);
				network[targetGene]->touch();
//				BooleanTree * formula = network[targetGene];
//				if (typeid(*formula->getRoot()) == typeid(BooleanOperator))
//					((BooleanOperator *)formula->getRoot())->recalculateNumberOfElements();
			}

		}
		for (unsigned int i = 0; i < varNames.size(); ++i)
		{
			network[i]->simplify();
		}
	}
	else
	{
		for (unsigned int i = 0; i < numMutations; )
		{
			BooleanTree * formula = network[intrand((unsigned int)network.size())];
			i += mutateFormula(formula, (unsigned int)varNames.size(), 1, insertionProb, completeRemovalAllowed);
		}
	}

	std::fstream output(mutatedFile.c_str(),std::ios::out);
	output << "targets, factors\n";

	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		output << varNames[i] << ", " << network[i]->toString(&varNames[0]) << std::endl;
		delete network[i];
	}
	output.close();
}

/**
 * Removes all dependencies of <targetGene> on <dependency> in the network <originalFile>,
 * and writes it to <targetFile>.
 */
void removeDependency(std::string originalFile, std::string mutatedFile,
					  std::string targetGene, std::string dependency)
{
	std::vector<BooleanTree *> network;
	std::vector<std::string> varNames;

	std::transform(targetGene.begin(), targetGene.end(), targetGene.begin(), (int(*)(int))toupper);
	std::transform(dependency.begin(), dependency.end(), dependency.begin(), (int(*)(int))toupper);

	loadNetwork(originalFile,varNames,network,true);

	int targetGeneIdx = -1, dependencyIdx = -1;
	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		std::string gene = varNames[i];
		std::transform(gene.begin(), gene.end(), gene.begin(), (int(*)(int))toupper);
		if (gene.compare(targetGene) == 0)
			targetGeneIdx = i;
		if (gene.compare(dependency) == 0)
			dependencyIdx = i;
	}

	if (targetGeneIdx == -1)
		throw std::string("Gene "+targetGene+" not known!");
	else
	if (dependencyIdx == -1)
		throw std::string("Gene "+dependency+" not known!");

	if (network[targetGeneIdx]->getRoot() != NULL)
	{
		internal_removeDependency(network, targetGeneIdx, dependencyIdx);
	}

	std::fstream output(mutatedFile.c_str(),std::ios::out);
	output << "targets, factors\n";
	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		output << varNames[i] << ", " << network[i]->toString(&varNames[0]) << std::endl;
		delete network[i];
	}
	output.close();

//	std::fstream output(mutatedFile.c_str(),std::ios::out);
//	output << "targets, factors\n";
//
//	for (unsigned int i = 0; i < varNames.size(); ++i)
//	{
//		std::string gene = varNames[i];
//		std::transform(gene.begin(), gene.end(), gene.begin(), (int(*)(int))toupper);
//		if (gene == targetGene && network[i]->getRoot() != NULL)
//		{
//			std::list<TreePosition> queue;
//			queue.push_back(TreePosition(NULL,0,network[i]->getRoot()));
//			while(queue.size() != 0)
//			{
//				TreePosition pos = *queue.begin();
//				if (typeid(*pos.element) == typeid(BooleanOperator))
//				{
//					BooleanOperator * op = ((BooleanOperator *)pos.element);
//					unsigned int i = 0;
//					for (std::vector<BooleanFormula*>::iterator it = op->operands.begin();
//										 it != op->operands.end(); ++it, ++i)
//					{
//						queue.push_back(TreePosition(op,i,*it));
//					}
//				}
//				else
//				{
//					Literal * lit = ((Literal *)pos.element);
//					if (lit->literalIndex >= 0)
//					{
//						std::string dep = varNames[lit->literalIndex];
//						std::transform(dep.begin(), dep.end(), dep.begin(), (int(*)(int))toupper);
//						if (dep == dependency)
//						{
//							//std::cout << targetGene << " " << dep << " " << dependency << std::endl;
//							//if (pos.parent == NULL)
//							//{
//							//	delete network[i];
//							//	network[i] = NULL;
//							//}
//							//else
//							if (pos.parent != NULL)
//								pos.parent->deleteOperand(pos.index);
//						}
//					}
//				}
//				queue.pop_front();
//				//network[i]->simplify();
//			}
//		}
//		output << varNames[i] << ", " << network[i]->toString(&varNames[0]) << std::endl;
//		delete network[i];
//	}
//	output.close();
}
