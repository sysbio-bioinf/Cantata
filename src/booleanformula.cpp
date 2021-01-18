/*
 * booleanformula.cpp
 * Provides classes for the management of
 * Boolean expressions
 *
 *  Created on: Apr 20, 2010
 *      Author: muessel
 * 	Edited 2020 by schwab
 */
#include "booleanformula.h"
#include <algorithm>

using namespace std;

/**
 * Comparison operator for Boolean formula trees
 * to obtain a unique order of subtrees
 */
struct sortFormula
{
	bool operator() (const BooleanFormula * f1, const BooleanFormula * f2) const
	{
		if (typeid(*f1) == typeid(Literal) && typeid(*f2) == typeid(Literal))
		// compare literals by literal index and negation
		{
			Literal * l1 = ((Literal *)f1);
			Literal * l2 = ((Literal *)f2);
			return (l1->literalIndex * 2 + (unsigned int)l1->negated < l2->literalIndex * 2 + (unsigned int)l2->negated);
		}
		else
		// literals are "smaller" than operators
		if (typeid(*f1) == typeid(Literal))
		{
			return true;
		}
		else
		if (typeid(*f2) == typeid(Literal))
		{
			return false;
		}
		else
		// compare two operator subtrees
		{
			BooleanOperator * o1 = ((BooleanOperator *)f1);
			BooleanOperator * o2 = ((BooleanOperator *)f2);
			if (o1->operands.size() == 0)
				return true;
			if (o2->operands.size() == 0)
				return false;

			unsigned int i;
			for (i = 0; i < o1->operands.size() && i < o2->operands.size(); ++i)
			{
				if (!(o1->operands[i] == o2->operands[i]))
					break;
			}
			if (i == o1->operands.size() || i == o2->operands.size())
				// all operands are equal
				return false;

			// recursive comparison of subtrees
			return sortFormula::operator()(o1->operands[i],o2->operands[i]);
		}
	}
};

bool BooleanOperator::evaluate(const bool * variableSet) const
{
	bool res;
	switch (operatorType)
	{
		case OP_AND:
			res = true;
			for (std::vector<BooleanFormula*>::const_iterator it = operands.begin();
				 it != operands.end(); ++it)
			{
				// evaluate operand
				res &= (*it)->evaluate(variableSet);

				if (!res)
				// break if any operand is false
					return (negated);
			}
			return (!negated);
		case OP_OR:
			res = false;
			for (std::vector<BooleanFormula*>::const_iterator it = operands.begin();
				 it != operands.end(); ++it)
			{
				// evaluate operand
				res |= (*it)->evaluate(variableSet);
				if (res)
					// break if any operand is true
					return (!negated);
			}
			return (negated);
		case OP_XOR:
			if (operands.size() == 0)
				return false;
			res = operands[0]->evaluate(variableSet);
			for (unsigned int i = 1; i < operands.size(); ++i)
			{
				bool next = operands[i]->evaluate(variableSet);
				res = (res && !next) || (!res && next);
			}
			if (negated)
				return !res;
			else
				return res;
	}
	return false;
}

std::string BooleanOperator::toString(const std::string * variableNames) const
{
	if (operands.size() == 0)
		return "UNKNOWN";

	std::ostringstream res;

	if (negated)
		res << "!";

	res << "(";
	for (unsigned int i = 0; i < operands.size(); ++i)
	{
		if (i != 0)
		{
			switch (operatorType)
			{
				case OP_AND:
					res << " & ";

					break;
				case OP_OR:
					res << " | ";

					break;
				case OP_XOR:
					res << " XOR ";
					break;
			}
		}
		res << operands[i]->toString(variableNames);
	}
	res << ")";
	return res.str();
}

BooleanFormula * BooleanOperator::copy(BooleanOperator * parent)
{
	// create empty operator
	BooleanOperator * res = new BooleanOperator(this->operatorType, parent, negated);

	//copy operands
	for (std::vector<BooleanFormula*>::iterator it = operands.begin();
		 it != operands.end(); ++it)
	{
		(*it)->copy(res);
	}
	return res;
}

BooleanFormula * BooleanOperator::getSimplifiedCopy()
{
	// temporary operand list
	vector<BooleanFormula *> resList;
	// specifies whether the operand is included in the result
	vector<bool> useList;
	unsigned int useCount = 0;

	// obtain simplified versions of the operands
	for (std::vector<BooleanFormula*>::iterator it1 = operands.begin();
		 it1 != operands.end(); ++it1)
	{
		resList.push_back((*it1)->getSimplifiedCopy());
		useList.push_back(true);
		++useCount;
	}

	// determine which operands can be replaced or deleted
	for (unsigned int i = 0; i < resList.size(); ++i)
	{
		if (!useList[i])
		// this operand has already been excluded
			continue;

		if (typeid(*resList[i]) == typeid(BooleanOperator))
		{
			BooleanOperator * obj = (BooleanOperator *)(resList[i]);
			if (obj->operands.size() == 1)
			// "unlist" operators with only one element
			{
				resList[i] = obj->operands[0]->copy(NULL);
				resList[i]->negated = (!resList[i]->negated && obj->negated) || (resList[i]->negated && !obj->negated);
				delete obj;
				continue;
			}
			else
			if (obj->operatorType == this->operatorType && obj->negated == false)
			// "unlist" operators with the same type
			{
				useList[i] = false;
				--useCount;

				for (std::vector<BooleanFormula*>::iterator it2 = obj->operands.begin();
					 it2 != obj->operands.end(); ++it2)
				{
					resList.push_back((*it2)->copy(NULL));
					useList.push_back(true);
					++useCount;
				}
			}
		}
		else
		{
			Literal * lit1 = (Literal *)(resList[i]);
			for (unsigned int j = 0; j < resList.size(); ++j)
			{
				if (useList[j] && i != j
					&& typeid(*resList[j]) == typeid(Literal))
				{
					Literal * lit2 = (Literal *)(resList[j]);
					if (lit1->literalIndex == lit2->literalIndex)
					{
						if (lit1->negated == lit2->negated)
						// remove duplicate literals
						{
							useList[j] = false;
							--useCount;
						}
						else
						// found "Literal &/| !Literal)
						{
							if (this->operatorType == OP_AND)
							// replace "Literal & !Literal" by 0
							{
								lit1->literalIndex = -1;
								lit1->negated = true;
							}
							else
							if (this->operatorType == OP_OR)
							// replace "Literal | !Literal" by 1
							{
								lit1->literalIndex = -1;
								lit1->negated = false;
							}
							if (this->operatorType == OP_AND || this->operatorType == OP_OR)
							// remove all other occurrences of the literal
							{
								for (unsigned int j = 0; j < resList.size(); ++j)
								{
									if (i != j)
									{
										useList[j] = false;
										--useCount;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// remove duplicate 0 and 1 from the simplified list
	for (unsigned int i = 0; i < resList.size(); ++i)
	{
		if (useList[i] && typeid(*resList[i]) == typeid(Literal))
		{
			Literal * lit1 = (Literal *)(resList[i]);
			if (lit1->literalIndex == -1)
			{
				if (((lit1->negated && this->operatorType == OP_OR) ||
					(!lit1->negated && this->operatorType == OP_AND))
					&& useCount > 1)
				{
					useList[i] = false;
					--useCount;
				}
				else
				if ((lit1->negated && this->operatorType == OP_AND) ||
					(!lit1->negated && this->operatorType == OP_OR))
				{
					for (unsigned int j = 0; j < resList.size(); ++j)
					{
						if (i != j)
						{
							useList[j] = false;
							--useCount;
						}
					}
					break;
				}
			}
		}
	}

	// create result object
	BooleanOperator * res = new BooleanOperator(this->operatorType, NULL, this->negated);
	for (unsigned int i = 0; i < resList.size(); ++i)
	{
		if (useList[i])
		{
			res->addNewOperand(resList[i]);
			resList[i]->setParent(res);
		}
		else
		{
			delete resList[i];
		}
	}
	return res;
}

unsigned int BooleanOperator::recalculateNumberOfElements()
{
	numberOfSubtreeElements = 0;
	for (std::vector<BooleanFormula*>::iterator it = operands.begin();
			 it != operands.end(); ++it)
	{
		// recursive calculation of element count
		if (typeid(**it) == typeid(BooleanOperator))
		{
			numberOfSubtreeElements += ((BooleanOperator *)(*it))->recalculateNumberOfElements() + 1;
		}
		else
			++numberOfSubtreeElements;
	}
	return numberOfSubtreeElements;
}

unsigned int BooleanOperator::getNumberOfSubtreeElements() const
{
	return numberOfSubtreeElements;
}

void BooleanOperator::addNewOperand(BooleanFormula * operand)
{
	updateNumberOfElements(operand->getNumberOfSubtreeElements() + 1);
	operands.push_back(operand);
}

void BooleanOperator::deleteOperand(unsigned int index)
{
	updateNumberOfElements(-((int)operands[index]->getNumberOfSubtreeElements()) - 1);
	delete operands[index];
	operands.erase(operands.begin() + index);
}

void BooleanOperator::updateNumberOfElements(int additionalElements)
{
	numberOfSubtreeElements += additionalElements;
	if (parent != NULL)
		parent->updateNumberOfElements(additionalElements);
}

void BooleanOperator::setParent(BooleanOperator * parent)
{
	this->parent = parent;
}

TreePosition BooleanOperator::getElementByIndex(unsigned int index, unsigned int &runningIndex)
{
	for (unsigned int i = 0; i < operands.size(); ++i)
	{
		if (runningIndex == index)
		// the current element is the desired one
			return TreePosition(this,i,operands[i]);
		++runningIndex;

		if (typeid(*operands[i]) == typeid(BooleanOperator))
		{
			BooleanOperator * obj = (BooleanOperator *)(operands[i]);
			if (runningIndex + obj->getNumberOfSubtreeElements() - 1 >= index)
			// the element is in the current subtree
			{
				return obj->getElementByIndex(index,runningIndex);
			}
			else
			// the element is not in the current subtree
				runningIndex += obj->getNumberOfSubtreeElements();
		}
	}
	// the index was too large
	return TreePosition(NULL,0,NULL);
}

void BooleanOperator::replaceElementByIndex(unsigned int index, unsigned int &runningIndex, BooleanFormula * newElement)
{

	for (unsigned int i = 0; i < operands.size(); ++i)
	{
		if (runningIndex == index)
		// we found the element to replace
		{
			// recalculate element count
			updateNumberOfElements(-((int)operands[i]->getNumberOfSubtreeElements())
								   + ((int)newElement->getNumberOfSubtreeElements()));
			//replace element
			delete operands[i];
			operands[i] = newElement;
			newElement->setParent(this);
			return;
		}
		++runningIndex;
		if (typeid(*operands[i]) == typeid(BooleanOperator))
		{
			BooleanOperator * obj = (BooleanOperator *)(operands[i]);
			if (runningIndex + obj->getNumberOfSubtreeElements() - 1 >= index)
			// the element is in the current subtree
			{
				obj->replaceElementByIndex(index,runningIndex,newElement);
				return;
			}
			else
			// the element is not in the current subtree
				runningIndex += obj->getNumberOfSubtreeElements();
		}
	}
}

void BooleanOperator::getLiterals(LiteralMap &literals)
{
	for (std::vector<BooleanFormula *>::iterator it = operands.begin();
				 it != operands.end(); ++it )
	{
		if (typeid(**it) == typeid(BooleanOperator))
		// recursive call in a subtree
			((BooleanOperator *)*it)->getLiterals(literals);
		else
		{
			Literal * lit = (Literal *)*it;
			if (literals.count(lit) == 0)
			// initialize previously unseen literal
				literals[lit] = 1;
			else
			// increment counter of a literal
				++literals[lit];
		}
	}
}

BooleanOperator::BooleanOperator(OperatorType operatorType, BooleanOperator * parent, bool negated)
: BooleanFormula()
{
	this->negated = negated;
	this->operatorType = operatorType;
	this->numberOfSubtreeElements = 0;
	this->parent = parent;

	if (parent != 0)
	{
		parent->addNewOperand(this);
	}
}

BooleanOperator::~BooleanOperator()
{
	for (std::vector<BooleanFormula *>::iterator it = operands.begin();
		 it != operands.end(); ++it )
	{
		delete *it;
	}
}





bool Literal::evaluate(const bool * variableSet) const
	{
		if (literalIndex < 0)
			return !negated;

		if (negated)
			return !variableSet[literalIndex];
		else
			return variableSet[literalIndex];
	}

std::string Literal::toString(const std::string * variableNames) const
{
	if (literalIndex == -1)
	{
		if (negated)
			return "0";
		else
			return "1";
	}

	if (literalIndex == -2)
		return "UNKNOWN";

	if (negated)
		return "!" + variableNames[literalIndex];
	else
		return variableNames[literalIndex];
}

BooleanFormula * Literal::copy(BooleanOperator * parent)
{
	return new Literal(literalIndex, parent, negated);
}

unsigned int Literal::getNumberOfSubtreeElements() const
{
	return 0;
}

void Literal::setParent(BooleanOperator * parent)
{

}

Literal::Literal(int literalIndex, BooleanOperator * parent, bool negated)
: BooleanFormula()
{
	this->negated = negated;
	this->literalIndex = literalIndex;

	if (parent != 0)
	{
		parent->addNewOperand(this);
	}
}

bool operator==(BooleanFormula const& lhs, BooleanFormula const& rhs)
{
	if (typeid(lhs) != typeid(rhs))
	// comparison of an operator and a literal
		return false;

	if (typeid(lhs) == typeid(Literal))
	// comparison of two literals
	{
		Literal * l1 = ((Literal *) &lhs);
		Literal * l2 = ((Literal *) &rhs);
		return (l1->literalIndex == l2->literalIndex && l1->negated == l2->negated);
	}
	else
	// comparison of two operators
	{
		BooleanOperator * op1 = ((BooleanOperator *) &lhs);
		BooleanOperator * op2 = ((BooleanOperator *) &rhs);

		if (op1->operands.size() != op2->operands.size() || op1->operatorType != op2->operatorType
			|| op1->negated != op2->negated)
		// type or number of operands to not agree
			return false;

		// compare all operands
		for (std::vector<BooleanFormula *>::const_iterator it1 = op1->operands.begin(), it2 = op2->operands.begin();
			it1 != op1->operands.end(); ++it1, ++it2)
		{
			if (!(**it1 == **it2))
				return false;
		}
		return true;
	}
}




TreePosition BooleanTree::getElementByIndex(unsigned int index) const
{
	if (index == 0)
	// the root element is requested
		return TreePosition(NULL,0,root);
	else
	{
		unsigned int runningIndex = 0;
		if (typeid(*root) == typeid(BooleanOperator))
		// look up element in subtree
			return ((BooleanOperator *)root)->getElementByIndex(index-1, runningIndex);
		else
		// the index was too large
			return TreePosition(NULL,0,NULL);
	}
}

TreePosition BooleanTree::getNextElement(const TreePosition &previous) const
{
	if (typeid(*previous.element) == typeid(BooleanOperator))
	// the previous element was an operator => return its first element
	{
		BooleanOperator * op = (BooleanOperator *)previous.element;
		if (op->operands.size() > 0)
			return TreePosition(op,0,op->operands[0]);
	}

	if (previous.parent == NULL)
	// no next element if root is a literal
		return TreePosition(NULL,0,NULL);

	BooleanOperator * parent = previous.parent;
	int index = previous.index;

	while ((int)parent->operands.size() < index + 2)
	{
		BooleanOperator * old = parent;
		// go one level up in the tree
		parent = parent->parent;

		if (parent == NULL)
		// no next element, as root is reached
			return TreePosition(NULL,0,NULL);

		// determine index of parent in "grandparent"
		for (index = 0; index < (int)parent->operands.size() && parent->operands[index] != old; ++index);

		if (index == (int)parent->operands.size())
		// parent was not found (this should not occur)
			return TreePosition(NULL,0,NULL);

	}
	return TreePosition(parent,index+1,parent->operands[index+1]);
}

TreePosition BooleanTree::getNextOperator(const TreePosition &previous) const
{
	if (typeid(*previous.element) == typeid(BooleanOperator))
	{
		BooleanOperator * op = (BooleanOperator *)previous.element;

		unsigned int nextElem;

		// look for an operator among the operands of the current operator
		for (nextElem = 0; nextElem < op->operands.size(); ++nextElem)
		{
			if (typeid(*op->operands[nextElem]) == typeid(BooleanOperator))
				break;
		}

		if (nextElem != op->operands.size())
			return TreePosition(op,0,op->operands[nextElem]);
	}

	if (previous.parent == NULL)
		return TreePosition(NULL,0,NULL);
	BooleanOperator * parent = previous.parent;
	int index = previous.index, nextIndex;

	while (true)
	{
		// look for operators among the children of the parent operator
		for (nextIndex = index + 1; nextIndex < (int)parent->operands.size(); ++nextIndex)
		{
			if (typeid(*parent->operands[nextIndex]) == typeid(BooleanOperator))
				break;
		}

		if (nextIndex != (int)parent->operands.size())
			break;

		BooleanOperator * old = parent;
		// go one level up in the tree
		parent = parent->parent;

		if (parent == NULL)
		// the root is reached => there is no next operator
			return TreePosition(NULL,0,NULL);


		for (index = 0; index < (int)parent->operands.size() && parent->operands[index] != old; ++index);

		if (index == (int)parent->operands.size())
			return TreePosition(NULL,0,NULL);

	}
	return TreePosition(parent,nextIndex,parent->operands[nextIndex]);
}

void BooleanTree::replaceElement(TreePosition position, BooleanFormula * newElement, bool freeOld)
{
	if (position.parent == NULL)
	// exchange the root of the tree
	{
		if (freeOld)
			delete root;
		root = newElement;
	}
	else
	// exchange an element at a lower level
	{
		// recalculate number of elements in the subtree
		position.parent->updateNumberOfElements(-((int)position.element->getNumberOfSubtreeElements())
												+ ((int)newElement->getNumberOfSubtreeElements()));
		if (freeOld)
			delete position.parent->operands[position.index];

		// replace element
		position.parent->operands[position.index] = newElement;
		newElement->setParent(position.parent);
	}
	touch();
}

unsigned int BooleanTree::getNumberOfElements() const
{
	if (root == NULL)
		return 0;
	else
		return root->getNumberOfSubtreeElements() + 1;
}

LiteralMap * BooleanTree::getLiterals()
{
	if (modified)
	// recalculate the cache
	{
		if (root != NULL)
		{
			cachedLiterals.clear();
			if (typeid(*root) == typeid(BooleanOperator))
				((BooleanOperator *)root)->getLiterals(cachedLiterals);
			else
			{
				Literal * lit = (Literal *)root;
				cachedLiterals[lit] = 1;
			}
		}
		modified = false;
	}

	// return the cache
	return &cachedLiterals;
}

std::string BooleanTree::toString(const std::string * variableNames) const
{
	if (root != NULL)
		return root->toString(variableNames);
	else
		return "UNKNOWN";
}

bool BooleanTree::evaluate(const bool * variableSet) const
{
	if (root == NULL)
	// empty tree
		return false;
	else
		return root->evaluate(variableSet);
}


BooleanTree * BooleanTree::copy()
{
	if (root == NULL)
		return new BooleanTree(NULL, this->uncertainty);
	else
		return new BooleanTree(root->copy(NULL), this->uncertainty);
}

void BooleanTree::simplify()
{
	if (root != NULL && typeid(*root) == typeid(BooleanOperator))
	{
		BooleanOperator * r = ((BooleanOperator *)root);

		// calculate a simplified copy of the root operator
		r = (BooleanOperator *)r->getSimplifiedCopy();

		// replace the root
		delete root;
		root = r;

		if (r->operands.size() == 1)
		// if the root operator consists of only one element,
		// replace it with this element
		{
			root = r->operands[0]->copy(NULL);
			root->negated = (!r->negated && root->negated) || (r->negated && !root->negated);
			root->setParent(NULL);
			delete r;
		}
	}
	touch();
}


BooleanTree::BooleanTree(BooleanFormula * root, double uncertainty)
{
	this->root = root;
	this->modified = true;
	this->uncertainty = uncertainty;
}

BooleanTree::~BooleanTree()
{
	if (root != NULL)
		delete root;
}

void BooleanTree::normalizeOperator(BooleanOperator ** op)
{
	BooleanOperator * opToNormalize = *op;
	if (opToNormalize->operatorType == OP_XOR)
	// convert to AND/OR
	{

		BooleanOperator * opOr = new BooleanOperator(OP_OR);
		//BooleanOperator * root = opOr;
		BooleanFormula * op1 = opToNormalize->operands[0]->copy();
		BooleanFormula * op2 = opToNormalize->operands[1]->copy();

		for (unsigned int i = 1; i < opToNormalize->operands.size(); ++i)
		{
			BooleanOperator * opAnd1 = new BooleanOperator(OP_AND,opOr);
			BooleanOperator * opAnd2 = new BooleanOperator(OP_AND,opOr);

			opAnd1->addNewOperand(op1);
			op1->setParent(opAnd1);
			opAnd1->addNewOperand(op2);
			op2->setParent(opAnd1);

			BooleanFormula * op12 = op1->copy(opAnd2);
			op2->copy(opAnd2);

			op2->negated = !op2->negated;
			op12->negated = !op12->negated;
			if (i < opToNormalize->operands.size() -1)
			{
				op1 = opOr;
				op2 = opToNormalize->operands[i+1]->copy();
				opOr = new BooleanOperator(OP_OR);
			}
		}
		delete *op;
		*op = opOr;
		opToNormalize = opOr;
	}
	if (opToNormalize->negated)
	// apply De Morgan's rule
	{
		opToNormalize->negated = false;
		if (opToNormalize->operatorType == OP_AND)
			opToNormalize->operatorType = OP_OR;
		else
			opToNormalize->operatorType = OP_AND;

		for (std::vector<BooleanFormula*>::iterator it = opToNormalize->operands.begin();
					 it != opToNormalize->operands.end(); ++it)
		{
			(*it)->negated = !(*it)->negated;
		}
	}

	// recursively normalize subtrees
	for (std::vector<BooleanFormula*>::iterator it = opToNormalize->operands.begin();
						 it != opToNormalize->operands.end(); ++it)
	{
				if (typeid(**it) == typeid(BooleanOperator))
				{
					BooleanOperator * operand = (BooleanOperator *)*it;
					normalizeOperator(&operand);
					*it = operand;
				}
	}

	// obtain a unique order by sorting the subtrees
	sort(opToNormalize->operands.begin(),opToNormalize->operands.end(),sortFormula());
}



BooleanTree * BooleanTree::getNormalizedCopy()
{
	if (root == NULL)
		return new BooleanTree(NULL, this->uncertainty);
	if (typeid(*root) == typeid(BooleanOperator))
	{
		// first, simplify the tree
		BooleanOperator * res = (BooleanOperator *)(root->getSimplifiedCopy());

		// now apply a normalization to the simplified tree
		normalizeOperator(&res);

		return new BooleanTree(res, this->uncertainty);
	}
	else
		return new BooleanTree(root->copy(), this->uncertainty);
}

void BooleanTree::touch()
{
	this->modified = true;
}

bool operator==(BooleanTree const& lhs, BooleanTree const& rhs)
{
	if (lhs.root == NULL && rhs.root == NULL)
		return true;
	if (lhs.root == NULL || rhs.root == NULL)
		return false;

	return (*lhs.root == *rhs.root);
}
