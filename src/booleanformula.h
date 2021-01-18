/*
 * booleanformula.h
 * Provides classes for the management of
 * Boolean expressions
 *
 *  Created on: Apr 14, 2010
 *      Author: muessel
 */

#ifndef BOOLEANFORMULA_H_
#define BOOLEANFORMULA_H_

#include <string>
#include <sstream>
#include <vector>
#include <typeinfo>
#include <iostream>
#include <map>

class BooleanOperator;
class Literal;

struct CompareLiterals;

/**
 * Base class for Boolean expressions
 */
class BooleanFormula
{
public:
	/**
	 * Is the expression negated?
	 */
	bool negated;

	/**
	 * Evaluates the expression on a given set of variable values
	 * <variableSet>.
	 * Returns the value of the expression.
	 */
	virtual bool evaluate(const bool * variableSet) const = 0;

	/**
	 * Returns a string representation of the formula based on the variable names
	 * in <variableNames>.
	 */
	virtual std::string toString(const std::string * variableNames) const = 0;

	/**
	 * Creates an exact copy of this object, and adds it to <parent> if supplied.
	 */
	virtual BooleanFormula * copy(BooleanOperator * parent = NULL) = 0;

	/**
	 * Changes the parent expression of the object to <parent>
	 */
	virtual void setParent(BooleanOperator * parent) = 0;

	/**
	 * Returns the number of elements in the subtree denoted by the expression.
	 * Here, each operator and each literal is counted as one element.
	 */
	virtual unsigned int getNumberOfSubtreeElements() const = 0;

	/**
	 * Simplifies an expression by removing duplicates, resolving tautologies, etc.,
	 * and returns the simplified version as a copy.
	 */
	virtual BooleanFormula * getSimplifiedCopy()
	{
		return this->copy(NULL);
	}

	virtual ~BooleanFormula()
	{

	}
};

/**
 * Type of Boolean operators
 * -- XOR is currently not used
 */
enum OperatorType
{
	OP_AND,
	OP_OR,
	OP_XOR
};

/**
 * Structure defining a position in an expression tree
 * by its parent, the index of the sub-node and
 * the corresponding element itself
 */
class TreePosition
{
public:
	// the parent expression
	BooleanOperator * parent;

	// the index of the element in the parent operator
	unsigned int index;

	// the element itself
	BooleanFormula * element;

	/**
	 * Creates a new TreePosition with the following parameters:
	 * <parent>: the parent expression
	 * <index>: the index of the element in the parent operator
	 * <element>: the element itself
	 */
	TreePosition(BooleanOperator * parent, unsigned int index,	BooleanFormula * element)
	{
		this->parent = parent;
		this->index = index;
		this->element = element;
	}
};

/**
 * A class representing a single literal
 * in a formula
 */
class Literal : public BooleanFormula
{
public:
	// the index of the literal, -1 = constant
	int literalIndex;

	virtual bool evaluate(const bool * variableSet) const;

	virtual std::string toString(const std::string * variableNames) const;

	virtual BooleanFormula * copy(BooleanOperator * parent = NULL);

	virtual unsigned int getNumberOfSubtreeElements() const;

	virtual void setParent(BooleanOperator * parent);

	/**
	 * Creates a new Literal with the following parameters:
	 * <literalIndex>: the index of the literal, use -1 for constant TRUE/FALSE depending on <negated>
	 * <parent>: the expression in which the literal is included
	 * <negated>: is this literal negated?
	 */
	Literal(int literalIndex, BooleanOperator * parent = NULL, bool negated = false);

};

/**
 * Structure needed for sorting literals
 */
struct CompareLiterals
{

	bool operator() (const Literal * l1, const Literal * l2) const
	{
		if (l1 == NULL)
			return true;
		if (l2 == NULL)
			return false;

		return (l1->literalIndex * 2 + l1->negated < l2->literalIndex * 2 + l2->negated);
	}

};

/**
 * A map holding the unique literals and their frequencies in a formula,
 * as returned by getLiterals()
 */
typedef std::map<Literal *, unsigned int, CompareLiterals> LiteralMap;

/**
 * Class representing a Boolean operation on other Boolean expressions
 */
class BooleanOperator : public BooleanFormula
{
private:
	unsigned int numberOfSubtreeElements;
public:
	// The parent expression
	BooleanOperator * parent;

	// The type of the operator, i.e. AND/OR/XOR
	OperatorType operatorType;

	// The operands of the operator (usually >=2)
	std::vector<BooleanFormula *> operands;

	virtual bool evaluate(const bool * variableSet) const;

	virtual std::string toString(const std::string * variableNames) const;

	virtual BooleanFormula * copy(BooleanOperator * parent = NULL);

	virtual BooleanFormula * getSimplifiedCopy();

	virtual unsigned int getNumberOfSubtreeElements() const;

	/**
	 * Adds a new operand <operand> to the expression
	 */
	void addNewOperand(BooleanFormula * operand);

	/**
	 * Removes the <index>-th operand from the expression
	 */
	void deleteOperand(unsigned int index);

	/**
	 * Recalculate the internal counter of sub-elements by
	 * adding <additionalElements>
	 */
	void updateNumberOfElements(int additionalElements);

	virtual void setParent(BooleanOperator * parent);

	/**
	 * Retrieves the <index>th element of the expression by traversing
	 * the tree in a breadth-first search.
	 * <runningIndex> is an internal counter that is increased during the recursion
	 * to determine the position.
	 * Returns the element as a TreePosition object.
	 */
	TreePosition getElementByIndex(unsigned int index, unsigned int &runningIndex);

	/**
	 * Replaces the <index>th element by a new expression <newElement>. The element is determined
	 * as in getElementByIndex().
	 * <runningIndex> is an internal counter that is increased during the recursion
	 * to determine the position.
	 */
	void replaceElementByIndex(unsigned int index, unsigned int &runningIndex, BooleanFormula * newElement);

	/**
	 * Retrieves a map of all literals in the expression and their frequencies.
	 * Here, the keys of the map denote the literals, and the values denote the frequencies.
	 */
	void getLiterals(LiteralMap &literals);

	/**
	 * Recalculates the internal counter of sub-elements in the tree.
	 */
	unsigned int recalculateNumberOfElements();

	/**
	 * Creates a new BooleanOperator with parameters
	 * <operatorType>: The type of the operator, i.e. AND/OR/XOR
	 * <parent>: The parent expression
	 * <negated>: Is the operator negated or not?
	 */
	BooleanOperator(OperatorType operatorType, BooleanOperator * parent = NULL, bool negated = false);

	virtual ~BooleanOperator();

};

/**
 * Test for equality of two Boolean expressions
 * -- this usually requires normalization of both expressions
 */
bool operator==(BooleanFormula const& lhs, BooleanFormula const& rhs);

/**
 * Container class for a Boolean expression that provides
 * additional top-level management and analysis functions
 * and interfaces to recursive calls of formula objects
 */
class BooleanTree
{
private:
	friend bool operator==(BooleanTree const& lhs, BooleanTree const& rhs);

	// the root element of the container
	BooleanFormula * root;

	// has the formula been modified,
	// such that the literal map must be recalculated?
	bool modified;

	// A cached map of the literals in the tree,
	// which is only recalculated if necessary
	LiteralMap cachedLiterals;

	// the uncertainty of the formula
	// (0 = no uncertainty, 1 = high uncertainty)
	double uncertainty;

	/**
	 * Internal recursive method to return a normalized version
	 * if <op> in negation normal form.
	 */
	void normalizeOperator(BooleanOperator ** op);
public:

	/**
	 * Returns the root element of the container
	 */
	inline BooleanFormula * getRoot()
	{
		return root;
	}

	/**
	 * Retrieves the uncertainty of the tree,
	 * which determines how frequently it is mutated
	 */
	inline double getUncertainty()
	{
		return uncertainty;
	}

	/**
	 * Sets the uncertainty of the tree,
	 * which determines how frequently it is mutated,
	 * to a new value <uncertainty>
	 * (0 = no uncertainty, 1 = high uncertainty)
	 */
	inline void setUncertainty(double uncertainty)
	{
		this->uncertainty = uncertainty;
	}

	/**
	 * Retrieves the <index>th element of the expression by traversing
	 * the tree in a breadth-first search.
	 * Returns the element as a TreePosition object.
	 */
	TreePosition getElementByIndex(unsigned int index) const;

	/**
	 * Retrieves the element following <previous> in a breadt-first search
	 * on the tree as a TreePosition object
	 */
	TreePosition getNextElement(const TreePosition &previous) const;

	/**
	 * Retrieves the operator following <previous> in a breadt-first search
	 * on the tree as a TreePosition object
	 */
	TreePosition getNextOperator(const TreePosition &previous) const;

	/**
	 * Replaces the element specified by <position> by a new expression <newElement>.
	 * If <freeOld> is true, the old element is deleted.
	 */
	void replaceElement(TreePosition position, BooleanFormula * newElement, bool freeOld = true);

	/**
	 * Returns the number of elements in the expression tree
	 */
	unsigned int getNumberOfElements() const;

	/**
	 * Returns a map of literals and their frequencies in the expression.
	 * Here, the keys denote the literals, and the values denote their frequencies
	 */
	LiteralMap * getLiterals();

	/**
	 * Returns a string representation of the formula based on the variable names
	 * in <variableNames>.
	 */
	std::string toString(const std::string * variableNames) const;

	/**
	 * Evaluates the expression on a given set of variable values
	 * <variableSet>.
	 * Returns the value of the expression.
	 */
	bool evaluate(const bool * variableSet) const;

	/**
	 * Returns an exact copy of the container and its expression.
	 */
	BooleanTree * copy();

	/**
	 * Replaces the expression in the container by a simplified copy
	 * by eliminating duplicates, resolving tautologies, etc.
	 */
	void simplify();

	/**
	 * Returns a normalized version of the expression by transforming it
	 * into the negation normal form.
	 */
	BooleanTree * getNormalizedCopy();

	/**
	 * Marks this object as modified, which means that the cached literal
	 * map has to be recalculated.
	 */
	void touch();

	/**
	 * Creates a new BooleanTree with the expression specified by <root>
	 * and the uncertainty value <uncertainty>
	 */
	BooleanTree(BooleanFormula * root, double uncertainty=1.0);

	virtual ~BooleanTree();

};
#endif /* BOOLEANFORMULA_H_ */

/**
 * Test for equality of two Boolean expression containers
 * -- this usually requires normalization of both expressions
 */
bool operator==(BooleanTree const& lhs, BooleanTree const& rhs);
