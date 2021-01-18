/*
 * formulaparser.cpp
 * Defines a parser for Boolean expressions
 *
 *  Created on: Apr 27, 2010
 *      Author: muessel
 * 	Edited 2020 by schwab
 */

#include "formulaparser.h"
#include "helpers.h"

#include <string>
#include <algorithm>
#include <cctype>

#include <iostream>

using namespace std;

void FormulaParser::scan()
{
	symbolIndex = 0;
	clearSymbols();
	unsigned int pos = 0;
	while (pos < currentInput.size())
	{
		// first, check for one-character symbols
		if (currentInput[pos] == ' ')
		// ignore whitespace
			++pos;
		else
		if (currentInput[pos] == '(')
		// recognize opening bracket
		{
			symbols.push_back(new Symbol(ST_OPENBRACKET,0,pos));
			++pos;
		}
		else
		if (currentInput[pos] == ')')
		// recognize closing bracket
		{
			symbols.push_back(new Symbol(ST_CLOSEBRACKET,0,pos));
			++pos;
		}
		else
		if (currentInput[pos] == '&')
		// recognize AND operator
		{
			symbols.push_back(new Symbol(ST_OPERATOR,OP_AND,pos));
			++pos;
		}
		else
		if (currentInput[pos] == '|')
		// recognize OR operator
		{
			symbols.push_back(new Symbol(ST_OPERATOR,OP_OR,pos));
			++pos;
		}
		else
		if (currentInput[pos] == '!')
		// recognize NOT operator
		{
			symbols.push_back(new Symbol(ST_NEGATION,0,pos));
			++pos;
		}
		else
		if (currentInput[pos] == '1')
		// recognize TRUE
		{
			symbols.push_back(new Symbol(ST_CONSTANT,-1,pos));
			++pos;
		}
		else
		if (currentInput[pos] == '0')
		// recognize FALSE
		{
			symbols.push_back(new Symbol(ST_CONSTANT,-2,pos));
			++pos;
		}
		else
		// no one-character symbol
		{
			// find next whitespace or delimiter
			unsigned int nextPos = pos + 1;
			while (nextPos < currentInput.size() && currentInput[nextPos] != ' '
				   && currentInput[nextPos] != ')' && currentInput[nextPos] != '('
				   && currentInput[nextPos] != '[')
				++nextPos;

			// extract symbol and convert it to upper case
			string nextSym = currentInput.substr(pos,nextPos - pos);
			std::transform(nextSym.begin(), nextSym.end(), nextSym.begin(), (int(*)(int))toupper);

			if (nextSym == "AND")
			// recognize alternative AND
			{
				symbols.push_back(new Symbol(ST_OPERATOR,OP_AND,pos));
			}
			else
			if (nextSym == "OR")
			// recognize alternative OR
			{
				symbols.push_back(new Symbol(ST_OPERATOR,OP_OR,pos));
			}
			else
			if (nextSym == "XOR")
			// recognize XOR
			{
				symbols.push_back(new Symbol(ST_OPERATOR,OP_XOR,pos));
			}
			else
			if (nextSym == "NOT")
			// recognize alternative NOT
			{
				symbols.push_back(new Symbol(ST_NEGATION,0,pos));
			}
			else
			if (nextSym == "TRUE")
			// recognize alternative TRUE
			{
				symbols.push_back(new Symbol(ST_CONSTANT,-1,pos));
			}
			else
			if (nextSym == "FALSE")
			// recognize alternative FALSE
			{
				symbols.push_back(new Symbol(ST_CONSTANT,-2,pos));
			}
			else
			// this must be a variable name
			{
				unsigned int litIdx;

				// lookup symbol in variable name list
				for (litIdx = 0; litIdx < varNames.size(); ++litIdx)
				{
					string geneName = varNames[litIdx];
					transform(geneName.begin(), geneName.end(), geneName.begin(), (int(*)(int))toupper);
					if (geneName == nextSym)
						break;
				}

				if (litIdx == varNames.size())
				// variable was not found
					throw ParseException("Unknown symbol " + nextSym + "!",currentInput,pos);
				else
				// add variable index as a symbol
					symbols.push_back(new Symbol(ST_CONSTANT,(int)litIdx,pos));
			}
			pos = nextPos;
		}
	}
}

BooleanOperator * FormulaParser::op(bool negated)
{
	// start with current symbol
	Symbol * current = popSym();
	bool childNegated = false;

	BooleanOperator * res = new BooleanOperator(OP_AND,NULL,negated);
	BooleanFormula * child;

	int opType = -1;
	Symbol * next;
	do
	{
		switch (current->type)
		{
			case ST_OPENBRACKET:
				// a new operator starts here
				child = op(childNegated);
				res->addNewOperand(child);
				child->setParent(res);
				childNegated = false;
				break;
			case ST_CONSTANT:
				// this is a literal
				// check for an operator or closing bracket
				next = previewSym();
				if (next->type != ST_OPERATOR && next->type != ST_CLOSEBRACKET)
					throw ParseException("Operator or \")\" expected!",
										 currentInput,next->start);

				child = literal(childNegated);
				res->addNewOperand(child);
				child->setParent(res);
				childNegated = false;
				break;
			case ST_NEGATION:
				// this is the negation of an operand
				childNegated = !childNegated;
				// check for operand
				next = previewSym();
				if (next->type != ST_CONSTANT && next->type != ST_NEGATION
					&& next->type != ST_OPENBRACKET)
					throw ParseException("Operand expected after negation!",currentInput,next->start);
				break;
			case ST_OPERATOR:
				// check whether the operand fits with other operands in the same bracket
				if (opType != -1 && opType != current->additionalInfo)
					throw ParseException("Multiple different operators in an expression without proper bracketing!",
								 currentInput,current->start);

				// check for operand
				next = previewSym();
				if (next->type != ST_CONSTANT && next->type != ST_NEGATION
					&& next->type != ST_OPENBRACKET)
					throw ParseException("Operand expected after operator!",currentInput,next->start);
				opType = current->additionalInfo;
				break;
			default:
				throw ParseException("Unexpected symbol!",currentInput,current->start);
		}
		// proceed to next symbol
		current = popSym();
	}
	while (current->type != ST_CLOSEBRACKET);

	if (opType != -1)
		res->operatorType = (OperatorType)opType;

	return res;
}

Literal * FormulaParser::literal(bool negated)
{
	Symbol * current = currentSym();
	Literal * res;
	if (current->additionalInfo == -1)
	// create a TRUE constant
	{
		res = new Literal(-1, NULL, negated);
	}
	else
	if (current->additionalInfo == -2)
	// create a FALSE constant
	{
		res = new Literal(-1, NULL, !negated);
	}
	else
	// create a variable literal
	{
		res = new Literal(current->additionalInfo, NULL, negated);
	}

	return res;
}

BooleanFormula * FormulaParser::parse(std::string input)
{
	// remove leading and trailing blanks
	trim(input);
	this->currentInput = input;

	if (input.size() == 0)
		return NULL;

	// convert to upper case
	string tmp = input;
	transform(tmp.begin(), tmp.end(), tmp.begin(), (int(*)(int))toupper);

	if (tmp == "UNKNOWN" || tmp == "(UNKNOWN)")
	// special case: unknown formula
		return new Literal(-2,NULL,0);

	// check whether the expression is surrounded by brackets
	bool brackets = (input[0] == '(' && input[input.size()-1] == ')');
	if (brackets && input.size() >= 4)
	{
		unsigned int bracketCount = 0;
		for (unsigned int i = 1; i < input.size() - 1; ++i)
		{
			if (input[i] == '(')
				++bracketCount;
			else
			if (input[i] == ')' && bracketCount > 0)
				--bracketCount;
		}
		if (bracketCount != 0)
			brackets = false;
	}
	// add brackets if necessary
	if (!brackets)
		this->currentInput = '(' + this->currentInput + ')';

	// scan symbols
	scan();
	if (symbols.size() == 0)
		return NULL;

	Symbol * current = currentSym();
	bool negated = false;
	BooleanFormula * res = NULL;
	Symbol * next;
	do
	{
		switch (current->type)
		{
			case ST_OPENBRACKET:
				// this is the beginning of an operator
				res = op(negated);
				negated = false;
				break;
			case ST_CONSTANT:
				// this is the beginning of a literal
				res = literal(negated);
				negated = false;
				break;
			case ST_NEGATION:
				// this is the beginning of a negated operand
				negated = !negated;
				next = previewSym();
				if (next->type != ST_CONSTANT && next->type != ST_NEGATION
					&& next->type != ST_OPENBRACKET)
					throw ParseException("Operand expected after negation!",currentInput,next->start);
				break;
			default:
				throw ParseException("Unexpected symbol!",currentInput,current->start);
		}
		if (!eof())
			current = popSym();
	}
	while (res == NULL);

	if (!eof())
	// the expression is longer than expected
		throw ParseException("Unexpected symbol!",currentInput,previewSym()->start);

	if (typeid(*res) == typeid(BooleanOperator))
	{
		BooleanOperator * tmp = (BooleanOperator *) res;
		if (tmp->operands.size() == 1 && (typeid(*tmp->operands[0]) == typeid(BooleanOperator)))
		// if operator only consists of one operand, replace it by operand
		{
			res = tmp->operands[0]->copy(NULL);
			delete tmp;
		}
	}
	return res;
}

FormulaParser::FormulaParser(std::vector<std::string> &varNames)
{
	this->varNames.assign(varNames.begin(),varNames.end());
}

FormulaParser::~FormulaParser()
{
	for (std::vector<Symbol *>::iterator it = symbols.begin(); it != symbols.end(); ++it)
		delete *it;
}
