/*
 * formulaparser.h
 * Defines a parser for Boolean expressions
 *
 *  Created on: Apr 27, 2010
 *      Author: muessel
 * Edited 2020 by schwab
 */

#ifndef FORMULAPARSER_H_
#define FORMULAPARSER_H_

#include "booleanformula.h"
#include <vector>
#include <exception>
#include <sstream>

/**
 * The types of symbols recognized by the scanner:
 * Opening and closing brackets, negation operators,
 * binary operators and string constants (e.g. variables)
 */
enum SymbolType
{
	ST_OPENBRACKET,
	ST_CLOSEBRACKET,
	ST_NEGATION,
	ST_OPERATOR,
	ST_CONSTANT
};

/**
 * A class representing a scanner symbol
 */
class Symbol
{
public:
	// the type of the symbol
	SymbolType type;

	// additional information on a symbol,
	// e.g. the index of the associated variable
	int additionalInfo;

	// the start position of the symbol in the input string
	unsigned int start;

	/**
	 * Creates a new symbol with the following parameters:
	 * <type>: the type of the symbol
	 * <additionalInfo: additional information on a symbol,
	 * e.g. the index of the associated variable
	 * <start>: the start position of the symbol in the input string
	 */
	Symbol(SymbolType type, int additionalInfo,
		   unsigned int start)
	{
		this->type = type;
		this->additionalInfo = additionalInfo;
		this->start = start;
	}

};

/**
 * An exception class for parser exceptions
 */
class ParseException: public std::exception
{
private:
	// the error message
	std::string message;
public:

	virtual const char* what() const throw()
	{
		return message.c_str();
	}

	/**
	 * Creates a new parser exception highlighting
	 * the position of the error in the input.
	 * <error> is the error message.
	 * <input> is the input string.
	 * <position> denotes the position of the error in the input string.
	 */
	ParseException(std::string error, std::string input, unsigned int position)
	{
		std::ostringstream out;
		out << "Parsing error: " << error << "\n";
		out << input << "\n";
		for (unsigned int i = 0; i < position; ++i)
			out << " ";
		out << "^";
		this->message = out.str();
	}

	virtual ~ParseException() throw ()
	{

	}
};

/**
 * A parser class that parses Boolean expressions
 */
class FormulaParser
{
private:
	// The names of known variables
	std::vector<std::string> varNames;

	// The currently processed input string
	std::string currentInput;

	// A list of symbols identified in the input
	std::vector<Symbol *> symbols;

	// The index of the currently processed symbol in the list
	unsigned int symbolIndex;

	/**
	 * Removes all symbols
	 */
	inline void clearSymbols()
	{
		for (std::vector<Symbol *>::iterator it = symbols.begin(); it != symbols.end(); ++it)
			delete *it;
		symbols.clear();
	}

	/**
	 * Previews the next symbol to process without
	 * moving the symbol pointer
	 */
	inline Symbol * previewSym()
	{
		if (symbolIndex == symbols.size() - 1)
			throw ParseException("Unexpected end!", currentInput, (unsigned int)currentInput.size() - 1);
		return symbols[symbolIndex+1];
	}

	/**
	 * Reads the next symbol in the symbol list
	 * and move the symbol pointer.
	 */
	inline Symbol * popSym()
	{
		if (symbolIndex == symbols.size() - 1)
			throw ParseException("Unexpected end!", currentInput, (unsigned int)currentInput.size() - 1);
		return symbols[++symbolIndex];
	}

	/**
	 * Returns the last symbol read.
	 */
	inline Symbol * currentSym()
	{
		if (symbolIndex == symbols.size())
			throw ParseException("Unexpected end!", currentInput, (unsigned int)currentInput.size() - 1);
		return symbols[symbolIndex];
	}

	/**
	 * Determines whether all symbols have been read.
	 */
	inline bool eof()
	{
		return (symbolIndex == symbols.size() - 1);
	}

	/**
	 * Scans the input string for symbols, and fills
	 * the symbol list.
	 */
	void scan();

	/**
	 * Parses a new operator/expression
	 */
	BooleanOperator * op(bool negated);

	/**
	 * Parses a new literal
	 */
	Literal * literal(bool negated);

public:
	/**
	 * Creates a new formula parser that
	 * recognizes the variables specified in <varNames>.
	 */
	FormulaParser(std::vector<std::string> &varNames);

	/**
	 * Parses the input string <input>.
	 */
	BooleanFormula * parse(std::string input);

	virtual ~FormulaParser();
};

#endif /* FORMULAPARSER_H_ */
