#ifndef PARAMETERMANAGER_H_
#define PARAMETERMANAGER_H_

/**
 * This file defines a class that processes
 * the command line arguments of the program
 */

#include <map>
#include <sstream>
#include <vector>
#include <algorithm>

/**
 * A class that parses and validates the command line arguments
 * of a program and stores them in a typed form.
 */
class ParameterManager
{
protected:

	/**
	 * The data types supported by the manager
	 */
	enum paramType
	{
		STRING,
		INT,
		DOUBLE,
		VOID
	};

	/**
	 * A single parameter in the command line
	 */
	class Parameter
	{
	public:
		// a counter for the running parameter numbers
		static unsigned int counter;

		// a running number used to sort the parameters in std::map
		unsigned int number;

		// the parameter name
		std::string name;

		// a description of the parameter for the usage message
		std::string description;

		// indicates whether the parameter is optional or not
		bool optional;

		// indicates whether this parameter is satisfied in a
		// parsed command line, i.e. it is present or it is optional
		bool fulfilled;

		// indicates whether the parameter was found in a parsed
		// command line or not
		bool isSet;

		// the type of the parameter
		paramType type;

		// the pointers to the data,
		// depending on <type>
		union
	    {
			std::string * stringVal;
			int * intVal;
			double * doubleVal;
	    } value;

	    /**
	     * Creates a new string parameter with the option flag <name>. <optional> defines
	     * whether this flag must be supplied or not, and <description> is a short description
	     * of the parameter value for the usage message. <val> is the default value of the parameter.
	     */
	    Parameter(std::string name, std::string description, std::string &val, bool optional)
	    {
	    	this->number = counter++;
	    	this->name = name;
	    	this->description = description;
	    	this->type = STRING;
	    	this->value.stringVal = new std::string(val);
	    	this->optional = optional;
	    	this->isSet = false;
	    }

	    /**
	     * Creates a new int parameter with the option flag <name>. <optional> defines
	     * whether this flag must be supplied or not, and <description> is a short description
	     * of the parameter value for the usage message. <val> is the default value of the parameter.
	     */
	    Parameter(std::string name, std::string description, int &val, bool optional)
	    {
	    	this->number = counter++;
	    	this->name = name;
	    	this->description = description;
	    	this->type = INT;
	    	this->value.intVal = new int[1];
	    	*this->value.intVal = val;
	    	this->optional = optional;
	    	this->isSet = false;
	    }

	    /**
	     * Creates a new double parameter with the option flag <name>. <optional> defines
	     * whether this flag must be supplied or not, and <description> is a short description
	     * of the parameter value for the usage message. <val> is the default value of the parameter.
	     */
	    Parameter(std::string name, std::string description, double &val, bool optional)
	    {
	    	this->number = counter++;
	    	this->name = name;
	    	this->description = description;
	    	this->type = DOUBLE;
	    	this->value.doubleVal = new double[1];
	    	*this->value.doubleVal = val;
	    	this->optional = optional;
	    	this->isSet = false;
	    }

	    /**
	     * Creates a new flag parameter with the option flag <name>. <optional> defines
	     * whether this flag must be supplied or not. As this defines a Boolean flag,
	     * such parameters are usually optional - if they are present, they represent
	     * true, and if they are not present, they represent false.
	     */
	    Parameter(std::string name, bool optional)
	    {
	    	this->number = counter++;
	    	this->name = name;
	    	this->type = VOID;
	    	this->value.doubleVal = NULL;
	    	this->optional = optional;
	    	this->isSet = false;
	    }

	    virtual ~Parameter()
	    {
	    	switch(type)
	    	{
	    		case STRING:
	    			delete value.stringVal;
	    			break;
	    		case INT:
	    			delete [] value.intVal;
	    			break;
	    		case DOUBLE:
	    			delete [] value.doubleVal;
	    			break;
	    		case VOID:
	    			break;
	    	}
	    }

	    /**
	     * Returns the string value of this parameter (only if type is STRING)
	     */
	    std::string getString()
	    {
	    	return *value.stringVal;
	    }

	    /**
	     * Returns the int value of this parameter (only if type is INT)
	     */
	    int getInt()
	    {
	    	return *value.intVal;
	    }

	    /**
	     * Returns the double value of this parameter (only if type is DOUBLE)
	     */
	    double getDouble()
	    {
	    	return *value.doubleVal;
	    }

	};

	/**
	 * A comparison operator to sort the parameters by their
	 * running numbers in std::map
	 */
	struct sortParameters
	{
		bool operator() (const Parameter * p1, const Parameter * p2) const
		{
			return (p1->number < p2->number);
		}
	};

	// a map of main options - a main option is an option that determines
	// the set of parameters that come afterwards. These options are stored
	// in the value map
	std::map<std::string,std::map<std::string, Parameter *> *> mainOptions;

	// the parameter map of the active main option
	std::map<std::string, Parameter *> * activeOptionMap;

	// the active main option
	std::string activeMainOption;
public:
	ParameterManager()
	{
	}

	/**
	 * Adds a new main option. A main option is an option that
	 * determines the set of following parameters. <option> is
	 * the flag of the new main option.
	 */
	inline void addMainOption(std::string option)
	{
		std::map<std::string, Parameter *> * optionMap = new std::map<std::string, Parameter *>();
		mainOptions[option] = optionMap;
	}

	/**
	 * Adds a new string parameter with name <option> for main option <mainOption>. <valueDescription> is a short
	 * description that is printed in the usage message, and <value> determines the default value. <optional>
	 * specifies whether this flag can be left out.
	 */
	inline void addStringOption(std::string mainOption, std::string option, std::string valueDescription, std::string value, bool optional)
	{
		(*(mainOptions[mainOption]))[option] = new Parameter(option,valueDescription,value,optional);
	}

	/**
	 * Adds a new int parameter with name <option> for main option <mainOption>. <valueDescription> is a short
	 * description that is printed in the usage message, and <value> determines the default value. <optional>
	 * specifies whether this flag can be left out.
	 */
	inline void addIntOption(std::string mainOption, std::string option, std::string valueDescription, int value, bool optional)
	{
		(*(mainOptions[mainOption]))[option] = new Parameter(option,valueDescription,value,optional);
	}

	/**
	 * Adds a new double parameter with name <option> for main option <mainOption>. <valueDescription> is a short
	 * description that is printed in the usage message, and <value> determines the default value. <optional>
	 * specifies whether this flag can be left out.
	 */
	inline void addDoubleOption(std::string mainOption, std::string option, std::string valueDescription, double value, bool optional)
	{
		(*(mainOptions[mainOption]))[option] = new Parameter(option,valueDescription,value,optional);
	}

	/**
	 * Adds a new Boolean flag parameter with name <option> for main option <mainOption>. <optional>
	 * specifies whether this flag can be left out.
	 */
	inline void addVoidOption(std::string mainOption, std::string option, bool optional)
	{
		(*(mainOptions[mainOption]))[option] = new Parameter(option,optional);
	}

	/**
	 * Processes the <argc> parameters in <argv> and stores them in the internal
	 * parameter structures.  Afterwards, the parameters are analyzed, and the function
	 * returns true if all specified requirements are fulfilled or false otherwise
	 */
	bool analyzeCommandLine(int argc, char ** argv)
	{
		if (argc == 1)
			return false;
		activeMainOption = std::string(argv[1]);

		if (mainOptions.count(activeMainOption) == 0)
			return false;

		activeOptionMap = mainOptions[activeMainOption];

		for (std::map<std::string, Parameter *>::iterator it = activeOptionMap->begin(); it != activeOptionMap->end(); it++)
		{
			it->second->fulfilled = it->second->optional;
		}

		if (argc > 2)
		{
			char ** ptr = &argv[2];
			int numLeft = argc - 2;
			do
			{
				std::string param = std::string(ptr[0]);

				// check whether the parameter is known
				if (activeOptionMap->count(param) == 0)
					return false;

				Parameter * p = (*activeOptionMap)[param];

				// check if the parameter has been specified more than once
				if (p->isSet)
					return false;

				// read in parameter values
				switch (p->type)
				{
					case STRING:
						if (numLeft == 1)
							return false;
						*p->value.stringVal = std::string((++ptr)[0]);
						numLeft--;
						break;
					case INT:
						if (numLeft == 1)
							return false;
						*p->value.intVal =  atoi((++ptr)[0]);
						numLeft--;
						break;
					case DOUBLE:
						if (numLeft == 1)
							return false;
						*p->value.doubleVal =  atof((++ptr)[0]);
						numLeft--;
						break;
					case VOID:
						break;
				}
				p->fulfilled = true;
				p->isSet = true;
				numLeft--;
				if (numLeft > 0)
					ptr++;
			}
			while (numLeft > 0);
		}
		std::cout << "huhu" << std::endl;
		// check whether all mandatory parameters have been specified
		for (std::map<std::string, Parameter *>::iterator it = activeOptionMap->begin(); it != activeOptionMap->end(); it++)
		{
			std::cout << it->second->fulfilled << std::endl;
			if (!it->second->fulfilled)
				return false;
		}
		return true;
	}

	/**
	 * Prints out a usage message using <progname> as the name of the program.
	 */
	std::string usage(std::string progname)
	{
		// determine the program name without a path
		size_t slashPos = progname.find_last_of("/");

		if (slashPos == std::string::npos)
			slashPos = progname.find_last_of("\\");

		if (slashPos != std::string::npos)
			progname = progname.substr(slashPos+1);

		std::ostringstream res;
		res << "Usage: \n";

		for (std::map<std::string,std::map<std::string, Parameter *> *>::iterator it1 = mainOptions.begin(); it1 != mainOptions.end(); it1++)
		// iterate over main options
		{
			res << progname << " " << it1->first << " ";
			std::map<std::string, Parameter *> * current = it1->second;
			std::vector<Parameter *> parameters;
			for (std::map<std::string, Parameter *>::iterator it2 = current->begin(); it2 != current->end(); it2++)
			{
				parameters.push_back(it2->second);
			}
			std::sort(parameters.begin(),parameters.end(),sortParameters());
			for (std::vector<Parameter *>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++)
			// list the parameters for the current main option
			{
				Parameter * p = *it2;
				if (p->optional)
					res << "[";
				res << p->name;
				switch (p->type)
				{
					case STRING:
					case INT:
					case DOUBLE:
						res << " " << p->description;
						break;
					case VOID:
						break;
				}
				if (p->optional)
					res << "] ";
				else
					res << " ";
			}
			res << "\n";
		}
		return res.str();
	}

	/**
	 * Returns the main option that is active in the processed command line
	 */
	std::string getActiveMainOption()
	{
		return activeMainOption;
	}

	/**
	 * Retrieves the value of the string option <name>
	 */
	std::string getStringOption(std::string name)
	{
		if ((*activeOptionMap).count(name) == 0)
			throw "Option " + name + " not available!";
		return (*activeOptionMap)[name]->getString();
	}

	/**
	 * Retrieves the value of the int option <name>
	 */
	int getIntOption(std::string name)
	{
		if ((*activeOptionMap).count(name) == 0)
			throw "Option " + name + " not available!";
		return (*activeOptionMap)[name]->getInt();
	}

	/**
	 * Retrieves the value of the double option <name>
	 */
	double getDoubleOption(std::string name)
	{
		if ((*activeOptionMap).count(name) == 0)
			throw "Option " + name + " not available!";
		return (*activeOptionMap)[name]->getDouble();
	}

	/**
	 * Specifies whether the option <name> is set in the
	 * supplied command line. Possible for all types of
	 * options
	 */
	bool isOptionSet(std::string name)
	{
		if ((*activeOptionMap).count(name) == 0)
			throw "Option " + name + " not available!";
		return (*activeOptionMap)[name]->isSet;
	}

	virtual ~ParameterManager()
	{
		for (std::map<std::string,std::map<std::string, Parameter *> *>::iterator it1 = mainOptions.begin(); it1 != mainOptions.end(); it1++)
		{
			std::map<std::string, Parameter *> * current = it1->second;
			for (std::map<std::string, Parameter *>::iterator it2 = current->begin(); it2 != current->end(); it2++)
			{
				delete it2->second;
			}
			delete it1->second;
		}
	}

};

unsigned int ParameterManager::Parameter::counter = 0;

#endif /*PARAMETERMANAGER_H_*/
