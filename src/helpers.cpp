#include "helpers.h"

#include <algorithm>
#include <fstream>
#include "formulaparser.h"
#include <cstring>
#include <sstream>

using namespace std;

set<string> warnings;

/**
 * Inserts a number <num> at a position marked by the format
 * string %d in <fmt>
 */
string insertNumber(string fmt, unsigned int num)
{
	size_t pos = fmt.find("%d");
	if (pos != string::npos)
	{
		ostringstream s;
		s << num;
		return fmt.replace(pos,2,s.str());
	}
	return fmt;
}

void trim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos)
  {
	// erase trailing whitespaces
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos)
    // erase leading whitespaces
    	str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}

vector<double> * decodeVector(string vec)
{
	stringstream str(vec);
	string val;
	vector<double> * result = new vector<double>();
	while(getline(str, val, '/'))
	{
		result->push_back(atof(val.c_str()));
	}
	return(result);
}

inline void finishConstraint(NetworkConstraint * constraint)
{
	if (constraint->initialValues == NULL)
		constraint->initialValues = new SimpleCondition(0,0);

	if (constraint->fixedGenes == NULL)
			constraint->fixedGenes = new SimpleCondition(0,0);

	for (std::unordered_map<unsigned int, bool >::iterator it1 = constraint->fixedGenes->geneValues.begin();
			it1 != constraint->fixedGenes->geneValues.end(); ++it1)
	{
		if (constraint->initialValues->geneValues.count(it1->first) != 0)
		// check whether fixed genes and initial values contradict each other
		{
			if (constraint->initialValues->geneValues[it1->first] != it1->second)
				throw string("A gene with a fixed value cannot take a different value in the initial state!");
		}
		else
		// if the initial state does not specify the fixed gene, add it
		{
			constraint->initialValues->geneValues[it1->first] = it1->second;
		}

		for (vector<SimpleCondition *>::iterator it2 = constraint->stateSpecifications.begin();
				it2 != constraint->stateSpecifications.end(); ++it2)
		{
			if ((*it2)->geneValues.count(it1->first) != 0 &&
				(*it2)->geneValues[it1->first] != it1->second)
				throw string("A gene with a fixed value cannot take a different value in a state specification!");
		}
	}

	constraint->groupIndices.push_back((unsigned int)constraint->stateSpecifications.size());
}

void loadNetworkConstraints(std::string file, std::vector<std::string> &varNames, std::vector<NetworkConstraint *> &res)
{
	fstream f(file.c_str(),std::ios::in);
	FormulaParser p(varNames);

	if (f.fail())
		throw string("File " + file + " could not be opened!");

	string buf;
	NetworkConstraint * constraint = NULL;

	const int NUM_SECTIONS = 5;
	bool sections[NUM_SECTIONS];
	memset(sections, 0, sizeof(bool) * NUM_SECTIONS);


	enum Mode
	{
		SEC_INITIAL_STATE = 0,
		SEC_PRECONDITIONS = 1,
		SEC_FIXED_GENES = 2,
		SEC_SPECIFICATIONS = 3,
		SEC_IMPORTANCE = 4,
		OUTSIDE_SEC = 5,
	};

	Mode currentMode = OUTSIDE_SEC;
	bool isAttractor = false;
	double importanceSum = 0.0;

	while(true)
	{
		if (f.eof())
			break;

		// find next line that is not empty or a comment
		do
		{
			getline(f,buf);
			trim(buf);
		}
		while (!f.eof() &&  (buf.size() == 0 || buf[0] == '#'));

		if (f.eof())
			break;

		if (buf == "Chain:" || buf == "Attractor:")
		// start of a new rule
		{
			isAttractor = (buf == "Attractor:");
			if (constraint != NULL)
			// finish the last constraint
			{
				if (!sections[SEC_INITIAL_STATE] || !sections[SEC_SPECIFICATIONS])
				// the last read constraint was incomplete
					throw string("You must specify at least the \"Initial condition:\" and the \"State specifications:\" section!");

				finishConstraint(constraint);

				importanceSum += constraint->importance;

				// add the constraint
				res.push_back(constraint);
			}

			// initialize new constraint
			constraint = new NetworkConstraint();
			constraint->initialValues = NULL;
			constraint->isAttractor = isAttractor;
			constraint->importance = 1.0;

			// reset flags to read in a new constraint
			memset(sections, 0, sizeof(bool) * NUM_SECTIONS);
			currentMode = OUTSIDE_SEC;
		}
		else
		if (buf == "Initial condition:")
		{
			if (sections[SEC_INITIAL_STATE])
				throw string("Double section \"Initial condition:\"!");
			if (constraint == NULL)
				throw string("\"Attractor:\" or \"Chain:\" expected!");
			sections[SEC_INITIAL_STATE] = true;
			currentMode = SEC_INITIAL_STATE;
			constraint->groupIndices.push_back(0);
		}
		else
		if (buf == "Fixed genes:")
		{
			if (sections[SEC_FIXED_GENES])
				throw string("Double section \"Fixed genes:\"!");
			if (constraint == NULL)
				throw string("\"Attractor:\" or \"Chain:\" expected!");
			sections[SEC_FIXED_GENES] = true;
			currentMode = SEC_FIXED_GENES;
		}
		else
		if (buf == "Importance:")
		{
			if (sections[SEC_IMPORTANCE])
				throw string("Double section \"Importance:\"!");
			if (constraint == NULL)
				throw string("\"Attractor:\" or \"Chain:\" expected!");
			sections[SEC_IMPORTANCE] = true;
			currentMode = SEC_IMPORTANCE;
		}
		else
		if (buf == "Preconditions:" || buf == "State specifications:")
		{
			if (buf == "Preconditions:")
			{
				if (sections[SEC_PRECONDITIONS])
					throw string("Double section \"Preconditions:\"!");
				sections[SEC_PRECONDITIONS] = true;
				currentMode = SEC_PRECONDITIONS;
			}
			else
			{
				if (sections[SEC_SPECIFICATIONS])
					throw string("Double section \"State specifications:\"!");
				sections[SEC_SPECIFICATIONS] = true;
				currentMode = SEC_SPECIFICATIONS;
			}
			if (constraint == NULL)
				throw string("\"Attractor:\" or \"Chain:\" expected!");
		}
		else
		if (currentMode == OUTSIDE_SEC)
		// something is wrong with the file format
		{
			throw string("Expected a valid section header instead of \"" + buf + "\"!");
		}
		else
		if (buf=="Or:")
		// a new alternative is started
		{
			if (!sections[SEC_SPECIFICATIONS])
				throw string("Expected section \"State specifications:\" before section \"Or:\"!");
			constraint->groupIndices.push_back((unsigned int)constraint->stateSpecifications.size());
		}
		else
		{
			if (currentMode == SEC_INITIAL_STATE)
			// the current line is an initial condition
			{
				constraint->initialValues = new SimpleCondition(buf, varNames, 0, 0);
				currentMode = OUTSIDE_SEC;
			}
			else
			if (currentMode == SEC_FIXED_GENES)
			// the current line is a fixed gene specification
			{
				constraint->fixedGenes = new SimpleCondition(buf, varNames, 0, 0);
				currentMode = OUTSIDE_SEC;
			}
			else
			if (currentMode == SEC_IMPORTANCE)
			// the current line is the importance of the objective
			{
				constraint->importance = atof(buf.c_str());
				currentMode = OUTSIDE_SEC;
			}
			else
			// this is a precondition or state specification list
			{
				unsigned int index;
				if (currentMode == SEC_PRECONDITIONS)
				// the current line is a precondition
					index = (unsigned int)constraint->preconditions.size();
				else
				// the current line is an entry of a state specification list
					index = ((unsigned int)constraint->stateSpecifications.size()) - constraint->groupIndices[constraint->groupIndices.size() - 1];
				SimpleCondition * cond = new SimpleCondition(buf, varNames, index, ((unsigned int)constraint->groupIndices.size()) - 1);

				//cout << "Constraint: " << t->toString(&varNames[0]) << endl;
				if (currentMode == SEC_PRECONDITIONS)
					constraint->preconditions.push_back(cond);
				else
					constraint->stateSpecifications.push_back(cond);
			}
		}
	}
	if (constraint != NULL)
	// create the last constraint
	{
		if (!sections[SEC_INITIAL_STATE] || !sections[SEC_SPECIFICATIONS])
			throw string("You must specify at least the \"Initial condition:\" and the \"State specifications:\" section!");

		finishConstraint(constraint);

		importanceSum += constraint->importance;
		res.push_back(constraint);
	}
	f.close();

	for (vector<NetworkConstraint *>::iterator it = res.begin(); it != res.end(); ++it)
	{
//		cout << "\n\nAttractor: " << (*it)->isAttractor << endl;
//		cout << "Initial condition:" << endl;
//		for (std::::unordered_map<unsigned int, bool >::iterator it2 = (*it)->initialValues->geneValues.begin();
//			it2 != (*it)->initialValues->geneValues.end(); ++it2)
//			cout << varNames[it2->first] << " " << it2->second << " ";
//		cout << endl;
//
//		cout << "Rules:" << endl;
//		for (unsigned int j = 0; j < (*it)->groupIndices.size(); ++j)
//			cout << (*it)->groupIndices[j] << endl;
//
//		for (unsigned int j = 0; j < (*it)->rules.size(); ++j)
//		{
//			cout << "Alternative " << (*it)->rules[j]->alternative << endl;
//			for (std::tr1::unordered_map<unsigned int, bool >::iterator it2 = (*it)->rules[j]->geneValues.begin();
//						it2 != (*it)->rules[j]->geneValues.end(); ++it2)
//				cout << varNames[it2->first] << " " << it2->second << " ";
//			cout << endl;
//		}
//
//		cout << "Importance: " << (*it)->importance << endl;

		// normalize importance
		(*it)->importance /= importanceSum;
	}
}

void loadGeneDependencies(std::string file, std::vector<std::string> &varNames, double ** pos, double ** neg)
{
	map<string, int> assignment;
	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		assignment[varNames[i]] = i;
	}

	fstream f(file.c_str(),std::ios::in);



	if (f.fail())
		throw string("File " + file + " could not be opened!");

	string buf;

	getline(f,buf);
	buf = buf.substr(buf.find(' '));
	double dflt_pos = atof(buf.c_str());
	getline(f,buf);
	buf = buf.substr(buf.find(' '));
	double dflt_neg = atof(buf.c_str());

	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		for (unsigned int j = 0; j < varNames.size(); ++j)
		{
			pos[i][j] = dflt_pos;
			neg[i][j] = dflt_neg;
		}
	}

	while(true)
	{
		if (f.eof())
			break;
		do
		{
			getline(f,buf);
		}
		while (!f.eof() && (buf.size() == 0 || buf[0] == '#'));

		if (f.eof() || buf.size() == 0)
			break;

		size_t firstGenePos  = buf.find('>');
		if (firstGenePos < 0)
			throw string("\"=>\" or \"/>\" expected!");

		size_t secondGenePos  = buf.find(':');
		if (secondGenePos < 0)
			throw string("\":\" expected!");

		string firstGene = buf.substr(0,firstGenePos - 1);
		trim(firstGene);
		string secondGene = buf.substr(firstGenePos + 1,secondGenePos - firstGenePos - 1);
		trim(secondGene);

		string prob = buf.substr(secondGenePos + 1);

		cout << firstGene << " " << secondGene << " " << atof(prob.c_str()) << endl;
		if (buf[firstGenePos - 1] == '=')
			pos[assignment[firstGene]][assignment[secondGene]] = atof(prob.c_str());
		else
			neg[assignment[firstGene]][assignment[secondGene]] = atof(prob.c_str());
	}
	f.close();
}

void loadNetwork(std::string file, std::vector<std::string> &varNames, std::vector<BooleanTree *> &res, bool simplify)
{
	fstream f(file.c_str(),std::ios::in);

	vector<string> lines;
	vector<double> uncertainties;

	varNames.clear();

	if (f.fail())
		throw string("File " + file + " could not be opened!");

	string buf;
	getline(f,buf);
	while(true)
	{
		if (f.eof())
			break;

		// find the next non-empty line
		do
		{
			getline(f,buf);
		}
		while (!f.eof() && buf.size() == 0);

		if (f.eof())
			break;

		// find the separator
		size_t breakPos = buf.find(',');
		if (breakPos < 0)
			throw string("Missing separator!");

		// find the uncertainty of the variables, if any
		size_t uncertaintyPos = buf.find('[',breakPos+1);
		string geneName = buf.substr(0,breakPos);
		trim(geneName);

		// read the transition function
		string logic;
		if (uncertaintyPos < 0)
		{
			logic = buf.substr(breakPos+1);
			uncertainties.push_back(1.0);
		}
		else
		{
			logic = buf.substr(breakPos+1,uncertaintyPos-breakPos-1);

			size_t uncertaintyEndPos = buf.find(']', uncertaintyPos);
			if (uncertaintyEndPos < 0)
				throw string("Missing ] symbol!");
			string uncertainty = buf.substr(uncertaintyPos+1,uncertaintyEndPos-uncertaintyPos-1);
			uncertainties.push_back(1.0 - atof(uncertainty.c_str()));
		}

		// add the variable names and functions to lists
		varNames.push_back(geneName);
		lines.push_back(logic);
	}

	FormulaParser p(varNames);
	res.resize(varNames.size());

	// parse the transition functions
	for (unsigned int i = 0; i < varNames.size(); ++i)
	{
		BooleanTree * rule;
		try
		{
			rule = new BooleanTree(p.parse(lines[i]), uncertainties[i]);
		}
		catch (ParseException &ex)
		{
			throw string(ex.what());
		}
		if (simplify)
			rule->simplify();

		res[i] = rule;
	}
	f.close();
}

//void removeDuplicates(std::vector<Individual *> &population, bool freeDuplicates)
//{
//	vector<vector<BooleanTree *> > normalizedIndividuals;
//	vector<vector<BooleanTree *> > uniqueList;
//
//	// create normalized copies of all individuals
//	for (vector<Individual *>::iterator it = population.begin();it != population.end(); ++it)
//	{
//		vector<BooleanTree *> net;
//		for (vector<BooleanTree *>::iterator f = (*it)->network.begin(); f != (*it)->network.end(); ++f)
//		{
//			net.push_back((*f)->getNormalizedCopy());
//		}
//
//		normalizedIndividuals.push_back(net);
//	}
//
//	vector<Individual *> oldPop;
//	oldPop.assign(population.begin(), population.end());
//	population.clear();
//
//	// check for equality
//	for (unsigned int i = 0; i < normalizedIndividuals.size(); ++i)
//	{
//		bool duplicate = false;
//
//		// iterate over the list of current unique individuals
//		for (vector<vector<BooleanTree *> >::iterator it2 = uniqueList.begin(); it2 != uniqueList.end(); ++it2)
//		{
//			if (normalizedIndividuals[i].size() != it2->size())
//			// the individuals have different numbers of genes => unequal
//				continue;
//
//			bool equal = true;
//
//			// compare the functions of the two individuals
//			for (unsigned int k = 0; k < it2->size(); ++k)
//			{
//				if (!(*((*it2)[k]) == *normalizedIndividuals[i][k]))
//				{
//					equal = false;
//					break;
//				}
//			}
//			if (equal)
//			{
//				duplicate = true;
//				break;
//			}
//		}
//
//		if (!duplicate)
//		// if this individual is not in the current result list,
//		// push it there, and add its normalized copy to uniqueList
//		{
//			population.push_back(oldPop[i]);
//			uniqueList.push_back(normalizedIndividuals[i]);
//		}
//		else
//		{
//			if (freeDuplicates)
//				delete oldPop[i];
//		}
//	}
//
//	// free normalized copies
//	for (vector<vector<BooleanTree *> >::iterator it1 = normalizedIndividuals.begin(); it1 != normalizedIndividuals.end(); ++it1)
//	{
//		for (vector<BooleanTree *>::iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
//		{
//			delete *it2;
//		}
//	}
//}

void removeDuplicates(std::vector<Individual *> &population, bool freeDuplicates)
{

	vector<Individual *> oldPop;
	oldPop.assign(population.begin(), population.end());
	population.clear();

	// check for equality
	for (unsigned int i = 0; i < oldPop.size(); ++i)
	{
		bool duplicate = false;

		// iterate over the list of current unique individuals
		for (vector<Individual * >::iterator it2 = population.begin(); it2 != population.end(); ++it2)
		{
			if (oldPop[i]->normalizedNetwork.size() != (*it2)->normalizedNetwork.size())
			// the individuals have different numbers of genes => unequal
				continue;

			bool equal = true;

			// compare the functions of the two individuals
			for (unsigned int k = 0; k < (*it2)->normalizedNetwork.size(); ++k)
			{
				if (!(*((*it2)->normalizedNetwork[k]) == *(oldPop[i]->normalizedNetwork[k])))
				{
					equal = false;
					break;
				}
			}
			if (equal)
			{
				duplicate = true;
				break;
			}
		}

		if (!duplicate)
		// if this individual is not in the current result list,
		// push it there, and add its normalized copy to uniqueList
		{
			population.push_back(oldPop[i]);
		}
		else
		{
			if (freeDuplicates)
				delete oldPop[i];
		}
	}

}
