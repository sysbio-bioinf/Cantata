/*
 * main.cpp
 * Main function and command line
 * argument processing
 *
 *  Created on: May 11, 2010
 *      Author: muessel
 * Edited 2020 by schwab
 */
#define VERSION "13/04"

#include "helpers.h"
#include "gp.h"
#include "random.h"
#include "formulaparser.h"
#include "randomnet.h"
#include "parametermanager.h"
#include "constraintviolation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <numeric>

using namespace std;

int main(int argc, char ** argv)
{
	ParameterManager pm;
	bool verbose = true;

	// determine the parameters for optimization
	pm.addMainOption("--optimize");
	pm.addStringOption("--optimize","-n","<network file>","",false);
	pm.addStringOption("--optimize","-r","<rule file>","",false);
	pm.addStringOption("--optimize","-o","<result file>","result.txt",true);
	pm.addStringOption("--optimize","-on","<output network files>","resultnet_%d.txt",true);
#ifdef TRACE_EVOLUTION
	pm.addStringOption("--optimize","-an","<ancestor network files>","ancestor_%d_%d.txt",true);
#endif
	pm.addDoubleOption("--optimize","-me","<max error>", 0.0, true);
	pm.addIntOption("--optimize","-ps","<population size>",100,true);
	pm.addIntOption("--optimize","-no","<number of offspring>",200,true);
	pm.addDoubleOption("--optimize","-if","<frac. injected original networks>",0.1,true);
	pm.addIntOption("--optimize","-nf","<allow negation every xth offspring>",50,true);
	pm.addIntOption("--optimize","-ni","<number of iterations>",1000,true);
	pm.addIntOption("--optimize","-ns","<number of restarts>",1,true);
	pm.addIntOption("--optimize","-ms","<max states>", 200, true);
	pm.addIntOption("--optimize","-mt","<max transitions>", 100, true);
	pm.addStringOption("--optimize","-tw","<topology fitness weight vector>", "0.25/0.25/0.5", true);
	pm.addIntOption("--optimize","-im","<initial mutations for new individuals>", 1, true);
	pm.addDoubleOption("--optimize","-eps","<epsilon>",0.0005,true);
	pm.addIntOption("--optimize","-rs","<random seed>",(int)time(0),true);
	pm.addVoidOption("--optimize","-q",true);
	pm.addVoidOption("--optimize","-greedy",true);

	// determine the parameters for random network generation
	pm.addMainOption("--randomnet");
	pm.addStringOption("--randomnet","-o","<output network file>","",false);
	pm.addIntOption("--randomnet","-ng","<number of genes>",10,true);
	pm.addIntOption("--randomnet","-ns","<number of generation steps>",10,true);
	pm.addIntOption("--randomnet","-rs","<random seed>",(int)time(0),true);

	// determine the parameters for network truncation
	pm.addMainOption("--truncate");
	pm.addStringOption("--truncate","-n","<original network file>","",false);
	pm.addStringOption("--truncate","-o","<output network file>","",false);
	pm.addIntOption("--truncate","-nd","<number of deletions>",2,true);
	pm.addIntOption("--truncate","-rs","<random seed>",(int)time(0),true);
	pm.addStringOption("--truncate","-tg","<target gene>","",true);
	pm.addStringOption("--truncate","-dg","<dependency>","",true);
	pm.addDoubleOption("--truncate","-ip","<insertion probability>",0.0,true);
	pm.addVoidOption("--truncate","-fd",true);
	pm.addVoidOption("--truncate","-ae",true);


	// determine the parameters for validation of an existing network
	pm.addMainOption("--validate");
	pm.addStringOption("--validate","-n","<network file>","",false);
	pm.addStringOption("--validate","-r","<rule file>","",false);
	pm.addStringOption("--validate","-o","<result file>","result.txt",true);
	pm.addIntOption("--validate","-ms","<max states>", 100000, true);
	pm.addIntOption("--validate","-mt","<max transitions>", 1000, true);
	pm.addIntOption("--validate","-c","<collapse>", 50, true);
	pm.addVoidOption("--validate","-ds", true);
	pm.addIntOption("--validate","-rs","<random seed>", (int)time(0), true);
	pm.addVoidOption("--validate","-q",true);
	pm.addVoidOption("--validate","-greedy",false);

	if (!pm.analyzeCommandLine(argc,argv))
	// incorrect command line => usage message
	{
		cerr << "CANTATA build " << VERSION << "\n";
		cerr << pm.usage(string(argv[0]));
		return 1;
	}

	setrandomseed(pm.getIntOption("-rs"));

	if (pm.getActiveMainOption() == "--optimize")
	{

		vector<NetworkConstraint *> constraints;
		vector<BooleanTree *> network;
		vector<string> varNames;
		string networkFile = pm.getStringOption("-n");
		string ruleFile = pm.getStringOption("-r");

		verbose = !pm.isOptionSet("-q");

		try
		{
			loadNetwork(networkFile, varNames, network,true);
		}
		catch (string &s)
		{
			cerr << "An error occurred while reading the network file:\n\t" << s << endl;
			return 1;
		}

		try
		{
			loadNetworkConstraints(ruleFile, varNames, constraints);
		}
		catch (string &s)
		{
			cerr << "An error occurred while reading the rule file:\n\t" << s << endl;
			return 1;
		}

		vector<double>  * topologyWeights = decodeVector(pm.getStringOption("-tw"));

		if (topologyWeights->size() != 3 ||
			abs(std::accumulate(topologyWeights->begin(), topologyWeights->end(), 0.0) - 1.0) > 1e-3)
		{
			cerr << "Parameter \"-tw\" must be a vector of exactly 3 floating point numbers, "
				 <<	"separated by \"/\" and summing up to 1!" << endl;
			return 1;
		}

		if (verbose)
		{
			cout << "Input network:\n";
			for (unsigned int i = 0; i < varNames.size(); ++i)
			{
				cout << varNames[i] << " = " <<  network[i]->toString(&varNames[0]) << endl;
			}
			cout << endl;
		}

		GP g(varNames, network, constraints, topologyWeights,
			 pm.getIntOption("-ms"),
			 pm.getIntOption("-mt"),
			 pm.getDoubleOption("-eps"),
			 pm.isOptionSet("-greedy"));

		delete topologyWeights;

		g.start(pm.getIntOption("-ps"),
				pm.getIntOption("-no"),
				pm.getIntOption("-ni"),
				pm.getIntOption("-ns"),
				pm.getIntOption("-im"),
				pm.getDoubleOption("-if"),
				pm.getIntOption("-nf"),
				verbose);

		// write results to file or stdout
		ostream * out;
		if (pm.isOptionSet("-o"))
			out = new fstream(pm.getStringOption("-o").c_str(),ios::out);
		else
			out = &cout;

		// write header with configuration
		(*out) << "CANTATA build " << VERSION << "\n";
		(*out) << setw(30) << left << "Input network file: " << networkFile << "\n";
		(*out) << setw(30) << left << "Rule file: " << ruleFile << "\n";
		(*out) << setw(30) << left << "Random seed: " << pm.getIntOption("-rs") << "\n";
		(*out) << setw(30) << left << "Population size: " << pm.getIntOption("-ps") << "\n";
		(*out) << setw(30) << left << "Number of offspring: " << pm.getIntOption("-no") << "\n";
		(*out) << setw(30) << left << "Fract. of injected nets: " << pm.getDoubleOption("-if") << "\n";
		(*out) << setw(30) << left << "Neg. every i-th offspring: " << pm.getIntOption("-nf") << "\n";
		(*out) << setw(30) << left << "Number of generations: " << pm.getIntOption("-ni") << "\n";
		(*out) << setw(30) << left << "Number of restarts: " << pm.getIntOption("-ns") << "\n";
		(*out) << setw(30) << left << "Initial mutations: " << pm.getIntOption("-im") << "\n";
		(*out) << setw(30) << left << "Epsilon: " << pm.getDoubleOption("-eps") << "\n";
		(*out) << setw(30) << left << "Weights of topology scores: " << pm.getStringOption("-tw") << "\n";
		(*out) << setw(30) << left << "Max. number of start states: " << pm.getIntOption("-ms") << "\n";
		(*out) << setw(30) << left << "Max. number of transitions: " << pm.getIntOption("-mt") << "\n\n";
		//(*out) << setw(30) << left << "Greedy search enabled: " << pm.getIntOption("-greedy") << "\n\n";
		(*out) << "Best candidate networks: " << "\n\n";

		bool writePerfectNetworks = pm.isOptionSet("-on");
		string networkFilePattern = pm.getStringOption("-on");
		double maxNetworkError = pm.getDoubleOption("-me");
		vector<Individual *> resultList = g.getResultList();
		unsigned int perfectCount = 0;

		for (vector<Individual *>::iterator it = resultList.begin(); it != resultList.end(); ++it)
		{
			(*out) << (*it)->toString() << "\n" << endl;

			if (writePerfectNetworks && (*it)->fitness[0] <= maxNetworkError)
			// write separate network file
			{
				++perfectCount;
				fstream net(insertNumber(networkFilePattern, perfectCount).c_str(),ios::out);
				net << "targets, factors\n";
				net << (*it)->toString(", ", false, false) << "\n" << endl;
				net.close();

#ifdef TRACE_EVOLUTION
				unsigned int anc = 0;
				if (pm.isOptionSet("-an"))
				{
					for (vector<Individual *>::iterator it2 = (*it)->ancestors.begin(); it2 != (*it)->ancestors.end(); ++it2)
					{
						++anc;
						fstream net2(insertNumber(insertNumber(pm.getStringOption("-an"), perfectCount),anc).c_str(),ios::out);
						net2 << "targets, factors\n";
						net2 << (*it2)->toString(", ", false, false) << "\n" << endl;
						net2.close();
					}
				}
#endif
			}
		}

		if (pm.isOptionSet("-o"))
		{
			((fstream *)out)->close();
			delete out;
		}


//		g.start(1,1,1,1,"result.txt",false,"",0.0);
	}
	else
	if (pm.getActiveMainOption() == "--randomnet")
	{
		try
		{
			writeRandomNetworks(pm.getStringOption("-o"),
								pm.getIntOption("-ng"),
								pm.getIntOption("-ns"));
		}
		catch (string &s)
		{
			cerr << s << endl;
			return 1;
		}
	}
	else
	if (pm.getActiveMainOption() == "--truncate")
	{
		try
		{
			if ((pm.isOptionSet("-tg") || pm.isOptionSet("-dg")) &&
				 (pm.isOptionSet("-nd") || pm.isOptionSet("-ae") || pm.isOptionSet("-ip")))
				 throw string("Please provide either a dependency to delete via -tg and -dg or a number of modifications via -nd");
			if (pm.isOptionSet("-tg"))
			{
				if (!pm.isOptionSet("-dg"))
					throw string("Please provide a dependency to remove in -dg!");
				removeDependency(pm.getStringOption("-n"),
								 pm.getStringOption("-o"),
								 pm.getStringOption("-tg"),
								 pm.getStringOption("-dg"));
			}
			else
			{
				mutateNetwork(pm.getStringOption("-n"),
							  pm.getStringOption("-o"),
							  pm.getIntOption("-nd"),
							  pm.getDoubleOption("-ip"),
							  pm.isOptionSet("-ae"),
							  pm.isOptionSet("-fd"));
			}
		}
		catch (string &s)
		{
			cerr << s << endl;
			return 1;
		}
	}
	else
	if (pm.getActiveMainOption() == "--validate")
	{
		vector<NetworkConstraint *> constraints;
				vector<BooleanTree *> network;
				std::vector<string> varNames;
		string networkFile = pm.getStringOption("-n");
		string ruleFile = pm.getStringOption("-r");
		verbose = !pm.isOptionSet("-q");

		try
		{
			loadNetwork(networkFile, varNames, network,true);
		}
		catch (string &s)
		{
			cerr << "An error occurred while reading the network file:\n\t" << s << endl;
			return 1;
		}

		try
		{
			loadNetworkConstraints(ruleFile, varNames, constraints);
		}
		catch (string &s)
		{
			cerr << "An error occurred while reading the rule file:\n\t" << s << endl;
			return 1;
		}

		try
		{

			if (verbose)
			{
				cout << "Input network:\n";
				for (unsigned int i = 0; i < varNames.size(); ++i)
				{
					cout << varNames[i] << " = " <<  network[i]->toString(&varNames[0]) << endl;
				}
				cout << endl;
			}

			GP g(varNames, network, constraints,
				 NULL,
				 pm.getIntOption("-ms"),
				 pm.getIntOption("-mt"),
				 0.0,
				 pm.isOptionSet("-greedy"));
			Individual ind(g, network, varNames, 0, 0, false, true);

			bool writeOutput = pm.isOptionSet("-o");
			ostream * out;
			if (writeOutput)
				out = new fstream(pm.getStringOption("-o").c_str(),ios::out);
			else
				out = &cout;

			(*out) << "CANTATA build " << VERSION << "\n";
			(*out) << setw(30) << left << "Input network file: " << networkFile << "\n";
			(*out) << setw(30) << left << "Rule file: " << ruleFile << "\n";
			(*out) << setw(30) << left << "Random seed: " << pm.getIntOption("-rs") << "\n";
			(*out) << setw(30) << left << "Max. number of start states: " << pm.getIntOption("-ms") << "\n";
			(*out) << setw(30) << left << "Max. number of transitions: " << pm.getIntOption("-mt") << "\n\n";

			(*out) << "Violations of the rule set:\n" << endl;

			for (unsigned int i = 0; i < constraints.size(); ++i)
			{
				(*out) << "Rule " << (i+1) << ":\n";
				ConstraintViolation * v;
				if (pm.isOptionSet("-greedy"))
					v = constraints[i]->getConstraintViolationsGreedy(ind,
																	  pm.getIntOption("-ms"),
																	  pm.getIntOption("-mt"),
																	  !pm.isOptionSet("-ds"));
				else
					v = constraints[i]->getConstraintViolationsDP(ind,
																  pm.getIntOption("-ms"),
																  pm.getIntOption("-mt"),
																  !pm.isOptionSet("-ds"));
				(*out) << v->toString(&varNames[0], pm.getIntOption("-c")) << endl;
				delete v;
			}

			if (writeOutput)
			{
				((fstream *)out)->close();
				delete out;
			}
		}
		catch (string &s)
		{
			cerr << "An error occured while analyzing the network:\n\t" << s << endl;
			return 1;
		}

	}

	if (verbose)
	{
		if (warnings.size() > 0)
		{
			cerr << "There were warnings:\n\n";
			for (set<string>::iterator it = warnings.begin(); it != warnings.end(); ++it)
				cerr << *it << endl;
			cerr << endl;
		}
		cout << "Finished" << endl;
	}
}
