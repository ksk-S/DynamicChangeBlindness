/*
This file is part of the decisional state reconstruction algorithm
technique exposed in "Decisional States", by Nicolas Brodu.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
    MA  02110-1301  USA

See http://nicolas.brodu.numerimoire.net/en/programmation/decisional_states/index.html
for more information and possibly updates.

Copyright holder: Nicolas Brodu <nicolas.brodu@numerimoire.net>
File release Date: February 09
*/

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

#include <stdint.h>

#include "decisional_states.hpp"
#include "decisional_states_samplers.hpp"
#include "decisional_states_optimisers.hpp"

#include "decisional_states_aggregate_clustering.hpp"


using namespace std;
using namespace boost;
using namespace decisional_states;

map<char,int> symbolsmap;
vector<char> symbols;
int PastSize = 10, FutureSize = 1;

uint64_t ipow(uint64_t base, uint64_t exp) {
    if (exp==0) return 1;
    if (exp==1) return base;
    uint64_t halfexp = exp/2;
    if (halfexp*2 == exp) return ipow(base*base, halfexp);
    return ipow(base*base, halfexp) * base;
}

struct DataSet : public vector<pair<uint64_t,uint64_t> > {
    
    // Decisional states API
    typedef uint64_t DataType;
    typedef uint64_t PredictionType;
    typedef char SymbolType;
    
    // Symbol emitted when passing from a to b
    bool getSymbol(DataType a, DataType b, SymbolType& symbol) {
        // by construction (see main) this is the latest entry in the vector b
        symbol = symbols[b % symbols.size()];
        // and only if both data match (multi-series file have no transitions between series)
        return a % discard_first_symbol_factor == b / symbols.size();
    }
    
    // our own stuff
    uint64_t discard_first_symbol_factor;
    DataSet() {
        assert(!symbols.empty());
        discard_first_symbol_factor = ipow(symbols.size(),PastSize-1);
    }
};

typedef DiscreteDistributionManager<DataSet> DistManager;
typedef DecisionalStatesAnalyser<DistManager, DataSet> DSA;

// probability that the reconstructed model generates the given series
float getPastSeriesInferredProbability(uint64_t series, DSA& analyser, int totalRecCount) {
    vector<char> symbolseries(PastSize);
    for (int i=PastSize-1; i>=0; --i) {
        symbolseries[i] = symbols[series % symbols.size()];
        series /= symbols.size();
    }
    float proba = 0;
    // sum of all probabilities from all starting states
    for (DSA::CausalStates::iterator cit = analyser.causalStates.begin(); cit != analyser.causalStates.end(); ++cit) {
        if (!(*cit)->is_recurrent()) continue;
        DSA::CausalState * currentState = cit->get();
        // probability of emitting the string starting from the current state
        // = proba of current state * proba of all transitions
        float chainedProba = currentState->count / (float)totalRecCount;
        // for each symbol in the string, look at the transitions emitting this symbol
        // and update the state as we go
        for (int i=0; i<PastSize; ++i) {
            bool transitionFound = false;
            // lookup which transition from the current state emits the series symbol
            for (DSA::CausalStateTransitionSet::iterator it = currentState->transitions->begin(); it != currentState->transitions->end(); ++it) {
                if (it->label==symbolseries[i]) {
                    chainedProba *= it->probability;
                    transitionFound = true;
                    currentState = it->state;
                    break;
                }
            }
            if (!transitionFound) {
                chainedProba = 0;
                break;
            }
        }
        proba += chainedProba;
    }
    return proba;
}

float getRelativeEntropy(DistManager& distManager, DSA& analyser, int totalCount, int totalRecCount) {
    float relativeEntropy = 0;
    for (int i=0; i<(int)distManager.size(); ++i) {
        float observedProba = distManager.getCount(distManager[i]) / (float)totalCount;
        uint64_t series = *distManager.getData(i);
        relativeEntropy += observedProba * (log2(observedProba/getPastSeriesInferredProbability(series,analyser,totalRecCount)));
    }
    return relativeEntropy;
}

float getEntropyRate(DistManager& distManager, DSA& analyser, int totalRecCount) {
    float entropyRate = 0;
    // sum of all probabilities from all starting states
    for (DSA::CausalStates::iterator cit = analyser.causalStates.begin(); cit != analyser.causalStates.end(); ++cit) {
        if (!(*cit)->is_recurrent()) continue;
        DSA::CausalState * currentState = cit->get();
        float stateProba = currentState->count / (float)totalRecCount;
        float distentropy = 0;
        int disttotalcount = 0;
        for (int i=0; i<distManager.nz; ++i) disttotalcount += currentState->distribution.counts[i];
        for (int i=0; i<distManager.nz; ++i) {
            float p = currentState->distribution.counts[i] / (float)disttotalcount;
            if (p) distentropy -= p * log2(p);
        }
        entropyRate += stateProba * distentropy;
    }
    return entropyRate;
}

int main (int argc, char** argv) {

    if (argc<4) {
        cout << "Process a symbolic series and compute the epsilon-machine" << endl;
        cout << "Arguments required: series_file past_size future_size [tag] [leakThreshold]" << endl;
        cout << "- series_file shall contain one series per line" << endl;
        cout << "- tag defaults to the series file name. The program generates four output files:"<< endl;
        cout << "  . tag.dot: graph of the reconstructed epsilon-machine" << endl;
        cout << "  . tag_transients.dot: graph of the reconstructed machine including transient states" << endl;
        cout << "  . tag_states.txt: the series of states for each original data" << endl;
        cout << "  . tag_members.txt: the strings that are in each state" << endl;
        cout << " - leakThreshold sets a lower bound to the probability for the transitions that leave the recurrent states. By default this is 0: the recurrent states really must be recurrent. But if you have very rare transitions (ex: happening with p=1e-7) that are artifacts of your data (ex: measurement noise), these could destroy the detection of recurrent states if these spurious transitions lead out of the recurrent connected component (ex: to an \"error state\" you never come back from). In that case just specify a higher threshold (ex: 1e-6) so the recurrent states ignore these transitions. You may always look at \"tag_transients.dot\" to see all states anyway." << endl;
        return 0;
    }
    ifstream seriesfile(argv[1]);
    PastSize = atoi(argv[2]);
    FutureSize = atoi(argv[3]);
    string tag = argv[1];
    if (argc>4) {
        tag = argv[4];
    }
    double leakThreshold = 0;
    if (argc>5) {
        leakThreshold = atof(argv[5]);
    }

    list<vector<char> > rawdata;
    string line;
    while (seriesfile && !seriesfile.eof()) {
        getline(seriesfile, line);
        trim(line);
        if (line.empty()) continue;
        rawdata.push_back(vector<char>());
        stringstream linereader(line);
        char value;
        while (linereader >> value) {
            rawdata.back().push_back(value);
            symbolsmap[value] = 0;
        }
    }
    
    int nsymbols = symbolsmap.size();
    cout << "File " << argv[1] << " loaded, found " << nsymbols << " symbols: ";
    symbols.resize(nsymbols);
    int sidx = 0;
    for (map<char,int>::iterator it = symbolsmap.begin(); it!=symbolsmap.end(); ++it) {
        cout << it->first;
        symbols[sidx] = it->first;
        it->second = sidx++;
    }
    cout << endl;

    // ensure 64-bit identifiers are enough to hold the series
    if (PastSize * log2(nsymbols) >= 64) {
        cout << "Unsupported past size for this version of the program: too many combinations" << endl;
        return 1;
    }
    if (FutureSize * log2(nsymbols) >= 64) {
        cout << "Unsupported future size for this version of the program: too many combinations" << endl;
        return 1;
    }
    
    int ndata = 0;
    for (list<vector<char> >::iterator it = rawdata.begin(); it!=rawdata.end(); ++it) {
        int ndata_this_series = it->size() - PastSize - FutureSize + 1;
        if (ndata_this_series>0) ndata += ndata_this_series;
    }
    
    // Generate data set -- use a sliding window on each series
    DataSet dataset;
    dataset.reserve(ndata);
    uint64_t largestOrderPast = ipow(nsymbols,PastSize-1);
    uint64_t largestOrderFuture = ipow(nsymbols,FutureSize-1);
    set<uint64_t> allFutures;
    for (list<vector<char> >::iterator it = rawdata.begin(); it!=rawdata.end(); ++it) {
        int ndata_this_series = it->size() - PastSize - FutureSize + 1;
        if (ndata_this_series<=0) continue;
        vector<char>& series = *it;
        vector<char> history(PastSize);
        vector<char> future(FutureSize);
        
        uint64_t curval_past = 0;
        for (int i=0; i<PastSize; ++i) {
            history[i] = series[i];
            curval_past = curval_past*nsymbols + symbolsmap[history[i]];
        }
        
        uint64_t curval_future = 0;
        for (int i=0; i<FutureSize; ++i) {
            future[i] = series[i+PastSize];
            curval_future = curval_future*nsymbols + symbolsmap[future[i]];
        }
        
        dataset.push_back(make_pair(curval_past, curval_future));
        allFutures.insert(curval_future);
        
        for (int i=PastSize+FutureSize; i<(int)it->size(); ++i) {
            curval_past -= largestOrderPast * symbolsmap[history[0]];
            for (int j=0; j<PastSize-1; ++j) history[j] = history[j+1];
            history[PastSize-1] = future[0];
            curval_past = curval_past*nsymbols + symbolsmap[future[0]];
            curval_future -= largestOrderFuture * symbolsmap[future[0]];
            for (int j=0; j<FutureSize-1; ++j) future[j] = future[j+1];
            future[FutureSize-1] = series[i];
            curval_future = curval_future*nsymbols + symbolsmap[series[i]];
            dataset.push_back(make_pair(curval_past, curval_future));
            allFutures.insert(curval_future);
        }
    }
        
    DistManager distManager(dataset);
    DSA analyser(distManager,dataset);

    analyser.computeCausalStates();
    int ncs = analyser.causalStates.size();
    cout << "Number of states (including transient ones): " << ncs;

    // build the epsilon-machine
    // some optional leak threshold for the recurrent states may be allowed by the user
    // => you could also call analyser.buildCausalStateGraph();
    analyser.buildCausalStateGraph(leakThreshold);
    
    cout << ". Number of recurrent states: " << analyser.causalStatesCounts.nrecurrent << endl;

    ofstream dot((tag+".dot").c_str());
    analyser.writeCausalStateGraph(dot);
    dot.close();
    
    dot.open((tag+"_transients.dot").c_str());
    analyser.writeCausalStateGraph(dot,false);
    dot.close();

    // series of states
    ofstream states_series((tag+"_states.txt").c_str());
    map<int, set<uint64_t> > state_members;
    int datasetidx = 0;
    for (list<vector<char> >::iterator it = rawdata.begin(); it!=rawdata.end(); ++it) {
        int ndata_this_series = it->size() - PastSize - FutureSize + 1;
        if (ndata_this_series<=0) continue;
        for (int i=0; i<PastSize-1; ++i) {
            if (i>0) states_series<<" ";
            states_series<<"?";
        }
        for (int i=PastSize-1; i<(int)it->size()-FutureSize; ++i) {
            boost::shared_ptr<DSA::CausalState> state = analyser.getCausalState(dataset[datasetidx].first);
            states_series<<" "<<state->index;
            state_members[state->index].insert(dataset[datasetidx].first);
            ++datasetidx;
        }
        // estimate the states from the original series for the last points
        for (int i=0; i<FutureSize; ++i) {
            uint64_t curval_past = 0;
            for (int j=0; j<PastSize; ++j) {
                curval_past = curval_past*nsymbols + symbolsmap[(*it)[ndata_this_series+i+j]];
            }
            boost::shared_ptr<DSA::CausalState> state = analyser.getCausalState(curval_past);
            states_series<<" "<<state->index;
            state_members[state->index].insert(dataset[datasetidx].first);
        }
        states_series << endl;
    }
    states_series.close();
    
    // members of each state
    ofstream members((tag+"_members.txt").c_str());
    for (map<int, set<uint64_t> >::iterator it=state_members.begin(); it!=state_members.end(); ++it) {
        members << "State " << it->first << ":" << endl;
        for (set<uint64_t>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
            uint64_t past_series = *sit;
            char string_series[PastSize+1];
            for (int i=PastSize-1; i>=0; --i) {
                string_series[i] = symbols[past_series % symbols.size()];
                past_series /= symbols.size();
            }
            string_series[PastSize] = 0;
            members << string_series << endl;
        }
    }
    members.close();

    // Extra information like what's done in CSSR
    cout << "Statistical complexity: " << analyser.statisticalComplexity() << endl;
    cout << "Relative entropy: " << getRelativeEntropy(distManager, analyser, dataset.size(), analyser.causalStatesCounts.recurrentcount) << endl;
    cout << "Entropy rate: " << getEntropyRate(distManager, analyser, analyser.causalStatesCounts.recurrentcount) << endl;
    
    return 0;
}

