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

#include "decisional_states.hpp"
#include "decisional_states_optimisers.hpp"
#include "decisional_states_connected_clustering.hpp"

#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>
#include <boost/random.hpp>

#include <fstream>
#include <sstream>

#include <set>

using namespace std;
using namespace boost;
using namespace decisional_states;

#include <stdlib.h>
#ifndef NO_SDL
#include <SDL/SDL.h>
#endif

// wrapping vector in each direction, up to limited distance
template<typename T>
struct WrappingVector : public vector<T> {
    WrappingVector(int size) : vector<T>(size) {}
    T& operator[](int idx) {
        int size = this->size();
        if (idx<0) idx += size;
        else if (idx>=size) idx -= size;
        assert(idx>=0); assert(idx<size);
        return vector<T>::operator[](idx);
    }
};

struct AutomataWorld : public WrappingVector<int> {
    AutomataWorld() : WrappingVector<int>(worldSize) {}
    static int worldSize;
};
int AutomataWorld::worldSize;


struct Rule {

    int nextState[8];

    Rule(unsigned int number) {
        // convert CA rule notation to state array
        for (int i=0; i<8; ++i) {
            nextState[i] = number & 1;
            number = number >> 1;
        }
    }

    // apply the rule to the 3 cells that are either 0 or 1
    int apply(int left, int center, int right) {
        return nextState[ (left << 2) | (center << 1) | right ];
    }
};

// to generate the initial conditions
struct RandFunctor {
    boost::mt19937 rng;
    RandFunctor(unsigned int seed) {
        rng.seed(seed);
    }
    int operator()() {
        return (rng() >> 17) & 1; // use any bit
    }
};

#ifndef NO_SDL
SDL_Surface *screen;
void drawGrayPixel(int x, int y, int gray) {
    *((int*)screen->pixels + y*screen->pitch/4 + x) = SDL_MapRGB(screen->format, gray, gray
, gray);
}
#else
void drawGrayPixel(int x, int y, int gray) {}
#endif

/*
1   1
3   4
5   9
7   16
9   25
11  36
13  49
...
64 bits enough for most depth
*/

typedef vector<pair<uint64_t,uint64_t> > DataSet;

// Feeds the real transitions to the causal states analyser from the data set
struct TransitionFeeder {
    DataSet& dataset;
    int world_size;
    int max_steps;
    int ntrans;
    typedef int SymbolType;
    
    int past_row_size;
    Rule rule;
    
    TransitionFeeder(DataSet& _dataset, int _world_size, int _max_steps, int past_depth, int rule_number) : dataset(_dataset), world_size(_world_size), max_steps(_max_steps), past_row_size(past_depth*2-1), rule(rule_number) {
        ntrans = world_size * (max_steps-1);
    }
    
    typedef int iterator;
    iterator begin() {return 0;}
    iterator end() {return ntrans;}
    
    uint64_t& getDataBeforeTransition(iterator idx) {
        int worldcell = idx / (max_steps-1);
        int time = idx % (max_steps-1);
        return dataset[time * world_size + worldcell].first;
    }
    
    uint64_t& getDataAfterTransition(iterator idx) {
        int worldcell = idx / (max_steps-1);
        int time = idx % (max_steps-1);
        return dataset[(time+1) * world_size + worldcell].first;
    }
    
    SymbolType getSymbol(iterator idx) {
        int worldcell = idx / (max_steps-1);
        int time = idx % (max_steps-1);
        
        // Symbol emitted when passing from a to b
        //uint64_t a = dataset[time * world_size + worldcell].first;
        uint64_t b = dataset[(time+1) * world_size + worldcell].first;
        
        // the new information that was not in the last cone is on the edge
        // and is equivalent to a 2-bit symbol
        uint64_t bit1 = b & 1;
        uint64_t bit2 = (b>>(past_row_size-1)) & 1;
        return bit1*2+bit2;
    }
    
};


struct Optimiser {
    int nbits;
    Optimiser(int _nbits) : nbits(_nbits) {}

    typedef BasicPredictionSet<uint64_t> PredictionSet;

    // Optimise takes the functor as double even if utility return ints
    // because we are optimising the expected utility
    template<class Functor> double optimise(Functor f, PredictionSet& res) {
        res = PredictionSet();
        double bestvalue = -numeric_limits<double>::max();
        // exhaustive search on nbits
        for (uint64_t i=0; i<(uint64_t)(1<<nbits); ++i) {
            double result = f(i);
            if (result > bestvalue) {
                res->clear();
                res->insert(i);
                bestvalue = result;
            }
            else if (result == bestvalue) {
                res->insert(i);
            }
        }
        return bestvalue;
    }
};


// Utility function: bit distance between predictions, 1 - num of errors in the cone, num of correct CA cells
struct Utility {
    int nbits;
    Utility(int _nbits) : nbits(_nbits) {}
    int operator()(uint64_t p1, uint64_t p2) {
        // xor: same bit = 0, different = 1 : what we want
        uint64_t x = p1 ^ p2;
        // Trick to count the number of bits that are one
        x = (x & 0x5555555555555555ULL) + ((x & 0xaaaaaaaaaaaaaaaaULL) >> 1ULL);
        x = (x & 0x3333333333333333ULL) + ((x & 0xccccccccccccccccULL) >> 2ULL);
        x = (x & 0x0f0f0f0f0f0f0f0fULL) + ((x & 0xf0f0f0f0f0f0f0f0ULL) >> 4ULL);
        x = (x & 0x00ff00ff00ff00ffULL) + ((x & 0xff00ff00ff00ff00ULL) >> 8ULL);
        x = (x & 0x0000ffff0000ffffULL) + ((x & 0xffff0000ffff0000ULL) >> 16ULL);
        x = (x & 0x00000000ffffffffULL) + ((x & 0xffffffff00000000ULL) >> 32ULL);
        // return num bits that match as the utility (1's mean discrepancies)
        return nbits-(int)x;
    }
};


int main(int argc, char** argv) {

    // User-defined variables
    int past_depth = 4;
    int future_depth = 3;
    int max_steps = 500;
    int skip_steps = 0;
    int rule_number = -1;
    AutomataWorld::worldSize = 500;
    unsigned int seed = 42;
    float siglevel = 0.04321;

#ifndef NO_SDL
    bool batchMode = false;
#else
    bool batchMode = true;
#endif
    bool quietMode = false;
    bool zMode = false;
    bool oneByOne = false;
    bool noninc = false;
    bool savegraphics = true;
    bool computeCausalStatesLast = false;
    int c; opterr = 0;
    while ((c=getopt(argc,argv,"lzdhbqk:p:f:t:w:s:r:g:"))!=-1) switch(c) {
        case 'd': savegraphics = false; break;
        case 'b': batchMode = true; break;
        case 'q': quietMode = true; break;
        case 'z': quietMode = true; zMode=true; break;
        case 'l': computeCausalStatesLast = true; break;
        case 'g': {
            float arg = (float)atof(optarg);
            if (arg!=0.0f) siglevel = arg;
            break;
        }
        case 'k': {
            skip_steps = atoi(optarg);
            break;
        }
        case 'p': {
            int arg = atoi(optarg);
            if (arg!=0) past_depth=arg;
            break;
        }
        case 'f': {
            int arg = atoi(optarg);
            if (arg!=0) future_depth=arg;
            break;
        }
        case 't': {
            int arg = atoi(optarg);
            if (arg!=0) max_steps=arg;
            break;
        }
        case 's': {
            int arg = atoi(optarg);
            if (arg!=0) seed=arg;
            break;
        }
        case 'w': {
            int arg = atoi(optarg);
            if (arg!=0) AutomataWorld::worldSize=arg;
            break;
        }
        case 'r': {
            rule_number=atoi(optarg); // 0 is valid
            break;
        }
        case '?':
        case 'h':
        default:
            cout << "Usage:\n-p<integer>: past depth\n-f<integer>: future depth\n-t<integer>: max time\n-w<integer>: world width\n-r<integer>: rule number\n-k<integer>: Skip this number of steps in addition to the minimum of (past_depth+future_depth-1)\n-g<float>: chi-square significance level for clustering distributions\n-b: batch mode\n-q: quiet mode\n-z: really quiet (no output)\n-s<integer>: random seed\n-d: don't save graphics\n-l: compute causal states last (by splitting decisional states)" << endl;
            return 0;
    }

    // derived variables
    const int timeSize = past_depth + future_depth;

//    const int num_past_bits = 2*past_depth-1; // only the last row
    const int num_future_bits = future_depth*future_depth; // all bits in the cone

    if (rule_number==-1) rule_number = 110; // default

    if (noninc && oneByOne) {
        cout << "Options -1 and -n are incompatible." << endl;
        return 0;
    }

    if (!quietMode) {
        std::ios::sync_with_stdio(false); // GNU specific?
#ifndef NO_SDL
        cout << (batchMode?"batch":"graphics") << " mode" << endl;
#else
        cout << "Batch mode forced (graphics not compiled in, please use SDL)." << endl;
#endif
        cout << "Processing causal states " << (computeCausalStatesLast?"last":"first") << endl;
        cout << "past depth = " << past_depth << endl;
        cout << "future depth = " << future_depth << endl;
        cout << "max time = " << max_steps << endl;
        cout << "world size = " << AutomataWorld::worldSize << endl;
        cout << "rule number = " << rule_number << endl;
        cout << "random seed = " << seed << endl;
        cout << "significance for clustering = " << siglevel << endl;
        cout << "Processing observations " << (oneByOne?"one by one.":(noninc?"non-incrementally":"row by row.")) << endl;
    }

#ifndef NO_SDL
    if (!batchMode) {
        // Initialize the SDL library
        if( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
            cerr << "Couldn't initialize SDL: " << SDL_GetError() << endl;
            return 1;
        }

        // Clean up on exit
        atexit(SDL_Quit);

        // Initialize the display in a 32-bit mode,
        screen = SDL_SetVideoMode(AutomataWorld::worldSize*2, max_steps*2, 32, SDL_HWSURFACE);
        if (!screen) {
            cerr << "Couldn't initialize SDL: " << SDL_GetError() << endl;
            return 1;
        }
    }
#endif

    ofstream rawImage, cfiltImage, uImage, pfiltImage;

    if (savegraphics) {
        // create image files
        stringstream ss; ss << "rule" << rule_number << "_raw_" << AutomataWorld::worldSize << "x" << max_steps << ".pgm";
        rawImage.open(ss.str().c_str(), std::ios::binary);
        
        if (!quietMode) cout << "Saving images in " << ss.str().c_str() << ", " << flush;

        stringstream ss2; ss2 << "rule" << rule_number << "_cfilt_p" << past_depth << "f" << future_depth << "_" << AutomataWorld::worldSize << "x" << max_steps << ".pgm";
        cfiltImage.open(ss2.str().c_str(), std::ios::binary);
        if (!quietMode) cout << ss2.str().c_str() << ", " << flush;

        stringstream ss3; ss3 << "rule" << rule_number << "_util_p" << past_depth << "f" << future_depth << "_" << AutomataWorld::worldSize << "x" << max_steps << ".pgm";
        uImage.open(ss3.str().c_str(), std::ios::binary);
        if (!quietMode) cout << ss3.str().c_str() << ", " << flush;

        stringstream ss4; ss4 << "rule" << rule_number << "_pfilt_p" << past_depth << "f" << future_depth << "_" << AutomataWorld::worldSize << "x" << max_steps << ".pgm";
        pfiltImage.open(ss4.str().c_str(), std::ios::binary);
        if (!quietMode) cout << ss4.str().c_str() << endl;

        stringstream pgmheader; pgmheader << "P5 " << AutomataWorld::worldSize << " " << max_steps << " 255\n";
        rawImage << pgmheader.str();
        cfiltImage << pgmheader.str();
        uImage << pgmheader.str();
        pfiltImage << pgmheader.str();
    }

    if (!quietMode) {
        if (!savegraphics) cout << "Don't saving graphics" << endl;
        cout << endl;
    }

    // Create the Cellular Automata
    WrappingVector<AutomataWorld> spaceTime(timeSize);
    Rule rule(rule_number);
    
    DataSet dataset;

    // random initial condition
    generate(spaceTime[0].begin(), spaceTime[0].end(), RandFunctor(seed));

    int timeSlice = 0;
    //int displayCount = 0;

    // fill the first steps where complexity can't be computed
    for (int t = 1; t < timeSize-1; ++t) {
        int lastSlice = timeSlice++;
        // compute new state
        for (int i=0; i<AutomataWorld::worldSize; ++i) {
            spaceTime[timeSlice][i] = rule.apply(
                spaceTime[lastSlice][i-1],
                spaceTime[lastSlice][i],
                spaceTime[lastSlice][i+1]);
        }
    }

    // skip extra steps at the user request
    for (int t = 0; t < skip_steps; ++t) {
        int lastSlice = timeSlice++;
        if (timeSlice >= timeSize) timeSlice = 0;
        // compute new state
        for (int i=0; i<AutomataWorld::worldSize; ++i) {
            spaceTime[timeSlice][i] = rule.apply(
                spaceTime[lastSlice][i-1],
                spaceTime[lastSlice][i],
                spaceTime[lastSlice][i+1]);
        }
    }

    // Now feed the light cones to the analyzer
    if (!quietMode) cout << "Computing light cones" << endl;

    for (int t = 0; t < max_steps; ++t) {
        // build past/future light cones for each points at future_depth in the past
        int present = timeSlice - future_depth + 1;
        // do this for each point
        for (int i=0; i<AutomataWorld::worldSize; ++i) {
            char b = spaceTime[present][i] * 255;
            if (savegraphics) rawImage.write(&b, 1);
            if (!batchMode) drawGrayPixel(i, t, b);

            // Past Light Cone = plc
            // Past is fully determined by the rule and the last line in past
            uint64_t plc = 0;
            int nb = 0;
            int time = past_depth-1;
            for (int space = i-time; space <= i+time; ++space) {
                plc |= spaceTime[present-time][space] << nb;
                ++nb;
            }
            assert(nb==num_past_bits);

            // Future Light Cone = flc
            uint64_t flc = 0;
            nb = 0;
            for (int time = 0; time < future_depth; ++time) {
                for (int space = i-time; space <= i+time; ++space) {
                    flc |= spaceTime[present+time][space] << nb;
                    ++nb;
                }
            }
            assert(nb==num_future_bits);

            // collect all data
            dataset.push_back(make_pair(plc, flc));
        }

#ifndef NO_SDL
        if (!batchMode) {
            SDL_Event event;
            SDL_PollEvent(&event);
            switch (event.type) {
                case SDL_KEYDOWN: if (event.key.keysym.sym!=SDLK_ESCAPE) break;
                case SDL_QUIT: return 0;
            }
        }
#endif
#ifndef NO_SDL
        if (!batchMode) SDL_UpdateRect(screen, 0, t, AutomataWorld::worldSize*2, 1);
#endif

        int lastSlice = timeSlice++;
        if (timeSlice >= timeSize) timeSlice = 0;
        // compute new state
        for (int i=0; i<AutomataWorld::worldSize; ++i) {
            spaceTime[timeSlice][i] = rule.apply(
                spaceTime[lastSlice][i-1],
                spaceTime[lastSlice][i],
                spaceTime[lastSlice][i+1]);
        }

    }

    if (!quietMode) cout << "Analysing..." << endl;

    typedef DiscreteDistributionManager<DataSet> DistManager;
    typedef DecisionalStatesAnalyser<DistManager, DataSet, Optimiser, TransitionFeeder> DSA;
    Optimiser optimiser(num_future_bits);
    DistManager distManager(dataset);
    TransitionFeeder transitionFeeder(dataset, AutomataWorld::worldSize, max_steps, past_depth, rule_number);
    DSA analyser(distManager, transitionFeeder, optimiser);
    bool consistent;

    if (!computeCausalStatesLast) {
        consistent = analyser.computeCausalStates(ConnectedClustering<>(), DistManager::DistributionMatcher(1.0f - siglevel));
        cout << "Num causal states: " << analyser.causalStates.size() << ", consistent flag=" << consistent << endl;
    }

    cout << "Applying utility..." << flush;
    analyser.applyUtility( Utility(num_future_bits) );
    cout << " Done!" << endl;

    consistent = analyser.computeIsoPredictionStates(ConnectedClustering<>(), std::equal_to<Optimiser::PredictionSet>());
    cout << "Num iso prediction states: " << analyser.isoPredictionStates.size() << ", consistent flag=" << consistent << endl;

    consistent = analyser.computeIsoUtilityStates(ConnectedClustering<>(), ApproximateUtilityMatcher<double>());
    cout << "Num iso utility states: " << analyser.isoUtilityStates.size() << ", consistent flag=" << consistent << endl;
    
    analyser.computeDecisionalStates();
    cout << "Num decisional states: " << analyser.decisionalStates.size() << endl;

    if (computeCausalStatesLast) {
        consistent = analyser.computeCausalStates(ConnectedClustering<>(), DiscreteDistributionManager<DataSet>::ChiSquareMatcher(1.0f - siglevel) );
        cout << "Num causal states: " << analyser.causalStates.size() << ", consistent flag=" << consistent << endl;
    }

    analyser.buildCausalStateGraph();
    stringstream ssdot; ssdot << "rule" << rule_number << "_graph.dot";
    ofstream dot(ssdot.str().c_str());
    analyser.writeCausalStateGraph(dot);
    dot.close();

    analyser.buildIsoPredictionStateGraph();
    analyser.buildIsoUtilityStateGraph();
    analyser.buildDecisionalStateGraph();

    // Plot by statistical complexity
    // complexity = -log2(p) = -log2(count/total) = -log2(count) + log2(total)
    {
    float compDenum = log2f(dataset.size());
    float minComp = compDenum - log2f(analyser.causalStatesCounts.maxcount);
    float maxComp = compDenum - log2f(analyser.causalStatesCounts.mincount);
cout << "Min statistical complexity:" << minComp << endl;
cout << "Max statistical complexity:" << maxComp << endl;
    float compDiff = maxComp - minComp;
    if (compDiff<=0) compDiff = 1; // only one state, numerical rounding

    DataSet::iterator it = dataset.begin();
    for (int t = 0; t < max_steps; ++t) for (int i=0; i<AutomataWorld::worldSize; ++i) {
        assert(it!=dataset.end());
        uint64_t plc = it->first;
        shared_ptr<DSA::CausalState> cstate = analyser.getCausalState(plc);
        assert(cstate != 0);
        float comp = compDenum - log2f(cstate->count);
        float scaledComp = (comp - minComp) / compDiff; // between 0 and 1
        char grayScale = 255 - (int)(scaledComp * 255.99f);
        if (!batchMode) drawGrayPixel(AutomataWorld::worldSize+i, t, grayScale);
        if (savegraphics) cfiltImage.write(&grayScale, 1);
        ++it;
    }
    }

    // Plot by utility
    {
    double minU = numeric_limits<double>::max();
    double maxU = -numeric_limits<double>::max();
    for (DSA::IsoUtilityStates::iterator it = analyser.isoUtilityStates.begin(); it != analyser.isoUtilityStates.end(); ++it) {
        if ((*it)->utility < minU) minU = (*it)->utility;
        if ((*it)->utility > maxU) maxU = (*it)->utility;
    }
    double deltaU = maxU - minU; if (deltaU<=0.0) deltaU = 1.0;
    DataSet::iterator it = dataset.begin();
    for (int t = 0; t < max_steps; ++t) for (int i=0; i<AutomataWorld::worldSize; ++i) {
        assert(it!=dataset.end());
        uint64_t plc = it->first;
        shared_ptr<DSA::IsoUtilityState> isoUtilityState = analyser.getIsoUtilityState(plc);
        assert(isoUtilityState != 0);
        float scaledU = (isoUtilityState->utility - minU) / deltaU; // between 0 and 1
        char grayScale = 255 - (int)(scaledU * 255.99f);
        if (!batchMode) drawGrayPixel(i, max_steps+t, grayScale);
        if (savegraphics) uImage.write(&grayScale, 1);
        ++it;
    }
    }

    // Plot by difficulty to get the best prediction
    {
    float compDenum = log2f(dataset.size());
    float minComp = compDenum - log2f(analyser.isoPredictionStatesCounts.maxcount);
    float maxComp = compDenum - log2f(analyser.isoPredictionStatesCounts.mincount);
cout << "Min complexity to get the best prediction:" << minComp << endl;
cout << "Max complexity to get the best prediction:" << maxComp << endl;
    float compDiff = maxComp - minComp;
    if (compDiff<=0) compDiff = 1; // only one state, numerical rounding

    DataSet::iterator it = dataset.begin();
    for (int t = 0; t < max_steps; ++t) for (int i=0; i<AutomataWorld::worldSize; ++i) {
        assert(it!=dataset.end());
        uint64_t plc = it->first;
        shared_ptr<DSA::IsoPredictionState> isoPredictionState = analyser.getIsoPredictionState(plc);
        assert(isoPredictionState != 0);
        float comp = compDenum - log2f(isoPredictionState->count);
        float scaledComp = (comp - minComp) / compDiff; // between 0 and 1
        char grayScale = 255 - (int)(scaledComp * 255.99f);
        if (!batchMode) drawGrayPixel(AutomataWorld::worldSize+i, max_steps+t, grayScale);
        if (savegraphics) pfiltImage.write(&grayScale, 1);
        ++it;
    }
    }

#ifndef NO_SDL
        if (!batchMode) SDL_UpdateRect(screen, 0, 0, AutomataWorld::worldSize*2, max_steps*2);
#endif

#ifndef NO_SDL
    if (!batchMode) while (1) {
        SDL_Event event;
        SDL_WaitEvent(&event);
        switch (event.type) {
            case SDL_KEYDOWN: if (event.key.keysym.sym!=SDLK_ESCAPE) break;
            case SDL_QUIT: return 0;
        }
    }
#endif

}
