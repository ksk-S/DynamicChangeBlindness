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

#include <iostream>
#include <fstream>

#include <boost/random.hpp>

#include <stdlib.h>

using namespace std;
using namespace boost;

int main (int argc, char** argv) {
    
    int ndata;

    if (argc>2) {
        ndata = atoi(argv[2]);
    } else {
        cout << "Generates a symbolic series according to the Even Process." << endl;
        cout << "Arguments are: output_file_name ndata [random_seed]" << endl;
        return 0;
    }
    unsigned int seed = 42;
    if (argc>3) {
        seed = atoi(argv[3]);
    }

    cout << "Generating " << ndata << " samples from the Even Process in the file " << argv[1] << " using the random seed " << seed << endl;

    ofstream epseries(argv[1]);

    mt19937 rng(seed);
    uniform_int<> randomRange(0,1);
    variate_generator<mt19937&, uniform_int<> > randomBinary(rng,randomRange);
    
    // Generate data set -- Even process
    
    int x = randomBinary();
    int state = (x==0) ? 0 : randomBinary();
    epseries << x;
    for (int i=1; i<ndata; ++i) {
        if (state == 1) {x = 1; state = 0;}
        else {x = randomBinary(); state = x;}
        epseries << x;
    }
    
    epseries << endl;
    epseries.close();

    return 0;
}
