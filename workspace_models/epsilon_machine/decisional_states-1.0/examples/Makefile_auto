# This file is part of the decisional state reconstruction algorithm
# technique exposed in "Decisional States", by Nicolas Brodu.
#
#     This library is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 2.1 of the License, or (at your option) any later version.
#
#     This library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with this library; if not, write to the Free
#     Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
#     MA  02110-1301  USA
#
# See http://nicolas.brodu.numerimoire.net/en/programmation/decisional_states/index.html
# for more information and possibly updates.
#
# Copyright holder: Nicolas Brodu <nicolas.brodu@numerimoire.net>
# File release Date: February 09


##################################################
# CONFIGURE THESE VARIABLES
# or set them on the command line, like:
#   make ImageFilterGray DEBUG=1
##################################################

# Uncomment this if you don't want to use OPENMP (OPENMP makes the program parallel)
#NO_OPENMP = 1

# Uncomment this if you don't have SDL installed (SDL is for graphics)
# SDL is used for the automaton example
#NO_SDL = 1

# Uncomment this if you wish to use the ANN library for neighborhood queries
# in the image filter example
# This is for demonstration purpose only, you won't gain much, if any, performance.
# Just uncompress the official ANN archive in the example directory
# This was tested with ANN versions 1.1.1 and 1.1.2
#USE_ANN=1.1.2

#DEBUG=1

##################################################

# You may copy this makefile in an object directory and
# point this variable to the source directory. This is especially useful for
# compiling with different targets/options in a shared filesystem.
SRCDIR=.
INCDIR=$(SRCDIR)/../include

# set this variable to another compiler, or g++ is the default
CXX = g++

# optimization options
ifdef DEBUG
CXXFLAGS += -pipe -g
else
CXXFLAGS += -g -pipe -O3 -march=native -DNDEBUG
endif

# Correctness options
CXXFLAGS += -Wall
CPPFLAGS += -I$(SRCDIR) -I$(INCDIR)

ifdef NO_SDL
CXXFLAGS += -DNO_SDL
else
CXXFLAGS += -lSDL
endif

ifndef NO_OPENMP
CXXFLAGS += -fopenmp
endif

ifdef USE_ANN
ANNOBJ = -L$(SRCDIR)/ann_$(USE_ANN)/lib -lANN
ANNINC = -I$(SRCDIR)/ann_$(USE_ANN)/include
ANNTARGET = ANN
else 
ANNOBJ = 
ANNINC = 
ANNTARGET = 
endif

all: EvenProcess SymbolicSeries CellularAutomaton TimeSeries ImageFilter

ANN:
	(cd $(SRCDIR)/ann_$(USE_ANN); make linux-g++)

CellularAutomaton: $(SRCDIR)/CellularAutomaton.cpp $(INCDIR)/
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SRCDIR)/CellularAutomaton.cpp -o CellularAutomaton

ImageFilter: $(ANNTARGET) $(SRCDIR)/ImageFilter.cpp $(INCDIR)/
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ANNINC) $(SRCDIR)/ImageFilter.cpp -ljpeg -lpng $(ANNOBJ) -o ImageFilter

EvenProcess: $(SRCDIR)/EvenProcess.cpp $(INCDIR)/
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SRCDIR)/EvenProcess.cpp -o EvenProcess

TimeSeries: $(SRCDIR)/TimeSeries.cpp $(INCDIR)/
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SRCDIR)/TimeSeries.cpp -o TimeSeries

SymbolicSeries: $(SRCDIR)/SymbolicSeries.cpp $(INCDIR)/
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SRCDIR)/SymbolicSeries.cpp -o SymbolicSeries

clean:
	rm -f EvenProcess SymbolicSeries CellularAutomaton TimeSeries ImageFilter
