# Sample Makefile to compile ProBT examples.
# Probayes, 2007
# 
# Targets:

PROBT_INCLUDE=/home/esau/updates/probt/include
PROBT_LIB=/home/esau/updates/probt/lib

CXX=g++ -Wall

# Environment variable containing the path to the dynamically loaded libraries
ifeq ($(shell uname -s),Darwin)
	# libraries for Mac OS X
	DLL_PATH_VAR := DYLD_LIBRARY_PATH
else
	# libraries for Linux
	DLL_PATH_VAR := LD_LIBRARY_PATH
endif

all: run

# Compiling a C++ source file: add the path to the ProBT includes, the path to the ProBT libraries,
# and link with the ProBT libraries.
redbayesiana: RedBayesiana.cpp
	$(DLL_PATH_VAR)=$(PROBT_LIB):${$(DLL_PATH_VAR)} $(CXX) -I$(PROBT_INCLUDE) RedBayesiana.cpp -L$(PROBT_LIB) -lspl -lm -o redbayesiana

# Running an example: adjust the LD_LIBRARY_PATH environment variable to point to the directory
# containing the ProBT libraries.
run: redbayesiana
	@echo ">>> Running the example. It will output bayesian_network.fig (figure in in xfig format)."
	$(DLL_PATH_VAR)=$(PROBT_LIB):${$(DLL_PATH_VAR)} ./redbayesiana

show: run
	xfig red_bayesiana.fig

.PHONY: run show

