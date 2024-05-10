#################################################
# Makefile for frodaN
# June 16, 2010
#
##################################################

# List external programs that need to be compiled
# independently.
EXT_PROGS := FIRST/src Seqaln/src ticpp-read-only frodaN

BINARY := FIRST/src/FIRST frodaN/main frodaN/util/projections2d frodaN/util/matchState1ToState2 frodaN/util/assignAtomTypes frodaN/util/measureRMSF frodaN/util/measureDistanceTraveled Seqaln/src/alignment_amino
prefix = $(shell pwd)/bin

# Define C++11 flag
CXX11_FLAG := -std=c++11

# The default make will depend on FIRST and any other
# INDEPENDENT executables. 
all:  external_progs

clean:
	@for i in $(EXT_PROGS); do \
        (cd $$i;  rm -f *.o *.d; $(MAKE) clean); done

############################################################
# Any independent executable programs 
# should have their own makefile and be compiled
# recursively. 
external_progs:
	@for i in $(EXT_PROGS); do \
		echo "Compiling in $$i..."; \
	        (cd $$i; $(MAKE) CXXFLAGS=$(CXX11_FLAG) $(MFLAGS)); done
	@for i in $(BINARY); do \
		echo "Moving $$i to $(prefix) ...";  \
		(mv $$i $(prefix);); done

