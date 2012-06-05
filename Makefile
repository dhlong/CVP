SYSTEM     = x86-64_rhel4.0_3.4
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /app1/ia32-64/iLog/cplex101
CONCERTDIR    = /app1/ia32-64/iLog/concert23
EIGENDIR      = /home/svu/isedhl/eigen-eigen-6e7488e20373
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++
CC  = gcc

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -O3 -fPIC -fexceptions -DIL_STD -g

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread 

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(EIGENDIR)

all: CVP

ALL_O = main.o cvp.o network.o function.o dijkstra.o

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------

CVP: $(ALL_O)
	$(CCC) $(CCFLAGS) $(ALL_O) -o CVP $(CCLNFLAGS) -lrt
%.o: %.cpp
	$(CCC) -c $(CCFLAGS) $< -o $@
