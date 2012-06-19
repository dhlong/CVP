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
GRBDIR        = /home/svu/isedhl/gurobi500/linux64/
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
GRBLIBLIB     = $(GRBDIR)/bin
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
GRBINCDIR     = $(GRBDIR)/include

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -L$(GRBLIBDIR) -lgurobi50 -lm -lpthread 

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(EIGENDIR) -I$(GRBINCDIR)

all: CVP matrix_test CVP_alter

CVP_O = main.o cvp.o network.o function.o dijkstra.o cvputility.o
MATRIX_TEST_O = matrix_test.o network.o dijkstra.o cvputility.o
CVP_ALTER_O = network.o dijkstra.o cvputility.o function.o cvp_alter.o my_sparse_vector.o solver.o

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf CVP matrix_test  CVP_alter
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp
	/bin/rm *~

# ------------------------------------------------------------

CVP: $(CVP_O)
	@echo Compiling CVP
	@$(CCC) $(CCFLAGS) $(CVP_O) -o CVP $(CCLNFLAGS) -lrt
matrix_test: $(MATRIX_TEST_O)
	@echo Compiling matrix_test
	@$(CCC) $(CCFLAGS) $(MATRIX_TEST_O) -o $@ $(CCLNFLAGS) -lrt
CVP_alter: $(CVP_ALTER_O)
	@echo Compiling CVP Alter
	@$(CCC) $(CCFLAGS) $(CVP_ALTER_O) -o $@ $(CCLNFLAGS) -lrt
%.o: %.cpp
	@echo Compiling $<
	@$(CCC) -c $(CCFLAGS) $< -o $@
