#################################################
#                                               #
#                   Makefile                    #
#                                               #
# Note: flags that can be modified for personal #
# convenience are indicated in the following by #
# the presence of the string CHANGE IF YOU WANT #
#                                               #
#################################################


## Usual variables

# Compilator
CC = gcc
# Flags for compilator
# CCOPTIONS = -O 

MAKE = make
RM   = rm -f
LN   = ln -s
MV   = mv
CP   = cp

OBJEXT     = o
SRCEXT = c
LIBEXT     = a
SHAREDEXT  = so

## Names of the methods

# Name of the method: CHANGE IF YOU WANT

MF=multifractal
MSM=msm
WTMM=wtmm

# Names of the directories where the methods are stored

PACKWTMM=wtmm
PACKMSM=msm
PACKMF=multifractal

# Anyway, be careful that if you change this variable, you 
# change some flags of the Makefile for compilation

## Variables for location

# Relative paths to the programs and source files

REL=..
CUR=.

# Location of the sources

DIRMSM=$(REL)/$(PACKMSM)
DIRWTMM=$(REL)/$(PACKWTMM)
DIRMF=$(REL)/$(PACKMF)

# Where some inputs/outputs are stored
NAMEBIN=bin
NAMEINC=include
NAMELIB=lib

DIRBIN=$(REL)/$(NAMEBIN)
DIRINC=$(REL)/$(NAMEINC)
DIRLIB=$(REL)/$(NAMELIB)

## Main variables, definitions and rules for the Makefile

# List of source files to compile
# In each papckage directory, the FILELIST file gives the list
# of files that need to be compiled.
FILELIST=FileList

include $(DIRMSM)/$(FILELIST)
# provides the variables MSMLIST, MSMSRCS, MSMLOBJS and MSMHEADS

include $(DIRWTMM)/$(FILELIST)
# provides the variables WTMMLIST, WTMMSRCS, WTMMLOBJS and WTMMHEADS

MFLIST = $(MSMLIST) $(WTMMLIST) 
MFLOBJS = $(MSMLOBJS) $(WTMMLOBJS)
MFHEADS = $(MSMHEADS) $(WTMMHEADS)

include $(CUR)/$(FILELIST) 
# provides the variables MFSRC and MFOBJ, MSMSRC and MSMOBJ, 
# WTMMSRC and WTMMOBJ

# Flag of shared library
# Name of extension for shared libraries: CHANGE IF YOU WANT
LIBMF = lib$(MF)
LFLAGMF = -l$(MF)
LDFLAGSMF = -L$(DIRLIB) -L$(CUR) $(LFLAGMF) -lm
# predominance of directory DIRLIB over CUR for the choice of the library

# Flag of distinguished libraries
LIBMSM = lib$(MSM)
LFLAGMSM = -l$(MSM)
LDFLAGSMSM = -L$(DIRLIB) -L$(CUR) $(LFLAGMSM) -lm

LIBWTMM = lib$(WTMM)
LFLAGWTMM = -l$(WTMM)
LDFLAGSWTMM = -L$(DIRLIB) -L$(CUR) $(LFLAGWTMM) -lm

# List of compiled libraries
LLIBS= $(LIBMF) $(LIBMSM) $(LIBWTMM)

## CRUCIAL FLAG FOR COMPILATION OF SHARED LIBRARY!!!
#
# Note that the objects that will be compiled for the stand alone 
# distributions and the combined distribution are different because
# of the presence of the flags FLAG_MSM and FLAG_WTMM in the code
#
# Flags for the shared library
FLAG_WTMM=FLAG_WTMM
FLAG_MSM=FLAG_MSM
DFLAGSMF= -D$(FLAG_WTMM) -D$(FLAG_MSM)
#
## !!!

# Name of the extension for the main program: CHANGE IF YOU WANT
PROCESSOR=_processor
MFPROG=$(MF)$(PROCESSOR)
MSMPROG=$(MSM)$(PROCESSOR)
WTMMPROG=$(WTMM)$(PROCESSOR)

LPROGS= $(MFPROG) $(MSMPROG) $(WTMMPROG)

# Options
INCLUDE= -I$(DIRINC) -I./
CFLAGS=  $(INCLUDE) 
# $(CCOPTIONS)

# MAIN RULES FOR THE MAKEFILE

# Generic target
all: clean libs progs
# not that we clean the directories after the creation
# of the stand-alone libraries

# fixed flag for make

lib:$(LIBMF)

libs: $(LLIBS)

progs: $(LPROGS)

msm:$(MSMPROG)

wtmm:$(WTMMPROG)

multifractal:$(MFPROG)

# Help for the Makefile
help:
	@echo "Main flags for the make command"
	@echo "   - make lib : create the library" $(LIBMF) "only"
	@echo "   - make libs : create all the libraries" $(LLIBS)
	@echo "   - make wtmm : create the program" $(WTMMPROG)
	@echo "   - make msm : create the program" $(MSMPROG)
	@echo "   - make multifractal : create program" $(MFPROG)
	@echo "   - make progs : create all the programs" $(LPROGS)
	@echo "   - make install : install the progs and the libraries"
	@echo "   - make clean : clean the object files"
	@echo "   - make cleanall : clean the object files, the libraries and the programs"
	@echo " You can type directly:"
	@echo "     make all"
	@echo " then:"
	@echo "     make install"
	@echo " to compile and install all programs and libraries."

# Stand-alone programs
$(WTMMPROG): $(LIBWTMM) $(WTMMMAINOBJ) 
	$(CC) $(CFLAGS) -c $(WTMMMAINSRC) -o $(WTMMMAINOBJ) 
	$(CC) $(CFLAGS) $(WTMMMAINOBJ)  -o $@ $(LDFLAGSWTMM)

$(MSMPROG): $(LIBMSM) $(MSMMAINOBJ) 
	$(CC) $(CFLAGS) -c $(MSMMAINSRC) -o $(MSMMAINOBJ) 
	$(CC) $(CFLAGS) $(MSMMAINOBJ)  -o $@ $(LDFLAGSMSM)

# Stand-alone libraries
$(LIBMSM): cleanmsm $(MSMHEADS)
	@for f in $(MSMLIST); do \
	echo $(CC) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	$(CC) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	done; 
	ar cr $@.$(LIBEXT) $(MSMLOBJS)
	ranlib $@.$(LIBEXT)

$(LIBWTMM): cleanwtmm $(WTMMHEADS)
	@for f in $(WTMMLIST); do \
	echo 	$(CC) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	$(CC) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	done; 
	ar cr $@.$(LIBEXT) $(WTMMLOBJS)
	ranlib $@.$(LIBEXT)

# Shared program
$(MFPROG): $(LIBMF) $(MFMAIN)
	$(CC) $(DFLAGSMF) $(CFLAGS) -c $(MFMAINSRC) -o $(MFMAINOBJ) 
	$(CC) $(DFLAGSMF) $(CFLAGS) $(MFMAINOBJ)  -o $@ $(LDFLAGSMF) 

# Shared library
$(LIBMF): cleanwtmm cleanmsm $(WTMMHEADS) $(MSMHEADS)
	@for f in $(MFLIST); do \
	echo $(CC) $(DFLAGSMF) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	$(CC) $(DFLAGSMF) $(CFLAGS) -c $$f.$(SRCEXT) -o $$f.$(OBJEXT); \
	done; 
	ar cr $@.$(LIBEXT) $(MFLOBJS)
	ranlib $@.$(LIBEXT)

install:
	@if [ -e $(LIBMF).$(LIBEXT) ]; then \
		$(MV) $(LIBMF).$(LIBEXT) $(DIRLIB); \
	fi 
	@if [ -e $(MFPROG) ]; then \
		$(MV) $(MFPROG) $(DIRBIN); \
	fi 
	@if [ -e $(LIBWTMM).$(LIBEXT) ]; then \
		$(MV) $(LIBWTMM).$(LIBEXT) $(DIRLIB); \
	fi 
	@if [ -e $(WTMMPROG) ]; then \
		$(MV) $(WTMMPROG) $(DIRBIN); \
	fi 
	@if [ -e $(LIBMSM).$(LIBEXT) ]; then \
		$(MV) $(LIBMSM).$(LIBEXT) $(DIRLIB); \
	fi 
	@if [ -e $(MSMPROG) ]; then \
		$(MV) $(MSMPROG) $(DIRBIN); \
	fi 


cleanmsm:
	$(RM) $(DIRMSM)/*.$(OBJEXT)

cleanwtmm:
	$(RM) $(DIRWTMM)/*.$(OBJEXT)

clean: cleanmsm cleanwtmm 
	$(RM) $(CUR)/*.$(OBJEXT) 

cleanall: clean
	@if [ -e $(DIRLIB)/$(LIBMF).$(LIBEXT) ]; then \
		echo $(RM) $(DIRLIB)/$(LIBMF).$(LIBEXT); \
		$(RM) $(DIRLIB)/$(LIBMF).$(LIBEXT); \
	fi 
	@if [ -e  $(DIRBIN)/$(MFPROG) ]; then \
		echo $(RM) $/$(MFPROG) ; \
		$(RM) $/$(MFPROG) ; \
	fi 
	@if [ -e $(DIRLIB)/$(LIBWTMM).$(LIBEXT) ]; then \
		echo $(RM) $(DIRLIB)/$(LIBWTMM).$(LIBEXT); \
		$(RM) $(DIRLIB)/$(LIBWTMM).$(LIBEXT); \
	fi 
	@if [ -e $(DIRBIN)/$(WTMMPROG) ]; then \
		echo $(RM) $(DIRBIN)/$(WTMMPROG) ; \
		$(RM) $(DIRBIN)/$(WTMMPROG) ; \
	fi 
	@if [ -e $(DIRLIB)/$(LIBMSM).$(LIBEXT) ]; then \
		echo $(RM) $(DIRLIB)/$(LIBMSM).$(LIBEXT) ; \
		$(RM) $(DIRLIB)/$(LIBMSM).$(LIBEXT) ; \
	fi 
	@if [ -e  $(DIRBIN)/$(MSMPROG) ]; then \
		echo $(RM)  $(DIRBIN)/$(MSMPROG); \
		$(RM)  $(DIRBIN)/$(MSMPROG); \
	fi 


