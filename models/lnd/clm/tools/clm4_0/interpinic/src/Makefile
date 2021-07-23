# Makefile for interpinic

EXENAME = ../interpinic_cheyenne

# Set shared memory to off by default
ifeq ($(SMP),$(null))
  SMP := FALSE
endif

# Set optimization on by default
ifeq ($(OPT),$(null))
  OPT := TRUE
endif

include Makefile.common
