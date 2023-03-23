#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Wed 30 Dec 2020 09:32:34 PM CST
#-------------------------------------------------------------------------------

FC := gfortran
CLFS := -fopenmp

CFLAGS :=
LFLAGS :=

SRCDIR := src
OBJDIR := obj
BINDIR := bin

STATIC :=
DEBUG :=
WITHQ :=
IEEE := ON
COLORPRINT :=
DFLAG_LIST := DEBUG WITHQ IEEE COLORPRINT

DERIV := _BETA # acceptable value: _BETA, _ALPHA, _RHO
VFLAG_LIST := DERIV

DFLAGS := $(foreach flag, $(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
VFLAGS := $(foreach flag, $(VFLAG_LIST),$(if $($(flag)),-D$(flag)=$($(flag)),))
DFLAGS := $(DFLAGS) $(VFLAGS)

ifdef STATIC
  CLFS += -fPIC
  ifeq ($(notdir $(FC)), gfortran)
    LFLAGS += -static-libgfortran
  else ifeq ($(notdir $(FC)), ifort)
    LFLAGS += -static-intel
  endif
endif

ifdef DEBUG
  CLFS += -O0
  CFLAGS += -g
  ifeq ($(notdir $(FC)), gfortran)
	CLFS += -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow
  else ifeq ($(notdir $(FC)), ifort)
	CLFS += -warn all -check all -fpe:0
  endif
else
  CLFS += -O2
endif

ifeq ($(notdir $(FC)), gfortran)
  CFLAGS += -J$(OBJDIR)
else ifeq ($(notdir $(FC)), ifort)
  CFLAGS += -module $(OBJDIR)
endif

CFLAGS += -I$(OBJDIR)
CFLAGS += $(DFLAGS)
CFLAGS := $(strip $(CFLAGS))
LFLAGS := $(strip $(LFLAGS))

SRC_MOD := constant.F90 math.F90 parameter.F90 grtc.F90
SRC_GRT := main.F90
EXE_GRT := grtcKer

OBJ_MOD := $(foreach f,$(SRC_MOD),$(OBJDIR)/$(f:.F90=.o))
OBJ_GRT := $(OBJDIR)/$(SRC_GRT:.F90=.o)

all : prep link

prep:
	@mkdir -p $(BINDIR) $(OBJDIR)

link: $(BINDIR)/$(EXE_GRT)

dispec:
	python setup.py build
	python setup.py install --user

$(BINDIR)/$(EXE_GRT): $(OBJ_MOD) $(OBJ_GRT)
	$(FC) $(LFLAGS) $(CLFS) $^ -o $@

clear: clean
	-rm -f $(BINDIR)/*
	-rm -rf build/*

clean:
	-rm -f $(OBJDIR)/*

vpath %.F90 $(SRCDIR)
.SUFFIXES: .F90 .o

$(OBJDIR)/%.o : %.F90
	$(FC) $(CFLAGS) $(CLFS) $(SRCDIR)/$(<F) -c -o $@

# vim:ft=make noet
