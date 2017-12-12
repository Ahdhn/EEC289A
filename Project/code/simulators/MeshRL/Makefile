# $Id$ 
# This should be a pretty standard makefile for all PG projects

# Uncomment for ATLAS neural net. Make sure similar lines are uncommented
# in pg library makefile

# ATLASCFLAGS=-DHAVE_ATLAS
# ATLASLIBSPATH=-L$(HOME)/src/ATLAS/lib/OSX_PPCG4AltiVec
# ATLASLIBS=-lcblas -latlas

TARGET=MeshRLApp
PG=../../
BOOSTINC=
LIBS=-lpg $(ATLASLIBS)
PGLIB=$(PG)/libpg.a
PGINC=-I$(PG)
OPTIMS=-O3 -DNDEBUG
DEBUG=-g
CFLAGS=-Wall $(DEBUG) $(OPTIMS) $(ATLASCFLAGS) $(PGINC) $(BOOSTINC)
LFLAGS=$(OPTIMS) -L$(PG) 
CCDEP=$(CXX) -MM

#
# Editing below this line not advised unless you know what's going on.
#

CXX=g++
OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))

.PHONY: clean distclean $(PGLIB)


$(TARGET) : $(OBJS) $(PGLIB)
	$(CXX) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

$(PGLIB) :
	make -C $(PG)

%.o : %.cc 
	$(CXX) $(CFLAGS) -c $<


clean:
	rm -fr *.o *~ *.bak *.d

distclean: clean
	rm -fr $(TARGET)

# All this stuff looks after automatic header dependencies

%.d: %.cc 
	@echo "Finding dependencies for $<"
	@$(SHELL) -ec '$(CCDEP) $(CFLAGS) $(INCLUDES) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o  $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(OBJS:%.o=%.d)
