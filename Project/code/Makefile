# PG $Id: Makefile.default 127 2007-09-10 17:20:04Z daa $
# Most recent svn revision number released
REVISION=126

# Uncomment below to define what you have.
# You need Boost anyway under any circumstances

# Boost sandbox is a set of user contributed bindings for CLAPACK 
# ATLAS. This is not needed for basic compilation, but is needed
# if you intend ot use ATLAS *OR* CLAPACK
#HAVESANDBOX=-DHAVE_BOOST_SANDBOX

# Atlas speeds up large neural nets
#HAVEATLAS=-DHAVE_ATLAS

# CLAPACK provides svd and other linear algebra
#HAVECLAPACK=-DHAVE_CLAPACK

# If not in standard places, point to Boost, Atlas, and CLAPACK
#INCLUDES

# These optimisations make a HUGE difference, but leave them off
# until your sure your code works
OPTIMS=-O3 -DNDEBUG
DEBUG=-g -O0

##################### DO NOT EDIT BELOW ############################
# Unless you know what you're doing.

TARGET=libpg.a

EXTRADEFINES=$(HAVEATLAS) $(HAVECLAPACK) $(HAVESANDBOX)
CFLAGS=-Wall $(DEFINES) $(EXTRADEFINES) $(DEBUG) $(OPTIMS) $(ATLASCFLAGS) $(INCLUDES)
OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
TESTS_CC=$(wildcard tests/*.cc)
TESTS=$(patsubst %.cc,%,$(TESTS_CC))
RM =rm -f
AR=ar -cr
RANLIB=ranlib
CCDEP=$(CXX) -MM

# Generate a new rule for each element of list TESTS

$(TARGET): $(OBJS)
	$(AR) $@ $(OBJS)  
	$(RANLIB) $@

tests: $(TESTS)

$(TESTS) : % : %.cc $(OBJS)  
	$(CXX) $(INCLUDES) -o $@ $^ 

%.o : %.cc
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean distclean publish
clean:
	$(RM) *.o *~ *.d tests/*.o tests/*~

distclean: clean
	$(RM) $(TARGET)

publish: distclean
	echo "Don't forget to run Doxygen!"
	find . -name 'Makefile' -execdir make clean \;
	tar -X exclude_from_tar.txt -czf ../libpg-$(REVISION).tgz ../libpg

# All this stuff looks after automatic header dependencies

%.d: %.cc 
	@echo "Finding dependencies for $<"
	@$(SHELL) -ec '$(CCDEP) $(CFLAGS) $(INCLUDES) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o  $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(OBJS:%.o=%.d)
