SWIG = swig -c++ -python
SFLAGS = -c++ -python
GCC = g++
CFLAGS += -O3 -fPIC -lm

# Here we set IDIR variable, which gives path to a directory with Python.h file. 
# The path will be different on different machines and can also be different for different 
# python distributions installed on the same computer. There are a few ways you can use to find 
# the appropriate Python.h file:
#  1) run:
#
#  python -c 'import sys; print(sys.executable)'
#  
#  and try finding going up and down in directory tree; in the case of popular anaconda 
#  distribution you will most probably go to ../../ and then to include/python3.5m/
#  
#  2) run:
#   
#  locate Python.h
#
#  and if many files are found, then choose the correct one
#
#  3) combine 1) and 2) to be sure you're right.
#
#
#  Then set the IDIR variable e.g.:
#  
#  IDIR = /usr/custom/Anaconda/include/python3.5m/

ifndef IDIR
$(error IDIR is not set. See the Makefile instructions on how to set it properly)
endif

.PHONY: all
all: _VBBinaryLensingLibrary.so

%_wrap.cxx: %.i
	$(SWIG) $^

%_wrap.o: %_wrap.cxx
	if [ ! -f $(IDIR)/Python.h ] ; then echo Problem IDIR variable: $(IDIR)/Python.h does not exist && exit 1  ; fi;
	$(GCC) $(CFLAGS) -c $^ -I $(IDIR)

%.o: %.cpp
	$(GCC) $(CFLAGS) -c $^

_%.so: %.o %_wrap.o
	$(GCC) -shared $^ -o $@

.PHONY: clean
clean:
	$(RM) _*.so *.pyc *.py *_wrap.cxx *.o

