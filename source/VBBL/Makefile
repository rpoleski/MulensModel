SWIG = swig -c++ -python
SFLAGS = -c++ -python
GCC = g++
CFLAGS += -O3 -fPIC -lm

# The following directory may be difference in difference machines. To find the correct one in your machine, do "locate Python.h" #
IDIR = /usr/custom/Anaconda3/include/python3.5m/
# /usr/custom/Anaconda-2.0.1/include/python2.7/

.PHONY: all
all: _VBBinaryLensingLibrary.so

%_wrap.cxx: %.i
	$(SWIG) $^

%_wrap.o: %_wrap.cxx
	$(GCC) $(CFLAGS) -c $^ -I $(IDIR)

%.o: %.cpp
	$(GCC) $(CFLAGS) -c $^

_%.so: %.o %_wrap.o
	$(GCC) -shared $^ -o $@

.PHONY: clean
clean:
	$(RM) _*.so *.pyc *_wrap.cxx *.o

