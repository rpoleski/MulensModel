GCC = gcc
CFLAGS += -lm -lstdc++ -fPIC

wrapper = VBBinaryLensingLibrary_wrapper.so
wrapper_source = VBBinaryLensingLibrary_wrapper.cpp
file_o = VBBinaryLensingLibrary.o

OS := $(shell uname)
ifeq ($(OS), Darwin)
  OPTIONS = -Wl,-install_name,$(wrapper)
else
  OPTIONS = -Wl,-soname,$(basename $wrapper)
endif
OPTIONS += -shared

.PHONY: all
all: $(wrapper) 

$(wrapper): $(file_o) $(wrapper_source)
	$(GCC) $(OPTIONS) -o $(wrapper) $(wrapper_source) $(CFLAGS) $(file_o)

%.o: %.cpp
	$(GCC) $(CFLAGS) -c $^ 

.PHONY: clean
clean:
	$(RM) *.so *.o

