GCC = gcc
LDFLAGS += -lm
CFLAGS += -O2 -fPIC -lstdc++

wrapper = VBBinaryLensingLibrary_wrapper.so
wrapper_source = VBBinaryLensingLibrary_wrapper.cpp
file_o = VBBinaryLensingLibrary.o

OS := $(shell uname)
ifeq ($(OS), Darwin)
  LDFLAGS += -Wl,-install_name,$(wrapper)
else
  LDFLAGS += -Wl,-soname,$(basename $(wrapper))
endif
LDFLAGS += -shared

.PHONY: all
all: $(wrapper) 

$(wrapper): $(file_o) $(wrapper_source)
	$(GCC) $(LDFLAGS) -o $(wrapper) $(wrapper_source) $(CFLAGS) $(file_o)

%.o: %.cpp
	$(GCC) $(CFLAGS) -c $^ 

.PHONY: clean
clean:
	$(RM) *.so *.o

