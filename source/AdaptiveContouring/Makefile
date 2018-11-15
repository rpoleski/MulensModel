.PHONY: all clean

GCC = gcc
LDFLAGS += -lm 
CFLAGS += -O2 -fPIC
# CFLAGS += -Wall

wrapper = AdaptiveContouring_wrapper.so
wrapper_source = AdaptiveContouring_wrapper.c
o_files = binext_adap.o adaptive_contour.o complex.o binary_pt.o integrate.o dyn_array.o erdlcaust.o

OS := $(shell uname)
ifeq ($(OS), Darwin)
  LDFLAGS += -Wl,-install_name,$(wrapper)
else
  LDFLAGS += -Wl,-soname,$(basename $(wrapper))
endif
LDFLAGS += -shared


all: $(wrapper)

%.o:	%.c
	$(GCC) $(CFLAGS) -c $^

$(wrapper): $(o_files) $(wrapper_source)
	$(GCC) $(LDFLAGS) -o $(wrapper) $(wrapper_source) $(CFLAGS) $(o_files)

.PHONY: clean
clean:
	$(RM) $(wrapper) $(o_files)

