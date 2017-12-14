.PHONY: all clean

wrapper = AdaptiveContouring_wrapper.so
wrapper_source = AdaptiveContouring_wrapper.c

all: $(wrapper)

$(wrapper): $(file_o) $(wrapper_source)
	gcc -lm -lstdc++ -fPIC -c binext_adap.c -Wall
	gcc -lm -lstdc++ -fPIC -c adaptive_contour.c -Wall
	gcc -lm -lstdc++ -fPIC -c complex.c -Wall
	gcc -lm -lstdc++ -fPIC -c binary_pt.c -Wall
	gcc -lm -lstdc++ -fPIC -c -Wall integrate.c
	gcc -lm -lstdc++ -fPIC -c -Wall dyn_array.c
	gcc -lm -lstdc++ -fPIC -c -Wall erdlcaust.c
	gcc -Wl,-soname,rapper -shared -o AdaptiveContouring_wrapper.so AdaptiveContouring_wrapper.c -lm -lstdc++ -fPIC binext_adap.o adaptive_contour.o -Wall complex.o binary_pt.o integrate.o dyn_array.o erdlcaust.o

.PHONY: clean
clean:
	$(RM) *.so *.o

