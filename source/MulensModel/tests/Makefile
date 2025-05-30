CC = g++ # define the C/C++ compiler to use
CFLAGS = -O3 -Wall -Wextra -pedantic -DNDEBUG  
DIR = /home/data/ceres/zbudzik/VBBinaryLensing/VBBinaryLensing/lib 
# define any directories containing header files other than /usr/include
INCLUDES = -I$(DIR)
# define library paths in addition to /usr/lib
LFLAGS = -L$(DIR)
# define any libraries to link into executable:
LIBS = -lVBB
# define the C/C++ source files
SRCS = test_VBB.cpp 
# define the C/C++ object files 
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
OBJS = $(SRCS:.c=.o)
# define the executable file 
TARGET = test_VBB

.PHONY: clean

all:    $(TARGET)
	@echo  Successfully compiled.

$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
clean:
	$(RM) *.o *~ $(TARGET)
