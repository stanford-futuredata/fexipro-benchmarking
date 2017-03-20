
CC := gcc
CXX := g++

CFLAGS := -g -O3 -std=c99
CXXFLAGS := -g -O3 -std=c++11

CXX_SRCS := main.cpp
CXX_OBJS := $(CXX_SRCS:.cpp=.o)

C_SRCS := 
C_OBJS := $(C_SRCS:.c=.o)

PKG_CONFIG_INCLUDE=$(shell pkg-config armadillo --cflags)
PKG_CONFIG_LIBRARIES=$(shell pkg-config armadillo --libs)

#USE_PKG_CONFIG := 0
#PKG_CONFIG_LIBRARIES = armadillo
#ifeq ($(USE_PKG_CONFIG), 1)
#	PKG_CONFIG_INCLUDE := $(foreach library,$(PKG_CONFIG_LIBRARIES),$(shell pkg-config $(library) --cflags))
#	PKG_CONFIG_LIBS := $(foreach library,$(PKG_CONFIG_LIBRARIES),$(shell pkg-config $(library) --libs))
#else
#	PKG_CONFIG_INCLUDE := 
#	PKG_CONFIG_LIBS := 
#endif

INCLUDE_DIRS := /usr/include/boost $(HOME)/packages/include
INCLUDES += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) $(PKG_CONFIG_INCLUDE)

LIBRARY_DIRS := /usr/lib/x86_64-linux-gnu/ $(HOME)/packages/lib
LIBRARIES :=  boost_program_options armadillo
LIBS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_CONFIG_LIBS) \
                $(foreach library,$(LIBRARIES),-l$(library))

# define the executable file 
MAIN := fexipro

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: clean

all: $(MAIN)

$(MAIN): $(CXX_OBJS) $(C_OBJS)
	$(CXX) $(CXXFLAGS) $(CXX_OBJS) $(C_OBJS) $(LIBS) $(LFLAGS) -o $(MAIN)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
clean:
	$(RM) *.o *~ $(MAIN)
