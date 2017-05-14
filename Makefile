
CC := gcc
CXX := g++

C_FLAGS := -g -O0 -std=c99

ifeq ($(DEBUG), 1)
  CXX_FLAGS := -g -O0 -D DEBUG
else
  CXX_FLAGS := -O3
endif
CXX_FLAGS += -std=c++11 -fopenmp

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
MAIN_DEBUG := fexipro_debug
MAIN_RELEASE := fexipro
ifeq ($(DEBUG), 1)
  MAIN := $(MAIN_DEBUG)
else
  MAIN := $(MAIN_RELEASE)
endif

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: clean

all: $(MAIN)

$(MAIN): $(CXX_OBJS) $(C_OBJS)
	$(CXX) $(CXX_FLAGS) $(CXX_OBJS) $(C_OBJS) $(LIBS) $(LFLAGS) -o $(MAIN)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CXX_FLAGS) $(INCLUDES) -c $<  -o $@

.c.o:
	$(CC) $(C_FLAGS) $(INCLUDES) -c $<  -o $@
clean:
	$(RM) *.o *~ $(MAIN_DEBUG) $(MAIN_RELEASE)
