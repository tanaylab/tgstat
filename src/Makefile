CC=g++

empty:=
space:= $(empty) $(empty)
GCCVERSIONSTRING := $(shell expr `$(CC) -dumpversion`)
#Create version number without "."
GCCVERSION := $(shell expr `echo $(GCCVERSIONSTRING)` | cut -f1 -d.)
GCCVERSION += $(shell expr `echo $(GCCVERSIONSTRING)` | cut -f2 -d.)
GCCVERSION += $(shell expr `echo $(GCCVERSIONSTRING)` | cut -f3 -d.)
# Make sure the version number has at least 3 decimals
GCCVERSION += 00
# Remove spaces from the version number
GCCVERSION := $(subst $(space),$(empty),$(GCCVERSION))
# Crop the version number to 3 decimals.
GCCVERSION := $(shell expr `echo $(GCCVERSION)` | cut -b1-3)
#$(warning $(GCCVERSION))

GCC430_OR_LATER := $(shell expr `echo $(GCCVERSION)` \>= 430)

CFLAGS=-O3 -DNDEBUG -fpic -I$(R_INCLUDE_DIR)

ifeq ($(GCC430_OR_LATER), 1)
	CFLAGS += -std=c++0x -Wno-deprecated
endif

TGT=tgstat.so
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
DEPS=$(wildcard *.d)

all: $(TGT)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $< -MP -MMD

$(TGT): $(OBJS)
	$(CC) -L$(R_HOME)/lib -lpthread -lrt -shared -o $@ $^ $(USR_LINK_FLAGS)

clean:
	rm -f $(OBJS) $(DEPS)

-include $(DEPS)

