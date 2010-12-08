TOPDIR :=.
SRCDIR :=$(TOPDIR)
OBJDIR :=g++

SRCS +=$(SRCDIR)/main.cpp
SRCS +=$(SRCDIR)/getacch.cpp
SRCS +=$(SRCDIR)/coord.cpp
SRCS +=$(SRCDIR)/io_input.cpp
SRCS +=$(SRCDIR)/structures.h

OBJS := $(SRCS)

# heh.
CXXFLAGS :=-pipe -g -fopenmp -O3 -march=native -ffast-math -mfpmath=sse -fomit-frame-pointer -fno-rtti -fno-exceptions
CXXFLAGS +=-DNDEBUG
#CXXFLAGS +=-fwhole-program
CXXFLAGS +=-Wall -Wextra -Wshadow -Wno-unused

CXXOLDFLAGS :=-pipe -g -fopenmp -O3 -march=native

## targets
scatr.exe: $(OBJS)
	g++ $(OBJS) -o scatr.exe $(CXXFLAGS)
	
old_scatr.exe: $(OBJS)
	g++ $(OBJS) -o old_scatr.exe $(CXXOLDFLAGS)
