# gtest required, see http://code.google.com/p/googletest/
# apt-get install libgtest-dev
TOPDIR :=.
SRCDIR :=$(TOPDIR)
OBJDIR :=$(TOPDIR)

# overkill, that's when you want different dirs for source, objects & binaries but you're a pig ;)
SRCS +=$(SRCDIR)/getacch.cpp
SRCS +=$(SRCDIR)/coord.cpp
SRCS +=$(SRCDIR)/io_input.cpp
SRCS +=$(SRCDIR)/getacch_sse.cpp

OBJS :=$(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

# check the OS.
OS:=$(shell uname)
BINSUFFIX := .exe
ifeq ($(OS), Linux)
	BINSUFFIX :=
endif
TARGETS:=$(OBJDIR)/scatr$(BINSUFFIX) # $(OBJDIR)/old_scatr$(BINSUFFIX)
GTEST_TARGETS := $(OBJDIR)/tests$(BINSUFFIX)

# heh.
# A nifty macro to make testing gcc features easier
# ie: $(call check_cxx,option-to-check-for,optional-fallback)
check_cxx=$(shell \
        if [ "$(1)" != "" ]; then \
                if $(CXX) $(1) -S -o /dev/null -xc /dev/null > /dev/null 2>&1; \
                then echo -n "$(1)"; else echo -n "$(2)"; fi \
        fi)

# override to set the C++ compiler, defaults to cxx (which is generally linked to g++).
# CXX :=g++-4.5
CXX := g++

# preprocessor flags, generate dependencies on the fly.
CPPFLAGS +=-Wp,-MT,$@,-MMD,$@.d
CPPFLAGS +=-DNDEBUG
CPPFLAGS +=-I.

CXXFLAGS :=-pipe -g
CXXFLAGS +=-Wall -Wextra -Wshadow -Wno-unused -Wno-format
CXXFLAGS +=-O3 -march=native -ffast-math -fopenmp -mfpmath=sse -fomit-frame-pointer -fno-rtti -fno-exceptions

# if LTO is supported, check once.
HAS_LTO:= # kludge for cygwin.
ifeq ($(OS), Linux) # meh
	HAS_LTO:=$(call check_cxx,-flto)
	ifneq ($(strip $(HAS_LTO)),)
		# compile & link with LTO.
		#CXXFLAGS +=-flto
		#LDFLAGS +=-flto -fwhole-program
	endif

	# GTest config
	# no pkg-config support. yay.
	GTEST_CPPFLAGS :=$(shell gtest-config --cppflags)
	GTEST_CXXFLAGS :=$(shell gtest-config --cxxflags)
	GTEST_LDFLAGS  :=$(shell gtest-config --ldflags --libs)
	# no way to query. just great.
	GTEST_LDFLAGS  +=$(shell gtest-config --libdir)/libgtest_main.a
else
	# rig it however you see fit.
	GTEST_CPPFLAGS +=-I.
	GTEST_CXXFLAGS +=
	GTEST_LDFLAGS  +=$(TOPDIR)/gtest_main.a
endif


## targets
.PHONY: all clean
all: $(TARGETS)
	@echo "done."

clean:
	-@rm *.d $(OBJS) $(TARGETS) $(GTEST_TARGETS) 2>/dev/null

dumpflags:
	@echo "CXXFLAGS $(CXXFLAGS)"
	@echo "LTO? $(HAS_LTO)"
	@echo "GTEST: cpp:$(GTEST_CPPFLAGS) cxx:$(GTEST_CXXFLAGS) ld:$(GTEST_LDFLAGS)"

disas: $(OBJDIR)/scatr$(BINSUFFIX)
	objdump -dC $(OBJDIR)/scatr$(BINSUFFIX)|less

bench: $(OBJDIR)/scatr$(BINSUFFIX)
	@echo "running..."
	./$(OBJDIR)/scatr$(BINSUFFIX)

gtest: $(GTEST_TARGETS)
	@echo "testing..."
	./$(GTEST_TARGETS)

# to specifically check stuff.
# $(OBJDIR)/getacch_sse.o: CXXFLAGS +=-ftree-vectorizer-verbose=1

## rules
# just so we can avoid being spammed
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "[CXX]  $<"
	@$(COMPILE.cc) $< -o $@

# main binary
$(OBJDIR)/scatr$(BINSUFFIX): $(OBJS)
	@echo "[LINK] $@"
	@$(LINK.cc) $(OBJS) main.cpp -o $@

# gtest binary
$(OBJDIR)/tests$(BINSUFFIX): CPPFLAGS +=$(GTEST_CPPFLAGS)
$(OBJDIR)/tests$(BINSUFFIX): CXXFLAGS +=$(GTEST_CXXFLAGS)
$(OBJDIR)/tests$(BINSUFFIX): LDFLAGS  +=$(GTEST_LDFLAGS)
$(OBJDIR)/tests$(BINSUFFIX): $(OBJS)
	@echo "[LINK] $@"
	@$(LINK.cc) $(OBJS) acceleration_tests.cpp -o $@
	
$(OBJDIR)/old_scatr$(BINSUFFIX): CXXFLAGS :=-pipe -g -fopenmp -O3 -march=native
$(OBJDIR)/old_scatr$(BINSUFFIX): $(OBJS)
	@echo "[LINK] $@"
	$(LINK.cc) $(OBJS) main.cpp -o $@

# include dependencies, if present.
-include $(OBJDIR)/*.d
