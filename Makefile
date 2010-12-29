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
BINSUFFIX := .exe
ifeq ($(shell uname), Linux)
	BINSUFFIX :=
endif
TARGETS:=$(OBJDIR)/scatr$(BINSUFFIX) # $(OBJDIR)/old_scatr$(BINSUFFIX)

# heh.
# A nifty macro to make testing gcc features easier
# ie: $(call check_cxx,option-to-check-for,optional-fallback)
check_cxx=$(shell \
        if [ "$(1)" != "" ]; then \
                if $(CXX) $(1) -S -o /dev/null -xc /dev/null > /dev/null 2>&1; \
                then echo -n "$(1)"; else echo -n "$(2)"; fi \
        fi)

# override to set the C++ compiler.
# CXX :=g++-4.5
CXX :=g++

# preprocessor flags, generate dependencies on the fly.
CPPFLAGS +=-Wp,-MT,$@,-MMD,$@.d
CPPFLAGS +=-DNDEBUG

CXXFLAGS :=-pipe -g
CXXFLAGS +=-Wall -Wextra -Wshadow -Wno-unused -Wno-format
CXXFLAGS +=-O3 -march=native -ffast-math -fopenmp -mfpmath=sse -fomit-frame-pointer -fno-rtti -fno-exceptions
#CXXFLAGS +=-ftree-vectorizer-verbose=1
# if LTO is supported, check once.
HAS_LTO:=$(call check_cxx,-flto)
ifneq ($(strip $(HAS_LTO)),)
	# compile & link with LTO.
	#CXXFLAGS +=-flto
	#LDFLAGS +=-flto -fwhole-program
endif


## targets
.PHONY: all clean
all: $(TARGETS)
	@echo "done."

clean:
	-@rm *.d $(OBJS) $(TARGETS) 2>/dev/null

dumpflags:
	@echo "CXXFLAGS $(CXXFLAGS)"
	@echo "LTO? $(HAS_LTO)"

disas: $(OBJDIR)/scatr$(BINSUFFIX)
	objdump -dC $(OBJDIR)/scatr$(BINSUFFIX)|less

bench: $(OBJDIR)/scatr$(BINSUFFIX)
	@echo "running..."
	./$(OBJDIR)/scatr$(BINSUFFIX)

# to specifically check stuff.
# $(OBJDIR)/getacch_sse.o: CXXFLAGS +=-ftree-vectorizer-verbose=1

## rules
# just so we can avoid being spammed
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "[CXX]  $<"
	@$(COMPILE.cc) $< -o $@

$(OBJDIR)/scatr$(BINSUFFIX): $(OBJS)
	@echo "[LINK] $@"
	@$(LINK.cc) $(OBJS) main.cpp -o $@

$(OBJDIR)/old_scatr$(BINSUFFIX): CXXFLAGS :=-pipe -g -fopenmp -O3 -march=native
$(OBJDIR)/old_scatr$(BINSUFFIX): $(OBJS)
	@echo "[LINK] $@"
	@$(LINK.cc) $(OBJS) main.cpp -o $@

# include dependencies, if present.
-include $(OBJDIR)/*.d
