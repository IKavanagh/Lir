# The C compiler to use
CC = icc

# User-friendly check for compiler
ifeq ($(shell which $(CC) > /dev/null 2>&1; echo $$?), 1)
	$(error The '$(CC)' compiler was not found. Make sure you have a GNU C99 compatible compiler installed, then set the CC make variable to point to the executable.)
endif

# -std=gnu99 needed to enable bessel functions
CFLAGS = -std=gnu99 \
         -O3 \
         -Iinclude \
         -Wall -Wextra \
         -pedantic \
         -Wfloat-equal -Wshadow -Wwrite-strings -Wstrict-prototypes \
         -fno-builtin \
         -Wundef -Wpointer-arith \
         -Wcast-qual -Wconversion -Wunreachable-code

# Link FFTW3 and Intel MKL
LDLIBS = -lfftw3

# TODO: Check if gcc bugs are present on Linux
ifeq ($(CC),icc)
	CFLAGS += -fopenmp
	LDLIBS += -lfftw3_omp -mkl=parallel -lm
else
	# Default to sequential when not using icc as there is a bug
	# in vefie.c which only arises when not using icc and
	# worsens when performing computations in parallel
	CFLAGS += -I${MKLROOT}/include -ftrapv -Wcast-align -Wswitch-default -Wswitch-enum
	LDLIBS += -L${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
endif

# Paths to search for C source and header file
vpath %.c src
vpath %.h include

DEPDIR = .deps
OBJDIR = objects

OUTDIR = output

# List of '.c' source files in all subdirectories recursively
SRCS := $(notdir $(wildcard **/*.c))

# Corresponding list of '.o' targets for '.c' source files
OBJS := $(addprefix $(OBJDIR)/, $(patsubst %.c,%.o, $(SRCS)))

2D: 2D.c $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(OBJS) $(LOADLIBES) $(LDLIBS) -o $@
	@echo
	@echo "Build finished. Created executable $@!"

all: docs 2D

.PHONY: clean docs

clean:
	-rm -rf 2D $(DEPDIR) $(OBJDIR) $(OUTDIR)/*

docs:
	@make -C doc html

# Order-only prerequisite, requires $(DEPDIR) to be created before any dependencies are
$(SRCS:%.c=$(DEPDIR)/%.d): | $(DEPDIR)

$(DEPDIR):
	mkdir $(DEPDIR)

# Order-only prerequisite, requires $(OBJDIR) to be created before any objects are
$(OBJS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

# Creates an object file and a dependencies file for each source file
# Object files are stored in $(OBJDIR) and dependencies in $(DEPDIR)
$(OBJDIR)/%.o: %.c
	$(CC) -Wp,-MMD,$(DEPDIR)/$*.d.part $(CPPFLAGS) $(CFLAGS) -c $< -o $@;
	@-rm -f $(DEPDIR)/$*.d;
	sed 's,\($*\)\.o[ :]*,\1.o $*.d: ,g' < $(DEPDIR)/$*.d.part > $(DEPDIR)/$*.d;
	@-rm -f $(DEPDIR)/$*.d.part
	@echo

# Only include dependency files if we need them
ifeq (,$(filter $(MAKECMDGOALS),clean docs tar))
-include $(SRCS:%.c=$(DEPDIR)/%.d)
endif
