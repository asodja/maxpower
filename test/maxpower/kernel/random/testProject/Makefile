USE_SLIC := 1
include $(MAXCOMPILERDIR)/lib/Makefile.include

SRC_DIR := ../../../../../src/maxpower/kernel/random/runtime

CFLAGS = -O3 -g -std=c99 -Wall -Werror $(MAXCOMPILER_INC) -I. -I$(SRC_DIR)
LDLIBS = $(MAXCOMPILER_LIBS)

MAXLIB_OBJS = $(SRC_DIR)/random_mt.o MersenneTwister.o slic_MersenneTwister.o

all: inform MersenneTwister

inform:
	@echo -e "\nBuild environment:"
	@echo -e "\tusing MAXELEROSDIR   = $(MAXELEROSDIR)"
	@echo -e "\tusing MAXCOMPILERDIR = $(MAXCOMPILERDIR)"
	@echo -e ""

MersenneTwister: $(MAXLIB_OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@
	
slic_%.o: %.max
	$(SLICCOMPILE) $< $@ 
	    
clean:
	$(RM) *.o MersenneTwister
