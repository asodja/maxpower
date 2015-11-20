NAME=superFifoTest
DESIGNNAME=SuperFifoTest_HW
CC=gcc

USE_SLIC = 1
include $(MAXCOMPILERDIR)/lib/Makefile.include

CFLAGS += -std=gnu99 -Wall -Wextra -DDESIGN_NAME=$(DESIGNNAME)
CFLAGS += $(MAXCOMPILER_INC) 
LDFLAGS += -ggdb $(MAXCOMPILER_LIBS)

$(NAME): $(NAME).o \
	$(DESIGNNAME).o 
	$(CC) $^ $(LDFLAGS) -o $@

$(DESIGNNAME).o: $(DESIGNNAME).max
	$(SLICCOMPILE) $^ $@ 

%.o : %.c $(DESIGNNAME).o
	gcc -c $(CFLAGS) $< -o $@

TARGETS := $(NAME) 

all: $(TARGETS)


clean:
	$(RM) *.o $(TARGETS)
