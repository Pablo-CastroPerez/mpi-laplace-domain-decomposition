MF      = Makefile

CC      = cc
CFLAGS  = -O3 -Wall
LFLAGS  = $(CFLAGS) -lm

EXE     = laplace
INC     = \
        config.h \
        solvers.h \
        utils.h

SRC     = \
        main.c \
        solvers.c \
        utils.c

.SUFFIXES:
.SUFFIXES: .c .o

OBJ     = $(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:    $(EXE)

$(OBJ): $(INC)

$(EXE): $(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ): $(MF)

clean:
	rm -f $(EXE) $(OBJ) core
