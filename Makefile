CC=gcc
EXEC=arithmetic

DEPS=$(wildcard *.h) 
SOURCE=$(wildcard *.c) 
OBJECTS=$(SOURCE:%.c=%.o)
FLAG=-Wall -O -lm -lgmp -lgmpxx -lmpfr -lmpc

$(EXEC): $(SOURCE) $(DEPS)
	$(CC) -g -o $@ $^ $(FLAG)

%.o: %.c $(DEPS)
	$(CC) -g -c -o $@ $< $(FLAG)

matrixPol: main.c test.c matrix_pol.c pol.c $(DEPS)
	$(CC) -g -o matrixPol $^ $(FLAG)


.PHONY: clean

clean:
	rm $(OBJECTS) $(EXEC)
