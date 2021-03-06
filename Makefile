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

exec:
	./${EXEC}
	
.PHONY: clean

clean:
	rm $(OBJECTS) $(EXEC)
