%.o: %.c main.h votrax_imp.h
	$(CC) -g -c -o $@ $< $(CFLAGS)

speak: main.o votrax_imp.o
	$(CC) -o speak main.o votrax_imp.o -lm -lpulse-simple

all: speak

.PHONY: clean all

clean:
	rm *.o speak
