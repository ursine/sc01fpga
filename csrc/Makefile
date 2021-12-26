%.o: %.c main.h votrax_imp.h
	$(CC) -c -o $@ $< $(CFLAGS)

speak: main.o votrax_imp.o
	$(CC) -o speak main.o votrax_imp.o -lm

all: speak

.PHONY: clean all

clean:
	rm *.o speak
