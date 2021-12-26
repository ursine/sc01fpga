%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

speak: main.o
	$(CC) -o speak main.o

all: speak

.PHONY: clean all

clean:
	rm *.o speak
