CC=gcc

CFLAGS= -O3 -Wall -Werror -o

all: lab5 dynamic

lab5: lab5.c
	$(CC) $(CFLAGS) lab5 lab5.c -lm -pthread -g
	
dynamic: lab5_d.c
	$(CC) $(CFLAGS) lab5d lab5_d.c -lm -pthread -g

clean:
	rm -f lab5 
