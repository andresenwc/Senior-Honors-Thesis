CC = g++
CFLAGS = -g -Wall -std=c++0x
OBJS = main.o ImpactFinder.o Hyperbola.o

.C.o:
	$(CC) -c $(CFLAGS) -o $@ $<

ImpactFinder: $(OBJS)
	$(CC) $(CFLAGS) -o ImpactFinder $(OBJS)

main.o: ImpactFinder.h Point.h config.h
    
ImpactFinder.o: ImpactFinder.h Hyperbola.h Point.h config.h

Hyperbola.o: Hyperbola.h Point.h config.h

clean:
	rm *.o ImpactFinder
