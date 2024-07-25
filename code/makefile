Compile=g++ -g -std=c++20 -O3 # improve run time e.g. 10s to 4s
semiflexible_polymer: main.o semiflexible_polymer.o
	$(Compile) -o semiflexible_polymer main.o semiflexible_polymer.o
main.o: main.cpp semiflexible_polymer.h
	$(Compile) -o main.o -c main.cpp
semiflexible_polymer.o: semiflexible_polymer.cpp semiflexible_polymer.h
	$(Compile) -o semiflexible_polymer.o -c semiflexible_polymer.cpp
clean:
	rm -f main.o semiflexible_polymer.o nohup.out