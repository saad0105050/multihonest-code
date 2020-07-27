
# Executables: 
#		1. mhestimates: 	calculate Pr[forkable]
#		2. mhtable: produce the table in the paper

CXXFLAGS = -g -Wall -std=c++11
LDFLAGS = -lm
CC = g++

all: mhestimates mhtable

mhestimates: prob_forkable.o ReachAndMargin.o MultiHonest.o 
	$(CC) -o $@ $^ $(LDFLAGS)

mhtable: forkability_table.o ReachAndMargin.o MultiHonest.o 
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile cpp files
prob_forkable.o: prob_forkable.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $<
forkability_table.o: forkability_table.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $<
%.o: %.cpp %.h 
	$(CC) $(CXXFLAGS) -o $@ -c $<


.PHONY: clean
clean:
	rm -rf *.o mhprob mhtable
