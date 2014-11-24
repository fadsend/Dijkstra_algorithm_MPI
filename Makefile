PROGRAM_NAME = dijkstra 
LIBS_DIR = /usr/lib/x86_64-linux-gnu/
LIBS =  
OPTIONS = -std=c++1y -Wall #-Werror
COMMAND = mpic++ $(PROGRAM_NAME).cpp -o  $(PROGRAM_NAME) -L$(LIBS_DIR) $(LIBS) $(OPTIONS)

all: dijkstra.o heap.o
	mpic++ dijkstra.o heap.o -o dijkstra 
heap.o: heap/heap.cpp
	g++ heap/heap.cpp -c $(OPTIONS) 

dijkstra.o: dijkstra.cpp
	mpic++ dijkstra.cpp -c $(OPTIONS)
compile: $(PROGRAM_NAME).cpp
	$(COMMAND)

debug: $(PROGRAM_NAME).cpp
	$(COMMAND) -g

	
