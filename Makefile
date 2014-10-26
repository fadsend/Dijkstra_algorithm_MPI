PROGRAM_NAME = dijkstra
LIBS_DIR = /usr/lib/x86_64-linux-gnu/
LIBS =  
OPTIONS = -std=c++1y -Wall #-Werror

COMMAND = mpic++ $(PROGRAM_NAME).cpp -o  $(PROGRAM_NAME) -L$(LIBS_DIR) $(LIBS) $(OPTIONS)

all: compile

compile: $(PROGRAM_NAME).cpp
	$(COMMAND)

debug: $(PROGRAM_NAME).cpp
	$(COMMAND) -g

	
