CC=g++
CFLAGS=-c -g -Wall -Wextra -pedantic
LDFLAGS=
SOURCES=main.cpp atom.cpp molecule.cpp ifeffit.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
		rm -rf *o main
