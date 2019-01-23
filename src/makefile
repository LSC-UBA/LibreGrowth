CC=g++

# Standard compilation
CFLAGS_STD=-c -Wall -Wextra -std=c++11
LDFLAGS_STD=-std=c++11

# Fast compilation
CFLAGS_FAST=-c -Wall -Wextra -std=c++11 -Ofast
LDFLAGS_FAST=-std=c++11

# March compilation
CFLAGS_MARCH=-c -Wall -Wextra -std=c++11 -march=native
LDFLAGS_MARCH=-std=c++11

# Debug compilation
CFLAGS_DBG=-c -g -Wall -Wextra -std=c++11 
LDFLAGS_DBG=-std=c++11

# With optimizations
CFLAGS_OPT=-c -O3 -Wall -Wextra -std=c++11 
LDFLAGS_OPT=-std=c++11

# With OpenMP
CFLAGS_OMP=-c -fopenmp -Wall -Wextra -std=c++11
LDFLAGS_OMP= -fopenmp  -std=c++11

# With optimizations and OpenMP
CFLAGS_OPT_OMP=-c -O3 -ffast-math -ftree-vectorize -funroll-loops -fopenmp -Wall -Wextra -std=c++11 
LDFLAGS_OPT_OMP= -fopenmp  -std=c++11


CFLAGS  = ${CFLAGS_${COMP_TYPE}}
LDFLAGS = ${LDFLAGS_${COMP_TYPE}}


SOURCES=Mesh.h ScalarField.h Params.h main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=LibreGrowth

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm $(EXECUTABLE) -f
	rm *.o *.gch -f
