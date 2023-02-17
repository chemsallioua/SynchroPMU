CC = gcc
CFLAGS = -I ../libs/iniparser -I ../src
LIBS = ../libs/iniparser/dictionary.c ../libs/iniparser/iniparser.c
SRC = ../src/pmu_estimator.c

SIMPLE_EXAMPLE_SRC = simple_example.c
SIMPLE_EXAMPLE1_OUT = build/simple_example_out

SIMPLE_EXAMPLE2_SRC = simple_example2.c
SIMPLE_EXAMPLE2_OUT = build/simple_example2_out

all: simple_example simple_example2

simple_example: $(SIMPLE_EXAMPLE_SRC)
	if [ ! -d "build" ]; then mkdir build; fi
	$(CC) $(CFLAGS) $(LIBS) $(SRC) $(SIMPLE_EXAMPLE_SRC) -o $(SIMPLE_EXAMPLE1_OUT)

simple_example2: $(SIMPLE_EXAMPLE2_SRC)
	if [ ! -d "build" ]; then mkdir build; fi
	$(CC) $(CFLAGS) $(LIBS) $(SRC) $(SIMPLE_EXAMPLE2_SRC) -o $(SIMPLE_EXAMPLE2_OUT)

clean:
	rm -rf build

