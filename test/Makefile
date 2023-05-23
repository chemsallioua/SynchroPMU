CC = gcc
CFLAGS = -pg -O3 -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -I../libs/iniparser -I../src
LIBS = ../libs/iniparser/dictionary.c ../libs/iniparser/iniparser.c ../src/pmu_estimator.c
SRC = test.c
SRC2 = test2.c
TARGET = test
TARGET2 = test2

# Check if the platform is Windows
ifeq ($(OS),Windows_NT)
	EXE_EXT = .exe
else
	EXE_EXT =
endif

all: $(TARGET) test2 profile

$(TARGET): $(LIBS) $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ -lm 

.PHONY: $(TARGET2)
$(TARGET2): $(LIBS) $(SRC2)
	$(CC) $(CFLAGS) -o $@ $^ -lm 

.PHONY: profile
profile: $(TARGET)$(EXE_EXT)
	./$(TARGET)$(EXE_EXT)
	gprof -b $(TARGET)$(EXE_EXT) gmon.out > profile.txt

.PHONY: clean
clean:
	rm -f $(TARGET)$(EXE_EXT) $(TARGET2)$(EXE_EXT) gmon.out *.csv *.txt 

.PHONY: benchmark
benchmark: $(TARGET2)$(EXE_EXT)
		  ./$(TARGET2)$(EXE_EXT)

