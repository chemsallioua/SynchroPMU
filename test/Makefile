CC = gcc
CFLAGS = -g -I../libs/iniparser -I../src 
LIBS = ../libs/iniparser/dictionary.c ../libs/iniparser/iniparser.c ../src/pmu_estimator.c
SRC = test.c
TARGET = test

# Check if the platform is Windows
ifeq ($(OS),Windows_NT)
	EXE_EXT = .exe
else
	EXE_EXT =
endif

all: $(TARGET)

$(TARGET): $(LIBS) $(SRC)
	$(CC) $(CFLAGS) -o $@ $^

# .PHONY: profile
# profile: $(TARGET)$(EXE_EXT)
# 	./$(TARGET)$(EXE_EXT)
# 	gprof $(TARGET)$(EXE_EXT) gmon.out > profile.txt

	
profile: $(TARGET)
	perf record -g ./$(TARGET)$(EXE_EXT)
	perf report

.PHONY: clean
clean:
	rm -f $(TARGET)$(EXE_EXT) gmon.out
