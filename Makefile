# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wno-conversion-null -Wno-format -Wno-unsequenced -Wno-write-strings

# Source files
SRCS = bitknit.cpp kraken.cpp lzna.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = ooz

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Test decompression of all test files
test: $(TARGET)
	@for f in testdata/*; do \
		./$(TARGET) -b "$$f" || exit 1; \
	done

.PHONY: all clean test