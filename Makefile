# Compiler
CXX = c++

# Compiler flags
CXXFLAGS =

# Source files
SRCS = bitknit.cpp kraken.cpp lzna.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = ooz

# .cpp source suffix
.SUFFIXES: .cpp .o

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile source files to object files
.cpp.o:
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
