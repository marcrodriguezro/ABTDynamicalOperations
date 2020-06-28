CXX = g++
CXXFLAGS = -g -Wall -Wextra -Wno-sign-compare -O3 -std=c++1z
SOURCES = src/abt_compression.h src/abt_compression.cc

LIB = -lgtest -lgtest_main -lpthread
GTEST_ROOT= lib/gtest-1.6.0


bin:
	mkdir -p bin
	mkdir -p bin/compress
	mkdir -p bin/test


bin/compress: src/compress_main.cc ${SOURCES}
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/tests: ../Test/test.cpp ${SOURCES} 
	$(CXX) $(CXXFLAGS) -o $@ $^ -I $(GTEST_ROOT)/include -L $(GTEST_ROOT)/lib $(LIB)
	./$@

.PHONY:	tests clean

clean:
	rm -rf bin
