CFLAGS = -O2 -fopenmp -std=c++0x -lz

%: %.cpp %.h $(DEPS)
	$(CXX) $(CFLAGS) $< -o $@
%: %.cpp $(DEPS)
	$(CXX) $(CFLAGS) $< -o $@