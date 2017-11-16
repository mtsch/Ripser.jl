CXXFLAGS := -std=gnu++11 -Ofast -D USE_COEFFICIENTS -D NDEBUG
LIBFLAGS := -fpic -shared -L.

all: libripser.so ripser.out
clean:
	rm -f libripser.so ripser.out

%.so: %.cpp
	$(CXX) $(LIBFLAGS) $(CXXFLAGS) libripser.cpp -o libripser.so

ripser.out:
	$(CXX) $(CXXFLAGS) ripser/ripser.cpp -o ripser.out
