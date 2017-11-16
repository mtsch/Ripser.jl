CXXFLAGS := -std=gnu++11 -Ofast -D USE_COEFFICIENTS -D NDEBUG
LIBFLAGS := -fPIC -shared -L.

ifeq ($(platform), Windows)
	EXT=dll
else
	EXT=so
endif

all: libripser.so ripser.out
clean:
	rm -f libripser.$(EXT) ripser.out

%.so: %.cpp
	$(CXX) $(LIBFLAGS) $(CXXFLAGS) libripser.cpp -o libripser.$(EXT)

ripser.out:
	$(CXX) $(CXXFLAGS) ripser/ripser.cpp -o ripser.out
