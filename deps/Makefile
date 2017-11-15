CXXFLAGS := -std=gnu++11 -Ofast -D USE_COEFFICIENTS -D NDEBUG
LIBFLAGS := -fpic -shared -L.

define cxx_compile
  $(CXX) $^ $(CXXFLAGS) -o $@
endef

all: libripser.so ripser.out
clean:
	rm libripser.so ripser.out

%.so: %.cpp
	$(CXX) $(CXXFLAGS) libripser.cpp -o libripser.so

ripser.out:
	$(CXX) $(LIBFLAGS) $(CXXFLAGS) ripser/ripser.cpp -o ripser.out
