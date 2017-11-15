CXXFLAGS := -fpic -shared -L. -Ofast -D USE_COEFFICIENTS -D NDEBUG

define cxx_compile
  $(CXX) $^ $(CXXFLAGS) -o $@
endef

all: libripser.so ripser.out
clean:
	rm libripser.so ripser.out

%.so: %.cpp
	$(call cxx_compile)

ripser.out:
	$(CXX) ripser/ripser.cpp -o ripser.out -Ofast -D NDEBUG -D USE_COEFFICIENTS
