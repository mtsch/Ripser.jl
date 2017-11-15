CXXFLAGS := -fpic -shared -L.

define cxx_compile
  $(CXX) $^ $(CXXFLAGS) -o $@
endef

all: libripser.so ripser.out

%.so: %.cpp
	$(call cxx_compile)

ripser.out:
	$(CXX) ripser/ripser.cpp -o ripser.out
