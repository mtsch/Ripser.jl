CXXFLAGS := -fpic -shared -L.

define cxx_compile
  $(CXX) $^ $(CXXFLAGS) -o $@
endef

all: ripser-wrapper.so

%.so: %.cpp
	$(call cxx_compile)
