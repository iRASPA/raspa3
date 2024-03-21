CXX=clang++
DBGCXXFLAGS = -g -O0 -DDEBUG -fsanitize=address
RELCXXFLAGS = -g -O3 -march=native -ffast-math 
MODULEPATHS = -fprebuilt-module-path=../mathkit -fprebuilt-module-path=../foundationkit -fprebuilt-module-path=../symmetrykit -fprebuilt-module-path=../raspakit
CXXFLAGS= -I/usr/local/include -fpic -std=c++2b -fmodules -fbuiltin-module-map -stdlib=libc++ -fopenmp -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wno-gnu-anonymous-struct -Werror -fomit-frame-pointer -ftree-vectorize -fno-stack-check -funroll-loops
PYTHONFLAGS = $(shell python3 -m pybind11 --includes)
PYTHONSITEPACKAGES = $(shell python3 -c 'import site; print(site.getsitepackages()[0])')

all: release;

release:
	$(MAKE) -C foundationkit CXX="$(CXX)" CXXFLAGS="$(RELCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C mathkit CXX="$(CXX)" CXXFLAGS="$(RELCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C symmetrykit CXX="$(CXX)" CXXFLAGS="$(RELCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C raspakit CXX="$(CXX)" CXXFLAGS="$(RELCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS) $(PYTHONFLAGS)"
	$(CXX) $(RELCXXFLAGS) $(CXXFLAGS) -fprebuilt-module-path=mathkit -fprebuilt-module-path=foundationkit -fprebuilt-module-path=raspakit -c raspa3/main/main.cpp -o raspa3/main/main.o
	$(CXX) -flto -stdlib=libc++ -o raspa3.exe raspa3/main/main.o raspakit/libraspakit.a symmetrykit/libsymmetrykit.a foundationkit/libfoundationkit.a mathkit/libmathkit.a -llapack -lblas -fopenmp
	$(CXX) -std=c++20 -shared -undefined dynamic_lookup `python3 -m pybind11 --includes` mathkit/*.o foundationkit/*.o symmetrykit/*.o raspakit/*.o -llapack -lblas -fopenmp `python3-config --ldflags` -o raspa3/raspa3.so

debug:
	$(MAKE) -C foundationkit CXX="$(CXX)" CXXFLAGS="$(DBGCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C mathkit CXX="$(CXX)" CXXFLAGS="$(DBGCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C symmetrykit CXX="$(CXX)" CXXFLAGS="$(DBGCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS)"
	$(MAKE) -C raspakit CXX="$(CXX)" CXXFLAGS="$(DBGCXXFLAGS) $(CXXFLAGS) $(MODULEPATHS) $(PYTHONFLAGS)"
	$(CXX) $(DBGCXXFLAGS) $(CXXFLAGS) -fprebuilt-module-path=mathkit -fprebuilt-module-path=foundationkit -fprebuilt-module-path=raspakit -c raspa3/main/main.cpp -o raspa3/main/main.o
	$(CXX) -flto -stdlib=libc++ -fsanitize=address -o raspa3.exe raspa3/main/main.o raspakit/libraspakit.a symmetrykit/libsymmetrykit.a foundationkit/libfoundationkit.a mathkit/libmathkit.a -llapack -lblas -fopenmp
	$(CXX) -std=c++20 -shared -undefined dynamic_lookup `python3 -m pybind11 --includes` mathkit/*.o foundationkit/*.o symmetrykit/*.o raspakit/*.o -llapack -lblas -fopenmp `python3-config --ldflags` -o raspa3/raspa3.so

clean:
	rm -f raspa3.exe
	cd foundationkit && $(MAKE) clean
	cd mathkit && $(MAKE) clean
	cd symmetrykit && $(MAKE) clean
	cd raspakit && $(MAKE) clean
	rm raspa3/*.so