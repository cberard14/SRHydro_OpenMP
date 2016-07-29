CXX=g++
CXXFLAGS=-O2 -g -std=c++11 -fopenmp 
LIBS=-lgsl -lgslcblas

all: Hydro RiemannTester ConservedToPrimitiveTester

Hydro: Hydro.o RiemannDriver.o ConservedToPrimitiveDriver.o
	${CXX} -L /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/lib -L /scinet/gpc/Libraries/boost_1_54_0-gcc4.8.1/lib -fopenmp -o Hydro Hydro.o RiemannDriver.o ConservedToPrimitiveDriver.o ${LIBS}

RiemannTester: RiemannTester.o RiemannDriver.o
	${CXX} -L /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/lib -L /scinet/gpc/Libraries/boost_1_54_0-gcc4.8.1/lib -fopenmp -o RiemannTester RiemannTester.o RiemannDriver.o ${LIBS}

RiemannTester.o: RiemannTester.cpp ./include/RiemannDriver.h
	${CXX} ${CXXFLAGS} -I /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/include  -c RiemannTester.cpp

ConservedToPrimitiveTester: ConservedToPrimitiveTester.o ./include/ConservedToPrimitiveDriver.h
	${CXX} -L /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/lib -L /scinet/gpc/Libraries/boost_1_54_0-gcc4.8.1/lib -fopenmp -o ConservedToPrimitiveTester ConservedToPrimitiveTester.o ConservedToPrimitiveDriver.o ${LIBS}

ConservedToPrimitiveTester.o: ConservedToPrimitiveTester.cpp ./include/ConservedToPrimitiveDriver.h
	${CXX} ${CXXFLAGS} -I /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/include  -c ConservedToPrimitiveTester.cpp

Hydro.o: Hydro.cpp ./include/ConservedToPrimitiveDriver.h ./include/RiemannDriver.h 
	${CXX} ${CXXFLAGS} -I /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/include  -c Hydro.cpp

ConservedToPrimitiveDriver.o: ./src/ConservedToPrimitiveDriver.cc
	${CXX} ${CXXFLAGS} -I /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/include  -c ./src/ConservedToPrimitiveDriver.cc

RiemannDriver.o: ./src/RiemannDriver.cc 
	${CXX} ${CXXFLAGS} -I /scinet/gpc/Libraries/gsl-1.16-gcc-4.8.1/include -I /scinet/gpc/Libraries/boost_1_54_0-gcc4.8.1/include -c ./src/RiemannDriver.cc

clean:
	rm Hydro RiemannTester ConservedToPrimitiveTester Hydro.o RiemannDriver.o ConservedToPrimitiveDriver.o
