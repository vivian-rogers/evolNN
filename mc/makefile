final: binarydriver.cpp crossbar.o trainingdata.o
	g++ -o test binarydriver.cpp -I ../eigen/eigen/ crossbar.o trainingdata.o 
#binarydriver.o: binarydriver.cpp crossbar.o
#	g++ -c binarydriver.cpp crossbar.o -I ./eigen/eigen/
trainingdata.o: trainingdata.cpp trainingdata.h
	g++ -c trainingdata.cpp -I ../eigen/eigen/
crossbar.o: crossbar.cpp crossbar.h
	g++ -c crossbar.cpp -I ../eigen/eigen/
