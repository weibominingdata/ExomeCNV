all: CNVinference.cpp HMModel.h HMModel.cpp MathTools.h MathTools.cpp
	g++ -O3 -o GENSENG CNVinference.cpp HMModel.h HMModel.cpp MathTools.h MathTools.cpp
clean:
	rm *.o
