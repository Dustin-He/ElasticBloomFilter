PP = g++
FLAG = -std=c++11 -fopenmp -O3 -mavx2

test:exp.cpp bloom.h param.h
	$(PP) $(FLAG) -o test exp.cpp 

.PHONY:clean
clean:
	rm test
