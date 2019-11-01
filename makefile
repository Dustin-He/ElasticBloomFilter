PP = g++
FLAG = -std=c++11

test:exp.cpp bloom.h param.h
	$(PP) $(FLAG) -O1 -o test exp.cpp

##mda:MakeData.cpp param.h
##	$(PP) $(FLAG) -o mda MakeData.cpp

.PHONY:clean
clean:
	rm test