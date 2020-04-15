#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>
#include <string>
#include "bloom.h"
#include "CBF.h"
#include "param.h"

#define START_FILE_NO 1
#define END_FILE_NO 1
#define szVal 1 << 16

#define load_rate_one_layer 1
#define fp_while_expanding 2
#define false_positive 3
#define speed_test 4
#define CBF_test 5

using namespace std;

struct FIVE_TUPLE{char key[14];};
typedef vector<FIVE_TUPLE> TRACE;
TRACE traces[END_FILE_NO - START_FILE_NO + 1];
uint32_t c[szVal];

extern uint32_t total_slot;
uint32_t experimentNo = 0;

void ReadInTraces(const char *trace_prefix) {
	for (int datafilCnt = START_FILE_NO; datafilCnt <= END_FILE_NO; ++datafilCnt) {
		char datafileName[100];
		sprintf(datafileName, "%s%dunique.dat", trace_prefix, datafilCnt - 1);
		FILE *funique = fopen(datafileName, "rb");

		FIVE_TUPLE tmp_five_tuple;
		traces[datafilCnt - 1].clear();
		char *tempstr = new char[13];
		printf("reading in file\n");

		while(fread(tempstr, 1, 13, funique) == 13) {
			for (int i = 0; i < 13; i++)
				tmp_five_tuple.key[i] = tempstr[i];
			traces[datafilCnt - 1].push_back(tmp_five_tuple);
		}
		fclose(funique);

		printf("Successfully read in %s, %ld packets\n", datafileName, traces[datafilCnt - 1].size());
	}
	printf("\n");
}

inline void panic() {
	cout << "Usage: ./test [-h] [-e Experiment Number]\n";
	exit(1);
}

inline void help() {
	cout << "Usage: ./test [-h] [-e Experiment Number]\n";
	cout << "-h help" << endl << "-e The experiment you want to run" << endl;
	cout << load_rate_one_layer << " To test the Load rate of the one-layer version." << endl;
	cout << fp_while_expanding << " To test the false positive rate while expanding." << endl;
	cout << false_positive << " To test the false positive rate." << endl;
	cout << speed_test << " To test the insertion and query speed of EBF." << endl;
	cout << CBF_test << " To test the CBF." << endl;
}

int main(int argc, char **argv) {
	expandOrNot = true;
	Ebloom_filter bf = Ebloom_filter(BLOOM_SIZE, HASH_NUM);
	int sz = szVal;
	int i;
	int stash = 0;

	if (argc < 2) {
		panic();
	}
	else {
		for (int i = 1; i < argc; ++i) {
			if (string(argv[i]) == "-h")
				help();
			else if (string(argv[i]) == "-e") {
				if ((i + 1) < argc) {
					experimentNo = atoi(argv[i + 1]);
					if (experimentNo > CBF_test || !experimentNo) {
						cout << "Wrong experiment number.\n";
						exit(1);
					}
					i++;
				}
				else
					panic();
			}
		}
	}

	/* To read the string data */
	//ReadInTraces("./");

	/* This is used for insert the string data. */
	/* You also need to modify the type of the parameters of the isnert and query function */
	/*
	for (i = 0; i < ending; ++i) {
		if (!bf.insert(string(traces[0][i].key))) {
			break;
		}
	}
	*/

	/* To test the Load rate of the one-layer version. */
	if (experimentNo == load_rate_one_layer) {
		expandOrNot = false;
		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);

		for (int i = 0; i < sz; ++i) {
			if (!bf.insert(c[i])) {
				stash++;
				if (stash > 16) {
					cout << "Items: " << i << endl;
					cout << "LoadRate: " << (double)(i) * HASH_NUM / (double)total_slot << endl;
					break;
				}
			}
		}
		cout << "_1_rate: " << bf.get_1_num() / (double)bf.get_size() << endl;
	}

	/* To test the false positive rate while expanding. */
	if (experimentNo == fp_while_expanding) {
		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);
	
		for (i = sz - 1; i >= 0; --i) {
			if (!bf.insert(c[i])) {
				cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
				break;
			}
			uint32_t cnt = 0;
			if ((double)bf.get_1_num() / (double)bf.get_size() >= EXPAND_THRESHOLD) {
				for (int j = sz + 1; j < sz * 10; ++j) {
					if (bf.query(j))
						cnt++;
				}
				cout << "False Positive: " << (double)cnt / (sz * 10 - sz -1) << endl;
				cout << "Expand times: " << -bf.get_compression() << endl;
				bf.expand();
			}
		}
	}
	
	/* To test the false positive rate. */
	if (experimentNo == false_positive) {
		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);

	
		for (i = sz - 1; i >= 0; --i) {
			if (!bf.insert(c[i])) {
				cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
				break;
			}
		}
		uint32_t cnt = 0;
		for (i = sz + 1; i < sz * 10; ++i) {
			if (bf.query(i))
				cnt++;
		}
		cout << "False Positive: " << (double)cnt / (sz * 10 - sz -1) << endl;
		cout << "Expand times: " << -bf.get_compression() << endl;
	}

	/*for (i = 0; i < sz; ++i) {
		int tmp = bf.insert(c[i]);
		if (!tmp) {
			printf("failed at i: %d\n", i);
			break;
		}
		else if (tmp == -1) {
			printf("expand at i: %d\n", i);
		}
	}*/

	/* To test the insertion and query speed of EBF. */
	if (experimentNo == speed_test) {
		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);
		timespec time1, time2;
		long long resns;
		clock_gettime(CLOCK_MONOTONIC, &time1);
		for (i = 0; i < sz; ++i) {
			bf.insert(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &time2);
		resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
		double th = (double)1000.0 * sz / resns;
		cout << "inTime: " << resns << endl;
		cout << "inThroughput: " << th << endl;

		clock_gettime(CLOCK_MONOTONIC, &time1);
		for (i = 0; i < sz; ++i) {
			bf.query(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &time2);
		resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
		th = (double)1000.0 * sz / resns;
		cout << "quTime: " << resns << endl;
		cout << "quThroughput: " << th << endl;
	}

	/* To test the CBF. */
	if (experimentNo == CBF_test) {
		CBF cbf = CBF(CBF_SIZE, HASH_NUM);

		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);

		for (i = 0; i < sz; ++i) {
			if (cbf.insert(c[i])) {
				printf("failed at i: %d\n", i);
				break;
			}
		}
	}

	return 0;
}