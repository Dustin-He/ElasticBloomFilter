#include <iostream>
#include <ctime>
#include <vector>
#include <random>
#include "bloom.h"
//#include "cbloom.h"
#include "param.h"


#define START_FILE_NO 1
#define END_FILE_NO 1
using namespace std;

struct FIVE_TUPLE{char key[14];};
typedef vector<FIVE_TUPLE> TRACE;
TRACE traces[END_FILE_NO - START_FILE_NO + 1];
uint32_t c[1 << 19];

extern uint32_t total_slot;

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

int main() {
	/* To read the string data */
	//ReadInTraces("./");
	bloom_filter bf = bloom_filter(BLOOM_SIZE, HASH_NUM);
	int i;
	int sz = 1 << 19;
	int stash = 0;
	/* This is used for insert the string data. */
	/* You also need to modify the type of the parameters of the isnert and query function */
	/*
	for (i = 0; i < ending; ++i) {
		if (!bf.insert(string(traces[0][i].key))) {
			break;
		}
	}
	*/

	/* To test the Load rate of the one-layer version. Comment the expand operation in the insert function.*/
	/*
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
	*/

	/* To test the false positive rate while expanding */
	/*
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
	*/
	
	/* To test the false positive rate. */
	
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
	
	

	/* To test the insertion and query speed of EBF */
	/*
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
	*/
	



	
	return 0;
}