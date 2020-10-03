#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>
#include <string>
#include "bloom.h"
#include "CBF.h"
#include "DynamicBF.h"
#include "param.h"
#include "ScalableBF.h"

#define START_FILE_NO 1
#define END_FILE_NO 1
#define szVal 1 << 18

#define load_rate_one_layer 1
#define fp_while_expanding 2
#define false_positive 3
#define speed_test 4
#define CBF_test 5
#define mem_copy_test 6

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
	cout << mem_copy_test << " To test the cost of memory copy." << endl;
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
					if (experimentNo > mem_copy_test || !experimentNo) {
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
//
//	/* This is used for insert the string data. */
//	/* You also need to modify the type of the parameters of the isnert and query function */
//	
//	for (i = 0; i < traces[0].size(); ++i) {
//		if (!bf.insert(traces[0][i].key)) {
//			break;
//		}
//	}
//	for (i = 0; i < traces[0].size(); ++i) {
//		if (!bf.query(traces[0][i].key)) {
//			break;
//		}
//	}

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
		expandOrNot = false;
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

		clock_gettime(CLOCK_MONOTONIC, &time1);
		for (i = 0; i < sz; ++i) {
			bf.deleteEle(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &time2);
		resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
		th = (double)1000.0 * sz / resns;
		cout << "delTime: " << resns << endl;
		cout << "delThroughput: " << th << endl;
		/****************************************************/
		DynamicBF dbf = DynamicBF(BLOOM_SIZE, EXPAND_THRESHOLD);
		timespec dtime1, dtime2;
		long long dresns;
		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		for (i = 0; i < sz; ++i) {
			dbf.insert(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		dresns = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		double dth = (double)1000.0 * sz / dresns;
		cout << "dbf inTime: " << dresns << endl;
		cout << "dbf inThroughput: " << dth << endl;

		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		for (i = 0; i < sz; ++i) {
			dbf.query(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		dresns = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		dth = (double)1000.0 * sz / dresns;
		cout << "dbf quTime: " << dresns << endl;
		cout << "dbf quThroughput: " << dth << endl;

		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		for (i = 0; i < sz; ++i) {
			dbf.deleteEle(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		dresns = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		dth = (double)1000.0 * sz / dresns;
		cout << "dbf delTime: " << dresns << endl;
		cout << "dbf delThroughput: " << dth << endl;
		/****************************************************/
		timespec stime1, stime2;
		long long sresns;
		clock_gettime(CLOCK_MONOTONIC, &stime1);
		ScalableBF sbf = ScalableBF(BLOOM_SIZE, EXPAND_THRESHOLD);
		for (i = 0; i < sz; ++i) {
			sbf.insert(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &stime2);
		sresns = (long long)(stime2.tv_sec - stime1.tv_sec) * 1000000000LL + (stime2.tv_nsec - stime1.tv_nsec);
		double sth = (double)1000.0 * sz / sresns;
		cout << "sbf inTime: " << sresns << endl;
		cout << "sbf inThroughput: " << sth << endl;

		clock_gettime(CLOCK_MONOTONIC, &stime1);
		for (i = 0; i < sz; ++i) {
			sbf.query(c[i]);
		}
		clock_gettime(CLOCK_MONOTONIC, &stime2);
		sresns = (long long)(stime2.tv_sec - stime1.tv_sec) * 1000000000LL + (stime2.tv_nsec - stime1.tv_nsec);
		sth = (double)1000.0 * sz / sresns;
		cout << "sbf quTime: " << sresns << endl;
		cout << "sbf quThroughput: " << sth << endl;
	}

	/* To test the CBF. */
	if (experimentNo == CBF_test) {
		CBF cbf = CBF(CBF_SIZE, HASH_NUM);
		ReadInTraces("./");
		expandOrNot = false;
		uint32_t cnt = 0;
		for (i = 0; i < (1 << 17); ++i) {
			if (!cbf.insert(traces[0][i].key)) {
				cout << "failed at i: " << i << endl;
				break;
			}
			if ((double)cbf._1_num / (double)cbf.get_size() >= EXPAND_THRESHOLD) {
				cout << (double)cbf._1_num / (double)cbf.get_size() << endl;
				cbf.expand();
				cnt++;
				cout << "CBF Expand times: " << cnt << " " << i << endl;
			}
		}
		//printf("percentage of 1: %lf\n", cbf._1_num / (double)(1 << CBF_SIZE));
		cout << endl;		
		for (i = 0; i < (1 << 17); ++i) {
			if (!bf.insert(traces[0][i].key)) {
				cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
				break;
			}
			if ((double)bf.get_1_num() / (double)bf.get_size() >= EXPAND_THRESHOLD) {
				cout << bf.get_1_num() / (double)bf.get_size() << endl;
				bf.expand();
				cout << "Expand times: " << -bf.get_compression() << " " << i << endl;
			}
		}
	}

	if (experimentNo == mem_copy_test) {
		#define vol 30000
		uint32_t sumtime = 0;
		int times;
		int tmp = 0;
		uint32_t sz = 1 << 17;
		uint16_t array1[sz][9];
		uint16_t array2[sz][9];
		uint32_t count = 0;
		vector<BOBHash32 *> hash;
		for (int i = 0; i < 3; ++i)
			hash.push_back(new BOBHash32(i + 750));
		ReadInTraces("./");
	  	//build
		times = 20;
		sumtime = 0;
		while (times--) {
			memset(array1, 0, sizeof(array1));
			memset(array2, 0, sizeof(array2));
			count = 0;
			timespec time1, time2;
			long long resns;
			clock_gettime(CLOCK_MONOTONIC, &time1);
			for (int i = 0; i < vol; ++i) {
				for (int j = 0; j < 3; ++j) {
					uint32_t hv = hash[j]->run(traces[0][i].key, KEY_LEN2);
					uint32_t pos = hv % sz;
					if (array1[pos][0] >= 8) {
						cout << "wrong: " << (double)count / (1 << 17) << endl;
						return 0;
					}
					else if (!array1[pos][0]) {
						count++;
					}
					array1[pos][++array1[pos][0]] = hv / sz;
				}
			}
			
			clock_gettime(CLOCK_MONOTONIC, &time2);
			resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
			sumtime += resns;
		}
		cout << "1_rate: " << (double)count / (1 << 17) << endl;
		cout << "Build time with an array of size " << sz * 8 << ": " << sumtime / 20 << endl;
		//lazy update
		sumtime = 0;
		times = 20;
		while (times--) {
			timespec time1, time2;
			long long resns;
			clock_gettime(CLOCK_MONOTONIC, &time1);

			memcpy(array2, array1, sizeof(array1));
			
			clock_gettime(CLOCK_MONOTONIC, &time2);
			resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
			sumtime += resns;
		}
		cout << "Lazy update time with an array of size " << sz * 8 << ": " << sumtime / 20 << endl;
		//normal update
		sumtime = 0;
		times = 20;
		while (times--) {
			memset(array2, 0, sizeof(array2));
			timespec time1, time2;
			long long resns;
			clock_gettime(CLOCK_MONOTONIC, &time1);
			for (int i = 0; i < sz; ++i) {
				uint32_t last_pos = 0;
				uint32_t cnt = 0;
				uint16_t ending = array1[i][0];
				for (int j = 1; j <= ending; ++j) {
					if (array1[i][j] % 2 == 1) {
						array2[i][++array2[i][0]] = array1[i][j];
						last_pos = j;
						cnt++;
					}
					else if (last_pos){
						array1[i][last_pos] = array1[i][j];
						last_pos = j;
					}
				}
				array1[i][0] -= cnt;
			}
			
			clock_gettime(CLOCK_MONOTONIC, &time2);
			resns = (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
			sumtime += resns;
		}
		cout << "Expand time with an array of size " << sz * 8 << ": " << sumtime / 20 << endl;
	}

	return 0;
}
