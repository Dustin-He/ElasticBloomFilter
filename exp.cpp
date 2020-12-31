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
#define szVal 1 << 24
#define SZLOG 20
#define DataREP 32
#define TraceSZVal (1 << SZLOG) * DataREP
#define QueryNUM 4

#define load_rate_one_layer 1
#define fp_while_expanding 2
#define false_positive 3
#define speed_test 4
#define CBF_test 5
#define mem_copy_test 6

using namespace std;

struct FIVE_TUPLE{char key[14];};
uint32_t c[szVal];

extern uint32_t total_slot;
uint32_t experimentNo = 0;

FIVE_TUPLE traces[szVal];
uint32_t trace_size = 0;

void ReadInTraces(const char *trace_prefix) {
	FILE *funique = fopen(trace_prefix, "rb");

	if (!funique) cout << "can not open the data file" << endl;

	FIVE_TUPLE tmp_five_tuple;
	//traces.clear();
	char *tempstr = new char[13];
	printf("reading in file\n");

	int i = 0;

	while(fread(tempstr, 1, 13, funique) == 13) {
		if (i == szVal) break;
		for (int i = 0; i < 13; i++)
			//tmp_five_tuple.key[i] = tempstr[i];
			traces[trace_size].key[i] = tempstr[i];
		trace_size++;
		i++;
		//traces.push_back(tmp_five_tuple);
	}
	fclose(funique);

	printf("Successfully read in %s, %d packets\n\n", trace_prefix, trace_size);
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

Ebloom_filter bf = Ebloom_filter(INI_BLOOM_SIZE, HASH_NUM);
int main(int argc, char **argv) {
	expandOrNot = true;
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

//	/* To read the string data */
////
////	/* This is used for insert the string data. */
////	/* You also need to modify the type of the parameters of the isnert and query function */
////	
////	for (i = 0; i < traces[0].size(); ++i) {
////		if (!bf.insert(traces[0][i].key)) {
////			break;
////		}
////	}
////	for (i = 0; i < traces[0].size(); ++i) {
////		if (!bf.query(traces[0][i].key)) {
////			break;
////		}
////	}
//
	/* To test the Load rate of the one-layer version. */
	if (experimentNo == load_rate_one_layer) {
		expandOrNot = false;
		for (i = 0; i < sz; ++i) {
			c[i] = i;
		}
		for(int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);
		cout << sz << endl;
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
//
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
				cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)bf.size << endl;
				break;
			}
			uint32_t cnt = 0;
			if ((double)bf.get_1_num() / (double)bf.get_size() >= EXPAND_THRESHOLD) {
				for (int j = sz + 1; j < sz * 10; ++j) {
					if (bf.query(j))
						cnt++;
				}
				cout << "False Positive: " << (double)cnt / (sz * 10 - sz -1) << endl;
				cout << "Expand times: " << -bf.compression << endl;
				bf.expand();
			}
		}
	}
	
	/* To test the false positive rate. */
	if (experimentNo == false_positive) {
		sz /= 16;
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
		cout << "Expand times: " << -bf.compression << endl;
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
		for (i = 0; i < sz; ++i)
		{
			c[i] = i;
		}
		for (int i = 0; i < sz; i++)
			swap(c[i], c[rand() % sz]);
		double start_time, end_time;
		double resns = 0;
		int falseCnt = 0;
		//sz = TraceSZVal;
		sz = 1 << SZLOG;
		bool queryRet = 0;
		timespec dtime1, dtime2;
		long long dresns;
		double dth, th;
		double cnttrue;
		int oldsz, datarep;
// ==============================EBF=============================================================================
		// 		// ------------------------------------insertion---------------------------
		// 		for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			// bf.print();
		// 			bf.sumtime = 0;
		// 			resns = 0;
		// 			datarep = (1 << 20) / sz;
		// 			for (int k = 1; k <= datarep; k++)
		// 			{
		// 				start_time = omp_get_wtime();
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{

		// // #pragma omp task firstprivate(i, k)
		// 					// bf.insert(c[i]);
		// 					bf.insert_with_duplicates(traces[i + (k & 31) * sz].key);
		// 					// bf.insert_with_duplicates(c[i]);
		// 				}
		// 				end_time = omp_get_wtime();
		// 				resns += (end_time - start_time) * 1e9;
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{

		// 					// #pragma omp task firstprivate(i, k)
		// 					// bf.deleteEle(traces[i + (k%DataREP) * sz].key);
		// 				}
		// 				bf.clear();
		// 				// cout<<"1num "<<bf._1_num<<endl;
		// 				while (bf.sizelog > INI_BLOOM_SIZE)
		// 					bf.compress();
		// 				// bf.clear();
		// 				// cout<<"1num "<<bf._1_num<<endl;
		// 			}

		// 			th = (double)1000.0 * sz * datarep / resns;
		// 			cout << "data size: " << sz << endl;
		// 			cout << "duplicate ratio: " << falseCnt / (double)sz << endl;
		// 			cout << "insertionTime: " << resns << endl;
		// 			cout << "insertionThroughput: \t\t\t" << th << endl;
		// 			cout << "expand time percent: " << bf.sumtime / (resns / 1e9) << endl;
		// 			cout << "1-bit-num: " << bf.get_1_num() << "\tEBF size:" << bf.get_size() << endl;
		// 		}
		// 		// ------------------------------------query-- -------------------------------------

		// 		for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			int cnttrue3=0;
		// 			cnttrue = 0;
		// 			resns = 0;
		// 			datarep = (1 << 24) / sz;
		// 			dresns=0;
		// 			for (int k = 1; k <= datarep; k++)
		// 			{
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{

		// // #pragma omp task firstprivate(i, k)
		// 					bf.insert_with_duplicates(traces[i + (k & 31) * sz].key);
		// 				}

		// 				for (int i = 0; i < sz; ++i)
		// 				{
		// 					cnttrue3 += bf.query(traces[i + (k & 31) * sz].key + 1);
		// 				}

		// 				start_time = omp_get_wtime();
		// 				clock_gettime(CLOCK_MONOTONIC, &dtime1);
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{
		// // #pragma omp task firstprivate(i, k)
		// 					cnttrue += bf.query(traces[i + (k & 31) * sz].key + 1);
		// 				}
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{
		// // #pragma omp task firstprivate(i, k)
		// 					cnttrue += bf.query(traces[i + (k & 31) * sz].key);
		// 				}

		// 				clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 				end_time = omp_get_wtime();
		// 				dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 				resns += (end_time - start_time) * 1e9;
		// 				bf.clear();
		// 				while (bf.sizelog > INI_BLOOM_SIZE)
		// 					bf.compress();
		// 			}
		// 			dth = (double)1000.0 * sz * datarep * 2 / dresns;
		// 			th = (double)1000.0 * sz * datarep * 2 / resns;
		// 			cout << "QueryTime: " << resns << endl;
		// 			cout << "QueryThroughput: \t\t\t" << th << endl;
		// 			// cout << "QueryThroughput: \t\t\t" << dth << endl;
		// 			cout << "size: \t" << bf.size << endl;
		// 			// cout << "True Rate: "<<cnttrue/(sz*datarep*2)<<endl;
		// 			printf("True Rate: %.8lf\n", cnttrue / (double)(sz * datarep * 2));
		// 			cout<<cnttrue3<<"_";
		// 		}
		// 		//-----------------------------deletion --------------------------------------
		// 		for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// 		{

		// 			resns = 0;
		// 			datarep = (1 << 22) / sz;
		// 			for (int k = 1; k <= datarep; k++)
		// 			{
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{

		// // #pragma omp task firstprivate(i, k)
		// 					bf.insert_with_duplicates(traces[i + (k & 31) * sz].key);
		// 				}

		// 				start_time = omp_get_wtime();
		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{

		// // #pragma omp task firstprivate(i, k)
		// 					bf.deleteEle(traces[i + (k & 31) * sz].key);
		// 				}
		// 				end_time = omp_get_wtime();
		// 				resns += (end_time - start_time) * 1e9;
		// 			}

		// 			th = (double)1000.0 * sz * datarep / resns;
		// 			cout << "DeletionTime: " << resns << endl;
		// 			cout << "DeletionThroughput: \t\t\t" << th << endl;
		// 		}
		// ------------------------------------insertion mix query-- -------------------------
		
		for (sz = (1 << 16); sz <= (1 << 20); sz <<= 4)
		{
			bf.sumtime=0;
			resns = 0;
			datarep = (1 << 24) / sz;
			for (int k = 1; k <= datarep; k++)
			{
				start_time = omp_get_wtime();
				#pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
				for (int i = 0; i < sz; ++i)
				{
					#pragma omp task firstprivate(i, k)
					{
						bf.insert_with_duplicates(traces[i + (k & 31) * sz].key);

						if (i & 1)
							for (int j = 0; j < QueryNUM; j++)
								cnttrue += bf.query(traces[c[(i + i * i * j) & 0xFFFFF] % i + (k & 31) * sz].key);
						if ((i & 1) == 0)
							for (int j = 0; j < QueryNUM; j++)
								cnttrue += bf.query(traces[i + ((k + j) & 31) * sz].key + 1 + j);
					}
				}
				end_time = omp_get_wtime();
				resns += (end_time - start_time) * 1e9;

				// #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
				for (int i = 0; i < sz; ++i)
				{
					// #pragma omp task firstprivate(i, k)
					bf.deleteEle(traces[i + (k & 31) * sz].key);
				}
				bf.clear();
				while (bf.sizelog > INI_BLOOM_SIZE)
					bf.compress();
			}

			th = (double)1000.0 * sz * datarep * (QueryNUM) / resns;
			cout << "data size: " << sz << endl;
			// cout << "true rate: " << cnttrue / (double)sz / (QueryNUM) << endl;
			cout << "insertionQueryTime: " << resns << endl;
			// cout << "insertionQueryThroughput: \t\t\t" << th << endl;
			// cout << "QueryNUM: " << QueryNUM << endl;
			cout << "expand time percent: " << bf.sumtime / (resns / 1e9) << endl;
			cout << th << endl;
		}
		//------------------------------------insertion mix query---------------------------

		// 		for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			resns = 0;
		// 			datarep = (1 << 20) / sz;
		// 			for (int k = 1; k <= datarep; k++)
		// 			{
		// 				start_time = omp_get_wtime();
		// #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{
		// #pragma omp task firstprivate(i, k)
		// 					{
		// 						bf.insert_with_duplicates(traces[i + (k%DataREP) * sz].key);

		// 						if (i & 1)
		// 							for (int j = 0; j < QueryNUM; j++)
		// 								cnttrue += bf.query(traces[c[((uint32_t)(i + i * i * j)) % sz] % i + (k%DataREP) * sz].key);
		// 						if ((i & 1) == 0)
		// 							for (int j = 0; j < QueryNUM; j++)
		// 								cnttrue += bf.query(traces[i + ((k + j) & 31) * sz].key + 1 + j);
		// 					}
		// 				}
		// 				end_time = omp_get_wtime();
		// 				resns += (end_time - start_time) * 1e9;

		// #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 				for (int i = 0; i < sz; ++i)
		// 				{
		// #pragma omp task firstprivate(i, k)
		// 					bf.deleteEle(traces[i + (k%DataREP) * sz].key);
		// 				}
		// 				bf.clear();
		// 				while (bf.sizelog > INI_BLOOM_SIZE)
		// 					bf.compress();
		// 			}

		// 			th = (double)1000.0 * sz * datarep * (QueryNUM) / resns;
		// 			// cout << "data size: " << sz << endl;
		// 			// cout << "true rate: " << cnttrue / (double)sz / (QueryNUM) << endl;
		// 			// cout << "EBF insertionQueryTime: " << resns << endl;
		// 			// cout << "EBF insertionQueryThroughput: \t\t\t" << th << endl;
		// 			// cout << "QueryNUM: " << QueryNUM << endl;
		// 			cout << th << endl;
		// 		}
		// 		sz = oldsz;

		//
		
		
		
		// ===============================================================parBF=================================================================================================================

for (sz = (1 << 16); sz <= (1 << 20); sz <<= 4)
		{
			// datarep = (1ull << 22) * (1ull << 18) / sz / (sz);
			datarep = 1;
			if (sz == (1 << 14))
				datarep <<= 4;
			if (sz == (1 << 16))
				datarep <<= 2;
			resns = 0;
			for (int k = 1; k <= datarep; k++)
			{
				DynamicBF parBF = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD, HASH_NUM) * (1 << INI_BLOOM_SIZE) / (HASH_NUM * sz)), (1. / (HASH_NUM - 1))));
				start_time = omp_get_wtime();
				#pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
				for (int i = 0; i < sz; ++i)
				{
					#pragma omp task firstprivate(i, k)
					{
						parBF.insert(traces[i + (k & 31) * sz].key);
						if (i & 1)
							for (int j = 0; j < QueryNUM; j++)
								cnttrue += parBF.query(traces[c[(i + i * i * j) & 0xFFFFF] % i + (k & 31) * sz].key);
						if ((i & 1) == 0)
							for (int j = 0; j < QueryNUM; j++)
								cnttrue += parBF.query(traces[i + ((k + j) & 31) * sz].key + 1 + j);
					}
				}
				end_time = omp_get_wtime();
				resns += (end_time - start_time) * 1e9;

				for (int i = 0; i < sz; ++i)
				{
					parBF.deleteEle(traces[i + (k & 31) * sz].key);
				}
			}
			th = (double)1000.0 * sz * datarep * (QueryNUM) / resns;
			// cout << "insertionQueryTime: " << resns << endl;
			// cout << "parBF insertionQueryThroughput: \t\t\t" << th << endl;
			cout << th << endl;
		}
		
		
		// ===============================================================DBF=================================================================================================================
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = (1 << 22) / sz;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD, HASH_NUM) * (1 << INI_BLOOM_SIZE) / (HASH_NUM * sz)), (1. / (HASH_NUM - 1))));
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			// if(i<10) cout<<(int)*(traces[i + (k & 31) * sz].key)<<"_";
		// 			dbf.insert(traces[i + (k & 31) * sz].key);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 		dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 		// if (k == 1)
		// 		// 	cout << "dbf w: " << dbf.w<<" 1num "<< dbf.b[0]->num << endl;
		// 	}
		// 	dth = (double)1000.0 * sz * datarep / dresns;
		// 	cout << "dbf inTime: " << dresns << endl;
		// 	cout << "dbf inThroughput: \t\t\t" << dth << endl;
		// }
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = (1ull << 22) * (1ull << 18) / sz / (sz);
		// 	double cnttrue2 = 0;
		// 	cout << datarep << endl;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD, HASH_NUM) * (1 << INI_BLOOM_SIZE) / (HASH_NUM * sz)), (1. / (HASH_NUM - 1))));
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			dbf.insert(traces[i + (k & 31) * sz].key);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			cnttrue += dbf.query(traces[i + (k & 31) * sz].key);
		// 		}

		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			cnttrue2 += dbf.query(traces[i + (k & 31) * sz].key + 1);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 		dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 	}
		// 	dth = (double)1000.0 * sz * 2 * datarep / dresns;
		// 	cout << "dbf queryTime: " << dresns << endl;
		// 	cout << "dbf queryThroughput: \t\t\t" << dth << endl;
		// 	cout << "FPR: \t\t\t" << cnttrue2 / sz / datarep << endl;
		// }
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = (1ull << 22) * (1ull << 18) / sz / (sz);
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD, HASH_NUM) * (1 << INI_BLOOM_SIZE) / (HASH_NUM * sz)), (1. / (HASH_NUM - 1))));
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			dbf.insert(traces[i + (k & 31) * sz].key);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			dbf.deleteEle(traces[i + (k & 31) * sz].key);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 		dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 	}
		// 	dth = (double)1000.0 * sz * datarep / dresns;
		// 	cout << "dbf deleteTime: " << dresns << endl;
		// 	cout << "dbf deleteThroughput: \t\t\t" << dth << endl;
		// }

		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	// datarep = (1ull << 22) * (1ull << 18) / sz / (sz);
		// 	datarep = 1;
		// 	if (sz == (1 << 14))
		// 		datarep <<= 4;
		// 	if (sz == (1 << 16))
		// 		datarep <<= 2;
		// 	resns = 0;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD, HASH_NUM) * (1 << INI_BLOOM_SIZE) / (HASH_NUM * sz)), (1. / (HASH_NUM - 1))));
		// 		start_time = omp_get_wtime();
		// 		for (int i = 0; i < sz; ++i)
		// 		{
		// 			{
		// 				dbf.insert(traces[i + (k & 31) * sz].key);
		// 				if (i & 1)
		// 					for (int j = 0; j < QueryNUM; j++)
		// 						cnttrue += dbf.query(traces[c[(i + i * i * j) & 0xFFFFF] % i + (k & 31) * sz].key);
		// 				if ((i & 1) == 0)
		// 					for (int j = 0; j < QueryNUM; j++)
		// 						cnttrue += dbf.query(traces[i + ((k + j) & 31) * sz].key + 1 + j);
		// 			}
		// 		}
		// 		end_time = omp_get_wtime();
		// 		resns += (end_time - start_time) * 1e9;

		// 		for (int i = 0; i < sz; ++i)
		// 		{
		// 			dbf.deleteEle(traces[i + (k & 31) * sz].key);
		// 		}
		// 	}
		// 	th = (double)1000.0 * sz * datarep * (QueryNUM) / resns;
		// 	// cout << "insertionQueryTime: " << resns << endl;
		// 	// cout << "DBF insertionQueryThroughput: \t\t\t" << th << endl;
		// 	cout << th << endl;
		// }

		/****************************************************/

		// ===============================================================sbf=================================================================================================================
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = (1 << 22) / sz;
		// 	dresns = 0;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD / pow(log(2*sz*HASH_NUM/(INI_BLOOM_SIZE*EXPAND_THRESHOLD)/log(2)), 1. / HASH_NUM));
		// 		// if(k==1)cout<<"EXPAND_THRESHOLD / pow(1, 1. / HASH_NUM)"<<EXPAND_THRESHOLD / pow(1, 1. / HASH_NUM)<<endl;
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime1);
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			// if(i<10) cout<<(int)*(traces[i + (k & 31) * sz].key)<<"_";
		// 			sbf.insert(traces[i + (k & 31) * sz].key);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 		dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 	}
		// 	dth = (double)1000.0 * sz * datarep / dresns;
		// 	cout << "sbf inTime: " << dresns << endl;
		// 	cout << "sbf inThroughput: \t\t\t" << dth << endl;
		// }
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = (1 << 24) / sz;
		// 	dresns = 0;
		// 	double cnttrue2 = 0;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD / pow(log(2*sz*HASH_NUM/(INI_BLOOM_SIZE*EXPAND_THRESHOLD)/log(2)), 1. / HASH_NUM));
		// 		// ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD / pow(1, 1. / HASH_NUM));
		// 		// cout<<"EXPAND_THRESHOLD / pow(10, 1. / HASH_NUM)"<<(EXPAND_THRESHOLD / pow(10, 1. / HASH_NUM))<<endl;
		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			sbf.insert(traces[i + (k & 31) * sz].key);
		// 		}

		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			cnttrue2 += sbf.query(traces[i + (k & 31) * sz].key + 2);
		// 		}
		// 		clock_gettime(CLOCK_MONOTONIC, &dtime1);

		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			cnttrue2 += sbf.query(traces[i + (k & 31) * sz].key + 1);
		// 		}

		// 		for (i = 0; i < sz; ++i)
		// 		{
		// 			cnttrue += sbf.query(traces[i + (k & 31) * sz].key);
		// 		}

		// 		clock_gettime(CLOCK_MONOTONIC, &dtime2);
		// 		dresns += (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
		// 		if (k == 1)
		// 			cout << "size: \t" << (1 << (sbf.m + sbf.w)) << endl;
		// 	}
		// 	dth = (double)1000.0 * sz * datarep * 2 / dresns;
		// 	cout << "sbf queryTime: " << dresns << endl;
		// 	cout << "sbf queryThroughput: \t\t\t" << dth << endl;
		// 	cout << "sbf FPR: \t\t\t" << cnttrue2 / sz / datarep << endl;
		// 	cout << cnttrue<<endl;
		// }
		// for (sz = (1 << 14); sz <= (1 << 20); sz <<= 2)
		// {
		// 	dresns = 0;
		// 	datarep = 1;
		// 	if (sz == (1 << 14))
		// 		datarep <<= 4;
		// 	if (sz == (1 << 16))
		// 		datarep <<= 2;
		// 	resns = 0;
		// 	for (int k = 1; k <= datarep; k++)
		// 	{
		// 		ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD / pow(log(2 * sz * HASH_NUM / (INI_BLOOM_SIZE * EXPAND_THRESHOLD) / log(2)), 1. / HASH_NUM));
		// 		start_time = omp_get_wtime();
		// 		for (int i = 0; i < sz; ++i)
		// 		{
		// 			{
		// 				sbf.insert(traces[i + (k & 31) * sz].key);
		// 				if (i & 1)
		// 					for (int j = 0; j < QueryNUM; j++)
		// 						cnttrue += sbf.query(traces[c[(i + i * i * j) & 0xFFFFF] % i + (k & 31) * sz].key);
		// 				if ((i & 1) == 0)
		// 					for (int j = 0; j < QueryNUM; j++)
		// 						cnttrue += sbf.query(traces[i + ((k + j) & 31) * sz].key + 1 + j);
		// 			}
		// 		}
		// 		end_time = omp_get_wtime();
		// 		resns += (end_time - start_time) * 1e9;
		// 	}

		// 	th = (double)1000.0 * sz * datarep * (QueryNUM) / resns;
		// 	// cout << "insertionQueryTime: " << resns << endl;
		// 	// cout << "sBF insertionQueryThroughput: \t\t\t" << th << endl;
		// 	cout << th << endl;
		// }

		cout << cnttrue << endl;

	}

	/* To test the CBF. */
	if (experimentNo == CBF_test) {
		CBF cbf = CBF(CBF_SIZE, HASH_NUM);
		ReadInTraces("./130000.dat");
		expandOrNot = false;
		uint32_t cnt = 0;
		for (i = 0; i < (1 << 17); ++i) {
			if (!cbf.insert(traces[i].key)) {
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
			if (!bf.insert(traces[i].key)) {
				cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)bf.size << endl;
				break;
			}
			if ((double)bf.get_1_num() / (double)bf.size >= EXPAND_THRESHOLD) {
				cout << bf.get_1_num() / (double)bf.size << endl;
				bf.expand();
				cout << "Expand times: " << -bf.compression << " " << i << endl;
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
		ReadInTraces("./130000.dat");
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
					uint32_t hv = hash[j]->run(traces[i].key, KEY_LEN2);
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
