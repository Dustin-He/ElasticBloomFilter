#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>
#include <string>
#include "bloom.h"
#include "CBF.h"
#include "DynamicBF.h"
#include "parBF.h"
#include "param.h"
#include "ScalableBF.h"

#define START_FILE_NO 1
#define END_FILE_NO 1
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

struct FIVE_TUPLE
{
	char key[14];
};
//typedef vector<FIVE_TUPLE> TRACE;
//TRACE traces;
FIVE_TUPLE traces[TraceSZVal];
uint32_t trace_size = 0;
uint32_t c[TraceSZVal];

extern uint32_t total_slot;
uint32_t experimentNo = 0;

void ReadInTraces(const char *trace_prefix)
{
	FILE *funique = fopen(trace_prefix, "rb");

	if (!funique)
		cout << "can not open the data file" << endl;

	FIVE_TUPLE tmp_five_tuple;
	//traces.clear();
	char *tempstr = new char[13];
	printf("reading in file\n");

	int i = 0;

	while (fread(tempstr, 1, 13, funique) == 13)
	{
		if (i == TraceSZVal)
			break;
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

inline void panic()
{
	cout << "Usage: ./test [-h] [-e Experiment Number]\n";
	exit(1);
}

inline void help()
{
	cout << "Usage: ./test [-h] [-e Experiment Number]\n";
	cout << "-h help" << endl
		 << "-e The experiment you want to run" << endl;
	cout << load_rate_one_layer << " To test the Load rate of the one-layer version." << endl;
	cout << fp_while_expanding << " To test the false positive rate while expanding." << endl;
	cout << false_positive << " To test the false positive rate." << endl;
	cout << speed_test << " To test the insertion and query speed of EBF." << endl;
	cout << CBF_test << " To test the CBF." << endl;
	cout << mem_copy_test << " To test the cost of memory copy." << endl;
}

Ebloom_filter bf = Ebloom_filter(INI_BLOOM_SIZE, HASH_NUM);

int main(int argc, char **argv)
{
	expandOrNot = true;
	int sz = 1 << SZLOG;
	int i;
	int stash = 0;

	if (argc < 2)
	{
		panic();
	}
	else
	{
		for (int i = 1; i < argc; ++i)
		{
			if (string(argv[i]) == "-h")
				help();
			else if (string(argv[i]) == "-e")
			{
				if ((i + 1) < argc)
				{
					experimentNo = atoi(argv[i + 1]);
					if (experimentNo > mem_copy_test || !experimentNo)
					{
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

	//const char *path = "./130000.dat";
	const char *path = "130000.dat";
	ReadInTraces(path);

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

	/* To test the false positive rate. */
	if (experimentNo == false_positive)
	{
		// 		sz = 1 << 20;

		// 		for (i = 0; i < sz * datarep; ++i)
		// 		{
		// 			c[i] = i;
		// 		}
		// 		for (int i = 0; i < sz * datarep; i++)
		// 			swap(c[i], c[rand() % (sz * datarep)]);

		// 		// for (sz = (1 << 12); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			for (i = sz - 1; i >= 0; --i)
		// 			{
		// 				bf.insert(c[i]);
		// 					// if (!bf.insert(c[i])) {
		// 					// 	cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
		// 					// 	break;
		// 					// }

		// 				double cnt = 0,
		// 						oldfpr = 1;
		// 				if (((sz - i) % (1 << 10) == 0))
		// 				{

		// 					int t = log(sz - i) / log(2) + 1e-10;
		// 					t = t - 2;
		// 					if ((((sz - i) >> t) << t) != sz - i)
		// 						continue;

		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 					for (int i = sz + 1; i < sz * 10; ++i)
		// 					{
		// 						if (bf.query(c[i]))
		// 						{
		// // #pragma omp atomic
		// 							cnt++;
		// 						}
		// 						if ((i & 0x3FFFF) == 0)
		// 							if (fabs((double)cnt / (i - sz) - oldfpr) < 1e-2 * oldfpr )//&&oldfpr>1e2
		// 								break;
		// 							else
		// 							{
		// 								oldfpr = (double)cnt / (i - sz);
		// 							}
		// 					}
		// 					// cout << <<"\t" << oldfpr << endl;
		// 					printf("%.5lf\t%.10lf\n", double(sz - i) / (1 << 14), oldfpr);
		// 				}
		// 			}
		// 		}

		// 		// x^5*(x^)=0.2^5   e^(-INS*5/BF)=1-x
		// 		// INS=-ln(1-x)*BF/5
		// 		// MAX/(-ln(1-x)*BF/5)*x^5=0.2^5
		// 		// x=(0.2^5*(BF)/(5*MAX))^(1/4)
		// 		// DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD);
		// 		// x^4*num=0.2^4
		// 		sz = 1 << 20;
		// 		// DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, 0.1);
		// 		DynamicBF dbf = DynamicBF(INI_BLOOM_SIZE, pow((pow(EXPAND_THRESHOLD,HASH_NUM)*(1<<INI_BLOOM_SIZE)/(HASH_NUM*sz)),(1./(HASH_NUM-1))));
		// 		// cout<<"DBF_EXPAND_THRESHOLD: "<<pow((pow(EXPAND_THRESHOLD,HASH_NUM)*(1<<INI_BLOOM_SIZE)/(HASH_NUM*sz)),(1./(HASH_NUM-1)))<<endl;

		// 		// for (sz = (1 << 12); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			for (i = sz - 1; i >= 0; --i)
		// 			{
		// 				dbf.insert(c[i]);
		// 					// if (!dbf.insert(c[i])) {
		// 					// 	cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
		// 					// 	break;
		// 					// }

		// 				double cnt = 0,
		// 						oldfpr = 1;
		// 				if (((sz - i) % (1 << 10) == 0))
		// 				{

		// 					int t = log(sz - i) / log(2) + 1e-10;
		// 					t = t - 2;
		// 					if ((((sz - i) >> t) << t) != sz - i)
		// 						continue;

		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 					for (int i = sz + 1; i < sz * 10; ++i)
		// 					{
		// 						if (dbf.query(c[i]))
		// 						{
		// // #pragma omp atomic
		// 							cnt++;
		// 						}
		// 						if ((i & 0x3FFFF) == 0)
		// 							if (fabs((double)cnt / (i - sz) - oldfpr) < 1e-2 * oldfpr )//&&oldfpr>1e2
		// 								break;
		// 							else
		// 							{
		// 								oldfpr = (double)cnt / (i - sz);
		// 							}
		// 					}
		// 					// cout << <<"\t" << oldfpr << endl;
		// 					printf("%.5lf\t%.10lf\n", double(sz - i) / (1 << 14), oldfpr);
		// 				}
		// 			}
		// 		}
		// 					cout<<"w"<<dbf.w<<endl;

		// 		// x^5*(x^)=0.2^5   e^(-INS*5/BF)=1-x
		// 		// INS=-ln(1-x)*BF/5
		// 		// ln(MAX/(-ln(1-x)*BF/5))/ln(2)*x^5=0.2^5
		// 		// 10*x^5=0.2^5
		// 		// (ln(MAX)-ln(x)-ln(BF)+ln(5))*x^5=0.2^5*ln2
		// 		// (ln(MAX)+4-ln(BF)+ln(5))*x^5=0.2^5*ln2
		// 		// x=(0.2^5*(BF)/(5*MAX))^(1/4)
		// 		// x^4*num=0.2^4
		// 		sz = 1 << 20;
		// 		// ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD);
		// 		ScalableBF sbf = ScalableBF(INI_BLOOM_SIZE, EXPAND_THRESHOLD/pow(10,1./HASH_NUM));
		// 		sz = 1 << 20;

		// 		// for (sz = (1 << 12); sz <= (1 << 20); sz <<= 2)
		// 		{
		// 			for (i = sz - 1; i >= 0; --i)
		// 			{
		// 				sbf.insert(c[i]);
		// 					// if (!sbf.insert(c[i])) {
		// 					// 	cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
		// 					// 	break;
		// 					// }

		// 				double cnt = 0,
		// 						oldfpr = 1;
		// 				if (((sz - i) % (1 << 10) == 0))
		// 				{

		// 					int t = log(sz - i) / log(2) + 1e-10;
		// 					t = t - 2;
		// 					if ((((sz - i) >> t) << t) != sz - i)
		// 						continue;

		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 					for (int i = sz + 1; i < sz * 10; ++i)
		// 					{
		// 						if (sbf.query(c[i]))
		// 						{
		// // #pragma omp atomic
		// 							cnt++;
		// 						}
		// 						if ((i & 0x3FFFF) == 0)
		// 							if (fabs((double)cnt / (i - sz) - oldfpr) < 1e-2 * oldfpr )//&&oldfpr>1e2
		// 								break;
		// 							else
		// 							{
		// 								oldfpr = (double)cnt / (i - sz);
		// 							}
		// 					}
		// 					// cout << <<"\t" << oldfpr << endl;
		// 					printf("%.5lf\t%.10lf\n", double(sz - i) / (1 << 14), oldfpr);
		// 				}
		// 			}
		// 		}
		// 					cout<<"w"<<sbf.w<<endl;

		// 		sz = 1 << 20;
		// 		{

		// 			bf.clear();
		// 			while (bf.sizelog > INI_BLOOM_SIZE)
		// 				bf.compress();
		// 			expandOrNot=0;
		// 			for (i = sz - 1; i >= 0; --i)
		// 			{
		// 				bf.insert(c[i]);
		// 					// if (!bf.insert(c[i])) {
		// 					// 	cout << "Load Rate: " << sz - i << ' ' << (double)(sz - i) * HASH_NUM / (double)total_slot << endl;
		// 					// 	break;
		// 					// }

		// 				double cnt = 0,
		// 						oldfpr = 1;
		// 				if (((sz - i) % (1 << 10) == 0))
		// 				{

		// 					int t = log(sz - i) / log(2) + 1e-10;
		// 					t = t - 2;
		// 					if ((((sz - i) >> t) << t) != sz - i)
		// 						continue;

		// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		// 					for (int i = sz + 1; i < sz * 10; ++i)
		// 					{
		// 						if (bf.query(c[i]))
		// 						{
		// // #pragma omp atomic
		// 							cnt++;
		// 						}
		// 						if ((i & 0x3FFFF) == 0)
		// 							if (fabs((double)cnt / (i - sz) - oldfpr) < 1e-2 * oldfpr )//&&oldfpr>1e2
		// 								break;
		// 							else
		// 							{
		// 								oldfpr = (double)cnt / (i - sz);
		// 							}
		// 					}
		// 					// cout << <<"\t" << oldfpr << endl;
		// 					printf("%.5lf\t%.10lf\n", double(sz - i) / (1 << 14), oldfpr);
		// 				}
		// 			}
		// 			bf.clear();
		// 			while (bf.sizelog > INI_BLOOM_SIZE)
		// 				bf.compress();
		// 		}
	}

	/* To test the insertion and query speed of EBF. */
	if (experimentNo == speed_test)
	{
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

	return 0;
}
