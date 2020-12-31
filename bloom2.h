#ifndef _BLOOM_H
#define _BLOOM_H

#include "param.h"
#include "BOBHash32.h"
#include <assert.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <immintrin.h>
#include <cstdint>
#include <bits/stdc++.h>
#include "murmur3.h"

#define FINGER_LENGTH finger_length
#define POSMASK posmask
#define SIZELOG sizelog
// #define FINGER_LENGTH 8
// #define POSMASK 0x00FFFFFF
// #define SIZELOG 24

using namespace std;

void fastMemcpy(void *pvDest, void *pvSrc, size_t nBytes)
{
	assert(nBytes % 32 == 0);
	assert((intptr_t(pvDest) & 31) == 0);
	assert((intptr_t(pvSrc) & 31) == 0);
	const __m256i *pSrc = reinterpret_cast<const __m256i *>(pvSrc);
	__m256i *pDest = reinterpret_cast<__m256i *>(pvDest);
	int64_t nVects = nBytes / sizeof(*pSrc);

	for (; nVects > 0; nVects--, pSrc++, pDest++)
	{
		const __m256i loaded = _mm256_stream_load_si256(pSrc);
		_mm256_stream_si256(pDest, loaded);
	}
	_mm_sfence();
}

extern uint32_t total_slot;

class Ebloom_filter
{
public:
	Ebloom_filter(uint32_t size, int hash_num);
	~Ebloom_filter();
	void lazy_update(int layer, uint32_t bid);
	void insert(uint32_t key);
	void insert(char *key);
	void insert_with_duplicates(uint32_t key);
	void insert_with_duplicates(char *key);
	bool query(uint32_t key);
	bool query(char *key);
	void expand();
	void compress();
	void clear();
	void deleteEle(uint32_t key);
	void deleteEle(char *key);
	inline uint32_t get_0_num();
	inline int get_1_num();
	inline uint32_t get_size();
	double sumtime = 0;

	// private:
	alignas(64) uint16_t finger_buckets[1 << MAX_BLOOM_SIZE][BUCKET_SIZE + 2]; //to simulate finger buckets in slow memory

	uint32_t size; //the length of the bit map
	int sizelog;
	uint32_t posmask;
	uint32_t ExpandBitNum;

	int compression; //to record the expanding times(= -compression)

	int _1_num = 0; //the number of those set bits

	omp_lock_t lock[LOCK_NUM];
	uint32_t finger_length; //the length of the fingerprint

	char *bloom_arr;				  //the standard Bloom filter
	BOBHash32 HashFunList[HASH_NUM];  //the hash function list
	inline bool getbit(uint32_t pos); //get a bit in the bit map
	inline void setbit(uint32_t pos); //set a bit in the bit map
	inline void resetbit(uint32_t pos);
	void print();
};
void Ebloom_filter::print()
{
	cout << "myprint::::::::::::::::::\n";
	cout << "size " << size << endl;
	cout << "sizelog " << sizelog << endl;
	cout << "posmask " << posmask << endl;
	cout << "ExpandBitNum " << ExpandBitNum << endl;
	cout << "compression " << compression << endl;
	cout << "_1_num " << _1_num << endl;
	cout << "myprint------------------\n";
}

/*build the EBF*/
Ebloom_filter::Ebloom_filter(uint32_t sz, int hash_num)
{
	assert(hash_num > 0);
	assert(sz < MAX_SIZE && sz > MIN_SIZE);

	finger_length = MAX_SIZE - sz;

	/*to allocate the bit map and the hash function*/
	size = 1 << sz;
	sizelog = sz;
	posmask = (1 << sizelog) - 1;
	ExpandBitNum = min(size, 1u << SAMPLEBITNUM) * EXPAND_THRESHOLD;
	cout << "hash_num\t" << hash_num << "\t"
		 << "ExpandBitNum\t" << ExpandBitNum << "\t" << endl;
	bloom_arr = new char[size >> 3];
	memset(bloom_arr, 0, size >> 3);
	for (int i = 0; i < hash_num; ++i)
		HashFunList[i].initialize(i + 1005);

	for (int i = 0; i < LOCK_NUM; i++)
		omp_init_lock(&(lock[i]));
}

/*clear the EBF*/
void Ebloom_filter::clear()
{
	_1_num = 0;
	memset(bloom_arr, 0, size >> 3);
	memset(&finger_buckets[0][0], 0, size * (BUCKET_SIZE + 2) * sizeof(uint16_t));
}

/*to free a class*/
Ebloom_filter::~Ebloom_filter()
{
	if (bloom_arr != NULL)
	{
		delete[] bloom_arr;
	}

	for (int i = 0; i < LOCK_NUM; i++)
		omp_destroy_lock(&(lock[i]));
}
/*Insertion*/

void Ebloom_filter::insert(uint32_t key)
{
	for (int i = 0; i < HASH_NUM; ++i)
	{
		// uint32_t pos = HashFunList[i].run((const char *)&key, KEY_LEN);
		uint32_t pos = MurmurHash3_x86_32(&key, KEY_LEN, i);
		uint32_t bid = pos & POSMASK;
		omp_set_lock(&(lock[pos & (LOCK_NUM - 1)]));

		lazy_update(0, bid);

		finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
		finger_buckets[bid][++finger_buckets[bid][0]] = pos >> SIZELOG;

		if (finger_buckets[bid][0] == BUCKET_SIZE + 1)
		{
			finger_buckets[bid][0]--;
			finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
			cout << "bucket overflows";
		}
		if ((bid >> SAMPLEBITNUM) == 0)
			if (!((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1))
				// #pragma opm atomic
				_1_num++;
		// _1_num ++;
		bloom_arr[bid >> SHIFT] |= (1 << (bid & MASK));

		omp_unset_lock(&(lock[pos & (LOCK_NUM - 1)]));

		if (((pos & 0x3FF) == 0) && i == 0)
		if (_1_num >= ExpandBitNum && expandOrNot)
		{
			expand();
		}
	}

	return;
}

void Ebloom_filter::insert(char *key)
{
	for (int i = 0; i < HASH_NUM; ++i)
	{
		// uint32_t pos = HashFunList[i].run((const char *)key, KEY_LEN2);
		uint32_t pos = MurmurHash3_x86_32(key, KEY_LEN2, i);
		uint32_t bid = pos & POSMASK;
		omp_set_lock(&(lock[pos & (LOCK_NUM - 1)]));

		lazy_update(0, bid);

		finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
		finger_buckets[bid][++finger_buckets[bid][0]] = pos >> SIZELOG;

		if (finger_buckets[bid][0] == BUCKET_SIZE + 1)
		{
			finger_buckets[bid][0]--;
			finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
			cout << "bucket overflows";
		}
		if ((bid >> SAMPLEBITNUM) == 0)
			if (!((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1))
#pragma opm atomic
				_1_num++;
		bloom_arr[bid >> SHIFT] |= (1 << (bid & MASK));

		omp_unset_lock(&(lock[pos & (LOCK_NUM - 1)]));

		if (((pos & 0x3FF) == 0) && i == 0)
			if (_1_num >= ExpandBitNum && expandOrNot)
			{
				expand();
			}
	}
	return;
}

inline void Ebloom_filter::lazy_update(const int layer, uint32_t bid)
{

	if (((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1) && finger_buckets[bid][0] && finger_buckets[bid][BUCKET_SIZE + 1] > FINGER_LENGTH)
	{

		uint32_t bucket_size = finger_buckets[bid][0],
				 dlen = (finger_buckets[bid][BUCKET_SIZE + 1] - FINGER_LENGTH),
				 t = (1 << dlen) - 1,
				 cnt = 0;

		for (int j = 1; j <= bucket_size; j++)
			if ((finger_buckets[bid][j] & t) == (bid >> (sizelog - dlen)))
				finger_buckets[bid][++cnt] = finger_buckets[bid][j] >> dlen;
		finger_buckets[bid][0] = cnt;
		finger_buckets[bid][BUCKET_SIZE + 1] = (cnt ? FINGER_LENGTH : 0);

		if (!cnt)
		{
			bloom_arr[bid >> SHIFT] &= ~(1u << (bid & MASK));
		}
		else
		{
			bloom_arr[bid >> SHIFT] |= (1u << (bid & MASK));
			if ((bid >> SAMPLEBITNUM) == 0)
#pragma opm atomic
				_1_num++;
		}
	}
}

void Ebloom_filter::insert_with_duplicates(uint32_t key)
{
	uint32_t pos[HASH_NUM];
	for (int i = 0; i < HASH_NUM; ++i)
		pos[i] = MurmurHash3_x86_32(&key, KEY_LEN, i);
		// pos[i] = HashFunList[i].run((const char *)&key, KEY_LEN);

	bool inBF = 1;
	for (int i = 0; i < HASH_NUM; ++i)
	{
		uint32_t bid = pos[i] & POSMASK;

		// omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));
		bool s = (bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1;

		// omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
		if (s == 0)
		{
			inBF = 0;
			break;
		}
	}

	if (inBF == true)
	{
		int i;
		for (i = min(32/(32-sizelog), HASH_NUM) - 1; i >= 0; --i)
		{
			uint32_t bid = pos[i] & POSMASK;
			omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));

			lazy_update(0, bid);

			int j;
			for (j = finger_buckets[bid][0]; j >= 1; --j)
				if (finger_buckets[bid][j] == (pos[i] >> SIZELOG))
					break;

			omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
			if (j == 0)
				break;
		}
		if (i == -1)
			return;
	}

	for (int i = 0; i < HASH_NUM; ++i)
	{
		uint32_t bid = pos[i] & POSMASK;

		omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));

		lazy_update(0, bid);

		finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
		finger_buckets[bid][++finger_buckets[bid][0]] = pos[i] >> SIZELOG;
		//cout<<"INS "<<bid<<endl;
		if (finger_buckets[bid][0] == BUCKET_SIZE + 1)
		{
			finger_buckets[bid][0]--;
			finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
			cout << "bucket overflows";
		}

		// set Bloom filter
		if ((bid >> SAMPLEBITNUM) == 0)
			if (!((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1))
#pragma opm atomic
				_1_num++;
		bloom_arr[bid >> SHIFT] |= (1 << (bid & MASK));
		omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
	}

	if (((pos[0] & 0x3FF) == 0))
		if (_1_num >= ExpandBitNum && expandOrNot)
		{
			expand();
		}
	return;
}


void Ebloom_filter::insert_with_duplicates(char *key)
{
	uint32_t pos[HASH_NUM];
	for (int i = 0; i < HASH_NUM; ++i)
		pos[i] = MurmurHash3_x86_32(key, KEY_LEN2, i);
		// pos[i] = HashFunList[i].run((const char *)key, KEY_LEN2);

	bool inBF = 1;
	for (int i = 0; i < HASH_NUM; ++i)
	{
		uint32_t bid = pos[i] & POSMASK;

		omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));
		lazy_update(0, bid);
		bool s = (bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1;

		omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
		if (s == 0)
		{
			inBF = 0;
			break;
		}
	}

	if (inBF == true)
	{
		int i;
		for (i = min(32/(32-sizelog), HASH_NUM) - 1; i >= 0; --i)
		{
			uint32_t bid = pos[i] & POSMASK;
			omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));

			lazy_update(0, bid);

			int j;
			for (j = finger_buckets[bid][0]; j >= 1; --j)
				if (finger_buckets[bid][j] == (pos[i] >> SIZELOG))
					break;

			omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
			if (j == 0)
				break;
		}
		if (i == -1)
			return;
	}

	for (int i = 0; i < HASH_NUM; ++i)
	{
		uint32_t bid = pos[i] & POSMASK;

		omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));

		lazy_update(0, bid);

		finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
		finger_buckets[bid][++finger_buckets[bid][0]] = pos[i] >> SIZELOG;
		//cout<<"INS "<<bid<<endl;
		if (finger_buckets[bid][0] == BUCKET_SIZE + 1)
		{
			finger_buckets[bid][0]--;
			finger_buckets[bid][BUCKET_SIZE + 1] = FINGER_LENGTH;
			cout << "bucket overflows";
		}

		// set Bloom filter
		if ((bid >> SAMPLEBITNUM) == 0)
			if (!((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1))
#pragma opm atomic
				_1_num++;
		bloom_arr[bid >> SHIFT] |= (1 << (bid & MASK));
		omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
	}

	if (((pos[0] & 0x3FF) == 0))
		if (_1_num >= ExpandBitNum && expandOrNot)
		{
			expand();
		}
	return;
}

void Ebloom_filter::deleteEle(uint32_t key)
{
	int i;
	for (int i = 0; i < HASH_NUM; i++)
	{
		// uint32_t pos = HashFunList[i].run((const char *)&key, KEY_LEN);
		uint32_t pos = MurmurHash3_x86_32(&key, KEY_LEN, i);
		uint32_t bid = pos & POSMASK;
		// omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));
		lazy_update(0, bid);
		int j;
		for (j = finger_buckets[bid][0]; j >= 1; --j)
			if (finger_buckets[bid][j] == (pos >> SIZELOG))
			{
				swap(finger_buckets[bid][j], finger_buckets[bid][finger_buckets[bid][0]]);
				finger_buckets[bid][0]--;
				break;
			}
		// if (j == 0)
		// 	printf("Deletion Error in bf!\n");
		if (j != 0)
			if (finger_buckets[bid][0] == 0)
			{
				if (((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1) == 0)
				{
					// printf("ERROR2\n");
				}
				else
				{
					bloom_arr[bid >> SHIFT] &= ~(1u << (bid & MASK));
					if ((bid >> SAMPLEBITNUM) == 0)
						_1_num--;
				}
			}
		// omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
	}
	return;
}

void Ebloom_filter::deleteEle(char *key)
{
	int i;
	for (int i = 0; i < HASH_NUM; i++)
	{
		// uint32_t pos = HashFunList[i].run((const char *)key, KEY_LEN2);
		uint32_t pos = MurmurHash3_x86_32(key, KEY_LEN2, i);
		uint32_t bid = pos & POSMASK;
		omp_set_lock(&(lock[bid & (LOCK_NUM - 1)]));
		lazy_update(0, bid);
		int j;
		for (j = finger_buckets[bid][0]; j >= 1; --j)
			if (finger_buckets[bid][j] == (pos >> SIZELOG))
			{
				swap(finger_buckets[bid][j], finger_buckets[bid][finger_buckets[bid][0]]);
				finger_buckets[bid][0]--;
				break;
			}
		// if (j == 0)
		// 	printf("Deletion Error in bf!\n");
		if (j != 0)
			if (finger_buckets[bid][0] == 0)
			{
				if (((bloom_arr[bid >> SHIFT] >> (bid & MASK)) & 1) == 0)
				{
					// cout<<ss<<" "<<POSMASK<<endl;
					// printf("ERROR2\n");
				}
				else
				{
					bloom_arr[bid >> SHIFT] &= ~(1u << (bid & MASK));
					if ((bid >> SAMPLEBITNUM) == 0)
						_1_num--;
				}
			}
		omp_unset_lock(&(lock[bid & (LOCK_NUM - 1)]));
	}
	return;
}

/*Query*/
bool Ebloom_filter::query(uint32_t key)
{
	for (int i = 0; i < HASH_NUM; ++i)
	{
		// uint32_t pos = HashFunList[i].run((const char *)&key, KEY_LEN);
		uint32_t pos = MurmurHash3_x86_32(&key, KEY_LEN, i);
		uint32_t bid = pos & POSMASK;

		// omp_set_lock(&(lock[pos&(LOCK_NUM-1)]));
		bool result = bloom_arr[bid >> SHIFT] & (1 << (bid & MASK));
		// omp_unset_lock(&(lock[pos&(LOCK_NUM-1)]));

		if (!result)
		{
			return false;
		}
	}
	return true;
}

// bool Ebloom_filter::query(char *key)
// {
// 	for (int i = 0; i < HASH_NUM; ++i)
// 	{
// 		// uint32_t pos = HashFunList[i].run((const char *)key, KEY_LEN2);
// 		uint32_t pos = MurmurHash3_x86_32(key, KEY_LEN2, i);
// 		uint32_t bid = pos & POSMASK;

// 		// omp_set_lock(&(lock[pos&(LOCK_NUM-1)]));
// 		bool result = bloom_arr[bid >> SHIFT] & (1 << (bid & MASK));
// 		// omp_unset_lock(&(lock[pos&(LOCK_NUM-1)]));

// 		if (!result)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }

bool Ebloom_filter::query(char *key)
{
	for (int i = 0; i < HASH_NUM; ++i)
	{
		// uint32_t pos = HashFunList[i].run((const char *)key, KEY_LEN2);
		uint32_t pos = MurmurHash3_x86_32(key, KEY_LEN2, i);
		omp_set_lock(&(lock[pos&(LOCK_NUM-1)]));
		lazy_update(0, pos & POSMASK);
		// if(!(bloom_arr[(pos & POSMASK) >> SHIFT] & (1 << (pos & MASK))))
		omp_unset_lock(&(lock[pos&(LOCK_NUM-1)]));
		if(!((bloom_arr[(pos & POSMASK) >> SHIFT] >> (pos & MASK))&1))
			return false;
	}
	return true;
}
/*Expand the buckets*/
void Ebloom_filter::expand()
{
#pragma omp critical
	{
		if (_1_num >= ExpandBitNum && expandOrNot)
		{
			timespec dtime1, dtime2;
			long long dresns;

			double rens=0;
			// cout << "EBF expansion." << endl;
			// cout << "\t Bloom sizelog" << sizelog << "\t _1_num:" << _1_num << " \tExpandBitNum:" << ExpandBitNum << endl;

			if (sizelog == MAX_BLOOM_SIZE)
			{
				cout << "EBF expansion ERROR. No enough bucket space." << endl;
			}
			assert((MAX_SIZE - finger_length + 1) > MIN_SIZE && (finger_length - 1) > 0);

			// #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
			rens-=omp_get_wtime();
			for (int i = LOCK_NUM - 1; i >= 0; i--)
				omp_set_lock(&(lock[i & (LOCK_NUM - 1)]));

			clock_gettime(CLOCK_MONOTONIC, &dtime1);
			uint32_t oldsize = size;
			char *old_bloom_arr = bloom_arr;

			char *new_arr = new char[(size << 1) >> 3];
			memcpy(new_arr, bloom_arr, size >> 3);
			memcpy(new_arr + (size >> 3), bloom_arr, size >> 3);
			fastMemcpy(&finger_buckets[size][0], &finger_buckets[0][0], size * (BUCKET_SIZE + 2) * sizeof(uint16_t));

			/*reallocate the bit map*/
			bloom_arr = new_arr;

			/*change the parameters*/
			compression--;
			finger_length--;
			_1_num = 0;

			/*expand the buckets*/
			size <<= 1;
			sizelog++;
			posmask = (1 << sizelog) - 1;
			ExpandBitNum = min(size, 1u << SAMPLEBITNUM) * EXPAND_THRESHOLD;

			for (int i = 0; i < min(size, 1u << SAMPLEBITNUM); i++)
			{
				lazy_update(0, i);
			}
			// cout << "expand end" << endl;
			for (int i = LOCK_NUM - 1; i >= 0; i--)
				omp_unset_lock(&(lock[i & (LOCK_NUM - 1)]));


			// for (int i = (1u << SAMPLEBITNUM); i < size; i++)
			// {
			// 	// omp_set_lock(&(lock[i & (LOCK_NUM - 1)]));
			// 	lazy_update(0, i);
			// 	// omp_unset_lock(&(lock[i & (LOCK_NUM - 1)]));
			// }
			rens+=omp_get_wtime();

			clock_gettime(CLOCK_MONOTONIC, &dtime2);
			dresns = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
			double dth = (double)dresns / 1e9;
			// cout << "time1----------" << dth << endl;
			// sumtime += dth;
			sumtime += rens;


			delete[] old_bloom_arr;

			// uint32_t local_1_num = 0;
			// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
			// for (int i = oldsize - 1; i >= 0; i--)
			// {
			// 	// //cout<<i<<endl;
			// 	// omp_set_lock(&(lock[i&(LOCK_NUM-1)]));

			// 	//if(sizelog==22)//cout<<i<<endl;
			// 	int BucketNum = finger_buckets[i][0];
			// 	// if(sizelog==22)//cout<<i<<" "<<BucketNum<<endl;
			// 	finger_buckets[i][0] = 0;
			// 	finger_buckets[i + oldsize][0] = 0;
			// 	// if(sizelog==22)//cout<<i<<" "<<BucketNum<<endl;
			// 	for (int j = 1; j <= BucketNum; j++){
			// 		// if(i==2097105)
			// 			// //cout<<j<<endl;
			// 		if (finger_buckets[i][j] & 1)
			// 		{
			// 			finger_buckets[i + oldsize][++finger_buckets[i + oldsize][0]] = finger_buckets[i][j]>>1;
			// 		}
			// 		else
			// 		{
			// 			finger_buckets[i][++finger_buckets[i][0]] = finger_buckets[i][j]>>1;
			// 		}
			// 	}
			// 	if (finger_buckets[i][0])
			// 	{
			// 		bloom_arr[i >> SHIFT] |= (1 << (i & MASK));
			// #pragma opm atomic// 		local
			// _1_num++;
			// // 	}
			// // 	if (finger_buckets[i + oldsize][0])
			// // 	{
			// // 		bloom_arr[(i  + oldsize) >> SHIFT] |= (1 << ((i  + oldsize) & MASK));
			// #pragma opm atomic// 		local
			// _1_num++;
			// 	}
			// 	// omp_unset_lock(&(lock[i&(LOCK_NUM-1)]));
			// }

			// _1_num=local_1_num;
			// // //cout<<"\t Bloom sizelog"<<sizelog<<endl;
			// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		}
	}
}

/*compress the buckets*/
void Ebloom_filter::compress()
{
#pragma omp critical
	{
		// if (_1_num <= ExpandBitNum/4)
		{
			timespec dtime1, dtime2;
			long long dresns;

			// cout << "EBF compress." << endl;
			// cout << "\t Bloom sizelog" << sizelog << "\t _1_num:" << _1_num << " \tExpandBitNum:" << ExpandBitNum << endl;

			// #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
			for (int i = LOCK_NUM - 1; i >= 0; i--)
				omp_set_lock(&(lock[i & (LOCK_NUM - 1)]));

			uint32_t oldsize = size;
			char *old_bloom_arr = bloom_arr;

			char *new_arr = new char[(size >> 1) >> 3];

			for (int pos = 0; pos < size; pos++)
			{
				lazy_update(0, pos);
			}

			for (int i = ((size >> 3) >> 1) - 1; i >= 0; i--)
			{
				new_arr[i] = old_bloom_arr[i] | old_bloom_arr[i + ((size >> 3) >> 1)];
			}
			_1_num = 0;
			for (int pos = 0; pos < (size >> 1); pos++)
			{
				if (finger_buckets[pos][0] + finger_buckets[pos + (size >> 1)][0] > BUCKET_SIZE)
				{
					printf("Compression Error!\n");
				}
				for (int j = finger_buckets[pos][0]; j >= 1; j--)
					finger_buckets[pos][j] <<= 1;
				for (int j = finger_buckets[pos + (size >> 1)][0]; j >= 1; j--)
					finger_buckets[pos][j + finger_buckets[pos][0]] = (finger_buckets[pos + (size >> 1)][j] << 1) | 1;
				finger_buckets[pos][0] += finger_buckets[pos + (size >> 1)][0];
				if (finger_buckets[pos][0])
				{
					if (finger_buckets[pos][BUCKET_SIZE + 1] == 0)
						cout << "ERROR3\n";
					finger_buckets[pos][BUCKET_SIZE + 1]++;
				}
				if ((pos >> SAMPLEBITNUM) == 0)
					if (finger_buckets[pos][0])
					{
						cout << pos << "AAAAAAAAAAAAAAAAA" << endl;
#pragma opm atomic
						_1_num++;
					}
			}
			/*reallocate the bit map*/
			bloom_arr = new_arr;

			/*change the parameters*/
			compression++;
			finger_length++;

			/*compress the buckets*/
			size >>= 1;
			sizelog--;
			posmask = (1 << sizelog) - 1;
			ExpandBitNum = min(size, 1u << SAMPLEBITNUM) * EXPAND_THRESHOLD;

			for (int i = 0; i < min(size, 1u << SAMPLEBITNUM); i++)
			{
				lazy_update(0, i);
			}
			for (int i = LOCK_NUM - 1; i >= 0; i--)
				omp_unset_lock(&(lock[i & (LOCK_NUM - 1)]));

			delete[] old_bloom_arr;
			for (int i = (1u << SAMPLEBITNUM); i < size; i++)
			{
				if(((bloom_arr[i >> SHIFT] >> (i & MASK)) & 1)){
					// omp_set_lock(&(lock[i & (LOCK_NUM - 1)]));
					lazy_update(0, i);
					// omp_unset_lock(&(lock[i & (LOCK_NUM - 1)]));
				}
			}
			// uint32_t local_1_num = 0;
			// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
			// for (int i = oldsize - 1; i >= 0; i--)
			// {
			// 	// //cout<<i<<endl;
			// 	// omp_set_lock(&(lock[i&(LOCK_NUM-1)]));

			// 	//if(sizelog==22)//cout<<i<<endl;
			// 	int BucketNum = finger_buckets[i][0];
			// 	// if(sizelog==22)//cout<<i<<" "<<BucketNum<<endl;
			// 	finger_buckets[i][0] = 0;
			// 	finger_buckets[i + oldsize][0] = 0;
			// 	// if(sizelog==22)//cout<<i<<" "<<BucketNum<<endl;
			// 	for (int j = 1; j <= BucketNum; j++){
			// 		// if(i==2097105)
			// 			// //cout<<j<<endl;
			// 		if (finger_buckets[i][j] & 1)
			// 		{
			// 			finger_buckets[i + oldsize][++finger_buckets[i + oldsize][0]] = finger_buckets[i][j]>>1;
			// 		}
			// 		else
			// 		{
			// 			finger_buckets[i][++finger_buckets[i][0]] = finger_buckets[i][j]>>1;
			// 		}
			// 	}
			// 	if (finger_buckets[i][0])
			// 	{
			// 		bloom_arr[i >> SHIFT] |= (1 << (i & MASK));
			// #pragma opm atomic// 		local
			// _1_num++;
			// // 	}
			// // 	if (finger_buckets[i + oldsize][0])
			// // 	{
			// // 		bloom_arr[(i  + oldsize) >> SHIFT] |= (1 << ((i  + oldsize) & MASK));
			// #pragma opm atomic// 		local
			// _1_num++;
			// 	}
			// 	// omp_unset_lock(&(lock[i&(LOCK_NUM-1)]));
			// }

			// _1_num=local_1_num;
			// // //cout<<"\t Bloom sizelog"<<sizelog<<endl;
			// // #pragma omp parallel for schedule(guided) num_threads(THREAD_NUM)
		}
	}
}

/*Set one bit*/
inline void Ebloom_filter::setbit(uint32_t pos)
{
	bloom_arr[pos >> SHIFT] |= (1 << (pos & MASK));
}

inline void Ebloom_filter::resetbit(uint32_t pos)
{
	bloom_arr[pos >> SHIFT] &= ~(1 << (pos & MASK));
}

/*Get the value of one bit*/
inline bool Ebloom_filter::getbit(uint32_t pos)
{
	if (bloom_arr[pos >> SHIFT] & (1 << (pos & MASK)))
		return true;
	return false;
}

/* get the number of 0 bits */
inline uint32_t Ebloom_filter::get_0_num()
{
	uint32_t ans = 0;
	for (int i = 0; i < size; ++i)
	{
		if (!getbit(i))
			++ans;
	}
	return ans;
}

/* get the number of 0 bits */
inline int Ebloom_filter::get_1_num()
{
	return _1_num;
}

/* get the size of bitmap */
inline uint32_t Ebloom_filter::get_size()
{
	return size;
}

#endif
