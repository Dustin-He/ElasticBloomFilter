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

using namespace std;

extern uint32_t total_slot;

class Ebloom_filter {
	public:
		Ebloom_filter(uint32_t size, int hash_num);
		~Ebloom_filter();
		int insert(uint32_t key);
		bool query(uint32_t key);
		void expand();
		inline int get_compression();
		inline uint32_t get_0_num();
		inline uint32_t get_size();
		inline uint32_t get_1_num();
	private:
		uint32_t size;			//the length of the bit map
		
		uint32_t finger_length;	//the length of the fingerprint
		int compression;		//to record the expanding times(= -compression)

		uint32_t _1_num;		//the number of those set bits
		
		char *bloom_arr;		//the bit map
		vector<BOBHash32*> HashFunList;						//the hash function list
		vector<vector<uint32_t>*> finger_buckets[LAYERS];	//to simulate finger buckets in slow memory
		inline bool getbit(uint32_t pos);					//get a bit in the bit map
		inline void setbit(uint32_t pos);					//set a bit in the bit map
		int insert_finger(uint32_t pos, uint32_t key);		//to be called when inserting a fingerprint
};


/*build the EBF*/
Ebloom_filter::Ebloom_filter(uint32_t sz, int hash_num) {
	assert(hash_num > 0);
	assert(sz < MAX_SIZE && sz > MIN_SIZE);

	finger_length = MAX_SIZE - sz;
	compression = _1_num = 0;

	/*to allocate the bit map and the hash function*/
	size = 1 << sz;
	bloom_arr = new char[size >> 3];
	memset(bloom_arr, 0, size >> 3);
	for (int i = 0; i < hash_num; ++i)
		HashFunList.push_back(new BOBHash32(i + 750));

	/*to allocate the buckets*/
	uint32_t sum = 0;
	uint32_t end = size;
	int j = 0;
	while (j < LAYERS) {
		sum += end;
		for (int i = 0; i < end; ++i) {
			finger_buckets[j].push_back(new vector<uint32_t>());
		}
		end = (end - 1) / BRANCH + 1;
		j++;
	}
	total_slot = sum * BUCKET_SIZE;
	//cout << "Total slot: " << total_slot << endl;
}

/*to free a class*/
Ebloom_filter::~Ebloom_filter() {
	if (bloom_arr != NULL) {
		delete []bloom_arr;
	}
	for (auto ele : HashFunList) {
		delete ele;
	}
	for (int j = 0; j < LAYERS; ++j) {
		for (auto ele : finger_buckets[j]) {
			delete ele;
		}
	}
}


int Ebloom_filter::insert_finger(uint32_t pos, uint32_t key) {
	int j = 0;
	while (j < LAYERS) {
		if (finger_buckets[j][pos]->size() < BUCKET_SIZE) {
			finger_buckets[j][pos]->push_back(key);
			return j;
		}
		/* We need one or several bits to record the bucket's path to the root */
		key = (key << BSHIFT) | (pos % BRANCH);
		pos /= BRANCH;
		j++;
	}
	if (j == LAYERS)
		printf("Insert Error in bf!\n");
	return LAYERS;
}

/*Insertion*/
int Ebloom_filter::insert(uint32_t key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)&key, KEY_LEN);
		if (!getbit(pos % size)) {
			setbit(pos % size);
			/* this is for the experiment, it can be commented when not needed. */
			++_1_num;
		}
		int ans;
		if ((ans = insert_finger(pos % size, pos / size)) == LAYERS) {
			return 0;
		}
		
	}
	/* This is the expanding operation. It can be commented when testing the insertion speed. */
	
	if ((double)_1_num / (double)size >= EXPAND_THRESHOLD && expandOrNot) {
		//cout << (double)_1_num / (double)size << endl;
		expand();
		return -1;
	}
	
	return 1;
}

/*Query*/
bool Ebloom_filter::query(uint32_t key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)&key, KEY_LEN);
		if (!getbit(pos % size)) {
			return false;
		}
	}
	return true;
}

/*Set one bit*/
inline void Ebloom_filter::setbit(uint32_t pos) {
	bloom_arr[pos >> SHIFT] |= (1 << (pos & MASK));
}

/*Get the value of one bit*/
inline bool Ebloom_filter::getbit(uint32_t pos) {
	if (bloom_arr[pos >> SHIFT] & (1 << (pos & MASK)))
		return true;
	return false;
}

/*Expand the buckets*/
void Ebloom_filter::expand() {
	assert((MAX_SIZE - finger_length + 1) > MIN_SIZE && (finger_length - 1) > 0);

	/*change the parameters*/
	compression--;
	finger_length--;

	/*expand the buckets*/
	uint32_t end = size;
	int h = 0;
	while (h < LAYERS) {
		for (int i = 0; i < end; ++i) {
			finger_buckets[h].push_back(new vector<uint32_t>());
		}
		end = (end - 1) / BRANCH + 1;
		h++;
	}

	uint32_t tmpsize = size << 1;

	/*change the fingerprints and insert into different buckets*/
	h = 0;
	end = size;
	uint32_t pos;
	uint32_t bias;

	while (h < LAYERS) {
		for (int i = 0; i < end; ++i) {
			vector<uint32_t> ins_ele;
			int sz = finger_buckets[h][i]->size();
			for (int k = 0; k < sz; ++k) {
				uint32_t ele = finger_buckets[h][i]->at(k);
				uint32_t tmp = (ele >> (BSHIFT * h)) & 1;
				if (tmp) {
					pos = i;
					for (int j = 0; j < h; ++j) {
						bias = (ele & (((1 << BSHIFT) - 1) << (j * BSHIFT))) >> (j * BSHIFT);
						pos = pos * BRANCH + bias;
						
					}
					insert_finger(pos + size, ele >> (BSHIFT * h + 1));
				}
				else {
					ins_ele.push_back(ele);
				}
			}
			if (sz) {
				finger_buckets[h][i]->clear();
				for (auto ele : ins_ele) {
					pos = i;
					for (int j = 0; j < h; ++j) {
						bias = (ele & (((1 << BSHIFT) - 1) << (j * BSHIFT))) >> (j * BSHIFT);
						pos = pos * BRANCH + bias;
					}
					insert_finger(pos, ele >> (BSHIFT * h + 1));
				}
			}
		}
		h++;
		end = (end - 1) / BRANCH + 1;
	}

	size = tmpsize;

	/*reallocate the bit map*/
	
	delete []bloom_arr;
	char *new_arr = new char[size >> 3];
	memset(new_arr, 0, size >> 3);
	bloom_arr = new_arr;
	_1_num = 0;
	for (int i = 0; i < size; ++i) {
		if (!finger_buckets[0][i]->empty()) {
			setbit(i);
			_1_num++;
		}
	}
	
}


/* get the compression */
inline int Ebloom_filter::get_compression() {
	return compression;
}

/* get the number of 0 bits */
inline uint32_t Ebloom_filter::get_0_num() {
	uint32_t ans = 0;
	for (int i = 0; i < size; ++i) {
		if (!getbit(i)) ++ans;
	}
	return ans;
}

/* get the number of set bits */
inline uint32_t Ebloom_filter::get_1_num() {
	return _1_num;
}

/* get the size of bitmap */
inline uint32_t Ebloom_filter::get_size() {
	return size;
}

#endif