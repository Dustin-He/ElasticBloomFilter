#ifndef _CBF_H
#define _CBF_H

#include "param.h"
#include "BOBHash32.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <string>
#include <assert.h>

class CBF {
	public:
		CBF(uint32_t sz, int hash_num);
		~CBF();
		uint32_t _1_num;
		int insert(char * key);
		bool query(char * key);
		int del(char * key);
		void expand();
		uint32_t get_size();
	private:
		uint32_t size;
		uint32_t csize;
		uint16_t *counter;
		char *bitmap;
		vector<BOBHash32*> HashFunList;
		inline bool getbit(uint32_t pos);
		inline void setbit(uint32_t pos);
		inline void resetbit(uint32_t pos);
};

CBF::CBF(uint32_t sz, int hash_num) {
	size = 1 << sz;
	csize = size << 10;
	counter = new uint16_t[csize];
	bitmap = new char[size >> 3];
	_1_num = 0;
	memset(bitmap, 0, size >> 3);
	memset(counter, 0, csize);

	for (int i = 0; i < hash_num; ++i)
		HashFunList.push_back(new BOBHash32(i + 750));
}

CBF::~CBF() {
	delete []counter;
	delete []bitmap;
	for (auto ele : HashFunList) {
		delete ele;
	}
}

int CBF::insert(char * key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t tmp = HashFunList[i]->run((const char *)key, KEY_LEN2);
		uint32_t pos = tmp % size;
		uint32_t pos2 = tmp % csize;
		if (counter[pos2] == (1 << 16) - 1)
			return 0;
		if (!getbit(pos)) _1_num++;
		setbit(pos);
		counter[pos2] += (counter[pos2] <= (1 << 16) - 1);
	}
	if ((double)_1_num / (double)size >= EXPAND_THRESHOLD && expandOrNot)
		expand();
	return 1;
}

void CBF::expand() {
	delete []bitmap;
	size = size << 1;
	bitmap = new char[size << 3];
	memset(bitmap, 0, size << 3);
	_1_num = 0;

	for (int i = 0; i < size; ++i) {
		for (int k = i; k < csize; k += size) {
			if (counter[k]) {
				_1_num++;
				setbit(i);
				break;
			}
		}
	}
}

bool CBF::query(char * key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)key, KEY_LEN2);
		if (!getbit(pos % size)) {
			return false;
		}
	}
	return true;
}

int CBF::del(char * key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t tmp = HashFunList[i]->run((const char *)key, KEY_LEN2);
		uint32_t pos =  tmp % size;
		uint32_t pos2 = tmp % csize;
		if (counter[pos2])
			counter[pos2]--;
		bool flag = false;
		for (int j = pos; j < csize; j += size)
			if (counter[j])
				flag = true;
		if (!flag)
			resetbit(pos);
	}
	return 0;
}

inline void CBF::setbit(uint32_t pos) {
	bitmap[pos >> SHIFT] |= (1 << (pos & MASK));
}

inline bool CBF::getbit(uint32_t pos) {
	if (bitmap[pos >> SHIFT] & (1 << (pos & MASK)))
		return true;
	return false;
}

inline void CBF::resetbit(uint32_t pos) {
	bitmap[pos >> SHIFT] &= ~(1 << (pos & MASK));
}

uint32_t CBF::get_size() {
	return size;
}

#endif
