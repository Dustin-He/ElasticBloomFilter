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
		int insert(uint32_t key);
		bool query(uint32_t key);
		int del(uint32_t key);
	private:
		uint32_t size;
		uint8_t *counter;
		char *bitmap;
		vector<BOBHash32*> HashFunList;
		inline bool getbit(uint32_t pos);
		inline void setbit(uint32_t pos);
		inline void resetbit(uint32_t pos);
};

CBF::CBF(uint32_t sz, int hash_num) {
	size = 1 << sz;
	counter = new uint8_t[size];
	bitmap = new char[size >> 3];

	memset(bitmap, 0, size >> 3);
	memset(counter, 0, size);

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

int CBF::insert(uint32_t key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)&key, KEY_LEN) % size;
		if (counter[pos] == 255) 
			return -1;
		setbit(pos);
		counter[pos] += (counter[pos] <= 255);
	}
	return 0;
}

bool CBF::query(uint32_t key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)&key, KEY_LEN);
		if (!getbit(pos % size)) {
			return false;
		}
	}
	return true;
}

int CBF::del(uint32_t key) {
	int hash_num = HashFunList.size();
	for (int i = 0; i < hash_num; ++i) {
		uint32_t pos = HashFunList[i]->run((const char *)&key, KEY_LEN) % size;
		if (counter[pos])
			counter[pos]--;
		if (!counter[pos])
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

#endif