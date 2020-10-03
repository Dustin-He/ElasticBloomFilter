#ifndef _DBF_H
#define _DBF_H

#include<iostream>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cmath>
#include "BOBHash32.h"
#include "param.h"

using namespace std;

class BFint
{
	public:
	int* c;
	int num;
	BFint() {};
	BFint(int m)
	{
		num = 0;
		c = new int[1<<m];
		memset(c, 0, sizeof(int) * (1<<m));
	}
};
class DynamicBF
{
	public:
	BOBHash32 * bobhash[HASH_NUM + 2];
	int w, m, maxn;
	double c;
	BFint* b[1010];
	DynamicBF(){};
	DynamicBF(int _m,double _c)
	{
		w = 0;
		m = _m;
		c = _c;
		maxn = c * (1<<m);
		memset(b, 0, sizeof(b));
		b[0] = new BFint(m);
		for(int i = 0; i < HASH_NUM; i++)
			bobhash[i] = new BOBHash32(i + 1005);
	}
	void insert(int x)
	{
		bool f = false;
		if(b[w]->num >= maxn)
			b[++w] = new BFint(m);
		for(int i = 0; i < HASH_NUM; i++)
		{
			int z = bobhash[i]->run((char *)&x,4) & ((1<<m)-1);
			if(!b[w]->c[z])
			{
				b[w]->num++;
			}
			b[w]->c[z]+=1;
		}
	}

	void insert(char *x)
	{
		bool f = false;
		if(b[w]->num >= maxn)
			b[++w] = new BFint(m);
		for(int i = 0; i < HASH_NUM; i++)
		{
			int z = bobhash[i]->run(x,13) & ((1<<m)-1);
			if(!b[w]->c[z])
			{
				b[w]->num++;
			}
			b[w]->c[z]+=1;
		}
	}

	bool query(int x)
	{
		int *z = new int[HASH_NUM];
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run((char *)&x,4)&((1<<m)-1);
		for(int i = 0; i <= w;i++)
		{
			bool f=1;
			for(int j = 0; j < HASH_NUM; j++)
			{
				if(!b[i]->c[z[j]])
				{
					f=0;
					break;
				}
			}
			if(f)
			{
				delete []z;
				return true;
			}
		}
		delete []z;
		return false;
	}
	
	bool query(char * x)
	{
		int *z = new int[HASH_NUM];
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run(x,13) & ((1<<m)-1);
		for(int i = 0; i <= w;i++)
		{
			bool f=1;
			for(int j = 0; j < HASH_NUM; j++)
			{
				if(!b[i]->c[z[j]])
				{
					f=0;
					break;
				}
			}
			if(f)
			{
				delete []z;
				return true;
			}
		}
		delete []z;
		return false;
	}

	void deleteEle(int x) {
		int cnt = 0;
		int deletePos = -1;
		int *z = new int[HASH_NUM];
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run((char *)&x,4) & ((1<<m)-1);
		for(int i = 0; i <= w;i++) {
			int f = 1;
			for(int j = 0; j < HASH_NUM; j++) {
				if(!b[i]->c[z[j]]) {
					f=0;
					break;
				}
			}
			if (f) {
				cnt += f;
				deletePos = i;
			}
		}
		if (cnt == 1) {
			for (int j = 0; j < HASH_NUM; ++j) {
				if (b[deletePos]->c[z[j]] > 0)
					b[deletePos]->c[z[j]]--;
			}
		}
	}

	void deleteEle(char * x) {
		int cnt = 0;
		int deletePos = -1;
		int *z = new int[HASH_NUM];
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run(x,13) & ((1<<m)-1);
		for(int i = 0; i <= w;i++) {
			int f = 1;
			for(int j = 0; j < HASH_NUM; j++) {
				if(!b[i]->c[z[j]]) {
					f=0;
					break;
				}
			}
			if (f) {
				cnt += f;
				deletePos = i;
			}
		}
		if (cnt == 1) {
			for (int j = 0; j < HASH_NUM; ++j) {
				if (b[deletePos]->c[z[j]] > 0)
					b[deletePos]->c[z[j]]--;
			}
		}
	}
};

#endif
