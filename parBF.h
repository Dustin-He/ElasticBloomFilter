#ifndef _PBF_H
#define _PBF_H

#include<iostream>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cmath>
#include "BOBHash32.h"
#include "param.h"

using namespace std;

class BFuint
{
	public:
	int* c;
	int num;
	BFuint() {};
	BFuint(int m)
	{
		num = 0;
		c = new int[1<<m];
		memset(c, 0, sizeof(int) * (1<<m));
	}
};
class parBF
{
	public:
	BOBHash32 * bobhash[HASH_NUM + 2];
	int w, m, maxn;
	double c;
	BFuint* b[1010];
	parBF(){};
	parBF(int _m,double _c)
	{
		w = 0;
		m = _m;
		c = _c;
		maxn = c * (1<<m);
		memset(b, 0, sizeof(b));
		b[0] = new BFuint(m);
		for(int i = 0; i < HASH_NUM; i++)
			bobhash[i] = new BOBHash32(i + 1005);
	}
	void insert(int x)
	{
		bool f = false;
		if(b[w]->num >= maxn)
			b[++w] = new BFuint(m);
		for(int i = 0; i < HASH_NUM; i++)
		{
			int z = bobhash[i]->run((char *)&x,4) & ((1<<m)-1);
			if(!b[w]->c[z])
			{
				#pragma opm atomic
				b[w]->num++;
			}
			#pragma opm atomic
			b[w]->c[z]++;
		}
	}

	void insert(char *x)
	{
		bool f = false;
		if(b[w]->num >= maxn)
			b[++w] = new BFuint(m);
		for(int i = 0; i < HASH_NUM; i++)
		{
			int z = bobhash[i]->run(x,13) & ((1<<m)-1);
			if(!b[w]->c[z])
			{
				#pragma opm atomic
				b[w]->num++;
			}
			#pragma opm atomic
			b[w]->c[z]++;
		}
	}

	
	bool query(int x)
	{
		int z[HASH_NUM];
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run((char *)&x,4)&((1<<m)-1);
		 //#pragma omp parallel for schedule(guided) num_threads(16) reduction(|:ret)
		for(int i = 0; i <= w;i++)
		{
			
			{
				bool f=1;
				for(int j = 0; j < HASH_NUM; j++)
				{
					if(!b[i]->c[z[j]])
					{
						f= 0;
						break;
					}
				}
				if(f)
				{
					return true;
				}
			}
		}
		return false;
	}
	
	bool query(char * x)
	{
		int z[HASH_NUM];
		 //#pragma omp parallel for schedule(guided) num_threads(4)
		for(int j = 0; j < HASH_NUM; j++)
			z[j] = bobhash[j]->run(x,13) & ((1<<m)-1);
		bool ret=0;
		for(int i = 0; i <= w;i++)
		{
			
			{
				bool f=1;
				for(int j = 0; j < HASH_NUM; j++)
				{
					if(!b[i]->c[z[j]])
					{
						f = 0;
						break;
					}
				}
				if(f)
				{
					ret=1;
				}
			}
		}
		return ret;
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
