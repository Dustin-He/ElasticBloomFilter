#ifndef _DBF_H
#define _DBF_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "BOBHash32.h"
#include "param.h"
#include "murmur3.h"

using namespace std;

class BFint
{
public:
	int *c;
	int num;
	BFint(){};
	BFint(int m)
	{
		num = 0;
		c = new int[1 << m];
		memset(c, 0, sizeof(int) * (1 << m));
	}
};
class DynamicBF
{
public:
	BOBHash32 *bobhash[HASH_NUM + 2];
	int w, m, maxn;
	double c;
	BFint *b[1010];
	DynamicBF(){};
	DynamicBF(int _m, double _c)
	{
		w = 0;
		m = _m;
		c = _c;
		maxn = c * (1 << m);
		memset(b, 0, sizeof(b));
		b[0] = new BFint(m);
		for (int i = 0; i < HASH_NUM; i++)
			bobhash[i] = new BOBHash32(i + 1005);
	}
	void insert(int x)
	{
		bool f = false;
#pragma omp critical
		{
			if (b[w]->num >= maxn)
			{
				{
#pragma opm atomic
					b[++w] = new BFint(m);
				}
			}
		}
		for (int i = 0; i < HASH_NUM; i++)
		{
			int z = MurmurHash3_x86_32((char *)&x, 4, i) & ((1 << m) - 1);
			// int z = bobhash[i]->run((char *)&x, 4) & ((1 << m) - 1);
			if (!b[w]->c[z])
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
#pragma omp critical
		{
			if (b[w]->num >= maxn)
			{
				{
#pragma opm atomic
					b[++w] = new BFint(m);
				}
			}
		}
		for (int i = 0; i < HASH_NUM; i++)
		{
			int z = MurmurHash3_x86_32((char *)x, 13, i) & ((1 << m) - 1);
			// int z = bobhash[i]->run(x, 13) & ((1 << m) - 1);
			if (!b[w]->c[z])
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
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)&x, 4, j) & ((1 << m) - 1);
			// z[j] = bobhash[j]->run((char *)&x, 4) & ((1 << m) - 1);
		// #pragma omp parallel for schedule(guided) num_threads(16) reduction(|:ret)
		for (int i = 0; i <= w; i++)
		{

			{
				bool f = 1;
				for (int j = 0; j < HASH_NUM; j++)
				{
#pragma opm atomic
					if (!b[i]->c[z[j]])
					{
						f = 0;
						break;
					}
				}
				if (f)
				{
					return true;
				}
			}
		}
		return false;
	}

	bool query(char *x)
	{
		int z[HASH_NUM];
		// #pragma omp parallel for schedule(guided) num_threads(4)
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)x, 13, j) & ((1 << m) - 1);
			// z[j] = bobhash[j]->run(x, 13) & ((1 << m) - 1);
		bool ret = 0;
		for (int i = 0; i <= w; i++)
		{

			{
				bool f = 1;
				for (int j = 0; j < HASH_NUM; j++)
				{
#pragma opm atomic
					if (!b[i]->c[z[j]])
					{
						f = 0;
						break;
					}
				}
				if (f)
				{
					ret = 1;
				}
			}
		}
		return ret;
	}

	void deleteEle(int x)
	{
		int cnt = 0;
		int deletePos = -1;
		int *z = new int[HASH_NUM];
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)&x, 4, j) & ((1 << m) - 1);
			// z[j] = bobhash[j]->run((char *)&x, 4) & ((1 << m) - 1);
		for (int i = 0; i <= w; i++)
		{
			int f = 1;
			for (int j = 0; j < HASH_NUM; j++)
			{
#pragma opm atomic
				if (!b[i]->c[z[j]])
				{
					f = 0;
					break;
				}
			}
			if (f)
			{
				cnt += f;
				deletePos = i;
			}
		}
		if (cnt == 1)
		{
			for (int j = 0; j < HASH_NUM; ++j)
			{
#pragma opm atomic
				if (b[deletePos]->c[z[j]] > 0)
#pragma opm atomic
					b[deletePos]->c[z[j]]--;
			}
		}
	}

	void deleteEle(char *x)
	{
		int cnt = 0;
		int deletePos = -1;
		int *z = new int[HASH_NUM];
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)x, 13, j) & ((1 << m) - 1);
			// z[j] = bobhash[j]->run(x, 13) & ((1 << m) - 1);
		for (int i = 0; i <= w; i++)
		{
			int f = 1;
			for (int j = 0; j < HASH_NUM; j++)
			{
#pragma opm atomic
				if (!b[i]->c[z[j]])
				{
					f = 0;
					break;
				}
			}
			if (f)
			{
				cnt += f;
				deletePos = i;
			}
		}
		if (cnt == 1)
		{
			for (int j = 0; j < HASH_NUM; ++j)
			{
#pragma opm atomic
				if (b[deletePos]->c[z[j]] > 0)
#pragma opm atomic
					b[deletePos]->c[z[j]]--;
			}
		}
	}
};

#endif
