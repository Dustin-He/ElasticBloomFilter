#ifndef _SBF_H
#define _SBF_H

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "BOBHash32.h"
#include "param.h"
#include "murmur3.h"

using namespace std;

class BF
{
public:
	char *c;
	int num;
	BF(){};
	BF(int m)
	{
		num = 0;
		c = new char[1 << m];
		memset(c, 0, 1 << m);
	}
};
class ScalableBF
{
public:
	BOBHash32 *bobhash[HASH_NUM + 2];
	int w, m;
	double c;
	BF *b[10100];
	ScalableBF(){};
	ScalableBF(int _m, double _c)
	{
		w = 0;
		m = _m;
		c = _c;
		memset(b, 0, sizeof(b));
		b[0] = new BF(m);
		for (int i = 0; i < HASH_NUM; i++)
			bobhash[i] = new BOBHash32(i + 1005);
	}
	void insert(int x)
	{
		if (b[w]->num >= c * (1 << (m + w)))
		{
			w++;
			b[w] = new BF(m + w);
			// c *= 0.9;
		}
		for (int i = 0; i < HASH_NUM; i++)
		{
			// int z = bobhash[i]->run((char *)&x, 4) & ((1 << (m + w)) - 1);
			int z = MurmurHash3_x86_32((char *)&x, 4, i) & ((1 << (m + w)) - 1);
			if (!((b[w]->c[z>>3]>>(z&7))&1))
			{
				b[w]->num++;
				b[w]->c[z>>3]|=1<<(z&7);
			}
		}
	}
	void insert(char *x)
	{
		if (b[w]->num >= c * (1 << (m + w)))
		{
			w++;
			b[w] = new BF(m + w);
			// c *= 0.9;
		}
		for (int i = 0; i < HASH_NUM; i++)
		{
			// int z = bobhash[i]->run((char *)x, 13) & ((1 << (m + w)) - 1);
			int z = MurmurHash3_x86_32((char *)x, 13, i) & ((1 << (m + w)) - 1);
			if (!((b[w]->c[z>>3]>>(z&7))&1))
			{
				b[w]->num++;
				b[w]->c[z>>3]|=1<<(z&7);
			}
		}
	}
	bool query(int x)
	{
		int z[HASH_NUM];
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)&x, 4, j) & ((1 << (m + w)) - 1);
			// z[j] = bobhash[j]->run((char *)&x, 4) & ((1 << (m + w)) - 1);
		for (int i = w; i >= 0; i--)
		{
			bool f = 1;
			for (int j = 0; j < HASH_NUM; j++)
			{
				if (!((b[i]->c[(z[j] & ((1 << (m + i)) - 1))>>3]>>(z[j]&7))&1) )
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
		
		return false;
	}
	bool query(char *x)
	{
		int z[HASH_NUM];
		for (int j = 0; j < HASH_NUM; j++)
			z[j] = MurmurHash3_x86_32((char *)x, 13, j) & ((1 << (m + w)) - 1);
			// z[j] = bobhash[j]->run((char *)x, 13) & ((1 << (m + w)) - 1);
		for (int i = w; i >= 0; i--)
		{
			bool f = 1;
			for (int j = 0; j < HASH_NUM; j++)
			{
				if (!((b[i]->c[(z[j] & ((1 << (m + i)) - 1))>>3]>>(z[j]&7))&1) )
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
		
		return false;
	}
};

#endif
