#ifndef _PARAM_H
#define _PARAM_H

#include <stdint.h>

#define BYTE_SIZE 8
#define SHIFT 3					//a shift
#define MASK 0x7				//a mask
#define KEY_LEN 4				//length of the key(bytes)
#define MAX_SIZE 32				//maximum length of the bits array(2^MAX_SIZE)
#define MIN_SIZE 3				//minimum length of the bits array(2^MIN_SIZE)
#define MAX_UINT32_T 0xFFFFFFFF	//maximum value of type uint32_t
#define MAX_UINT16_T 0xFFFF
#define BUCKET_SIZE 8			//number of elements in one bucket
#define SEARCH_RANGE 32			//range for serching a empty bucket
#define EXPAND_THRESHOLD 0.2	//expand when the _1_rate reaches it

#define HASH_NUM 5				//number of hash functions
#define BLOOM_SIZE 20			//length of the bits array
#define CBF_SIZE 22				//length of the CBF

#define LAYERS 2				//the number of layers of a bucket tree
#define BRANCH 3				//the number of branches of a bucket tree
#define BSHIFT 2				//a shift for calculating the position in a bucket tree

uint32_t total_slot;			//total number of element slots
bool expandOrNot;				//whether to expand

#endif