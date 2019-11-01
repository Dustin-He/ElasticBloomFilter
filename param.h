#ifndef _PARAM_H
#define _PARAM_H

#define BYTE_SIZE 8
#define SHIFT 3					//移位
#define MASK 0x7				//掩码
#define KEY_LEN 4				//key的字节数
#define MAX_SIZE 32				//bits最长长度
#define MIN_SIZE 3				//bits最短长度
#define MAX_UINT32_T 0xFFFFFFFF	//无符号数最大值
#define MAX_UINT16_T 0xFFFF
#define BUCKET_SIZE 8			//桶内元素个数
#define SEARCH_RANGE 32			//找空桶的范围
#define EXPAND_THRESHOLD 0.2	//1的个数达到此比例则扩大

#define HASH_NUM 5				//哈希函数个数
#define BLOOM_SIZE 18			//bits长度

//#define DataNum 20000			//数据规模

#define LAYERS 1
#define BRANCH 3
#define BSHIFT 2

uint32_t total_slot;

#endif