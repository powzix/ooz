// stdafx.h : includefile for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#ifdef _WIN32
#include "targetver.h"
#endif

#define _CRT_SECURE_NO_WARNINGS 1

#include <stdio.h>
#ifdef _WIN32
#include <tchar.h>
#endif
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <intrin.h>
#include <assert.h>
#include <Windows.h>

#pragma warning (disable: 4244)
#else
#include <assert.h>
#include <emmintrin.h>
#include <stdint.h>
#include <sys/stat.h>
#include <time.h>
#endif

#ifndef __int64
#define __int64 long long
#endif

#ifndef LARGE_INTEGER
#define LARGE_INTEGER int64
#endif

// TODO: reference additional headers your program requires here
typedef unsigned char byte;
typedef unsigned char uint8;
typedef unsigned int uint32;
typedef unsigned __int64 uint64;
typedef signed __int64 int64;
typedef signed int int32;
typedef unsigned short uint16;
typedef signed short int16;
typedef unsigned int uint;

#ifndef __forceinline
#define __forceinline inline
#endif

#ifndef WINAPI
#ifdef _WIN32
#define WINAPI __stdcall
#else
#define WINAPI
#endif
#endif

#ifndef _BitScanReverse
inline uint32 _BitScanReverse( unsigned long* a, uint32 b ) {
//  return _BitScanReverse( (uint32*)a, b );
  a[0] = 31-__builtin_clz(b);
  return 0;
}
#endif

#ifndef _BitScanForward
inline uint32 _BitScanForward( unsigned long* a, uint32 b ) {
//  return _BitScanForward( (uint32*)a, b );
  a[0] = __builtin_ctz(b);
  return 0;
}
#endif

#ifndef _rotl
inline uint32 _rotl( uint32 x, int n ) {
  return (((x) << (n)) | ((x) >> (32-(n))));
}
#endif

#ifndef QueryPerformanceFrequency
inline uint32 QueryPerformanceFrequency( LARGE_INTEGER* lpFrequency ) {
  lpFrequency[0] = 1000000000LL;  // nanoseconds per second
  return 1;
}
#endif

#ifndef QueryPerformanceCounter
inline uint32 QueryPerformanceCounter( LARGE_INTEGER* lpPerformanceCount ) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  lpPerformanceCount[0] = (int64)ts.tv_sec * 1000000000LL + ts.tv_nsec;
  return 1;
}
#endif

#ifndef _byteswap_ushort
inline unsigned short _byteswap_ushort ( unsigned short val ) {
  return __builtin_bswap16(val);
}
#endif

#ifndef _byteswap_ulong
inline unsigned long _byteswap_ulong ( unsigned long val ) {
  return __builtin_bswap32(val);
}
#endif

#ifndef _byteswap_uint64
inline uint64 _byteswap_uint64 ( uint64 val ) {
  return __builtin_bswap64(val);
}
#endif
