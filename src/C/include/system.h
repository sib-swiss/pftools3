/*******************************************************
                       PFTOOLS
 *******************************************************
  Sep 19, 2013 system.h
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#ifndef SYSTEM_H_
#define SYSTEM_H_
#include <sched.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#if defined(USE_KNC) && !(defined(__MIC__))
#include "/usr/include/intel-coi/source/COIEngine_source.h"
#endif

#define MM_NONE         0x00000000 /* Not SSE at all */
#define MM_MMX          0x00000001 /* standard MMX */
#define MM_MMXEXT       0x00000002 /* SSE integer instructions or AMD MMX ext */
#define MM_3DNOW        0x00000004 /* AMD 3DNOW */
#define MM_SSE          0x00000008 /* SSE instructions */
#define MM_SSE2         0x00000010 /* PIV SSE2 instructions */
#define MM_3DNOWEXT     0x00000020 /* AMD 3DNowExt */
#define MM_SSE3         0x00000040 /* Prescott SSE3 instructions */
#define MM_SSSE3        0x00000080 /* Conroe SSSE3 instructions */
#define MM_SSE41        0x00000100 /* Penryn SSE41 instructions */
#define MM_SSE42        0x00000200 /* Nehalem SSE42 instructions */
#define MM_POPCOUNT     0x00000400
#define MM_SSE4A        0x00000800 /* AMD Barcelona instructions */
#define MM_AVX          0x00001000 /* AVX technology 256 bits */
#define MM_XOP          0x00002000 /* AMD XOP instructions */
#define MM_FMA3         0x00004000 /* Float Multiply Accumulate functions with 3 operands */
#define MM_FMA4         0x00008000 /* Float Multiply Accumulate functions with 4 operands */
#define MM_AVX2         0x00010000 /* AVX2 technology 256 bits */
#define MM_AVX512       0x00020000 /* Intel MIC architecture 512 bit AVX */
#define MM_AVX512F      0x00040000 /* Intel Xeon AVC 512 F */
#define MM_AVX512BW     0x00080000 /* Intel Xeon AVC 512 Byte Word */
#define MM_AVX512CD     0x00100000 /* Intel Xeon AVC 512 Conflict Detection */

#define TLB_2MB         0x00000001
#define TLB_1GB         0x00000002

#if !defined(__MIC__)
/* TOPOLOGY structure is as follow
 *
 *  2 bits for hyperthreading threads b0, b1 : b1 is hyperthreading only, b0 is non hyperthreading
 * 54 bits for the cores
 *  8 bits for the sockets
 * ------------------------
 * 64 bits
 */
#define SOCKET_MASK_WIDTH                               8
#define SOCKET_MASK                    0xFF00000000000000
#define SOCKET_MASK_SHIFT                              56
#define CORE_MASK                      0x00FFFFFFFFFFFFFC
#define CORE_MASK_WIDTH                                54
#define HYPERTHREADING_MASK            0x0000000000000003
#define HYPERTHREADING_MASK_WIDTH                       2
typedef uint64_t topology_mask_t;

#else

/* TOPOLOGY structure is as follow
 *
 * 4  bits for hyperthreading threads number 0, 1, 2, 3
 * 60 bits for the cores
 * 0 bits for the sockets
 * ------------------------
 * 64 bits
 */
#define SOCKET_MASK_WIDTH 0
#define SOCKET_MASK 0x0
#define SOCKET_MASK_SHIFT 64
#define CORE_MASK   0xFFFFFFFFFFFFFFF0
#define CORE_MASK_WIDTH 60
#define HYPERTHREADING_MASK 0xF
#define HYPERTHREADING_MASK_WIDTH 4
typedef uint64_t topology_mask_t;

#endif

enum CacheType { Data=1, Instruction=2, Unified=3 };
enum CacheLevel { OneI=0, OneD=1, Two=2, Three=3 };

typedef struct Affinity_Mask { cpu_set_t data; } Affinity_Mask_t;

typedef struct CacheInfo {
  unsigned int Size;
  unsigned int Set;
  unsigned short int Line;
  unsigned short int Partition;
  unsigned short int Ways;
  unsigned char Level;
  enum CacheType Type;
} cache_t;


typedef struct SystemInfo {
    /* Operating system data */
    char Username[64];
    char Release[64];
    char Architecture[64];
    char Nodename[64];

    /* Architecture data */
    char CPU_Vendor[13];
    char CPU_Name[48];
    unsigned int Family;
    unsigned int Model;
    unsigned int Extensions;
    unsigned int TLBs;

    /* Cache */
    cache_t Cache[4];

    /* Topology */
    unsigned int nSockets;
    unsigned int nCores;
    unsigned int nOverallCores;
    _Bool HyperthreadingAvailable;
    _Bool HyperthreadingOn;
    unsigned int nComputeUnit;
    unsigned int nCorePerComputeUnit;

    unsigned int SocketSelectMask;
    unsigned int CoreSelectMask;
    unsigned int HyperthreadingMask;
    unsigned short int SocketSelectMaskShift;
    unsigned short int HyperthreadingMaskWidth;

    /* NUMA */
#ifdef __NUMA__
    unsigned int nNodes;
    unsigned int nCpusPerNode;
    _Bool NumaAble;
#endif /* __NUMA__ */

    /* KNC */
    unsigned int nKNCs;
#ifdef USE_KNC
#if defined(USE_AFFINITY) && !(defined(__MIC__))
    COI_ENGINE_INFO * MicEngines;
#endif
#endif

    topology_mask_t * TopologyMasks;
    unsigned int * APICID;

} __attribute__((aligned(16))) SystemInfo;

void printTopology(const SystemInfo * const info);
void printAPIC(const SystemInfo * const info);
void printSystemInfo(FILE * const fd, const SystemInfo * const info);
void getSystemInfo(SystemInfo * const info);

_Bool writeSystemInfo(FILE * const fd, const SystemInfo * const info);
_Bool readSystemInfo(FILE * const fd, SystemInfo * const info);

#ifdef USE_AFFINITY
void getCumulativeMask(const SystemInfo * const info, const topology_mask_t SocketId, const topology_mask_t CoreId,
                       const topology_mask_t ThreadId, Affinity_Mask_t * const Mask);
size_t getMasks(const SystemInfo * const info, const topology_mask_t SocketId, const topology_mask_t CoreId,
                const topology_mask_t ThreadId, Affinity_Mask_t * * const Masks);
size_t getMasksFromParent(const SystemInfo * const info, Affinity_Mask_t * * const Masks);
void printMasks(const Affinity_Mask_t * const Masks, const size_t N);
#endif

inline void __attribute__((always_inline)) freeSystemInfo(SystemInfo * const info) {
    if (info->APICID) {
        free(info->APICID);
    }
#ifdef USE_KNC
#if defined(USE_AFFINITY) && !(defined(__MIC__))
    if (info->MicEngines) free(info->MicEngines);
#endif
#endif
    /* WARNING: DO NOT FREE TOPOLOGY HAS ITS MEMORY IS LINKED TO APICID
    if (info->TopologyMasks) {
        free(info->TopologyMasks);
    }
    */
    char * ptr = (char*) info;
    for (size_t i=0; i < sizeof(SystemInfo); ++i) ptr[i] = '\0';
}

#endif /* SYSTEM_H_ */
