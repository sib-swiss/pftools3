/*******************************************************
                       FastEpistasis
 *******************************************************
  Sep 19, 2013 system.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#define _GNU_SOURCE
#include "config.h"
#include "system.h"
#ifdef USE_AFFINITY
#ifndef __USE_GNU
# error "GNU EXTENSION not activated."
#endif
#endif

#include <stdio.h>
#include <string.h>
#include <sys/utsname.h>
#ifndef NDEBUG
#include <sys/stat.h>
#endif
#include <sys/types.h>
#include <unistd.h>
#include <dlfcn.h>
#ifndef NO_USERNAME
#  include <pwd.h>
#endif

#ifdef __NUMA__
#  include <numa.h>
#endif

#if defined(USE_KNC) && !defined(__MIC__)
#include "/usr/include/intel-coi/common/COIEngine_common.h"
#endif

#ifdef _SYSTEM_TEST
#define __SYSTEM_DEBUG__
#endif

#ifndef __ALWAYS_INLINE
#define __ALWAYS_INLINE __attribute__((always_inline))
#endif

#define LINE_SIZE 70
#define LINE_BUFFER_SIZE 128
#define TCHAR(text) text
#define SINGLE_VERTICAL_BAR "|"
#define DOUBLE_LINE "===================================================================="
#define SINGLE_LINE "--------------------------------------------------------------------"
#define EMPTY_SPACE "                                                                    "
#define TOP_DOUBLE_LINE     "#" DOUBLE_LINE "#\n"
#define MIDDLE_DOUBLE_LINE  "#" DOUBLE_LINE "#\n"
#define MIDDLE_SINGLE_LINE  "+" SINGLE_LINE "+\n"
#define BOTTOM_DOUBLE_LINE  "#" DOUBLE_LINE "#\n"
#define BOTTOM_SINGLE_LINE  "#" SINGLE_LINE "#\n"
#define EMPTY_LINE          "|" EMPTY_SPACE  "|\n"

#define SGL_TEXT_LINE(text) SINGLE_VERTICAL_BAR text SINGLE_VERTICAL_BAR "\n"
#define PRINT_SGL_STRING_LINE(format, format_size, string, string_size) {\
    fprintf(fd,SINGLE_VERTICAL_BAR format, string);\
    for (int i=1+format_size+string_size; i<LINE_SIZE-1; i++) { fputc(' ', fd);}\
		fputs(SINGLE_VERTICAL_BAR "\n", fd);\
}
#define PRINT_SGL_TEXT_LINE(string, string_size) {\
    fputs(SINGLE_VERTICAL_BAR string, fd);\
    for (int i=1+string_size; i<LINE_SIZE-1; i++) { fputc(' ', fd);}\
		fputs(SINGLE_VERTICAL_BAR "\n", fd);\
}



extern inline unsigned long __ALWAYS_INLINE createMask(unsigned int numEntries, unsigned int *maskWidth)
{
    unsigned int i;
    unsigned long k = ((unsigned long) numEntries) * 2L - 1L;
    unsigned long index;
    unsigned int AllZero;

#ifndef __MIC__    
    __asm__ __volatile__ (
      "xorl   %0, %0; \n\t"
      "bsrq   %2,%1;  \n\t"
      "cmovzl %3,%0;  \n\t"
      : "=r"(AllZero), "=r"(index)
      : "r"(k), "rm"(1)
    );
#else
    __asm__ __volatile__ (
      "xorl   %0, %0; \n\t"
      "bsrq   %2,%1;  \n\t"
      "jnz 1f ;       \n\t"
      "movl   %3,%0;  \n\t"
      "1:"
      : "=r"(AllZero), "=r"(index)
      : "r"(k), "rm"(1)
    ); 
#endif
    if (AllZero) {
      *maskWidth = 0; 
      return 0L;
    }
    *maskWidth = (unsigned int) index;
    if (index == 63L) {
      return -1L;
    } else {
      return (1L << index) - 1L;
    }
}

extern inline void __ALWAYS_INLINE GetOSInfo(SystemInfo * const info)
{
  struct utsname uts_name;
  uname(&uts_name);
  size_t tmp_size = strlen(uts_name.release);
  if (tmp_size<63) {
    strncpy(info->Release, uts_name.release, 63);
  } else {
    memcpy(info->Release, uts_name.release, 60);
    info->Release[60] = '.';
    info->Release[61] = '.';
    info->Release[62] = '.';
    info->Nodename[63] = '\0';
  }
#ifndef NO_USERNAME
  uid_t uid = getuid();
  struct passwd *pass= getpwuid(uid);
  if ( pass != NULL) {
    if (strlen(pass->pw_gecos) > 0) {
      strncpy(info->Username,pass->pw_gecos,63);
    } else {
      strncpy(info->Username,pass->pw_name,63);
    }
  }
#else
  const char text[] = "Static linking prevents safe username query";
  strcpy(info->Username, text);
#endif
  
  tmp_size = strlen(uts_name.nodename);
  if (tmp_size<63) {
    strncpy(info->Nodename, uts_name.nodename, 63);
  } else {
    memcpy(info->Nodename, uts_name.nodename, 60);
    info->Nodename[60] = '.';
    info->Nodename[61] = '.';
    info->Nodename[62] = '.';
    info->Nodename[63] = '\0';
  }  
    
  tmp_size = strlen(uts_name.machine);
  if (tmp_size<63) {
    strncpy(info->Architecture, uts_name.machine, 63);
  } else {
    memcpy(info->Architecture, uts_name.machine, 60);
    info->Architecture[60] = '.';
    info->Architecture[61] = '.';
    info->Architecture[62] = '.';
    info->Nodename[63] = '\0';
  }
}

#ifdef USE_AFFINITY
void getCumulativeMask(const SystemInfo * const info, const topology_mask_t SocketId, const topology_mask_t CoreId,
                       const topology_mask_t ThreadId, Affinity_Mask_t * const Mask)
{
	size_t count = 0;
	topology_mask_t ValidMask;
	
	ValidMask  = ThreadId & HYPERTHREADING_MASK;  
	ValidMask |= (CoreId << HYPERTHREADING_MASK_WIDTH) & CORE_MASK ;
#if SOCKET_MASK_SHIFT < 64
	ValidMask |= (SocketId << SOCKET_MASK_SHIFT) & SOCKET_MASK;
#endif

	count = 0;
	CPU_ZERO_S(sizeof(Affinity_Mask_t), Mask);
	
	for (unsigned int thread=0; thread<info->nOverallCores; ++thread) {
		register const topology_mask_t tmp = info->TopologyMasks[thread] & ValidMask;
#if SOCKET_MASK_SHIFT < 64
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK) && (tmp & SOCKET_MASK))
#else
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK))
#endif
		{
			CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) Mask);
		}
	}
}

size_t getMasks(const SystemInfo * const info, const topology_mask_t SocketId, const topology_mask_t CoreId,
                const topology_mask_t ThreadId, Affinity_Mask_t * * const Masks)
{
	size_t count = 0;
	topology_mask_t ValidMask;
	
	if (SocketId == 0 || CoreId == 0 || ThreadId == 0) return 0;
	ValidMask  = ThreadId & HYPERTHREADING_MASK;  
	ValidMask |= (CoreId << HYPERTHREADING_MASK_WIDTH) & CORE_MASK ;
#if SOCKET_MASK_SHIFT < 64
	ValidMask |= (SocketId << SOCKET_MASK_SHIFT) & SOCKET_MASK;
#endif
#ifdef __SYSTEM_DEBUG__
	printf("Topology mask size is %lu\nThread mask is 0x%8.8x\nCore mask shift is %i\n"\
	       "Core mask is 0x%8.8x\nSocket mask shift is %i\nSocket mask is 0x%8.8x\n",
	       sizeof(topology_mask_t)*8,
	       info->HyperthreadingMask, info->HyperthreadingMaskWidth,
	       info->CoreSelectMask, info->SocketSelectMaskShift, info->SocketSelectMask);
	
	printf("Valid Mask  : 0x%16.16lx\n             : ", ValidMask);
	for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
		if (ValidMask & ((topology_mask_t)1 << ((8*sizeof(topology_mask_t)-1) - j)) ) 
			fputs("1",stdout);
		else 
			fputs("0", stdout);
		if ( ((8*sizeof(topology_mask_t)-1 - j) == (CORE_MASK_WIDTH+HYPERTHREADING_MASK_WIDTH)) 
				|| ((8*sizeof(topology_mask_t)-1 - j) == (HYPERTHREADING_MASK_WIDTH))) fputs(" ",stdout);
	}
	fputs("\n",stdout);
#endif  
	for (unsigned int thread=0; thread<info->nOverallCores; ++thread) {
		const topology_mask_t tmp = info->TopologyMasks[thread] & ValidMask;
		#if SOCKET_MASK_SHIFT < 64
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK) && (tmp & SOCKET_MASK)) ++count;
		#else
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK)) ++count; 
		#endif
	}
	  
  if (count == 0) return 0;
	
	*Masks = (Affinity_Mask_t*) malloc(count*sizeof(Affinity_Mask_t));
	if (*Masks == NULL) return 0;
	count = 0;
	for (unsigned int thread=0; thread<info->nOverallCores; ++thread) {
		register const topology_mask_t tmp = info->TopologyMasks[thread] & ValidMask;
#if SOCKET_MASK_SHIFT < 64
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK) && (tmp & SOCKET_MASK))
#else
		if ((tmp & HYPERTHREADING_MASK) && (tmp & CORE_MASK))
#endif
		{
			CPU_ZERO_S(sizeof(Affinity_Mask_t), *Masks + count);
			CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) (*Masks + count));
#ifdef __SYSTEM_DEBUG__
			printf("Thread %3u   : ", thread);
			for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
				if (info->TopologyMasks[thread] & (((topology_mask_t) 1) << ((8*sizeof(topology_mask_t)-1) - j)) ) 
					fputs("1",stdout);
				else 
					fputs("0", stdout);
				if ( ((8*sizeof(topology_mask_t)-1 - j) == (CORE_MASK_WIDTH+HYPERTHREADING_MASK_WIDTH)) 
					|| ((8*sizeof(topology_mask_t)-1 - j) == (HYPERTHREADING_MASK_WIDTH))) fputs(" ",stdout);
			}
      fputs("\n",stdout);
#endif
			++count;
		}
#ifdef __SYSTEM_DEBUG__
    else {
			printf("Thread %3u   : ", thread);
			for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
				if (info->TopologyMasks[thread] & (((topology_mask_t) 1) << ((8*sizeof(topology_mask_t)-1) - j)) ) 
					fputs("1",stdout);
				else 
					fputs("0", stdout);
				if ( ((8*sizeof(topology_mask_t)-1 - j) == (CORE_MASK_WIDTH+HYPERTHREADING_MASK_WIDTH)) 
					|| ((8*sizeof(topology_mask_t)-1 - j) == (HYPERTHREADING_MASK_WIDTH))) fputs(" ",stdout);
			}
      fputs(" FAILED\n",stdout);
		}
#endif
	}
	return count;
}

size_t getMasksFromParent(const SystemInfo * const info, Affinity_Mask_t * * const Masks)
{
  size_t count = 0;
  Affinity_Mask_t * const restrict CurrentMask = (Affinity_Mask_t*) malloc(sizeof(Affinity_Mask_t));
  if (CurrentMask == NULL) goto FIN;
  if (sched_getaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) CurrentMask) != 0) goto FIN1;
  
  /* Count potential cores */
  for (int i=0; i<sizeof(Affinity_Mask_t); i++) {
    if (CPU_ISSET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) CurrentMask)) ++count;
  }
  if (count == 0) return 0;
  
  *Masks = (Affinity_Mask_t*) malloc(count*sizeof(Affinity_Mask_t));
  if (*Masks == NULL) return 0;
  
  Affinity_Mask_t * restrict newMask = *Masks;
  for (int i=0; i<sizeof(Affinity_Mask_t); i++) {
    if (CPU_ISSET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) CurrentMask)) {
      CPU_ZERO_S(sizeof(Affinity_Mask_t), newMask);
      CPU_SET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) newMask);
      newMask++;
    }
  }
  
  
  free(CurrentMask);
  return count;
  
  FIN1:
    free(CurrentMask);
  FIN:
    return 0;
}
#endif

static void GetIntelTopology(SystemInfo * const info)
{
	unsigned int APIC_HyperthreadingMask, APIC_HyperthreadingMaskWidth;
	unsigned int APIC_SocketSelectMask, APIC_SocketSelectMaskShift;
	unsigned int APIC_CoreSelectMask;
	
	int eax_value, ebx_value, ecx_value, edx_value;
	
	/* Do we support Leaf 11 (0xB) */
	_Bool UseLeafB = false;
	__asm__ __volatile__ (
		"xorl %%eax, %%eax;"
		"xorl %%ecx, %%ecx;"
		"cpuid ;"
		: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	);
	const unsigned int MaxCPUID = eax_value;
	if (eax_value >= 0xB) {
		eax_value = 0xB;
		__asm__ __volatile__ (
			"xorl %%ecx, %%ecx;"
			"cpuid ;"
			: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
			: "0"(eax_value), "2"(ecx_value)
		);
		UseLeafB = (ebx_value != 0); 
	} 
  /* Use HWMT hyperthreading feature flag to treat different configurations */
	eax_value = 1;
	__asm__ __volatile__ (
		"xorl %%ecx, %%ecx;"
		"cpuid ;"
		: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
		: "0"(eax_value)
	);
		
	if (edx_value & (1<<28)) {
		int wasCoreReported = 0;
		int wasThreadReported = 0;
		unsigned int ThreadPerCore = 1;
		unsigned int MaxCore = 1;
		
		if (UseLeafB) {
			int subLeaf = 0, levelType, levelShift;
			unsigned int coreplusSMT_Mask;
			do {
				eax_value = 0xB;
				ecx_value = subLeaf;
				__asm__ __volatile__ (
					"cpuid ;"
					: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
					: "0"(eax_value), "2"(ecx_value)
				);
				if (ebx_value == 0) break;
				
				levelType  = (ecx_value >> 8) & 0xFF;
				levelShift = (eax_value & 0xF);
				switch (levelType)
				{
					case 1:
						//level type is SMT, so levelShift is the SMT_Mask_Width
						APIC_HyperthreadingMask = ~((-1) << levelShift);
						APIC_HyperthreadingMaskWidth = levelShift;
						wasThreadReported = 1;
						break;
					case 2: //level type is Core, so levelShift is the CorePlsuSMT_Mask_Width
						coreplusSMT_Mask = ~((-1) << levelShift);
						APIC_SocketSelectMaskShift = levelShift;
						APIC_SocketSelectMask = (-1) ^ coreplusSMT_Mask;
						wasCoreReported = 1;
						break;
					default:
						// handle in the future
						break;
				}
					++subLeaf;
			} while (1);
			
			if (wasThreadReported && wasCoreReported)
			{
				APIC_CoreSelectMask = coreplusSMT_Mask ^ APIC_HyperthreadingMask;
			}
			else if (!wasCoreReported && wasThreadReported)
			{
				APIC_CoreSelectMask = 0;
				APIC_SocketSelectMaskShift = APIC_HyperthreadingMaskWidth;
				APIC_SocketSelectMask = (-1) ^ APIC_HyperthreadingMask;
			}
			else //(case where !wasThreadReported)
			{
				// throw an error, this should not happen if hardware function normally
				fputs("Error in hardware info extraction\n", stderr);
				exit(1);
			}
			info->HyperthreadingAvailable = true;
			info->HyperthreadingOn = (wasThreadReported > 0); 
		}
		else {
			const unsigned int MaxCorePlusThread = (ebx_value >> 16) & 0xFF;
			if (MaxCPUID >= 4) {
				eax_value = 4;
				__asm__ __volatile__ (
					"xorl %%ecx, %%ecx;"
					"cpuid ;"
					: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
					: "0"(eax_value)
				);
				MaxCore = 1 + ( ( eax_value >> 26 ) & 0x3F );
				ThreadPerCore = MaxCorePlusThread / MaxCore;
				info->HyperthreadingAvailable = true;
				info->HyperthreadingOn = (MaxCorePlusThread > MaxCore ) ? true : false;
			} 
			else {
				MaxCore = 1;
				ThreadPerCore = MaxCorePlusThread;
				
				/* Check whether BIOS is preventing Hyperthreading */
				eax_value = 0x80000000;
				__asm__ __volatile__ (
					"xorl %%ecx, %%ecx;"
					"cpuid ;"
					: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
					: "0"(eax_value)
				);
				if ( MaxCPUID <= 4 && eax_value > 0x80000004) {
					info->HyperthreadingAvailable = true;
				} 
			}
      
		  /* Create the masks */
			APIC_HyperthreadingMask     = createMask(ThreadPerCore, &(APIC_HyperthreadingMaskWidth));
			APIC_CoreSelectMask         = createMask(MaxCore, &(APIC_SocketSelectMaskShift));
			APIC_SocketSelectMaskShift += APIC_HyperthreadingMaskWidth;
			APIC_CoreSelectMask       <<= APIC_HyperthreadingMaskWidth;
			APIC_SocketSelectMask       = (-1) ^ (APIC_CoreSelectMask | APIC_HyperthreadingMask); 
		}

		/* Get Cache topology */
		if (MaxCPUID >= 4) {
			int subLeaf = 0;
			do {
				eax_value = 0x4;
				ecx_value = subLeaf;
				__asm__ __volatile__ (
					"cpuid ;"
					: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
					: "0"(eax_value), "2"(ecx_value)
				);
				if ((eax_value & 0x1F) == 0) break;
				
				unsigned int cachetype = (eax_value & 0x1F);
				unsigned int cachelevel = (eax_value >> 5) & 0x7;
				unsigned int lineSize = 1 + (ebx_value & 0xFFF);
				unsigned int physicalLinePartition = 1+ ((ebx_value >> 12) & 0x3FF);
				unsigned int nsets = 1 + ecx_value;
				unsigned int ways = 1 + ((ebx_value >> 22) & 0x3FF);
				unsigned int cachesize = ways * lineSize * physicalLinePartition * nsets;
#ifdef __SYSTEM_DEBUG__
				printf("Cache Level %u - type %u - %u kB - ways %u - Line size %u - partition %u - set %u\n",
							 cachelevel, cachetype, cachesize >> 10, ways, lineSize, physicalLinePartition, nsets);
				if (eax_value & 0x100) printf("Previous is fully associative\n");
#endif
				int id;
				if (cachelevel == 1) {
					id = (cachetype == 1) ? 1 : 0;
				}
				else {
					id = cachelevel;
				}
				{
					cache_t * const restrict CachePtr = &(info->Cache[id]);  
					CachePtr->Size = cachesize;
					CachePtr->Set  = nsets;
					CachePtr->Line = (unsigned short int) lineSize;
					CachePtr->Partition = (unsigned short int) physicalLinePartition;
					CachePtr->Type = cachetype;
					CachePtr->Ways = (unsigned short int) ways;
					CachePtr->Level = cachelevel;
				}
				++subLeaf;
			} while (1);
		}
	}
	else {
		/* Prior to Hyperthreading Technology only one thread per core */
		APIC_SocketSelectMask = -1;
	}
    
#ifdef USE_AFFINITY
	/* Retrieve APIC data for each thread */
	Affinity_Mask_t Mask, BackupMask;
	if ( sched_getaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &BackupMask)) {
		fputs("Error getting affinity!\n",stderr);
		exit(1);
	}
						  
	info->APICID = (unsigned int*) malloc(info->nOverallCores*(sizeof(unsigned int) + sizeof(topology_mask_t)));
	info->TopologyMasks = (topology_mask_t*) &(info->APICID[info->nOverallCores]);

	if (UseLeafB) {
		for (unsigned int thread = 0; thread<info->nOverallCores; ++thread) {
			CPU_ZERO_S(sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask);
			CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask);
			if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask)) {
				fputs("Error setting affinity!\n",stderr);
				exit(1);
			}
			      __asm__ __volatile__ (
							"xorl %%ecx, %%ecx;"
							"cpuid;"
							: "=d"(edx_value)
							: "a"(0xB)
							: "memory", "ebx", "ecx"
						);
						info->APICID[thread] = edx_value;
		}
	} 
	else {
		for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
			CPU_ZERO_S(sizeof(Affinity_Mask_t), &Mask);
			CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask);
			if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask)) {
				fputs("Error setting affinity!\n",stderr);
				exit(1);
			}
			      __asm__ __volatile__ (
							"xorl %%ecx, %%ecx;"
							"cpuid;"
							: "=b"(ebx_value)
							: "a"(0x1)
							: "memory", "%ecx", "%edx"
						);
						unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
						info->APICID[thread] = tmp;
		}
	}

  if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &BackupMask)) {
		fputs("Error setting affinity!\n",stderr);
		exit(1);
	}
#ifdef __SYSTEM_DEBUG__          
	printf("HTT Mask (%2u)    : ",APIC_HyperthreadingMaskWidth);
	for (unsigned int j=0; j<32; ++j) {
		if ( APIC_HyperthreadingMask & (1 << (31 - j)) ) 
			fputs("1",stdout);
		else 
			fputs("0", stdout);
	}
	unsigned int CoreMaskWidth = 0;
	for (unsigned int j=0; j<32; ++j) {
		if ( APIC_CoreSelectMask & (1 << (31 - j)) ) ++CoreMaskWidth;
	}
  	printf("\nCore Mask (%2u)   : ",CoreMaskWidth);
	for (unsigned int j=0; j<32; ++j) {
		if ( APIC_CoreSelectMask & (1 << (31 - j)) ) 
			fputs("1",stdout);
		else 
			fputs("0", stdout);
	}
	printf("\nSocket Mask (%2u) : ", APIC_SocketSelectMaskShift);
	for (unsigned int j=0; j<32; ++j) {
		if ( APIC_SocketSelectMask & (1 << (31 - j)) ) 
			fputs("1",stdout);
		else 
			fputs("0", stdout);
	}
	fputs("\n",stdout);
#endif
	
	register unsigned int nSockets = 0;

	/* Set the topology structure to use based upon current ARCHITECTURE 
	 * that is:
	 *          T thread bits (minimum is 1)
	 *          C core bits   (minimum is 1)
	 *          S socket bits is the rest (32 - the above)
	 */
	info->HyperthreadingMaskWidth = APIC_HyperthreadingMaskWidth;
	info->HyperthreadingMask      = APIC_HyperthreadingMask;
	info->CoreSelectMask          = APIC_CoreSelectMask;
	info->SocketSelectMask        = APIC_SocketSelectMask;
	info->SocketSelectMaskShift   = APIC_SocketSelectMaskShift;
	
	unsigned char * const restrict CorePresent = alloca(info->nOverallCores*sizeof(unsigned char));
	memset(CorePresent, 0, info->nOverallCores*sizeof(unsigned char));
	
	for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
		const register unsigned int ThreadAPIC = info->APICID[thread];
		topology_mask_t TMask = ((topology_mask_t) 1) << (ThreadAPIC & APIC_HyperthreadingMask);
		unsigned int tmp = (ThreadAPIC & APIC_CoreSelectMask) >> APIC_HyperthreadingMaskWidth;
		CorePresent[tmp] = 1;
		TMask |= ((topology_mask_t) 1) << (HYPERTHREADING_MASK_WIDTH + tmp);
#if SOCKET_MASK_SHIFT < 64  
		tmp = (ThreadAPIC & APIC_SocketSelectMask) >> APIC_SocketSelectMaskShift;
		if (tmp > nSockets) nSockets = tmp;
		TMask |= ((topology_mask_t) 1) << ( SOCKET_MASK_SHIFT + tmp);
#endif
		info->TopologyMasks[thread] = TMask;
	}
	unsigned int nCores = 0;
	for (unsigned int i=0; i<info->nOverallCores; i++) if(CorePresent[i]) nCores++;
	info->nCores = nCores;
	info->nSockets = 1 + nSockets;
#endif
}
    
static void GetAMDTopology(SystemInfo * const info)
{
  unsigned int APIC_HyperthreadingMask, APIC_HyperthreadingMaskWidth;
  unsigned int APIC_SocketSelectMask, APIC_SocketSelectMaskShift;
  unsigned int APIC_CoreSelectMask;
    
  int eax_value, ebx_value, ecx_value, edx_value;
  
  /* Use HTT hyperthreading feature flag to test single or multicore architecture */
  eax_value = 1;
  __asm__ __volatile__ (
	  "xorl %%ecx, %%ecx;"
	  "cpuid ;"
	  : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	  : "0"(eax_value)
	  );
  
  if (edx_value & (1<<28)) {
    /* Get the number of logical core per processor */
    info->nCores = 1 + ((ebx_value >> 16) & 0xFF);
  } 
  
  /*Get the largest extension leaf supported */
  eax_value = 0x80000000;
  int MaxSupportedExtLeaf;
  __asm__ __volatile__ (
	"xorl %%ecx, %%ecx;"
	"cpuid ;"
	: "=a"(MaxSupportedExtLeaf)
	: "0"(eax_value)
	: "%ebx", "%ecx", "%edx"
  );
  
  /* Bulldozer */
  const _Bool IsBulldozer = eax_value >= 0x8000001E ? true : false;
  
  /* Get the number of physical core per processor */
  eax_value = 0x80000008;
  __asm__ __volatile__ (
	"xorl %%ecx, %%ecx;"
	"cpuid ;"
	: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	: "0"(eax_value)
  );
  
  info->nCores = 1 + (ecx_value & 0xFF );
  unsigned int coreplusSMT_MaskWidth = (ecx_value >> 12) & 0xF;
  
  /* Get Cache topology */
  if (MaxSupportedExtLeaf >= 0x80000006 ) {
    /* L1 cache */
    eax_value = 0x80000005;
    __asm__ __volatile__ (
	"xorl %%ecx, %%ecx;"
	"cpuid ;"
	: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	: "0"(eax_value)
    );
    
    /* Data */
    info->Cache[OneD].Level = 1;
    info->Cache[OneD].Type  = Data;
    info->Cache[OneD].Line  = ecx_value & 0xFF;
    info->Cache[OneD].Set   = (ecx_value >> 8) & 0xFF;
    info->Cache[OneD].Ways  = 1 << ((ecx_value >> 17) & 0xFF);
    info->Cache[OneD].Size  = ((ecx_value >> 24) & 0xFF) << 10;
  
    /* Instructions */
    info->Cache[OneI].Level = 1;
    info->Cache[OneI].Type  = Instruction;
    info->Cache[OneI].Line  = edx_value & 0xFF;
    info->Cache[OneI].Set   = (edx_value >> 8) & 0xFF;
    info->Cache[OneI].Ways  = 1 << ((edx_value >> 17) & 0xFF);
    info->Cache[OneI].Size  = ((edx_value >> 24) & 0xFF) << 10;
    
    /* L2 */
    eax_value = 0x80000006;
    __asm__ __volatile__ (
	"xorl %%ecx, %%ecx;"
	"cpuid ;"
	: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	: "0"(eax_value)
    );
    info->Cache[Two].Level = 2;
    info->Cache[Two].Type  = Data;
    info->Cache[Two].Line  = ecx_value & 0xFF;
    info->Cache[Two].Set   = (ecx_value >> 8) & 0xF;
    info->Cache[Two].Ways  = 1 << ((ecx_value >> 13) & 0x7);
    info->Cache[Two].Size  = ((ecx_value >> 16) & 0xFFFF) << 10;
    
    /* L3 */
    info->Cache[Three].Level = 3;
    info->Cache[Three].Type  = Data;
    info->Cache[Three].Line  = edx_value & 0xFF;
    info->Cache[Three].Set   = (edx_value >> 8) & 0xF;
    info->Cache[Three].Ways  = 1 << ((edx_value >> 13) & 0x7);
    info->Cache[Three].Size  = ((edx_value >> 18)) << 19;
    
  }
  
#ifdef USE_AFFINITY
  /* Retrieve APIC data for each thread */
  Affinity_Mask_t Mask, BackupMask;
  if ( sched_getaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &BackupMask)) {
      fputs("Error getting affinity!\n",stderr);
      exit(1);
  }
  
  info->APICID = (unsigned int*) malloc(info->nOverallCores*(sizeof(unsigned int) + sizeof(topology_mask_t)));
  info->TopologyMasks = (topology_mask_t*) &(info->APICID[info->nOverallCores]);
  
  if (!IsBulldozer) {
    for (unsigned int thread=0U; thread<info->nOverallCores; ++thread) {
      CPU_ZERO_S(sizeof(Affinity_Mask_t), &Mask);
      CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask);
      if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask)) {
	  fputs("Error setting affinity!\n",stderr);
	  exit(1);
      }
      __asm__ __volatile__ (
	  "xorl %%ecx, %%ecx;"
	  "cpuid;"
	  : "=b"(ebx_value)
	  : "a"(0x1)
	  : "%ecx", "%edx"
	);
      unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
      info->APICID[thread] = tmp;
    }
  } 
  else {
    /* Get number of core per compute unit assuming all are identical */
    __asm__ __volatile__ (
	  "xorl %%ecx, %%ecx;"
	  "cpuid;"
	  : "=b"(ebx_value)
	  : "a"(0x8000001E)
	  : "%ecx", "%edx"
	);
    info->nCorePerComputeUnit = 1 + ( ( ebx_value >> 8) & 0xFF );
    unsigned int nCU = 0U;
    for (unsigned int thread=0U; thread<info->nOverallCores; ++thread) {
      CPU_ZERO_S(sizeof(Affinity_Mask_t), &Mask);
      CPU_SET_S(thread, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask);
      if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask)) {
	  fputs("Error setting affinity!\n",stderr);
	  exit(1);
      }
      
      __asm__ __volatile__ (
	  "xorl %%ecx, %%ecx;"
	  "cpuid;"
	  : "=b"(ebx_value)
	  : "a"(0x1)
	  : "%ecx", "%edx"
	);
      unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
      info->APICID[thread] = tmp;
#ifdef __SYSTEM_DEBUG__	  
      fprintf(stderr,"APIC DATA (%2u)   : ", thread);
      for (unsigned int j=0; j<8; ++j) {
	if ( tmp & (1 << (7 - j)) ) 
	  fputs("1",stderr);
	else 
	  fputs("0", stderr);
      }
#endif
      __asm__ __volatile__ (
	  "xorl %%ecx, %%ecx;"
	  "cpuid;"
	  : "=b"(ebx_value)
	  : "a"(0x8000001E)
	  : "%ecx", "%edx"
	);
      
      tmp = ebx_value & 0xFF;
      if (tmp > nCU) nCU = tmp;
#ifdef __SYSTEM_DEBUG__		  
      fputs(" CU ",stderr);
      for (unsigned int j=0; j<8; ++j) {
	if ( tmp & (1 << (7 - j)) ) 
	  fputs("1",stderr);
	else 
	  fputs("0", stderr);
      }
      fputs("\n", stderr);
#endif
    }
    info->nComputeUnit = 1 + nCU;
  }
	
  /* Create the masks */
  if ( info->nCorePerComputeUnit ) {
    APIC_HyperthreadingMask = createMask(info->nCorePerComputeUnit, &(APIC_HyperthreadingMaskWidth));
    coreplusSMT_MaskWidth -= APIC_HyperthreadingMaskWidth; 
  }
  else {
    APIC_HyperthreadingMask = 0;
  }
  APIC_SocketSelectMaskShift = coreplusSMT_MaskWidth;
  APIC_CoreSelectMask = ( 1 << coreplusSMT_MaskWidth ) - 1; //createMask(coreplusSMT_MaskWidth, &(APIC_SocketSelectMaskShift));
  APIC_SocketSelectMaskShift += APIC_HyperthreadingMaskWidth;
  APIC_CoreSelectMask <<= APIC_HyperthreadingMaskWidth;
  APIC_SocketSelectMask = (-1) ^ (APIC_CoreSelectMask | APIC_HyperthreadingMask);  
    
  if ( sched_setaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &BackupMask)) {
      fputs("Error setting affinity!\n",stderr);
      exit(1);
  }
#ifdef __SYSTEM_DEBUG__     
  printf("HTT Mask (%2u)    : ", 0);
  for (unsigned int j=0; j<32; ++j) {
    if ( APIC_HyperthreadingMask & (1 << (31 - j)) ) 
      fputs("1",stdout);
    else 
      fputs("0", stdout);
  }
  fprintf(stdout, "\nCore Mask (%2u)   : ", APIC_HyperthreadingMaskWidth);
  for (unsigned int j=0; j<32; ++j) {
    if ( APIC_CoreSelectMask & (1 << (31 - j)) ) 
      fputs("1",stdout);
    else 
      fputs("0", stdout);
  }
  printf("\nSocket Mask (%2u) : ", APIC_SocketSelectMaskShift);
  for (unsigned int j=0; j<32; ++j) {
    if ( APIC_SocketSelectMask & (1 << (31 - j)) ) 
      fputs("1",stdout);
    else 
      fputs("0", stdout);
  }
  fputs("\n",stdout);
#endif
  register unsigned int nSockets = 0;
  register unsigned int nCU = 0;
  for (unsigned int thread=0U; thread<info->nOverallCores; ++thread) {
    topology_mask_t TMask = ((topology_mask_t) 1) << (info->APICID[thread] & APIC_HyperthreadingMask);
    unsigned int tmp = (info->APICID[thread] & APIC_CoreSelectMask) >> APIC_HyperthreadingMaskWidth;
    if (tmp > nCU) nCU = tmp;
    TMask |= ((topology_mask_t) 1) << (HYPERTHREADING_MASK_WIDTH + tmp);  
    tmp = (info->APICID[thread] & APIC_SocketSelectMask) >> APIC_SocketSelectMaskShift;
    if (tmp > nSockets) nSockets = tmp;
    TMask |= ((topology_mask_t) 1) << ( SOCKET_MASK_SHIFT + tmp);
    info->TopologyMasks[thread] = TMask;
#ifdef __SYSTEM_DEBUG__
    fprintf(stderr,"APIC DATA (%2u)   : ", thread);
    for (unsigned int j=0; j<32; ++j) {
      if ( info->APICID[thread] & (1 << (31 - j)) ) 
	fputs("1",stderr);
      else 
	fputs("0", stderr);
    }
    fputs("\n",stderr);
#endif
  }
  info->nComputeUnit = 1 + nCU;
  info->nCores = (1 + nCU)*info->nCorePerComputeUnit;
  info->nSockets = 1 + nSockets;
  
  info->HyperthreadingMaskWidth = APIC_HyperthreadingMaskWidth;
  info->HyperthreadingMask      = APIC_HyperthreadingMask;
  info->CoreSelectMask          = APIC_CoreSelectMask;
  info->HyperthreadingMask      = APIC_HyperthreadingMask;
  info->SocketSelectMask        = APIC_SocketSelectMask;
  info->SocketSelectMaskShift   = APIC_SocketSelectMaskShift;
#endif
}


void getSystemInfo(SystemInfo * const info)
{
  time_t creation_time;

  /* clear all data */
  memset(info, 0, sizeof(SystemInfo));
  
  /*
   *  OPERATING SYSTEM INFORMATIONS ------------------------------------------------
   */
  GetOSInfo(info);
 
  /*
   * ARCHITECTURE INFORMATION -----------------------------------------------
   */
  
  for (size_t i=0; i<13; ++i) info->CPU_Vendor[i] = '\0';
  
  /* Available number of logical processor seen by OS */
  info->nOverallCores = (unsigned int) sysconf(_SC_NPROCESSORS_ONLN);
  const char NOBrand[] = "Not supported by this cpu"; 
  strcpy(info->CPU_Name, NOBrand);
  
  long a,c;
  __asm__ __volatile__ (
    /* See if CPUID instruction is supported ... */
    /* ... Get copies of EFLAGS into eax and ecx */
    "pushf\n\t"
    "pop %0\n\t"
    "mov %0, %1\n\t"

    /* ... Toggle the ID bit in one copy and store */
    /*     to the EFLAGS reg */
    "xor $0x200000, %0\n\t"
    "push %0\n\t"
    "popf\n\t"

    /* ... Get the (hopefully modified) EFLAGS */
    "pushf\n\t"
    "pop %0\n\t"
    : "=a" (a), "=c" (c)
    :
    : "cc"
    );
  
  unsigned int rval = 0;
  if (a != c) {
    int eax_value, ebx_value, ecx_value, edx_value;
    
    /* Get the CPU vendor string */
    __asm__ __volatile__ (
            "xorl %%eax, %%eax;"
	    "xorl %%ecx, %%ecx;"
	    "cpuid ;"
	    "movl  %%ebx,   (%2);"
	    "movl  %%edx,  4(%2);"
	    "movl  %%ecx,  8(%2);"
	    //"mov  %%ebx, 4(%0) ;" 
	    : "=&a"(eax_value)
	    : "0" (eax_value), "r"(info->CPU_Vendor)
	    : "memory", /*"%ebp",*/ "%ebx", "%ecx", "%edx");
    
    /* if available query SSE and MMX functionalities */
    if (eax_value >= 0x2) {
      __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx;"
	      "cpuid ;"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(0x1)
	      );
      
      if (edx_value & (1<<23))    rval |= MM_MMX;
      if (edx_value & (1<<25))    rval |= MM_MMXEXT | MM_SSE;
      if (edx_value & (1<<26))    rval |= MM_SSE2;
      if (ecx_value & 1)          rval |= MM_SSE3;
      if (ecx_value & (1<<9) )    rval |= MM_SSSE3;
      if (ecx_value & (1<<19))    rval |= MM_SSE41;
      if (ecx_value & (1<<20))    rval |= MM_SSE42;
      if (ecx_value & (1<<23))    rval |= MM_POPCOUNT;
      if (ecx_value & (1<<28))    rval |= MM_AVX;
      if (ecx_value & (1<<12))    rval |= MM_FMA3;
      
      /* Test OSXSAVE for MIC */
      const int OSXSAVE = (ecx_value & (1<<27));
      
      info->Family = (( eax_value >> 8 ) & 0xF) + (( eax_value >> 20 ) & 0xFF);
      info->Model  = (( eax_value >> 4 ) & 0xF) + ((( eax_value >> 16 ) & 0xF) << 4);

//       const unsigned int BaseFamily = (( eax_value >> 8 ) & 0xF);
//       if (BaseFamily == 0xF) {
// 	info->Family = (( eax_value >> 20 ) & 0xFF) << 4;
// 	info->Model  = ((( eax_value >> 16 ) & 0xF) << 4 ) | (( eax_value >> 4 ) & 0xF);
//       }
//       else {
// 	info->Family = BaseFamily;
// 	info->Model  = (( eax_value >> 4 ) & 0xF);
//       }
      
      /* Test for MIC architecture */
#ifndef __APPLE__
      if (OSXSAVE > 0) {
				__asm__ __volatile__ (
							"xorl %%ecx, %%ecx ; "
							"xgetbv ;"
							: "=a"(eax_value)
							:
							: "%ecx", "%edx"
						);
				if ( ((eax_value & 0x3) == 0x3) && (((eax_value >> 5) & 0x7) == 0x7) ) {
						__asm__ __volatile__ (
							"xorl %%ecx, %%ecx;"
							"cpuid ;"
							: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
							: "0"(7)
						);
					if (eax_value == ebx_value == ecx_value == edx_value == 0) {
#if _SYSTEM_TEST
						printf("Leaf 0x7H is not supported\n");
#endif
					} else {
						if (ebx_value & (1<<16) ) rval |= MM_AVX512F;
						if (ebx_value & (1<<30) ) rval |= MM_AVX512BW;
						if (ebx_value & (1<<28) ) rval |= MM_AVX512CD;
						if (ebx_value & (1<<5) )  rval |= MM_AVX2;
					}
				}
      }
#ifdef _SYSTEM_TEST
      else {
				printf("OSXSAVE is not active");
				if (ecx_value & (1<<26)) {
						printf("but XSAVE is!\n");
				}
				else {
						printf(", neither is XSAVE\n");
				}
				__asm__ __volatile__ (
							"xorl %%ecx, %%ecx;"
							"cpuid ;"
							: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
							: "0"(7));
				if (eax_value == ebx_value == ecx_value == edx_value == 0) {
					printf("Lead 0x7H is not supported\n");

				} 
				else {
					if (ebx_value & (1<<16) ) rval |= MM_AVX512;
					if (ebx_value & (1<<5) )  rval |= MM_AVX2;
				}
      }
#endif
/* 
 * MIC OS does not report MIC featues !!!
 * Enforce when __MIC__ is active 
 */
#ifdef __MIC__
			rval |= MM_AVX512;
#endif
#endif
    }
     __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx; "
	      "cpuid ;"
	      : "=a"(eax_value)
	      : "0"(0x80000000)
	      : "%ebx", "%ecx", "%edx"
	      );
    const unsigned int MaxLeaf = eax_value;
#ifdef __SYSTEM_DEBUG__
    fprintf(stderr, "Maximum extended leaf : 0x%8xh\n", MaxLeaf);
#endif
    if (MaxLeaf >= 0x80000001) {
      __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx;"
	      "cpuid ;"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(0x80000001)
	      );
      
      if (edx_value & (1<<31))  rval |= MM_3DNOW;
      if (edx_value & (1<<30))  rval |= MM_3DNOWEXT;
      if (edx_value & (1<<23))  rval |= MM_MMX;
      if (edx_value & (1<<22))  rval |= MM_MMXEXT;
			if (edx_value & (1<<26))  info->TLBs |= TLB_1GB;
      if (ecx_value & (1<<6))   rval |= MM_SSE4A;
      if (ecx_value & (1<<16))  rval |= MM_FMA4;
      if (ecx_value & (1<<11))  rval |= MM_XOP;
      if (ecx_value & (1<<12))  rval |= MM_FMA3;
    } 
    
    if (MaxLeaf >= 0x80000004) {
      /* Get the CPU brand string */
      __asm__ __volatile__ (
	    "xorl  %%ecx, %%ecx;"
	    "movl  $0x80000002, %%eax;"
	    "cpuid ;"
	    "movl  %%eax,   (%0);"
	    "movl  %%ebx,  4(%0);"
	    "movl  %%ecx,  8(%0);"
	    "movl  %%edx, 12(%0);"
	    "movl  $0x80000003, %%eax;"
	    "cpuid ;"
	    "movl  %%eax, 16(%0);"
	    "movl  %%ebx, 20(%0);"
	    "movl  %%ecx, 24(%0);"
	    "movl  %%edx, 28(%0);"
	    "movl  $0x80000004, %%eax;"
	    "cpuid ;"
	    "movl  %%eax, 32(%0);"
	    "movl  %%ebx, 36(%0);"
	    "movl  %%ecx, 40(%0);"
	    "movl  %%edx, 44(%0);"
	    : 
	    : "r"(info->CPU_Name)
	    : "memory", "%eax", "%ebx", "%ecx", "%edx");
      {
				/* Correct space at the beginning */
				char * ptr = info->CPU_Name;
				size_t i = 0;
				while (ptr[i] == ' ' && i < 48 ) ++i;
				while (i < 48 ) { *ptr++ = info->CPU_Name[i++]; } 
      }
    } 
  
    info->Extensions = rval;
    
    /*
     * TOPOLOGY INFORMATIONS ---------------------------------------------------------
     */
    /* Intel Topology*/
    if (info->CPU_Vendor[0] == 'G' && info->CPU_Vendor[1] == 'e') 
      GetIntelTopology(info);
    /* AMD Topology*/
    else if (info->CPU_Vendor[0] == 'A' && info->CPU_Vendor[1] == 'u')
      GetAMDTopology(info);
  }

#ifdef __NUMA__
  info->NumaAble = (numa_available() < 0) ? false : true;
  if (info->NumaAble) {
    info->nNodes       = (unsigned int) numa_num_configured_nodes();
    info->nCpusPerNode = numa_num_configured_cpus() / info->nNodes;
  }
#endif

#if defined(USE_KNC) && !(defined(__MIC__))
	uint32_t num_engines = 0U;
#ifndef NDEBUG
	{
		struct stat st;
		if (stat("/usr/lib64/libscif.so.0", &st) != 0) {
                        perror("stat");
                        fputs("Unable to find /usr/lib64/libscif.so.0\n", stdout);
                }
		
		if (stat("/usr/lib64/libcoi_host.so.0", &st) != 0) {
			perror("stat");
			fputs("Unable to find /usr/lib64/libcoi_host.so.0\n", stdout);			
		}
	}
#endif
	void * const SCIF = dlopen("/usr/lib64/libscif.so.0", RTLD_LAZY | RTLD_GLOBAL);
	if (SCIF == NULL) {
		fprintf(stderr, "Unable to load SCIF library: %s\n", dlerror());
	}
	void * const COI = dlopen("/usr/lib64/libcoi_host.so.0", RTLD_LAZY);
	if (COI != NULL) {
			COIRESULT result;
			COIACCESSAPI COIRESULT (*COIEngineGetCountPtr)(
				COI_ISA_TYPE in_DeviceType,
				uint32_t *out_pNumEngines);
			COIEngineGetCountPtr = dlsym(COI, "COIEngineGetCount");
			if (COIEngineGetCountPtr)	{
			result = COIEngineGetCountPtr(COI_ISA_KNC, &num_engines);
#ifdef USE_AFFINITY
			if (result == COI_SUCCESS) {
				info->MicEngines = (COI_ENGINE_INFO*) malloc(num_engines*sizeof(COI_ENGINE_INFO));
				if (info->MicEngines) {
					for (uint32_t iMic=0; iMic<num_engines; iMic++) {
						COIENGINE Engine;
						COIACCESSAPI COIRESULT (*COIEngineGetHandlePtr)(
							COI_ISA_TYPE in_DeviceType,
							uint32_t in_EngineIndex,
							COIENGINE      *out_pEngineHandle);
						COIACCESSAPI COIRESULT (*COIEngineGetInfoPtr)(
							COIENGINE in_EngineHandle,
							uint32_t  in_EngineInfoSize,
							COI_ENGINE_INFO *out_pEngineInfo);
						COIACCESSAPI const char * (*COIResultGetNamePtr) (COIRESULT in_ResultCode);
						COIEngineGetHandlePtr = dlsym(COI, "COIEngineGetHandle");
						COIResultGetNamePtr = dlsym(COI, "COIResultGetName");
						COIEngineGetInfoPtr = dlsym(COI, "COIEngineGetInfo");
						if (COIEngineGetHandlePtr && COIResultGetNamePtr && COIEngineGetInfoPtr) {
							result = COIEngineGetHandlePtr(COI_ISA_KNC, iMic, &Engine);
							if (result != COI_SUCCESS) 
								fprintf(stderr,"COIEngineGetHandle result %s\n", COIResultGetNamePtr(result));
							else {
								result = COIEngineGetInfoPtr(Engine, sizeof(COI_ENGINE_INFO), &(info->MicEngines[iMic]));
								if (result != COI_SUCCESS) 
									fprintf(stderr,"COIEngineGetInfo result %s\n", COIResultGetNamePtr(result));
							}
						}
					}
				}
				else {
					fputs("WARNING: getSystemInfo cannot allocate memory for MIC Engines information\n", stderr);
					return;
				}
			}
		}
#endif
	}
#ifndef NDEBUG
	else {
		fprintf(stderr, "Unable to load /usr/lib64/libcoi_host.so.0: %s\n", dlerror());
	}
#endif
	info->nKNCs = num_engines;
#endif
}

void printSystemInfo(FILE * const fd, const SystemInfo * const info)
{
  const static char CacheTypeText[4][12] = { {"None"}, {"Data"}, {"Instruction"}, {"Unified"} };
  char Buffer[LINE_BUFFER_SIZE] __attribute__((aligned(16)));
	
#ifdef MPI_VERSION
  fputs(  TOP_DOUBLE_LINE
          SGL_TEXT_LINE("  Master node system information                                    ")
          MIDDLE_SINGLE_LINE, fd);
#else 
  fputs(  TOP_DOUBLE_LINE
          SGL_TEXT_LINE("  System informations                                               ")
          MIDDLE_SINGLE_LINE, fd);
#endif

  PRINT_SGL_STRING_LINE("     User name     : %s", 21, info->Username, strlen(info->Username));
  PRINT_SGL_STRING_LINE("     Node name     : %s", 21, info->Nodename, strlen(info->Nodename));
  PRINT_SGL_TEXT_LINE("", 0);
  
  PRINT_SGL_STRING_LINE("     Kernel        : %s", 21, info->Release, strlen(info->Release));
  PRINT_SGL_STRING_LINE("     Architecture  : %s", 21, info->Architecture, strlen(info->Architecture));
  PRINT_SGL_TEXT_LINE("", 0);
  
  PRINT_SGL_STRING_LINE("     CPU vendor    : %s", 21, info->CPU_Vendor, strlen(info->CPU_Vendor));
  PRINT_SGL_STRING_LINE("     CPU brand     : %s", 21, info->CPU_Name, strlen(info->CPU_Name));
  PRINT_SGL_STRING_LINE("     Family        : 0x%2.2xh", 24, info->Family, 2);
  PRINT_SGL_STRING_LINE("     Model         : 0x%2.2xh", 24, info->Model, 2);
  
  char * ptr = Buffer;
  *ptr = '\0';
  if (info->Extensions & MM_MMXEXT)   ptr += sprintf(ptr,"MMXEXT ");
  if (info->Extensions & MM_SSE)      ptr += sprintf(ptr,"SSE ");
  if (info->Extensions & MM_SSE2)     ptr += sprintf(ptr,"SSE2 ");
  if (info->Extensions & MM_SSE3)     ptr += sprintf(ptr,"SSE3 ");
  if (info->Extensions & MM_SSSE3)    ptr += sprintf(ptr,"SSSE3 ");
  if (info->Extensions & MM_SSE41)    ptr += sprintf(ptr,"SSE41 ");
  if (info->Extensions & MM_SSE42)    ptr += sprintf(ptr,"SSE42 ");
  if (info->Extensions & MM_3DNOW)    ptr += sprintf(ptr,"3DNow ");
  if (info->Extensions & MM_3DNOWEXT) ptr += sprintf(ptr,"3DNowExt ");
  if (info->Extensions & MM_SSE4A)    ptr += sprintf(ptr,"SSE4a ");
  if (info->Extensions & MM_POPCOUNT) ptr += sprintf(ptr,"POPCOUNT ");
  if (info->Extensions & MM_AVX)      ptr += sprintf(ptr,"AVX ");
  if (info->Extensions & MM_XOP)      ptr += sprintf(ptr,"XOP ");
  if (info->Extensions & MM_FMA3)     ptr += sprintf(ptr,"FMA3 ");
  if (info->Extensions & MM_FMA4)     ptr += sprintf(ptr,"FMA4");
  if (info->Extensions & MM_AVX2)     ptr += sprintf(ptr,"AVX2 ");
  if (info->Extensions & MM_AVX512)   ptr += sprintf(ptr,"AVX512 ");
	if (info->Extensions & MM_AVX512F)   ptr += sprintf(ptr,"AVX512F ");
	if (info->Extensions & MM_AVX512BW)   ptr += sprintf(ptr,"AVX512BW ");
  
  if (ptr >= &Buffer[LINE_SIZE-2-21]) {
    ptr = &Buffer[LINE_SIZE-2-21];
    while ( ptr > &Buffer[0] && *ptr != ' ') --ptr;
    *ptr = '\0';
    PRINT_SGL_STRING_LINE("     CPU extensions: %s", 21, Buffer, strlen(Buffer));
    PRINT_SGL_STRING_LINE("                     %s", 21, ptr+1, strlen(ptr+1));
  } else {
    PRINT_SGL_STRING_LINE("     CPU extensions: %s", 21, Buffer, strlen(Buffer));
  }
  ptr = Buffer;
	if (info->TLBs & TLB_1GB) ptr += sprintf(ptr, "1GB ");
	*ptr = '\0';
	PRINT_SGL_STRING_LINE(  "     TLBs          : %s", 21, Buffer, strlen(Buffer));
  PRINT_SGL_TEXT_LINE("", 0);
  for (int i=0; i<4; i++) {
   snprintf(Buffer, LINE_BUFFER_SIZE, "     Cache Level %1u : %5u kB - %3u ways - %2u B - %-11s",
	    (unsigned int) info->Cache[i].Level,
	    info->Cache[i].Size >> 10,
	    info->Cache[i].Ways,
	    info->Cache[i].Line,
	    &CacheTypeText[info->Cache[i].Type][0]
 	  );
   PRINT_SGL_STRING_LINE("%s",0, Buffer, strlen(Buffer)); 
  }
  PRINT_SGL_TEXT_LINE("", 0);
  if (info->HyperthreadingAvailable) {
    if ( info->HyperthreadingOn) {
      PRINT_SGL_TEXT_LINE("     Hyperthreading: Available", 30);
    } else {
      PRINT_SGL_TEXT_LINE("     Hyperthreading: Available but BIOS disabled", 48);
    }
  }
  
#ifdef USE_AFFINITY
  PRINT_SGL_STRING_LINE("     Sockets       : %3u", 22, info->nSockets, 2);
  if (info->nComputeUnit) {
    PRINT_SGL_STRING_LINE("     Compute unit  : %3u", 22, info->nComputeUnit, 2);
    PRINT_SGL_STRING_LINE("     Core per unit : %3u", 22, info->nCorePerComputeUnit, 2);
  }
  PRINT_SGL_STRING_LINE("     Cores         : %3u", 22, info->nCores, 2);
  if (info->HyperthreadingAvailable) {
     PRINT_SGL_STRING_LINE("     Threads       : %3u", 22, info->nOverallCores/(info->nSockets*info->nCores), 2);
  }
#endif
  PRINT_SGL_STRING_LINE("     Overall cores : %3u", 21, info->nOverallCores, 3);
#ifdef __NUMA__
  if (info->NumaAble) {
    fputs("     NUMA          : Available\n",fd);
    fprintf(fd,"     NUMA nodes    : %u\n", info->nNodes);
  } else {
    fputs("              NUMA : Not available\n",fd);
  }
#endif
#if defined(USE_KNC) && !(defined(__MIC__))
	PRINT_SGL_TEXT_LINE("", 0);
	PRINT_SGL_STRING_LINE("   MIC KNC engines : %3u", 21, info->nKNCs, 3);
	for (unsigned int iMic=0; iMic<info->nKNCs; iMic++) {
		snprintf(Buffer, LINE_BUFFER_SIZE, "%1.1u: %2.2u cores, %1.1u threads @ %u MHz", 
						 iMic, info->MicEngines[iMic].NumCores, info->MicEngines[iMic].NumThreads,
						 info->MicEngines[iMic].CoreMaxFrequency);
		PRINT_SGL_STRING_LINE("     MIC %s", 9, Buffer, strlen(Buffer));
	}
#endif
	PRINT_SGL_TEXT_LINE("", 0);
  fputs( MIDDLE_DOUBLE_LINE, fd);
}

_Bool writeSystemInfo(FILE * const fd, const SystemInfo * const info)
{
    _Bool res = true;
   
    /* Write down a character ID for future modification */
    {	
#ifdef __NUMA__
      const unsigned char c = 1;
#else
      const unsigned char c = 0;
#endif
      res &= (fwrite(&c, sizeof(unsigned char), 1, fd) == 1);
    }
    
    /* Operating system data */
    res &= (fwrite(info->Username, sizeof(char), 64, fd) == 64);
    res &= (fwrite(info->Release, sizeof(char), 64, fd) == 64);
    res &= (fwrite(info->Architecture, sizeof(char), 64, fd) == 64);
    res &= (fwrite(info->Nodename, sizeof(char), 64, fd) == 64);

    /* Architecture data */
    res &= (fwrite(info->CPU_Vendor, sizeof(char), 13, fd) == 13);
    res &= (fwrite(info->CPU_Name, sizeof(char), 48, fd) == 48);
    res &= (fwrite(&(info->Family), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->Family), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->Model), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->Extensions), sizeof(unsigned int), 1, fd) == 1);

    /* Topology */
    res &= (fwrite(&(info->nSockets), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->nCores), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->nOverallCores), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->nComputeUnit), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->nCorePerComputeUnit), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->SocketSelectMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->CoreSelectMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->HyperthreadingMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->SocketSelectMaskShift), sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(&(info->HyperthreadingMaskWidth), sizeof(unsigned int), 1, fd) == 1);
 
#ifdef __NUMA__
    res &= (fwrite(info->nNodes, sizeof(unsigned int), 1, fd) == 1);
    res &= (fwrite(info->nCpusPerNode, sizeof(unsigned int), 1, fd) == 1);
#endif
    {
      const unsigned char c = (info->HyperthreadingAvailable ? 1 : 0 ) | (info->HyperthreadingOn ? 2 : 0 );
      res &= (fwrite(&c, sizeof(unsigned char), 1, fd) == 1);
    }
    
#ifdef __NUMA__
    {
      const unsigned char c = (info->NumaAble) ? 1 : 0;
      res &= (fwrite(&c, sizeof(unsigned char), 1, fd) == 1);
    }
#endif /* __NUMA__ */

//     unsigned int * TopologyMasks;
//     unsigned int * APICID;
    
  return res;
}

_Bool readSystemInfo(FILE * const fd, SystemInfo * const info)
{
    _Bool res = true;
    _Bool PotentialNuma = false;

    /* Read character ID */
    {	
      unsigned char c;
      res &= (fread(&c, sizeof(unsigned char), 1, fd) == 1);
      PotentialNuma = ( c == 1 ) ? true : false;
    }
    
    /* Operating system data */
    res &= (fread(info->Username, sizeof(char), 64, fd) == 64);
    res &= (fread(info->Release, sizeof(char), 64, fd) == 64);
    res &= (fread(info->Architecture, sizeof(char), 64, fd) == 64);
    res &= (fread(info->Nodename, sizeof(char), 64, fd) == 64);

    /* Architecture data */
    res &= (fread(info->CPU_Vendor, sizeof(char), 13, fd) == 13);
    res &= (fread(info->CPU_Name, sizeof(char), 48, fd) == 48);
    res &= (fread(&(info->Family), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->Family), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->Model), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->Extensions), sizeof(unsigned int), 1, fd) == 1);

    /* Topology */
    res &= (fread(&(info->nSockets), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->nCores), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->nOverallCores), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->nComputeUnit), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->nCorePerComputeUnit), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->SocketSelectMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->CoreSelectMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->HyperthreadingMask), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->SocketSelectMaskShift), sizeof(unsigned int), 1, fd) == 1);
    res &= (fread(&(info->HyperthreadingMaskWidth), sizeof(unsigned int), 1, fd) == 1);
 
    if (PotentialNuma) {
#ifdef __NUMA__
      res &= (fread(info->nNodes, sizeof(unsigned int), 1, fd) == 1);
      res &= (fread(info->nCpusPerNode, sizeof(unsigned int), 1, fd) == 1);
#else
      unsigned int i,j;
      res &= (fread(&i, sizeof(unsigned int), 1, fd) == 1);
      res &= (fread(&j, sizeof(unsigned int), 1, fd) == 1);
#endif
    }
    
    {
      unsigned char c;
      res &= (fread(&c, sizeof(unsigned char), 1, fd) == 1);
      info->HyperthreadingAvailable = (c & 0x1) ? true : false;
      info->HyperthreadingOn        = (c & 0x2) ? true : false;
    }
    
    if (PotentialNuma) {
      unsigned char c;
      res &= (fread(&c, sizeof(unsigned char), 1, fd) == 1);
#ifdef __NUMA__
      info->NumaAble = (c & 0x1) ? true : false;
#endif
    }

   info->TopologyMasks = NULL;
   info->APICID        = NULL;
  
   return res;
}

#ifdef USE_AFFINITY
void printMasks(const Affinity_Mask_t * const Masks, const size_t N)
{
  //char Mask[32] __attribute__((aligned(16)));
  fputs(SGL_TEXT_LINE("     CPUSET Masks                                                   ")
	EMPTY_LINE, stdout);
  
  for (size_t iMask = 0; iMask<N; ++iMask) {
    fprintf(stdout, SINGLE_VERTICAL_BAR "%3lu:", iMask+1);
    for (int i=63; i>=0; --i) {
	const int digit = (int) '0' + CPU_ISSET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) &Masks[iMask]);
	fputc(digit, stdout);
    }
    fputs(SINGLE_VERTICAL_BAR "\n", stdout);
  }
  fputs(EMPTY_LINE, stdout);
}
#endif

void printAPIC(const SystemInfo * const info)
{
   fputs(SGL_TEXT_LINE("     CPU APIC Id                                                    ")
	MIDDLE_SINGLE_LINE, stdout);
  if (info->APICID) {
    for (size_t iThread = 0; iThread<info->nOverallCores; ++iThread) {
      fprintf(stdout, SINGLE_VERTICAL_BAR " %3lu :                             ", iThread+1);
      for (int i=31; i>=0; --i) {
	  const int digit = info->APICID[iThread] & (1 << i) ? (int) '1' : (int) '0'; 
	  fputc(digit, stdout);
      }
      fputs(" " SINGLE_VERTICAL_BAR "\n", stdout);
    }
  } else {
    fputs(SGL_TEXT_LINE(" NO APIC ID FOUND                                                   "),stdout);
  }
  fputs(BOTTOM_SINGLE_LINE, stdout);
}

void printTopology(const SystemInfo * const info)
{
  fputs(SGL_TEXT_LINE("     System topology                                                ")
	MIDDLE_SINGLE_LINE
	SGL_TEXT_LINE("           Socket      Core/Compute unit         Thread             "), stdout);
  if (info->APICID) {
    for (size_t iThread = 0; iThread<info->nOverallCores; ++iThread) {
      fprintf(stdout, SINGLE_VERTICAL_BAR " %3lu :  ", iThread+1);
      const unsigned int APIC = info->APICID[iThread];
      unsigned int value = ( APIC & info->SocketSelectMask ) >> (info->SocketSelectMaskShift);
      fprintf(stdout,"%9u ", value+1);
      value = ( APIC & info->CoreSelectMask ) >> info->HyperthreadingMaskWidth;
      fprintf(stdout,"    %18u ", value+1);
      value = APIC & info->HyperthreadingMask;
      fprintf(stdout,"        %6u            ", value+1);
     
      fputs(" " SINGLE_VERTICAL_BAR "\n", stdout);
    }
  } else {
    fputs(SGL_TEXT_LINE(" NO APIC ID FOUND                                                   "),stdout);
  }
  fputs(BOTTOM_SINGLE_LINE, stdout);
}

#ifdef _SYSTEM_TEST
#include <getopt.h>
#include "gtlVersion.h"
/* Usage function */
void usage()
{
   fputs(   SGL_TEXT_LINE(" Command line arguments                                             ")
            MIDDLE_SINGLE_LINE
            SGL_TEXT_LINE("      --apic                show APIC informations                  ")
            SGL_TEXT_LINE("    MODE                                                            ")
            SGL_TEXT_LINE("      --test                                                        ")
            SGL_TEXT_LINE("      --query                                                       ")
            SGL_TEXT_LINE("    INPUT                                                           ")
            SGL_TEXT_LINE("      --socket              number of threads to use                ")
	          SGL_TEXT_LINE("      --core                number of threads to use                ")
	          SGL_TEXT_LINE("      --thread              number of threads to use                ")	    
            BOTTOM_DOUBLE_LINE
    ,stderr);
    exit(1);
}

int main (int argc, char * argv[]) 
{
  int c;
  const char opt_to_test[] = "s:c:t:qTh";
  unsigned int SocketId=-1, CoreId=-1, ThreadId=-1;
  _Bool isQuery=false, isTest=false, ShowAPIC=false;
  const struct option long_options[] = {
    /* These options set a flag. */
    {"help",		no_argument,       	0,	'h'},
    /* These options don't set a flag. We distinguish them by their indices. */
    {"socket",	required_argument,	0,	's'},
    {"core",		required_argument,	0,	'c'},
    {"thread",	required_argument,	0,	't'},
    {"query",		no_argument,    		0,	'q'},
    {"test",		no_argument,				0,	'T'},
    {"apic",		no_argument,				0,	'a'},
    {0, 			0, 		0,  	0  }
  };
  SystemInfo Info;
  
fputs(
  TOP_DOUBLE_LINE\
  SINGLE_VERTICAL_BAR "       System tool      v" GTL_VERSION "          built on " __DATE__ "  " SINGLE_VERTICAL_BAR "\n",
	stdout);

    
  while (1) {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c) {
        /* -----------------------------------------------------*/
        /*              STD ARGUMENTS                           */
        /* -----------------------------------------------------*/
        case 0:
            /* If this option set a flag, do nothing else now. */
            //if (long_options[option_index].flag != 0)
            break;
	case 'a':
	    ShowAPIC = true;
	    break;
	case 's':
	    SocketId = (unsigned int) atoi(optarg);
	    break;
	case 'c':
	    CoreId = (unsigned int) atoi(optarg);
	    break;
	case 't':
	    ThreadId = (unsigned int) atoi(optarg);
	    break;
	case 'T':
	    isTest = true;
	    break;
	case 'q':
	    isQuery = true;
	    break;
        /* -------------------------------------------------*/
        /*              OTHER ARGUMENTS                     */
        /* -------------------------------------------------*/
        case 'h':
            usage();
            break;
        default:
            usage();
        }
    }
  
//   if (optind == argc) usage();
  
  if (!isQuery && !isTest) isQuery=true;
  if (isQuery && isTest) {
    fputs("You may not have both query and test together\n", stderr);
    exit(1);
  }

  getSystemInfo(&Info);
  printSystemInfo(stdout, &Info);
  if (ShowAPIC) { printAPIC(&Info); printTopology(&Info); }
  
  if (isQuery) {
#ifdef USE_AFFINITY  
    fprintf(stderr,"Input %u %u %u\n", SocketId, CoreId, ThreadId);
    Affinity_Mask_t * Masks;
    const unsigned int count = getMasks(&Info, SocketId, CoreId, ThreadId, &Masks); 
    printf("%3u masks satisfy the criteria\n"
           "-------------------------------\n", count);
    unsigned int nTotalThreads = Info.nOverallCores; 
    for (unsigned int i=0; i<count; ++i) {
      for (unsigned int j=0; j<nTotalThreads; ++j) {
	if ( CPU_ISSET_S(nTotalThreads-1-j, sizeof(Affinity_Mask_t), (cpu_set_t*) &Masks[i]) ) 
		fputs("1",stdout);
	else 
		fputs("0", stdout);
      }
      fputs("\n",stdout);
    }
#else
   fputs("Query only possible when affinity was build in!\n", stderr);
#endif /* USE_AFFINITY */
  }
  
  if (isTest) {
#ifdef USE_AFFINITY
      Affinity_Mask_t Mask, BackupMask;
      if ( sched_getaffinity(0, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask)) {
				fputs("Error getting affinity!\n",stderr);
				exit(1);
      }
      puts("Current mask:"); 
      unsigned int nTotalThreads = Info.nOverallCores; 
      for (unsigned int j=0U; j<nTotalThreads; ++j) {
	if ( CPU_ISSET_S(nTotalThreads-1-j, sizeof(Affinity_Mask_t), (cpu_set_t*) &Mask) )
		fputs("1",stdout);
	else
		fputs("0", stdout);
      }
      fputs("\n",stdout);
#else
      fputs("test option not available without affinity enabled.", stdout);
#endif
  }
  
  freeSystemInfo(&Info);
  
  return 1; 
}
#endif /* _SYSTEM_TEST */
