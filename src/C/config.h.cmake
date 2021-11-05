#ifndef _CONFIG_H
#define _CONFIG_H
/*
 ********************************************************************************
 * SPECIFIC ARCHITECTURE
 ********************************************************************************
 */
#cmakedefine USE_WINAPI
#ifdef USE_WINAPI
#define __USE_WINAPI__
#endif

/* NOTE Compilation fails if _GNU_SOURCE is not defined with affinity off */
#if defined(__GNUC__) && !defined(USE_WINAPI)
# define _GNU_SOURCE
#endif

/*
 ********************************************************************************
 * Math Library
 ********************************************************************************
 */
#ifdef __INTEL_COMPILER
# include <mathimf.h>
#elif __IBMC__
# include <math.h>
# include <mass.h>
#else
# include <math.h>
#endif

/*
 ********************************************************************************
 * ALLOCA header file
 ********************************************************************************
 */
#cmakedefine HAVE_ALLOCA_H
#if !defined(HAVE_ALLOCA_H) && defined(__GNUC__)
/* This is not defined in MinGW */
#define alloca __builtin_alloca
#endif

/*
 ********************************************************************************
 * MM_MALLOC header file
 ********************************************************************************
 */
#cmakedefine HAVE_MM_MALLOC_H 1
#ifdef HAVE_MM_MALLOC_H
# include <mm_malloc.h>
#else
#include <stdlib.h>
/* We can't depend on <stdlib.h> since the prototype of posix_memalign
   may not be visible.  */
#ifndef __cplusplus
extern int posix_memalign (void **, size_t, size_t);
#else
extern "C" int posix_memalign (void **, size_t, size_t) throw ();
#endif

static inline void * __ALWAYS_INLINE
_mm_malloc (size_t size, size_t alignment)
{
  void *ptr;
  if (alignment == 1)
    return malloc (size);
  if (alignment == 2 || (sizeof (void *) == 8 && alignment == 4))
    alignment = sizeof (void *);
  if (posix_memalign (&ptr, alignment, size) == 0)
    return ptr;
  else
    return NULL;
}

static inline void __ALWAYS_INLINE
_mm_free (void * ptr)
{
  free (ptr);
}
#endif

/*
 ********************************************************************************
 * Haru PDF Library
 ********************************************************************************
 */
#cmakedefine USE_PDF

/*
 ********************************************************************************
 * GD Library
 ********************************************************************************
 */
#cmakedefine USE_GD

/*
 ********************************************************************************
 * PLplot Library
 ********************************************************************************
 */
#cmakedefine USE_PLPLOT

/*
 ********************************************************************************
 * Mapping database to memory or opening it through libC
 ********************************************************************************
 */
#cmakedefine USE_MMAP
#ifndef USE_MMAP
# undef __USE_MMAP__
# define SETUP_DATABASE_ACCESS(fileoraddress) FILE* inSequence = fopen((fileoraddress), "r")
# define GET_DATABASE_SEQUENCE(dest, offsets, id) ReadSequenceIndex((dest), (size_t) (id), inSequence, (offsets))
# define UNSET_DATABASE_ACCESS() fclose(inSequence)
#else
# include <sys/mman.h>
# define __USE_MMAP__
# define SETUP_DATABASE_ACCESS(fileoraddress) const char * const restrict SequenceFileMap = (fileoraddress);
# ifndef MMAP_DEBUG
#  define GET_DATABASE_SEQUENCE(dest, offsets) MMAP_ReadSequenceIndex((dest), SequenceFileMap, (offsets), 0)
# else
#  define GET_DATABASE_SEQUENCE(dest, offsets) MMAP_ReadSequenceIndex((dest), SequenceFileMap, (offsets), 0\
                                                   ,((struct ThreadData*) _Data)->threadId, 0, *(((struct ThreadData*) _Data)->maplength))
# endif
# define UNSET_DATABASE_ACCESS()
#endif

/*
 ********************************************************************************
 * PfTools versioning
 ********************************************************************************
 */
#define PF_VERSION "${VERSION}"
#define PF_MAJOR_VERSION ${MAJOR_VERSION}
#define PF_MINOR_VERSION ${MINOR_VERSION}

/*
 ********************************************************************************
 * Executable build using inline functions
 ********************************************************************************
 */
#ifndef BUILD_LIBRARY
#define __USE_INLINE_FUNCTIONS__
#endif
#endif /* _CONFIG_H */
