/*
 *  Thierry Schuepbach
 *  Vital-IT, SIB Swiss Institute of Bioinformatics
 *
 *  Copyright (C) 2011-2019 Thierry Schuepbach
 *  Copyright (C) 1995-2010 Philip Bucher
 *
 *  This file is part of PFTools.
 *
 *  PFTools is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as published
 *  by the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  PFTools is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with PFTools; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
#ifndef _PFCONFIG_H
#define _PFCONFIG_H

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
 * PfTools integer format
 ********************************************************************************
 */
#cmakedefine USE_32BIT_INTEGER
#ifdef USE_32BIT_INTEGER
# define __USE_32BIT_INTEGER__
#endif
 /*
 ********************************************************************************
 * Enforce inline function
 ********************************************************************************
 */
#ifndef __ALWAYS_INLINE
# ifdef __GNUC__
#  define  __ALWAYS_INLINE __attribute__((__gnu_inline__, __always_inline__, __artificial__))
#  define __inline inline
# else
#  define  __ALWAYS_INLINE __attribute__((__always_inline__))
# endif
#endif

/*
 ********************************************************************************
 * Regular Expression library
 ********************************************************************************
 */
#cmakedefine USE_PCRE2
#cmakedefine USE_PCRE

/*
 ********************************************************************************
 * Import/Export definition
 ********************************************************************************
 */
#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define PFEXPORT    __declspec( dllexport )
      #define PFIMPORT    __declspec( dllimport )
    #else
      #define PFEXPORT
      #define PFIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define PFEXPORT    __declspec( dllexport )
    #define PFIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// PFEXPORT to be visible.
    #define PFEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define PFIMPORT
  #endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef PFEXPORT
  #define PFEXPORT
  #define PFIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case), if
// cmake didn't set the corresponding macro xxxx_EXPORTS when the
// corresponding library is built (DLLIMPEXP is set to DLLEXPORT
// then)
#if defined ( pftools_EXPORTS )
  #define PFIMPEXP    PFEXPORT
  #define PFIMPEXP_DATA( type )    PFEXPORT type
#else
  #define PFIMPEXP    PFIMPORT
  #define PFIMPEXP_DATA( type )    PFIMPORT type
#endif


#endif /* _PFCONFIG_H */
