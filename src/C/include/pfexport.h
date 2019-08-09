#ifndef _PF_EXPORT_H
#define _PF_EXPORT_H

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

#endif // _PF_EXPORT_H
