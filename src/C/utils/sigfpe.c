/*   Title:  Floating-point exception handling example
    Author:  David N. Williams
      File:  fe-handlng-example.c
   License:  Public Domain
   Version:  0.5.0
   Started:  21-Sep-09
   Revised:  22-Sep-09
   Revised:  30-Sep-09 (comment typo)

This code is an example of alternate, nondefault handling of
IEEE 754 floating-point exceptions in OS X and Linux, based on
the GNU functions feenableexcept(), fedisableeexcept(), and
fegetexcept() [in libm], plus POSIX sigaction().

The GNU functions above are not implemented in OS X Leopard,
gcc 4.x, but are present in Linux.  We implement them here for
OS X, at least until the underlying mechanism is no longer
supported by Apple.

The mechanism is to use the POSIX functions fegetenv() and
fesetenv(), which *are* present in OS X, to manipulate the ppc
and intel floating-point control registers, after changing bits
in fields corresponding to those registers in the fenv_t data
type.

Assembly language code to directly access the floating-point
status and control registers for ppc and intel is also included.

This example grew out of an update to legacy code for Apple
ppc's.  The original legacy code is in Listing 7-1 in "PowerPC
Numerics", 2004:

http://lists.apple.com/archives/unix-porting/2003/May/msg00026.html

Another version of the ppc legacy code is here:

http://developer.apple.com/documentation/Performance/Conceptual/Mac_OSX_Numerics/Mac_OSX_Numerics.pdf

Terry Lambert pointed out that our naive update of the legacy
example to Mac OS X Leopard made egregious unsupported use of
system context structures in the handler.  See his reply to

http://lists.apple.com/archives/Darwin-dev/2009/Sep/msg00091.html

The example in this file is more plain vanilla, and aims at
alternate handling that does not return to the application, but
rather aborts with a diagnostic message.

To compile it under Mac OS X, execute:

  cc -o fe-handling fe-handling-example.c

To compile it under Linux, execute:

  cc -DLINUX -lm -o fe-handling fe-handling-example.c
*/

#ifdef LINUX
/* BEGIN quote
http://graphviz.sourcearchive.com/documentation/2.16/gvrender__pango_8c-source.html
*/
/* _GNU_SOURCE is needed (supposedly) for the feenableexcept
 * prototype to be defined in fenv.h on GNU systems.
 * Presumably it will do no harm on other systems.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/* We are not supposed to need __USE_GNU, but I can't see
 * how to get the prototype for fedisableexcept from
 * /usr/include/fenv.h without it.
 */
#ifndef __USE_GNU
#define __USE_GNU
#endif
/* END quote */
#endif // LINUX

#include <fenv.h>

#define DEFINED_PPC      (defined(__ppc__) || defined(__ppc64__))
#define DEFINED_INTEL    (defined(__i386__) || defined(__x86_64__))

#ifndef LINUX
#if DEFINED_PPC

#define FE_EXCEPT_SHIFT 22  // shift flags right to get masks
#define FM_ALL_EXCEPT    FE_ALL_EXCEPT >> FE_EXCEPT_SHIFT

/* GNU C Library:
http://www.gnu.org/software/libc/manual/html_node/Control-Functions.html

     - Function: int fegetexcept (int excepts)

       The function returns a bitmask of all currently enabled
       exceptions.  It returns -1 in case of failure.

   The excepts argument appears in other functions in fenv.h,
   and corresponds to the FE_xxx exception flag constants.  It
   is unclear whether the bitmask is for the flags or the masks.
   We return that for the flags, which corresponds to the
   excepts argument in feenableexcept(excepts) and
   fedisableexcept(excepts).  In GNU/Linux the argument is void,
   and that's what we implement.  Linux "man fegetenv" appears
   to suggest that it's the mask corresponding to bits in
   excepts that is returned.
*/
static int
fegetexcept (void)
{
  static fenv_t fenv;

  return ( fegetenv (&fenv) ? -1 :
    (
      ( fenv & (FM_ALL_EXCEPT) ) << FE_EXCEPT_SHIFT )
    );
}

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;

  fenv = (fenv & ~new_excepts) | new_excepts;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

static int
fedisableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int still_on = ~( (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT ),
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;

  fenv &= still_on;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#elif DEFINED_INTEL

static int
fegetexcept (void)
{
  static fenv_t fenv;

  return fegetenv (&fenv) ? -1 : (fenv.__control & FE_ALL_EXCEPT);
}

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

static int
fedisableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // mask
  fenv.__control |= new_excepts;
  fenv.__mxcsr   |= new_excepts << 7;

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#endif  // PPC or INTEL enabling
#endif  // not LINUX

#if DEFINED_PPC

#define getfpscr(x)    asm volatile ("mffs %0" : "=f" (x));
#define setfpscr(x)    asm volatile ("mtfsf 255,%0" : : "f" (x));

typedef union {
    struct {
        unsigned long hi;
        unsigned long lo;
    } i;
    double d;
} hexdouble;

#endif  // DEFINED_PPC

#if DEFINED_INTEL

// x87 fpu
#define getx87cr(x)    __asm__ ("fnstcw %0" : "=m" (x));
#define setx87cr(x)    __asm__ ("fldcw %0"  : "=m" (x));
#define getx87sr(x)    __asm__ ("fnstsw %0" : "=m" (x));

// SIMD, gcc with Intel Core 2 Duo uses SSE2(4)
#define getmxcsr(x)    __asm__ ("stmxcsr %0" : "=m" (x));
#define setmxcsr(x)    __asm__ ("ldmxcsr %0" : "=m" (x));

#endif  // DEFINED_INTEL

#include <signal.h>
#include <stdio.h>   // printf()
#include <stdlib.h>  // abort(), exit()
#include <limits.h>
#include <float.h>

static char *fe_code_name[] = {
  "FPE_NOOP",
  "FPE_FLTDIV", "FPE_FLTINV", "FPE_FLTOVF", "FPE_FLTUND",
  "FPE_FLTRES", "FPE_FLTSUB", "FPE_INTDIV", "FPE_INTOVF"
  "FPE_UNKNOWN"
};

/* SAMPLE ALTERNATE FP EXCEPTION HANDLER

   The sample handler just reports information about the
   exception that invoked it, and aborts.  It makes no attempt
   to restore state and return to the application.

   More sophisticated handling would have to confront at least
   these issues:

     * interface to the system context for restoring state
     * imprecision of interrupts from hardware for the intel x87
       fpu (but not the SIMD unit, nor the ppc)
     * imprecision of interrupts from system software
*/
static void
fhdl ( int sig, siginfo_t *sip, ucontext_t *scp )
{
  int fe_code = sip->si_code;
  unsigned int excepts = fetestexcept (/*FE_ALL_EXCEPT*/-1);

  switch (fe_code)
  {
#ifdef FPE_NOOP  // occurs in OS X
    case FPE_NOOP:   fe_code = 0; break;
#endif
    case FPE_FLTDIV: fe_code = 1; break; // divideByZero
    case FPE_FLTINV: fe_code = 2; break; // invalid
    case FPE_FLTOVF: fe_code = 3; break; // overflow
    case FPE_FLTUND: fe_code = 4; break; // underflow
    case FPE_FLTRES: fe_code = 5; break; // inexact
    case FPE_FLTSUB: fe_code = 6; break; // invalid
    case FPE_INTDIV: fe_code = 7; break; // overflow
    case FPE_INTOVF: fe_code = 8; break; // underflow
            default: fe_code = 9;
   }

  if ( sig == SIGFPE )
  {
#if DEFINED_INTEL
    unsigned short x87cr,x87sr;
    unsigned int mxcsr;

    getx87cr (x87cr);
    getx87sr (x87sr);
    getmxcsr (mxcsr);
    printf ("X87CR:   0x%04X\n", x87cr);
    printf ("X87SR:   0x%04X\n", x87sr);
    printf ("MXCSR:   0x%08X\n", mxcsr);
#endif

#if DEFINED_PPC
   hexdouble t;

   getfpscr (t.d);
   printf ("FPSCR:   0x%08X\n", t.i.lo);
#endif

    printf ("signal:  SIGFPE with code %s\n", fe_code_name[fe_code]);
    fputs(  "except mask : ", stdout);
    for (int i=0; i<32; i++) {
	unsigned char c = ( excepts & (1 << i) ) ? '1' : '0';
	fputc(c, stdout);
    }

    printf ("\ninvalid flag        :  0x%04X\n", excepts & FE_INVALID);
    printf ("divByZero flag      :  0x%04X\n", excepts & FE_DIVBYZERO);
    printf ("float overflow flag :  0x%04X\n", excepts & FPE_FLTOVF);
    printf ("float precision flag:  0x%04X\n", excepts & FPE_FLTRES);
    printf ("float underflow flag:  0x%04X\n", excepts & FPE_FLTUND);
    printf ("int underflow flag  :  0x%04X\n", excepts & FPE_INTOVF);
  }
  else printf ("Signal is not SIGFPE, it's %i.\n", sig);

  abort();
}

int main (int argc, char **argv)
{
    double s;
    int a = INT_MIN;
    struct sigaction act;

    act.sa_sigaction = (void(*))fhdl;
    sigemptyset (&act.sa_mask);
    act.sa_flags = SA_SIGINFO;

    printf ("Old int underflow exception: 0x%08X\n", feenableexcept (FE_ALL_EXCEPT));
//    printf ("Old divByZero exception    : 0x%08X\n", feenableexcept (FE_DIVBYZERO));
//  printf ("Old invalid exception      : 0x%08X\n", feenableexcept (FE_INVALID));
    printf ("New fp exception           : 0x%08X\n", fegetexcept ());

    // set handler
    if (sigaction(SIGFPE, &act, (struct sigaction *)0) != 0)
    {
        perror("Yikes");
        exit(-1);
    }
    printf("Before %i\n", a);
    a = a  - 1;
    printf("After  %i\n", a);
//      s = 1.0 / 0.0;  // FE_DIVBYZERO
//     s = 0.0 / 0.0;  // FE_INVALID
//      a = 0;
//      while (1) a-=1;
    return 0;
}

/*

#define _GNU_SOURCE
#include <signal.h>
#include <ucontext.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <time.h>
#include <execinfo.h>
#include <limits.h>
#include <fenv.h>
#include <sys/types.h>

#define arrsizeof(x) (sizeof(x)/sizeof(x[0]))

//
//
//
char *getProcessName (void)
{
  char buffer [128];
  FILE *fp;

  sprintf (buffer, "/proc/%d/stat", getpid ());

  if ((fp = fopen (buffer, "r")))
  {
    if (fgets (buffer, sizeof (buffer), fp))
    {
      char *s;

      if ((s = index (buffer, ')')))
      {
        *s = '\0';

        if ((s = index (buffer, '(')))
          return ++s;
      }
    }
  }

  return NULL;
}

//
//
//
static void sighandlerPrint (FILE *fp, int signo, int code, ucontext_t *context, void *bt [], int bt_size)
{
  char *processName;
  time_t ttime = time (NULL);

  fprintf (fp, "%s", ctime (&ttime));
  fprintf (fp, "PID=%d (%s)\n", getpid (), (processName = getProcessName ()) ? processName : "unknown");
  fprintf (fp, "signo=%d/%s\n", signo, strsignal (signo));
  fprintf (fp, "code=%d (not always applicable)\n", code);
  fprintf (fp, "\nContext: 0x%08lx\n", (unsigned long) context);
  fprintf (fp, "    gs: 0x%08x   fs: 0x%08x   es: 0x%08x   ds: 0x%08x\n"
               "   edi: 0x%08x  esi: 0x%08x  ebp: 0x%08x  esp: 0x%08x\n"
               "   ebx: 0x%08x  edx: 0x%08x  ecx: 0x%08x  eax: 0x%08x\n"
               "  trap:   %8u  err: 0x%08x  eip: 0x%08x   cs: 0x%08x\n"
               "  flag: 0x%08x   sp: 0x%08x   ss: 0x%08x  cr2: 0x%08lx\n",
               context->uc_mcontext.gregs [REG_GS],     context->uc_mcontext.gregs [REG_FS],   context->uc_mcontext.gregs [REG_ES],  context->uc_mcontext.gregs [REG_DS],
               context->uc_mcontext.gregs [REG_EDI],    context->uc_mcontext.gregs [REG_ESI],  context->uc_mcontext.gregs [REG_EBP], context->uc_mcontext.gregs [REG_ESP],
               context->uc_mcontext.gregs [REG_EBX],    context->uc_mcontext.gregs [REG_EDX],  context->uc_mcontext.gregs [REG_ECX], context->uc_mcontext.gregs [REG_EAX],
               context->uc_mcontext.gregs [REG_TRAPNO], context->uc_mcontext.gregs [REG_ERR],  context->uc_mcontext.gregs [REG_EIP], context->uc_mcontext.gregs [REG_CS],
               context->uc_mcontext.gregs [REG_EFL],    context->uc_mcontext.gregs [REG_UESP], context->uc_mcontext.gregs [REG_SS],  context->uc_mcontext.cr2
  );

  fprintf (fp, "\n%d elements in backtrace\n", bt_size);
  fflush (fp);

  backtrace_symbols_fd (bt, bt_size, fileno (fp));
}

//
//
//
static void sighandlerABRT (int signo, struct siginfo *si, void *ctx)
{
  void *bt [128];
  int bt_size;

  bt_size = backtrace (bt, arrsizeof (bt));

  sighandlerPrint (stderr, signo, si->si_code, (ucontext_t *) ctx, bt, bt_size);

  exit (1);
}

//
//
//
int installSignalHandlers (void)
{
  struct sigaction sa;

  sigemptyset (&sa.sa_mask);
  sigaddset (&sa.sa_mask, SIGABRT);

  sa.sa_flags = SA_ONESHOT | SA_SIGINFO;
  sa.sa_sigaction = sighandlerABRT;

  if (sigaction (SIGFPE, &sa, NULL))
  {
    fprintf (stderr, "sigaction failed, line %d, %d/%s\n", __LINE__, errno, strerror (errno));
    exit (1);
  }

  return 1;
}

//
//
//
int main (void)
{
  #pragma STDC FENV_ACCESS ON
  int a = INT_MIN;
  int b = 2;
  int c = 0;

  feenableexcept (FE_NOMASK_ENV);
  installSignalHandlers ();

  c = a - b;

  printf ("c=%d\n", c);

  exit (0);
}*/