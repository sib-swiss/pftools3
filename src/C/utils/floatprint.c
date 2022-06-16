#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[])
{
  size_t i;
  int input;
  float output;
  unsigned int mask = 0x80000000;
  int pbit;
  const int bit = (int) '0';
  int exponent;

  if (argc >= 2) {
    input  = atoi(argv[1]);
    output = (float) input;
    printf("Integer %6i         : ", input);
    for (i=0; i<32; i++) {
	pbit = (*((unsigned int*) &input) & mask) ? bit + 1 : bit;
	fputc(pbit,stdout);
	mask >>= 1;
    }

    fputc('\n', stdout);
    fputs("is 32 bit float format : ",stdout);

    mask = 0x80000000;
    /* Sign */
    pbit = (*((unsigned int*) &output) & mask) ? bit + 1 : bit;
    fputc(pbit,stdout);fputc(' ',stdout);
    mask >>= 1;
    /* Exponent */
    for (i=0; i<8; i++) {
	pbit = (*((unsigned int*) &output) & mask) ? bit + 1 : bit;
	fputc(pbit,stdout);
	mask >>= 1;
    }
    exponent = (int) ((*((unsigned int*) &output) >> 23 ) & 0xFF);
    printf("(%4i) ",exponent);
    /* Mantissa */
    for (i=0; i<23; i++) {
	pbit = (*((unsigned int*) &output) & mask) ? bit + 1 : bit;
	fputc(pbit,stdout);
	mask >>= 1;
    }

    fputc('\n', stdout);
  }
  return 0;
}