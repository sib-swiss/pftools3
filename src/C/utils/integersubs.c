#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[])
{
  size_t i;
  int input1, input2;
  unsigned int mask = 0x80000000;
  int pbit;
  const int bit = (int) '0';
  int exponent;
  
  if (argc >= 3) {
    input1  = atoi(argv[1]);
    input2  = atoi(argv[2]);
    printf("Integer A %6i         : ", input1);
    for (i=0; i<32; i++) {
	pbit = (*((unsigned int*) &input1) & mask) ? bit + 1 : bit; 
	fputc(pbit,stdout);
	mask >>= 1;
    }
    mask = 0x80000000;
    printf("\nInteger B %6i         : ", input2);
    for (i=0; i<32; i++) {
	pbit = (*((unsigned int*) &input2) & mask) ? bit + 1 : bit; 
	fputc(pbit,stdout);
	mask >>= 1;
    }
    input1 -= input2;
    printf("\nInteger A-B %6i       : ", input1);
    mask = 0x80000000;
    for (i=0; i<32; i++) {
	pbit = (*((unsigned int*) &input1) & mask) ? bit + 1 : bit; 
	fputc(pbit,stdout);
	mask >>= 1;
    }
    
    fputc('\n', stdout);
  }
  return 0;
}