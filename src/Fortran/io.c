#include <stdio.h>

int getc_(ch, ch_len)
char *ch;
int *ch_len;
{
    int ret_val;
    *ch = getc(stdin);
    if(*(unsigned char *)ch == 255) { 
	ret_val = -1;
    } else {
	ret_val = -0;
    }
    return ret_val;
} 
