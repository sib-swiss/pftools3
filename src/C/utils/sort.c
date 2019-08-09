void iqsort(int * const restrict data, const size_t N)
{
  int i,j;
  int itmp;
  int t,v;

  if (N<=1) return;
  v = data[0];
  i = 0;
  j = N;
  for (;;)
    {
      while(data[++i] < v && i <  N) {}
      while(data[--j] > v) {}
      if (i>=j)
	{ break; }
      else
	{
	  t = data[i];
	  data[i] = data[j];
	  data[j] = t; 
	}
    }
  t = data[i-1];
  data[i-1] = data[0];
  data[0] = t;
  iqsort(data,i-1);
  iqsort(data+i,N-i);
}