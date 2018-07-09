#include "global.h"
#include <stdio.h>


void evenIndices(int *index_array, int *target_array)
{

}

void oddIndices(int *index_array, int *target_array)
{

}

/*
 * i(x) = x1 + N1*x2 + N1*N2*x3 + ...
 * We work under the assumption that N1 = N2 = .. = NN.
 * Might work always but I am too lazy too check.
 * x3 = (int) i(x)/(N1*N2)
 * x2 = (int) (i(x) - N1*N2*x3)/N1
 * x1 = (int) (i(x) - N1*x2 - N1*N2*x3)
 */
void indexArray(
                int index,
                int *target_array /* has to be of dim ndim */
                )
{
  int *subFact = (int* ) calloc(ndim, sizeof(int));
  subFact[0] = 0;               /* not used, just for consistency */
  subFact[1] = 1;
  for (int i=2; i<=ndim; i++)
    {
      subFact[i] = subFact[i-1] * lsize[i-1];
    }

  int calcIndex;
  double cand = 0.0;

  for (int i=ndim; i>0; i--)
    {
      cand = index;

      for (int j=(ndim-i); j>0; j--)
        {
          cand = cand - subFact[j+1]*target_array[j+1];

        }

      if (cand < 0.0) {
        cand = 0.0;
      }

      target_array[i] = (int) (cand / subFact[i]);
    }
}
