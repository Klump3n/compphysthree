<<<<<<< HEAD
#include "added_stuff.h"
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
#include "global.h"
#include <stdio.h>


<<<<<<< HEAD
void evenOddIndices(int *evenArray, int *oddArray)
{
  int j=0,
    k=0;

  int *target_array = (int *) calloc(ndim, sizeof(int));
  int evenOdd = 0;

  for (int i=0; i<nvol; i++)
    {

      int *target_array = (int *) malloc(ndim*sizeof(int));
      memset(target_array, 0, ndim*sizeof(int));
      indexArray(i, target_array);

      evenOdd = indexEvenOdd(target_array);

      if (evenOdd == 0)
        {
          evenArray[j] = i;
          j++;
        }
      else
        {
          oddArray[k] = i;
          k++;
        }
    }
}

int indexEvenOdd(int *index_array)
{
  int sum = 0;
  for (int i=1; i<=ndim; i++)
    {
      sum += index_array[i];
    }

  return sum % 2;               /* 0 if even, 1 if odd */
=======
void evenIndices(int *index_array, int *target_array)
{

}

void oddIndices(int *index_array, int *target_array)
{

>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
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
<<<<<<< HEAD
                int *target_array /* has to be of dim ndim, index 0 will not be used */
=======
                int *target_array /* has to be of dim ndim */
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
                )
{
  int *subFact = (int* ) calloc(ndim, sizeof(int));
  subFact[0] = 0;               /* not used, just for consistency */
  subFact[1] = 1;
  for (int i=2; i<=ndim; i++)
    {
      subFact[i] = subFact[i-1] * lsize[i-1];
    }

<<<<<<< HEAD
=======
  int calcIndex;
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
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
