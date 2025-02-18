/* gb_flip.c */
/* This file is part of the Stanford GraphBase (c) Stanford University 1993 */

#include <stdio.h>
#include "gb_flip.h"

/* Private declarations */
static long A[56] = {-1}; /* pseudo-random values */

/* External declarations */
long *gb_fptr=A; /* the next |A| value to be exported */

#define mod_diff(x,y) (((x)-(y))&0x7fffffff) /* difference modulo $2^{31}$ */

/* External functions */
long gb_flip_cycle(void)
{
  register long *ii, *jj;
  for (ii=&A[1],jj=&A[32];jj<=&A[55];ii++,jj++)
    *ii=mod_diff(*ii,*jj);
  for (jj=&A[1];ii<=&A[55];ii++,jj++)
    *ii=mod_diff(*ii,*jj);
  gb_fptr=&A[54];
  return A[55];
}

void gb_init_rand(long seed)
{
  register long i;
  register long prev=seed, next=1;
  seed=prev=mod_diff(prev,0); /* strip off the sign */
  A[55]=prev;
  for (i=21; i; i=(i+21)%55) {
    A[i]=next;
    next=mod_diff(prev,next);
    if (seed&1) seed=0x40000000+(seed>>1);
    else seed>>=1; /* cyclic shift right 1 */
    next=mod_diff(next,seed);
    prev=A[i];
  }
  (void) gb_flip_cycle();
  (void) gb_flip_cycle();
  (void) gb_flip_cycle();
  (void) gb_flip_cycle();
  (void) gb_flip_cycle();
}

#define two_to_the_31 ((unsigned long)0x80000000)

long gb_unif_rand(long m)
{
  register unsigned long t=two_to_the_31-(two_to_the_31 % m);
  register long r;
  do{
    r=gb_next_rand();
  }while (t<=(unsigned long)r);
  return r%m;
}
