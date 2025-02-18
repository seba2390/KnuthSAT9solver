/* gb_flip.h */
/* This file is part of the Stanford GraphBase (c) Stanford University 1993 */

#ifndef GB_FLIP_H
#define GB_FLIP_H

#define gb_next_rand() (*gb_fptr>=0?*gb_fptr--:gb_flip_cycle())
extern long *gb_fptr; /* the next |A| value to be used */
extern long gb_flip_cycle(void); /* compute 55 more pseudo-random numbers */
extern void gb_init_rand(long);
extern long gb_unif_rand(long);

#endif /* GB_FLIP_H */
