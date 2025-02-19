/* sat9.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gb_flip.h"

typedef unsigned int uint; /* a convenient abbreviation */
typedef unsigned long long ullng; /* ditto */

/* Type definitions */
typedef union {
  char ch8[8];
  uint u2[2];
  long long lng;
} octa;
typedef struct tmp_var_struct {
  octa name; /* the name (one to seven ASCII characters) */
  uint serial; /* 0 for the first variable, 1 for the second, etc. */
  int stamp; /* |m| if positively in clause |m|; |-m| if negatively there */
  struct tmp_var_struct *next; /* pointer for hash list */
} tmp_var;

typedef struct vchunk_struct {
  struct vchunk_struct *prev; /* previous chunk allocated (if any) */
  tmp_var var[341]; /* vars_per_vchunk = 341 */
} vchunk;

typedef struct chunk_struct {
  struct chunk_struct *prev; /* previous chunk allocated (if any) */
  tmp_var *cell[511]; /* cells_per_chunk = 511 */
} chunk;

typedef struct {
  double eta; /* the external force on this literal */
  double pi; /* this literal's current $\pi$ value */
  uint zf; /* the number of suppressed zero factors in |pi| */
  uint link; /* first occurrence of the literal in |mem|, plus 1 */
  int rating; /* $+1$ positive, $-1$ negative, $0$ wishy-washy or wild */
} literal; /* would it go faster if I added four more bytes of padding? */
typedef struct {
  uint start; /* where the literal list starts in |mem| */
  uint size; /* number of remaining literals in clause postprocessing phase */
} clause;
typedef struct {
  union {double d; ullng u;} eta; /* $\eta$ message for a literal */
  uint lit; /* number of that literal */
  uint next; /* where that literal next appears in |mem|, plus 1 */
} mem_item;

/* Global variables */
int random_seed=0; /* seed for the random words of |gb_rand| */
int verbose=1; /* show_basics = 1; level of verbosity */
int hbits=8; /* logarithm of the number of the hash lists */
int buf_size=1024; /* must exceed the length of the longest input line */
int max_iter=1000; /* maximum iterations */
int min_iter=5; /* minimum iterations before reinforcement kicks in */
int confidence=50; /* lower limit for confidence of setting a variable */
double damper=0.99; /* the damping factor for reinforcement */
double threshold=0.01; /* upper limit for convergence check */
ullng imems,mems; /* mem counts */
ullng thresh=0; /* report when |mems| exceeds this, if |delta!=0| */
ullng delta=0; /* report every |delta| or so mems */
ullng bytes; /* memory used by main data structures */
char *buf; /* buffer for reading the lines (clauses) of |stdin| */
tmp_var **hash; /* heads of the hash lists */
uint hash_bits[93][8]; /* random bits for universal hash function */
vchunk *cur_vchunk; /* the vchunk currently being filled */
tmp_var *cur_tmp_var; /* current place to create new |tmp_var| entries */
tmp_var *bad_tmp_var; /* the |cur_tmp_var| when we need a new |vchunk| */
chunk *cur_chunk; /* the chunk currently being filled */
tmp_var **cur_cell; /* current place to create new elements of a clause */
tmp_var **bad_cell; /* the |cur_cell| when we need a new |chunk| */
ullng vars; /* how many distinct variables have we seen? */
ullng clauses; /* how many clauses have we seen? */
ullng nullclauses; /* how many of them were null? */
ullng cells; /* how many occurrences of literals in clauses? */
clause *cmem; /* the master array of clauses */
literal *lmem; /* the master array of literals */
mem_item *mem; /* the master array of literals in clauses */
mem_item *cur_mcell; /* the current cell of interest in |mem| */
octa *nmem; /* the master array of symbolic variable names */
double *gam; /* temporary array to hold gamma ratios */
int iter; /* number of the current iteration */
double acc,etabar,pi0,pi1,old_eta,new_eta,new_gam,factor,rein,diff;
  /* intermediate registers for floating-point calculations */
double max_diff; /* biggest change from |old_eta| to |new_eta| */
int azf; /* number of zero factors suppressed from |acc| */
int bucket[101],unit;
int fixcount,unitcount;
char name_buf[32];
FILE *out_file;


/* Subroutines */
void print_clause(uint c) { /* the first clause is called clause 1, not 0 */
  register uint l,ll;
  fprintf(stderr,"%d:\n",c); /* show the clause number */
  for (l=cmem[c-1].start;l<cmem[c].start;l++) {
    ll=mem[l].lit;
    fprintf(stderr," %s%.8s(%d), eta=%.15g\n",
      ll&1? "~": "",nmem[ll>>1].ch8,ll>>1,mem[l].eta.d);
  }
}

void print_var(uint k) {
  register uint l=k<<1;
  fprintf(stderr,"pi(%.8s)=%.15g(%d), eta(%.8s)=%.15g, ",
      nmem[k].ch8,lmem[l].pi,lmem[l].zf,nmem[k].ch8,lmem[l].eta);
  fprintf(stderr,"pi(~%.8s)=%.15g(%d), eta(~%.8s)=%.15g\n",
      nmem[k].ch8,lmem[l+1].pi,lmem[l+1].zf,nmem[k].ch8,lmem[l+1].eta);
}

#define hack_in(q,t) (tmp_var*)(t|(ullng)q)
#define hack_out(q) (((ullng)q)&0x3)
#define hack_clean(q) ((tmp_var*)((ullng)q&-4))
#define cl(p) mem[p].eta.u /* new use for the |eta| field */
#define removed (uint)(-1)
#define o mems++
#define oo mems+=2
#define ooo mems+=3


int fixlist(register int k, int b);


int main (int argc, char *argv[]) {
  register uint c,g,h,i,j,k,l,p,q,r,ii,kk,ll,fcount;

  /* Process the command line */
  for (j=argc-1,k=0;j;j--) switch (argv[j][0]) {
  case 'v': k|=(sscanf(argv[j]+1,"%d",&verbose)-1);break;
  case 'h': k|=(sscanf(argv[j]+1,"%d",&hbits)-1);break;
  case 'b': k|=(sscanf(argv[j]+1,"%d",&buf_size)-1);break;
  case 's': k|=(sscanf(argv[j]+1,"%d",&random_seed)-1);break;
  case 'd': k|=(sscanf(argv[j]+1,"%lld",&delta)-1);thresh=delta;break;
  case 't': k|=(sscanf(argv[j]+1,"%d",&max_iter)-1);break;
  case 'l': k|=(sscanf(argv[j]+1,"%d",&min_iter)-1);break;
  case 'c': k|=(sscanf(argv[j]+1,"%d",&confidence)-1);break;
  case 'p': k|=(sscanf(argv[j]+1,"%lf",&damper)-1);break;
  case 'e': k|=(sscanf(argv[j]+1,"%lf",&threshold)-1);break;
  default: k=1; /* unrecognized command-line option */
  }
  if (k || hbits<0 || hbits>30 || buf_size<=0) {
    fprintf(stderr,
       "Usage: %s [v<n>] [h<n>] [b<n>] [s<n>] [d<n>] [t<n>] [l<n>] [c<n>] [p<f>] [e<f>]\n",
           argv[0]);
    exit(-1);
  }
  if (damper<0.0 || damper>1.0) {
    fprintf(stderr,"Parameter p should be between 0.0 and 1.0!\n");
    exit(-666);
  }
  if (confidence<0 || confidence>100) {
    fprintf(stderr,"Parameter c should be between 0 and 100!\n");
    exit(-667);
  }

  /* Initialize everything */
  gb_init_rand(random_seed);
  buf=(char*)malloc(buf_size*sizeof(char));
  if (!buf) {
    fprintf(stderr,"Couldn't allocate the input buffer (buf_size=%d)!\n",
              buf_size);
    exit(-2);
  }
  hash=(tmp_var**)malloc(sizeof(tmp_var)<<hbits);
  if (!hash) {
    fprintf(stderr,"Couldn't allocate %d hash list heads (hbits=%d)!\n",
             1<<hbits,hbits);
    exit(-3);
  }
  for (h=0;h<1<<hbits;h++) hash[h]=NULL;
  for (j=92;j;j--) for (k=0;k<8;k++)
    hash_bits[j][k]=gb_next_rand();

  /* Input the clauses */
    vars = 0;
    clauses = 0;
    nullclauses = 0;
    cells = 0;
    cur_vchunk = NULL;
    cur_tmp_var = NULL;
    bad_tmp_var = NULL;

    cur_chunk = NULL;
    cur_cell = NULL;
    bad_cell = NULL;


  while (1) {
    if (!fgets(buf,buf_size,stdin)) break;
    clauses++;
    if (buf[strlen(buf)-1]!='\n') {
      fprintf(stderr,
        "The clause on line %d (%.20s...) is too long for me;\n",clauses,buf);
      fprintf(stderr," my buf_size is only %d!\n",buf_size);
      fprintf(stderr,"Please use the command-line option b<newsize>.\n");
      exit(-4);
    }
    /* Input the clause in |buf| */
    for (j=k=0;;) {
      while (buf[j]==' ') j++; /* scan to nonblank */
      if (buf[j]=='\n') break;
      if (buf[j]<' ' || buf[j]>'~') {
        fprintf(stderr,"Illegal character (code #%x) in the clause on line %d!\n",
          buf[j],clauses);
        exit(-5);
      }
      if (buf[j]=='~') i=1,j++;
      else i=0;
      /* Scan and record a variable; negate it if |i==1| */
      {
        register tmp_var *p;
        if (cur_tmp_var==bad_tmp_var) /* Install a new |vchunk| */
        {
          register vchunk *new_vchunk;
          new_vchunk=(vchunk*)malloc(sizeof(vchunk));
          if (!new_vchunk) {
            fprintf(stderr,"Can't allocate a new vchunk!\n");
            exit(-6);
          }
          new_vchunk->prev=cur_vchunk, cur_vchunk=new_vchunk;
          cur_tmp_var=&new_vchunk->var[0];
          bad_tmp_var=&new_vchunk->var[341];
        }
        /* Put the variable name beginning at |buf[j]| in |cur_tmp_var->name|
           and compute its hash code |h| */
        cur_tmp_var->name.lng=0;
        for (h=l=0;buf[j+l]>' '&&buf[j+l]<='~';l++) {
          if (l>7) {
            fprintf(stderr,
                "Variable name %.9s... in the clause on line %d is too long!\n",
                buf+j,clauses);
            exit(-8);
          }
          h^=hash_bits[buf[j+l]-'!'][l];
          cur_tmp_var->name.ch8[l]=buf[j+l];
        }
        if (l==0) goto empty_clause; /* `\.\~' by itself is like `true' */
        j+=l;
        h&=(1<<hbits)-1;
        /* Find |cur_tmp_var->name| in the hash table at |p| */
        for (p=hash[h];p;p=p->next)
          if (p->name.lng==cur_tmp_var->name.lng) break;
        if (!p) { /* new variable found */
          p=cur_tmp_var++;
          p->next=hash[h], hash[h]=p;
          p->serial=vars++;
          p->stamp=0;
        }
        if (p->stamp==clauses || p->stamp==-clauses) /* Handle a duplicate literal */
        {
          if ((p->stamp>0)==(i>0)) goto empty_clause;
        }
        else {
          p->stamp=(i? -clauses: clauses);
          if (cur_cell==bad_cell) /* Install a new |chunk| */
          {
            register chunk *new_chunk;
            new_chunk=(chunk*)malloc(sizeof(chunk));
            if (!new_chunk) {
              fprintf(stderr,"Can't allocate a new chunk!\n");
              exit(-7);
            }
            new_chunk->prev=cur_chunk, cur_chunk=new_chunk;
            cur_cell=&new_chunk->cell[0];
            bad_cell=&new_chunk->cell[511];
          }
          *cur_cell=p;
          if (i==1) *cur_cell=hack_in(*cur_cell,1);
          if (k==0) *cur_cell=hack_in(*cur_cell,2);
          cur_cell++,k++;
        }
      }
    }
    if (k==0) {
      fprintf(stderr,"(Empty line %d is being ignored)\n",clauses);
      nullclauses++; /* strictly speaking it would be unsatisfiable */
    }
    goto clause_done;
    empty_clause: /* Remove all variables of the current clause */
    while (k) {
      /* Move |cur_cell| backward to the previous cell */
      if (cur_cell>&cur_chunk->cell[0]) cur_cell--;
      else {
        register chunk *old_chunk=cur_chunk;
        cur_chunk=old_chunk->prev;free(old_chunk);
        bad_cell=&cur_chunk->cell[511];
        cur_cell=bad_cell-1;
      }
      k--;
    }
    if ((buf[0]!='~')||(buf[1]!=' '))
      fprintf(stderr,"(The clause on line %d is always satisfied)\n",clauses);
    nullclauses++;
    clause_done: cells+=k;
  }
  if ((vars>>hbits)>=10) {
    fprintf(stderr,"There are %d variables but only %d hash tables;\n",
       vars,1<<hbits);
    while ((vars>>hbits)>=10) hbits++;
    fprintf(stderr," maybe you should use command-line option h%d?\n",hbits);
  }
  clauses-=nullclauses;
  if (clauses==0) {
    fprintf(stderr,"No clauses were input!\n");
    exit(-77);
  }
  if (vars>=0x80000000) {
    fprintf(stderr,"Whoa, the input had %llu variables!\n",vars);
    exit(-664);
  }
  if (clauses>=0x80000000) {
    fprintf(stderr,"Whoa, the input had %llu clauses!\n",clauses);
    exit(-665);
  }
  if (cells>=0x100000000) {
    fprintf(stderr,"Whoa, the input had %llu occurrences of literals!\n",cells);
    exit(-666);
  }

  if (verbose&1) /* show_basics */
    /* Report the successful completion of the input phase */
    fprintf(stderr,"(%d variables, %d clauses, %llu literals successfully read)\n",
                           vars,clauses,cells);
  /* Set up the main data structures */
  /* Allocate the main arrays */
  free(buf);free(hash); /* a tiny gesture to make a little room */
  lmem=(literal*)malloc((vars+vars+1)*sizeof(literal));
  if (!lmem) {
    fprintf(stderr,"Oops, I can't allocate the lmem array!\n");
    exit(-12);
  }
  bytes=(vars+vars+1)*sizeof(literal);
  nmem=(octa*)malloc(vars*sizeof(octa));
  if (!nmem) {
    fprintf(stderr,"Oops, I can't allocate the nmem array!\n");
    exit(-13);
  }
  bytes+=vars*sizeof(octa);
  mem=(mem_item*)malloc(cells*sizeof(mem_item));
  if (!mem) {
    fprintf(stderr,"Oops, I can't allocate the big mem array!\n");
    exit(-10);
  }
  bytes+=cells*sizeof(mem_item);
  cmem=(clause*)malloc((clauses+1)*sizeof(clause));
  if (!cmem) {
    fprintf(stderr,"Oops, I can't allocate the cmem array!\n");
    exit(-11);
  }
  bytes+=(clauses+1)*sizeof(clause);
  /* Zero the links */
  for (l=vars+vars; l; l--) o,lmem[l-1].link=0;
  /* Copy all the temporary cells to the |mem| and |cmem| arrays
     in proper format */
  for (c=clauses,cur_mcell=mem+cells,kk=0; c; c--) {
    o,cmem[c].start=cur_mcell-mem;
    k=0;
    /* Insert the cells for the literals of clause |c| */
    for (i=0;i<2;k++) {
      /* Move |cur_cell| back... */
      if (cur_cell>&cur_chunk->cell[0]) cur_cell--;
      else {
        register chunk *old_chunk=cur_chunk;
        cur_chunk=old_chunk->prev;free(old_chunk);
        bad_cell=&cur_chunk->cell[511];
        cur_cell=bad_cell-1;
      }
      i=hack_out(*cur_cell);
      p=hack_clean(*cur_cell)->serial;
      cur_mcell--;
      o,cur_mcell->lit=l=p+p+(i&1);
      oo,cur_mcell->next=lmem[l].link;
      o,lmem[l].link=cur_mcell-mem+1;
    }
    if (k>kk) kk=k; /* maximum clause size seen so far */
  }
  if (cur_mcell!=mem) {
    fprintf(stderr,"Confusion about the number of cells!\n");
    exit(-99);
  }
  o,cmem[0].start=0;
  gam=(double*)malloc(kk*sizeof(double));
  if (!gam) {
    fprintf(stderr,"Oops, I can't allocate the gamma array!\n");
    exit(-16);
  }
  bytes+=kk*sizeof(double);
  /* Copy all the temporary variable nodes to the |nmem| array in proper format */
  for (c=vars; c; c--) {
    /* Move |cur_tmp_var| back... */
    if (cur_tmp_var>&cur_vchunk->var[0]) cur_tmp_var--;
    else {
      register vchunk *old_vchunk=cur_vchunk;
      cur_vchunk=old_vchunk->prev;free(old_vchunk);
      bad_tmp_var=&cur_vchunk->var[341];
      cur_tmp_var=bad_tmp_var-1;
    }
    o,nmem[c-1].lng=cur_tmp_var->name.lng;
  }
  /* Check consistency */
  if (cur_cell!=&cur_chunk->cell[0] ||
       cur_chunk->prev!=NULL ||
       cur_tmp_var!=&cur_vchunk->var[0] ||
       cur_vchunk->prev!=NULL) {
    fprintf(stderr,"This can't happen (consistency check failure)!\n");
    exit(-14);
  }
  free(cur_chunk);free(cur_vchunk);
  imems=mems, mems=0;

  /* Solve the problem */
  factor=1.0;
  /* Initialize all $\eta$'s to random fractions */
  for (k=0;k<cells;k++)
    mems+=5,mem[k].eta.d=((double)(gb_next_rand()))/2147483647.0;
  for (k=0;k<vars+vars;k+=2)
    ooo,lmem[k].eta=0.0,lmem[k+1].eta=0.0;
  for (iter=0;iter<max_iter;iter++) {
    if ((iter&1) && iter>=min_iter) {
      /* Adjust the reinforcement fields */
      {
        factor*=damper;
        rein=1.0-factor;
        if (verbose&4) /* show_details */
          fprintf(stderr," (rein=%.15g)\n",rein);
        for (l=0;l<vars+vars;l+=2) {
          if (o,lmem[l].zf) pi0=0.0;
          else o,pi0=lmem[l].pi;
          if (o,lmem[l+1].zf) pi1=0.0;
          else o,pi1=lmem[l+1].pi;
          if (pi0+pi1==0.0) {
            if (verbose&1) /* show_basics */
              fprintf(stderr,
                "Sorry, a contradiction was found after iteration %d!\n",iter);
            goto contradiction;
          }
          if (pi1>pi0) {
            o,lmem[l].rating=(pi0>=0.5? 0: 1);
            if ((verbose&8) && lmem[l+1].eta) /* show_gory_details */
              fprintf(stderr," eta(~%.8s) reset\n",nmem[l>>1].ch8);
            oo,lmem[l].eta=rein*(pi1-pi0)/(pi0+pi1-pi0*pi1),lmem[l+1].eta=0.0;
          }else {
            o,lmem[l].rating=(pi1>=0.5? 0: -1);
            if ((verbose&8) && lmem[l].eta)  /* show_gory_details */
              fprintf(stderr," eta(%.8s) reset\n",nmem[l>>1].ch8);
            oo,lmem[l+1].eta=rein*(pi0-pi1)/(pi0+pi1-pi0*pi1),lmem[l].eta=0.0;
          }
        }
      }
      /* Exit if the clauses are pseudo-satisfied */
      for (k=c=0;c<clauses;c++) {
        for (o;k<cmem[c+1].start;k++) {
          oo,l=mem[k].lit, p=lmem[l&-2].rating;
          if (p==0) goto ok;
          if (((int)p<0)==(l&1)) goto ok;
        }
        goto not_ok; /* clause not pseudo-satisfied */
      ok: k=cmem[c+1].start;
        continue;
      }
      if (verbose&4) /* show_details */
        fprintf(stderr,"Clauses pseudo-satisfied on iteration %d\n",iter+1);
      break; /* yes, we made it through all of them */
      not_ok:;
    }
    if (verbose&2) /* show_choices */
      fprintf(stderr,"beginning iteration %d\n",iter+1);
    /* Compute the $\pi$'s */
    for (l=0;l<vars+vars;l++) {
      if (o,lmem[l].eta==1.0) acc=1.0,azf=1;
      else acc=1.0-lmem[l].eta,azf=0;
      for (j=lmem[l].link;j;j=mem[j-1].next) {
        o,etabar=1.0-mem[j-1].eta.d;
        if (etabar==0.0) azf++;
        else acc*=etabar;
      }
      oo,lmem[l].zf=azf, lmem[l].pi=acc;
    }
    /* Update the $\eta$'s */
    max_diff=0.0;
    for (k=c=0;c<clauses;c++) {
      acc=1.0, azf=0;
      for (o,j=0;k<cmem[c+1].start;j++,k++) {
        o,l=mem[k].lit;
        if (o,lmem[l^1].zf) pi0=0.0;
        else o,pi0=lmem[l^1].pi;
        o,old_eta=mem[k].eta.d;
        if (old_eta==1.0) {
          if (o,lmem[l].zf>1) pi1=0.0;
          else o,pi1=lmem[l].pi;
        }else if (o,lmem[l].zf) pi1=0.0;
          else o,pi1=lmem[l].pi/(1.0-old_eta);
        pi1=pi1*(1.0-pi0);
        if (pi1==0.0) azf++,o,gam[j]=0.0;
        else {
          new_gam=pi1/(pi1+pi0);
          o,gam[j]=new_gam;
          acc*=new_gam;
        }
      }
      for (i=j;i;i--) {
        if (o,gam[j-i]==0.0) {
          if (azf>1) new_eta=0.0;
          else new_eta=acc;
        }else if (azf) new_eta=0.0;
        else new_eta=acc/gam[j-i];
        o,diff=new_eta-mem[k-i].eta.d;
        if (diff>0) {
          if (diff>max_diff) max_diff=diff;
        }else if (-diff>max_diff) max_diff=-diff;
        o,mem[k-i].eta.d=new_eta;
      }
    }
    if (verbose&4)  /* show_details */
      fprintf(stderr," (max diff %.15g, %lld mems)\n",max_diff,mems);
    if (delta && (mems>=thresh)) {
      thresh+=delta;
      fprintf(stderr," after %lld mems, iteration %d had max diff %g\n",
             mems,iter+1,max_diff);
    }
    if (max_diff<threshold && iter>=min_iter) break;
  }
  /* Output a reduced problem */
  if (iter==max_iter) {
    if (verbose&1) /* show_basics */
      fprintf(stderr,"The messages didn't converge.\n");
    goto contradiction;
  }
  if (verbose&32) /* show_pis */
    /* Print all the $\pi$'s */
    {
      if (iter<max_iter)
        fprintf(stderr,"converged after %d iterations.\n",iter+1);
      else fprintf(stderr,"no convergence (diff %g) after %d iterations.\n",
                 max_diff,max_iter);
      fprintf(stderr,"variable      pi(v)        pi(~v)         1    0    *\n");
      for (k=0;k<vars;k++) {
        double den;
        fprintf(stderr,"%8.8s %10.7f(%d) %10.7f(%d)",
             nmem[k].ch8,lmem[k+k].pi,lmem[k+k].zf,lmem[k+k+1].pi,lmem[k+k+1].zf);
        pi0=lmem[k+k].pi;
        if (lmem[k+k].zf) pi0=0.0;
        pi1=lmem[k+k+1].pi;
        if (lmem[k+k+1].zf) pi1=0.0;
        den=pi0+pi1-pi0*pi1;
        fprintf(stderr,"    %4.2f %4.2f %4.2f\n",
              pi1*(1-pi0)/den,pi0*(1-pi1)/den,pi0*pi1/den);
      }
    };
  if (verbose&16) /* show_histogram */
    /* Print a two-dimension histogram of $\pi_v$ versus $\pi_{\bar v}$ */
    {
      uint hist[10][10];
      for (j=0;j<10;j++) for (k=0;k<10;k++) hist[j][k]=0;
      for (k=0;k<vars;k++) {
        i=(int)(10*lmem[k+k].pi), j=(int)(10*lmem[k+k+1].pi);
        if (lmem[k+k].zf) i=0;
        if (lmem[k+k+1].zf) j=0;
        if (i==10) i=9;
        if (j==10) j=9;
        hist[i][j]++;
      }
      fprintf(stderr,"Histogram of the pi's, after %d iterations:\n",iter+1);
      for (j=10;j;j--) {
        for (i=0;i<10;i++) fprintf(stderr,"%7d",hist[i][j-1]);
        fprintf(stderr,"\n");
      }
    };
  /* Decide which variables to fix */
  for (k=confidence;k<=100;k++) o,bucket[k]=2;
  unit=2;
    fixcount = 0;
    unitcount = 0;
  for (l=0;l<vars+vars;l+=2) {
    if (o,lmem[l].zf) pi0=0.0;
    else o,pi0=lmem[l].pi;
    if (o,lmem[l+1].zf) pi1=0.0;
    else o,pi1=lmem[l+1].pi;
    if (pi0+pi1==0.0) {
      if (verbose&1)/*show_basics*/
        fprintf(stderr,"Sorry, a contradiction was found!\n");
      goto contradiction;
    }
    acc=(pi1-pi0)/(pi0+pi1-pi0*pi1);
    o,lmem[l].rating=acc>0?+1:acc<0?-1:0;
    if (acc<0) acc=-acc;
    j=(int)(100.0*acc);
    if (j>=confidence) {
      oo,lmem[l+1].rating=bucket[j];
      o,bucket[j]=l+1;
      fixcount++;
    }
  }
  if (verbose&1)/*show_basics*/
    fprintf(stderr,"(fixing %d variables after %d iterations, e=%g)\n",
      fixcount,iter+1,max_diff);
    /* Preprocess the clauses for reduction */
  for (k=c=0;c<clauses;c++) {
    for (;k<cmem[c+1].start;k++) o,cl(k)=c;
    oo,cmem[c].size=k-cmem[c].start;
    if (cmem[c].size==1) {
      /* Enforce the unit literal |mem[k-1].lit| */
        ll=mem[k-1].lit;
        if (ll&1) {
          if (o,lmem[ll].rating) {
            if (o,lmem[ll-1].rating>0) goto contra;
          }else {
            o,lmem[ll-1].rating=-1;
            o,lmem[ll].rating=unit, unit=ll, unitcount++;
          }
        }else {
          if (o,lmem[ll+1].rating) {
            if (o,lmem[ll].rating<0) {
        contra: printf("\n");
            fprintf(stderr,"Oops, clause %d is contradicted!\n",c);
            goto contradiction;  // Corrected: Jump to main contradiction handling
            }
          }else {
            o,lmem[ll].rating=+1;
            o,lmem[ll+1].rating=unit, unit=ll+1, unitcount++;
          }
        }
    }
  }

  /* Reduce the problem */
  for (k=100;k>=confidence;k--)
    if (ooo,fixlist(bucket[k],k)==0) goto contradiction;
  while (unit&1) {
    p=unit, unit=2;
    if (oo,fixlist(p,-1)==0) goto contradiction;
  }
  printf("\n");
  if (unitcount && (verbose&1))/*show_basics*/
    fprintf(stderr,"(unit propagation fixed %d more variable%s)\n",
      unitcount,unitcount==1?"":"s");
  /* Output the reduced problem */
  system("mkdir -p output");
  sprintf(name_buf,"output/sat9-%d.dat",random_seed);
  out_file=fopen(name_buf,"w");
  if (!out_file) {
    fprintf(stderr,"I can't open `%s' for writing!\n", name_buf);
    exit(-668);
  }
  for (kk=k=p=c=0;c<clauses;c++) {
    o,i=cmem[c].size;
    if (i==0) {
      o,k=cmem[c+1].start;
      continue;
    }
    p++;
    while (i>kk) gam[kk++]=0;
    gam[i-1]+=1;
    for (o;k<cmem[c+1].start;k++) if (o,mem[k].next!=removed) {
      l=mem[k].lit;
      fprintf(out_file," %s%.8s",l&1?"~":"",nmem[l>>1].ch8);
    }
    fprintf(out_file,"\n");
  }
  fclose(out_file);
  fprintf(stderr,"Reduced problem of %d clauses written on file %s\n",
      p,name_buf);
  for (i=0;i<kk;i++) if (gam[i])
    fprintf(stderr," (%g %d-clauses)\n",gam[i],i+1);
  goto done;

  contradiction: printf("~~?\n");
  done:

  if (verbose&1)/*show_basics*/
    fprintf(stderr,"Altogether %llu+%llu mems, %llu bytes.\n",
               imems,mems,bytes);
  return 0;
}


int fixlist(register int k, int b) {
  register int c,j,l,ll,p,q;
  for (;k&1;o,k=lmem[k].rating) {
    if (o,lmem[k-1].rating<0) l=k;
    else l=k-1;
    printf(" %s%.8s",l&1?"~":"",nmem[l>>1].ch8);
    /* Mark the clauses that contain |l| satisfied */
    for (o,p=lmem[l].link; p; o,p=mem[p-1].next) {
      oo,c=cl(p-1), j=cmem[c].size;
      if (j) o,cmem[c].size=0;
    }
    /* Remove $\bar l$ from all clauses */
    for (o,p=lmem[l^1].link; p; p=q) {
      o,q=mem[p-1].next;
      oo,c=cl(p-1), j=cmem[c].size;
      if (j==0) continue; /* clause already satisfied */
      oo,mem[p-1].next=removed, cmem[c].size=j-1;
      if (j==2) {
        for (o,p=cmem[c].start;o,mem[p].next==removed;p++);
        /* Enforce the unit literal |mem[p].lit| */
        ll=mem[p].lit;
        if (ll&1) {
          if (o,lmem[ll].rating) {
            if (o,lmem[ll-1].rating>0) goto contra;
          }else {
            o,lmem[ll-1].rating=-1;
            o,lmem[ll].rating=unit, unit=ll, unitcount++;
          }
        }else {
          if (o,lmem[ll+1].rating) {
            if (o,lmem[ll].rating<0) {
        contra: printf("\n");
            fprintf(stderr,"Oops, clause %d is contradicted",c);
            if (b>=0) fprintf(stderr," in bucket %d!\n",b);
            else fprintf(stderr," while propagating unit literals!\n");
            return 0; // Correct return from fixlist
            }
          }else {
            o,lmem[ll].rating=+1;
            o,lmem[ll+1].rating=unit, unit=ll+1, unitcount++;
          }
        }
      }
    }
  }
  return 1;
}
