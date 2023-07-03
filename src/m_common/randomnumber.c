/***************************************************************************
* Lagged Fibonacci pseudo random-number generator

This code implements a lagged Fibonnacci pseudo-random number generator
with lags p and q in the notation of "A current view of random number generators"
published by G.A. Marsaglia in Computational Science and Statistics: The Interface,
ed. L. Balliard (Elsevier, Amsterdam, 1985).

The implementation also relies on the following references:

* P.D. Coddington, "Monte Carlo tests of random number generators", NPAC technical report SCCS-526.
* P.D. Coddington, "Random number generators for parallel computers." (1997), Syracuse University, SURFACE.

The followings are comments by Paul Coddington in October 1993.

Background
--------------------

This program implements a lagged Fibonnacci pseudo-random number generator
using multiplication, with lags p and q, i.e. F(p,q,*) in the notation of
[Marsaglia].

A sequence of odd random integers Xn is obtained from the formula

Xn= Xn-p * Xn-q mod M

where M is a large integer less than the machine precision (M is taken
to be 2^{31} here).

Unlike lagged Fibonnacci generators which use subtraction, addition, or
XOR, this generator seems to give "good" random numbers for ANY lag, not
necessarily large (see [Coddington]). However larger lags give larger
periods and empirically seem to give "better" random numbers (see
[Marsaglia], [Knuth], [Coddington]).

The period is given by (2^{p}-1)*2^{b-3}, where p is the largest lag
and b is the number of bits in the integers X_{n}, i.e. b = log_{2}(M).
The lag needs to be chosen large enough to give a period much longer
than the number of random numbers to be generated.

Only certain lags (p,q) will give the maximal period of the generator
(see [Marsaglia], [Knuth]). Here are some of those lags, with the
corresponding periods for 32-bit integers. As a comparison, note that
linear congruential generators have a period of about 10^{9} (which is
too short for most purposes), or 10^{18} for 64-bit integers (which is
fine for most purposes). A Teraflop-year is about 10^{21} floating point
operations.

P     Q        period
1279   418       10^{393}
607   273       10^{191}
521   168       10^{166}
250   103       10^{84}
127    63       10^{46}
97    33       10^{38}
89    38       10^{35}
55    24       10^{25}
43    22       10^{21}
31    13       10^{18}
24    10       10^{16}
17     5       10^{14}
7     3       10^{11}
5     2       10^{10}

I would recommend using at least (43,22), and ideally (607,273). This
will only add an extra 607x4 = 2428 bytes to the memory required for
the user program, which should be completely negligible.

Note that very little is known theoretically about this generator,
however it performs very well in all empirical statistical tests (see
[Marsaglia], [Coddington]).
THIS DOES NOT NECESSARILY MEAN IT WILL PERFORM WELL FOR YOUR APPLICATION.
It is a good idea to check by using another good generator such as RANF
or DRAND48 or a good 64-bit generator (see [Coddington]).

Program
--------------------

This program is based on the implementation of RANMAR in [James].
It is the same as RANMAR except it uses multiplication instead of
subtraction, odd integers instead of reals, and does not include an
extra Weyl generator.

The program accesses the "lag table", i.e. the array storing the previous
p values in the sequence, using a "circular list", i.e. the elements of
the array are not moved with every new addition, only the pointers pt0
and pt1 to the (n-p)th and (n-q)th value are changed.

NOTE THAT YOU *MUST* CALL THE INITIALIZATION SUBROUTINE (rand_init) BEFORE
USING THE GENERATOR.

In order for the sequence to be random, the initial values of the lag
table have to be "random enough". There have not to my knowledge been any
quantitive tests of this point, however initializing the lag table using
a good linear congruential generator (the Unix/C generator RAND is used
here) seems empirically to work fine.

The computation required to generate a new random number, apart from
updating 2 pointers, is 1 integer multiplication, which multiplies two
32-bit integers and returns the answer modulo 2^{32} (this is the standard
for unsigned ints in C), and 1 floating point multiplication, to convert
the integers into a floating point number in [0,1).

For simplicity and clarity the lag table is indexed from 1 to P rather
than 0 to P-1.

References
--------------------
P.D. Coddington, "Monte Carlo tests of random number generators",
NPAC technical report SCCS-526.

F. James, "A Review of Pseudo-random Number Generators",
Comput. Phys. Comm. 60 (1990) 329.

D.E. Knuth, "The Art of Computer Programming Vol. 2: Seminumerical
Methods", (Addison-Wesley, Reading, Mass., 1981).

P. L'Ecuyer, "Random numbers for simulation", Comm. ACM 33:10 (1990) 85.

G.A. Marsaglia, "A current view of random number generators",
in Computational Science and Statistics: The Interface,
ed. L. Balliard (Elsevier, Amsterdam, 1985).

****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define TWOTO31 2147483648.0
#define TWOTO30 1073741824 
double  TWOTONEG31=(1.0/TWOTO31);

/* Lags, P and Q */
#define P 1279
#define Q 418

/* Variables for the lagged Fibonacci generator */
unsigned int u[P+1];
int pointer0, pointer1;

/* Variables for the linear congruential generator */
#define M ( 0x7fffffff )     /* 2^{31}-1 -- mask by M to get modulo 2^{31} */
#define	A 1103515245
#define	K 12345

/***************************************************************************
   Initialize the lagged Fibonacci random number generator.
   
   Every bit of all the seeds in the lag table is initialized using the
   linear congruential generator.
   
   s = (A * s + K) mod M, where A = 1103515245, K = 12345, and M = 2^{31}-1
   
   This routine must be called before the random number generator is called.
***************************************************************************/

/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
  void rand_init(seed)
#else
  void rand_init_(seed)
#endif
int * seed;
{
   int i,j;
   int s,t,y;

   y = *seed;

   /* Initialize every bit of the lag table using RAND */
   for (i=1;i<=P;i++) {
      s = 0;
      t = 1;
      for (j=0;j<32;j++) {
      y = (A * y + K) & M;
         if ( y < TWOTO30 ) {
            s = s + t;
         }
      t = t << 1;
      }
      /* Seeds must be odd integers because the last bit should be 1 */
      s = s | 1;
      u[i] = s;
   }

   pointer0 = P;
   pointer1 = Q;

   return;
}

/***************************************************************************
   Return a random **double** in [0.0, 1.0).
   
   !!!Caution!!!
   Only the most significant 31 bits of drand1 will actually
   be random, since a random 32-bit integer is converted into a double.

   To return a float instead of a double, simply change the type of rand
   and drand1 to float.
***************************************************************************/

/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC  
  double drand1()
#else
  double drand1_()
#endif
{
   unsigned int unsi;
   double rand;

   unsi = u[pointer0] * u[pointer1];
   u[pointer0] = unsi;
   pointer0 --;
   if (pointer0 == 0) pointer0 = P;
   pointer1 --;
   if (pointer1 == 0) pointer1 = P;

   unsi = unsi >> 1;

   rand = (double)unsi * TWOTONEG31;
   return(rand);
}

/***************************************************************************
   Store the lag table and pointers.
   They are needed for restarting MCMC.
***************************************************************************/
/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
int write_seeds(seedfilename)
#else
int write_seeds_(seedfilename)
#endif 

char seedfilename[128];
{
   FILE *f1;
   int i,pointer[2];

   pointer[0] = pointer0;
   pointer[1] = pointer1;

   /* Open a file to store the lag table and pointers */
   f1 = fopen(seedfilename,"w") ;
   for (i=1;i<=P;i++) {
      fprintf(f1,"%u\n",u[i]);
   }
   fflush(f1);
   fprintf(f1,"%d %d\n",pointer[0],pointer[1]);
   fclose(f1);

   return 0;
}

/***************************************************************************
   Read in stored lag table and pointers.
   They are needed for restarting MCMC.
***************************************************************************/

/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
int read_seeds(seedfilename)
#else
int read_seeds_(seedfilename)
#endif
char seedfilename[128];
{
  FILE *f1;
  int i,pointer[2];

  /* Open the file containing the lag table and pointers */
  f1 = fopen(seedfilename,"r");
  if( f1 == NULL ) {
      fprintf( stderr, "WARNING error opening file: %s\n",seedfilename);
      return 1;
  } else {
    for(i=1;i<=P;i++)
      if(fscanf(f1,"%u",&u[i])<0) return 1;
      fprintf( stdout, "WAR1");
    if(fscanf(f1,"%d %d",&pointer[0],&pointer[1])<0) return 1;
      fprintf( stdout, "WAR2");
    pointer0 = pointer[0];
    pointer1 = pointer[1];
  }

  fclose(f1);

  return 0;
}

/***************************************************************************
   Store the lag table and pointers with MPI-IO.
   They are needed for restarting MCMC.
***************************************************************************/
#ifdef PARALLEL
#include <mpi.h>
/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
    int write_seeds_mpiio(seedfilename, rank_address)
#else
    int write_seeds_mpiio_(seedfilename, rank_address)
#endif 
char seedfilename[128];
int* rank_address;
{
  MPI_File fh;
  MPI_Offset offset1, offset2;
  MPI_Info info;
  MPI_Status status;
  int rank=*rank_address;
  int root=(rank==0);
  int result;
  unsigned int pointer[2]={pointer0,pointer1};

  if(root) fprintf(stdout, " Writing random seeds to %s with MPI-IO ...\n", seedfilename);

  /* Open a file to store the lag table and pointers */
  info = MPI_INFO_NULL;
  result = MPI_File_open(MPI_COMM_WORLD, seedfilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fh);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to open file %s for writing random seeds with MPI-IO!\n", seedfilename);
    return 1;
  }

  MPI_File_set_size(fh, 0);
  result = MPI_File_set_view(fh, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", info);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail in MPI_File_set_view!\n");
    return 1;
  }

  offset1 = (MPI_Offset) rank*(P+2);
  offset2 = (MPI_Offset) rank*(P+2)+P;
  // file size is rank*(P+2)*4 bytes.

  //fprintf(stdout,"rank %d, offset %d %d\n",rank,offset1,offset2);
  result = MPI_File_write_at_all(fh, offset1, &u[1]  , P, MPI_UNSIGNED, &status);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to write vector u to %s!\n", seedfilename);
    return 1;
  }

  result = MPI_File_write_at_all(fh, offset2, pointer, 2, MPI_UNSIGNED, &status);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to write pointer0 and pointer1 to %s!\n", seedfilename);
    return 1;
  }

  result = MPI_File_close(&fh);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to close file %s!\n", seedfilename);
    return 1;
  }

  return 0;
}

/***************************************************************************
  Read in stored lag table and pointers with MPI-IO.
  They are needed for restarting MCMC.
***************************************************************************/
/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
   int read_seeds_mpiio(seedfilename,rank_address)
#else
   int read_seeds_mpiio_(seedfilename,rank_address)
#endif
char seedfilename[128];
int* rank_address;
{
  MPI_File fh;
  MPI_Offset offset1, offset2;
  MPI_Info info;
  MPI_Status status;
  int rank=*rank_address;
  int root=(rank==0);
  int result;
  unsigned int pointer[2];

  if(root) fprintf(stdout, " Reading random seeds from %s with MPI-IO ...\n", seedfilename);

  /* Open a file to store the lag table and pointers */
  info = MPI_INFO_NULL;
  result = MPI_File_open(MPI_COMM_WORLD, seedfilename, MPI_MODE_RDONLY, info, &fh);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to open file %s for reading random seeds with MPI-IO!\n", seedfilename);
    return 1;
  }

  result = MPI_File_set_view(fh, 0, MPI_UNSIGNED, MPI_UNSIGNED, "native", info);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail in MPI_File_set_view!\n");
    return 1;
  }

  offset1 = (MPI_Offset) rank*(P+2);
  offset2 = (MPI_Offset) rank*(P+2)+P;

  //fprintf(stdout,"rank %d, offset %d %d\n",rank,offset1,offset2);
  result = MPI_File_read_at_all(fh, offset1, &u[1]  , P, MPI_UNSIGNED, &status);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to read vector u to %s!\n", seedfilename);
    return 1;
  }

  result = MPI_File_read_at_all(fh, offset2, pointer, 2, MPI_UNSIGNED, &status);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to read pointer0 and pointer1 to %s!\n", seedfilename);
    return 1;
  }
  pointer0=pointer[0]; pointer1=pointer[1];

  result = MPI_File_close(&fh);
  if(result != MPI_SUCCESS){
    if(root) fprintf(stderr, "ERROR: Fail to close file %s!\n", seedfilename);
    return 1;
  }

  return 0;
}

// end of PARALLEL for reading and writing MPI-IO random seeds.
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

/* avoid possible conflict, #ifdef RISC preprocessor */
#ifdef RISC
double cclock()
#else
double cclock_()
#endif
/* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970)
*/
{
    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}
