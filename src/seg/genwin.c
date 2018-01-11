
/*****************************************************************************/
/***   (genwin.c)                                                          ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"

/*---------------------------------------------------------------(defines)---*/

#define STRSIZE 100

/*----------------------------------------------------------------(protos)---*/

#ifndef MIN
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#endif

int nabets;
struct Alphabet **abets;
int ntvecs;
struct TransVector **tvecs;
int nsvecs;
struct ScoreVector **svecs;
int nsmats;
struct ScoreMatrix **smats;

int aaindex[128];
unsigned char	aaflag[128];
char aachar[20];

struct strlist
  {
   char string[STRSIZE];
   struct strlist *next;
  } *str, *curstr;

/*---------------------------------------------------------------(tmalloc)---*/

#include <stdint.h>

#define TESTMAX 1000
void *tmalloc();
//int record_ptrs[TESTMAX] = {0,0,0,0};
int64_t record_ptrs[TESTMAX] = {0,0,0,0};
int rptr = 0;

/*------------------------------------------------------------(genwininit)---*/

void genwininit()
{
	char	*cp, *cp0;
	int		i;
	char	c;

	for (i = 0; i < sizeof(aaindex)/sizeof(aaindex[0]); ++i) {
		aaindex[i] = 20;
		aaflag[i] = TRUE;
	}
		
	for (cp = cp0 = "ACDEFGHIKLMNPQRSTVWY"; (c = *cp) != '\0'; ++cp) {
		i = cp - cp0;
		aaindex[(int)c] = i;
		aaindex[tolower(c)] = i;
		aachar[i] = tolower(c);
		aaflag[(int)c] = FALSE;
		aaflag[tolower(c)] = FALSE;
	}
	return;
}
      
/*--------------------------------------------------------------(closeseq)---*/

extern void closeseq(struct Sequence *seq)
  {
   if (seq==NULL) return;

   if (seq->id!=NULL)          free(seq->id);
   if (seq->name!=NULL)        free(seq->name);
   if (seq->organism!=NULL)    free(seq->organism);
   if (seq->header!=NULL)      free(seq->header);
   if (seq->state!=NULL)       free(seq->state);
   if (seq->composition!=NULL) free(seq->composition);
   free(seq->seq);
   free(seq);
   return;
  }

/*---------------------------------------------------------------(openwin)---*/

extern struct Sequence *openwin(struct Sequence *parent, int start, int length)
  {struct Sequence *win;
   int i;

   if (start<0 || length<0 || start+length>parent->length)
     {
      return((struct Sequence *) NULL);
     }

   win = (struct Sequence *) malloc(sizeof(struct Sequence));

/*---                                          ---[set links, up and down]---*/

   win->parent = parent;
   if (parent->root==NULL)
     {win->root = parent;}
   else
     {win->root = parent->root;}
   win->children = (struct Sequence **) NULL;

/* parent->children = ***foo***                   ---[not yet implemented]---*/

   win->id = (char *) NULL;
   win->name = (char *) NULL;
   win->organism = (char *) NULL;
   win->header = (char *) NULL;

/*---                          ---[install the local copy of the sequence]---*/

   win->start = start;
   win->length = length;
#if 0
   win->seq = (char *) malloc(sizeof(char)*length + 1);
   memcpy(win->seq, (parent->seq)+start, length);
   win->seq[length] = '\0';
#else
	win->seq = parent->seq + start;
#endif

/*---                          ---[setup window implementation parameters]---*/

/*---                                                 ---[set local flags]---*/

	win->rubberwin = FALSE;
	win->floatwin = FALSE;
	win->punctuation = FALSE;

/*---                                   ---[initially unconfiguerd window]---*/

	win->entropy = -2.;
	win->state = (int *) NULL;
	win->composition = (int *) NULL;
	win->classvec = (char *) NULL;
	win->scorevec = (double *) NULL;

	stateon(win);

	return win;
}

/*---------------------------------------------------------------(nextwin)---*/

extern struct Sequence *nextwin(struct Sequence *win, int shift)
  {
   if ((win->start+shift)<0 ||
       (win->start+win->length+shift)>win->parent->length)
     {
      return((struct Sequence *) NULL);
     }
   else
     {
      return(openwin(win->parent, win->start+shift, win->length));
     }
  }

/*--------------------------------------------------------------(shiftwin1)---*/
static void	decrementsv(), incrementsv();

extern int shiftwin1(struct Sequence *win)
{
	register int	j, length;
	register int	*comp;

	length = win->length;
	comp = win->composition;

	if ((++win->start + length) > win->parent->length) {
		--win->start;
		return FALSE;
	}

	if (!aaflag[j = win->seq[0]])
		decrementsv(win->state, comp[aaindex[j]]--);

	j = win->seq[length];
	++win->seq;

	if (!aaflag[j])
		incrementsv(win->state, comp[aaindex[j]]++);

	if (win->entropy > -2.)
		win->entropy = entropy(win->state);

	return TRUE;
}

/*--------------------------------------------------------------(closewin)---*/

extern void closewin(struct Sequence *win)
  {
   if (win==NULL) return;

   if (win->state!=NULL)       free(win->state);
   if (win->composition!=NULL) free(win->composition);
   if (win->classvec!=NULL)    free(win->classvec);
   if (win->scorevec!=NULL)    free(win->scorevec);

   free(win);
   return;
  }

/*----------------------------------------------------------------(compon)---*/

extern void compon(struct Sequence *win)
{
	register int	*comp;
	register int	aa;
	register char	*seq, *seqmax;

	win->composition = comp = (int *) calloc(20*sizeof(*comp), 1);
	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		aa = *seq++;
		if (!aaflag[aa])
			comp[aaindex[aa]]++;
	}

	return;
}

/*---------------------------------------------------------------(stateon)---*/

static int state_cmp(const void *s1, const void *s2)
{
	return *((int*) s2) - *((int*) s1);
}

extern void stateon(struct Sequence *win)
{
	register int	aa, nel, c;

	if (win->composition == NULL)
		compon(win);

	win->state = (int *) malloc(21*sizeof(win->state[0]));

	for (aa = nel = 0; aa < 20; ++aa) {
		if ((c = win->composition[aa]) == 0)
			continue;
		win->state[nel++] = c;
	}
	for (aa = nel; aa < 21; ++aa)
		win->state[aa] = 0;

	qsort(win->state, nel, sizeof(win->state[0]), state_cmp);

	return;
}

/*-----------------------------------------------------------------(enton)---*/

extern void enton(struct Sequence *win)
  {
   if (win->state==NULL) {stateon(win);}

   win->entropy = entropy(win->state);

   return;
  }

/*---------------------------------------------------------------(entropy)---*/
static int		thewindow;
static double	*entray;

#define LN2	0.69314718055994530941723212145818

void
entropy_init(int window)
{
	int		i;
	double	x, xw;

	entray = (double *)malloc((window+1) * sizeof(*entray));
	xw = window;
	for (i = 1; i <= window; ++i) {
		x = i / xw;
		entray[i] = -x * log(x) / LN2;
	}

	thewindow = window;
}

extern double entropy(int *sv)
{
	int	*sv0 = sv;
	register double	ent;
	register int	i, total;
	register int	*svmax;
	register double	xtotrecip, xsv;

	for (total = 0; (i = *sv) != 0; ++sv)
		total += i;
	svmax = sv;
	ent = 0.0;
	if (total == thewindow) {
		for (sv = sv0; sv < svmax; ) {
			ent += entray[*sv++];
		}
		return ent;
	}
	if (total == 0)
		return 0.;

	xtotrecip = 1./(double)total;
	for (sv = sv0; sv < svmax; ) {
		xsv = *sv++;
		ent += xsv * log(xsv * xtotrecip);
	}
	return -ent * xtotrecip / LN2;
}

/*-----------------------------------------------------------(decrementsv)---*/

static void
decrementsv(int *sv, int class)
{
	register int	svi;

	while ((svi = *sv++) != 0) {
		if (svi == class && *sv < class) {
			sv[-1] = svi - 1;
			break;
		}
	}
}

/*-----------------------------------------------------------(incrementsv)---*/

static void
incrementsv(int *sv, int class)
{
	for (;;) {
		if (*sv++ == class) {
			sv[-1]++;
			break;
		}
	}
}

/*--------------------------------------------------------------(findchar)---*/

extern int findchar(char *str, char chr)
  {int i;

   for (i=0; ; i++)
     {
      if (str[i]==chr)
        {
         return(i);
        }
      if (str[i]=='\0')
        {
         return(-1);
        }
     }
   }

/*-----------------------------------------------------------------(upper)---*/

extern void upper(char *string, size_t len)
{
       register char   *stringmax, c;

       for (stringmax = string + len; string < stringmax; ++string)
               if (islower(c = *string))
                       *string = toupper(c);
}

/*-----------------------------------------------------------------(lower)---*/

extern void lower(char *string, size_t len)
{
       register char   *stringmax, c;

       for (stringmax = string + len; string < stringmax; ++string)
               if (isupper(c = *string))
                       *string = tolower(c);
}

/*-------------------------------------------------------------------(min)---*/

int min(int a, int b)
  {
   if (a<b) {return(a);}
   else {return(b);}
  }

/*-------------------------------------------------------------------(max)---*/

int max(int a, int b)
  {
   if (a<b) {return(b);}
   else {return(a);}
  }

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------(tmalloc)---*/

void *tmalloc(size_t size)
  {void *ptr;

   ptr = (void *) malloc(size);

   if (rptr>TESTMAX)
     {
      exit(2);
     }

   //record_ptrs[rptr] = (int) ptr;
   record_ptrs[rptr] = (int64_t) ptr;
   rptr++;

   return(ptr);
  }

/*-----------------------------------------------------------------(tfree)---*/

void tfree(void *ptr)
  {int i;

   for (i=0; i<rptr; i++)
     {
      //if (record_ptrs[i]==(int)ptr)
      if (record_ptrs[i]==(int64_t)ptr)
        {
         record_ptrs[i] = 0;
         break;
        }
      }

   free(ptr);
  }

/*---------------------------------------------------------------------------*/
