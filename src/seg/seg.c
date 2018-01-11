
/*****************************************************************************/
/***  (seg.c)                                                              ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"                   ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"
#include "lnfac.h"
#include "seg.h"

/*---------------------------------------------------------------(defines)---*/

#define LENCORRLIM 120
#define MIN(a,b)	((a) <= (b) ? (a) : (b))

/*---------------------------------------------------------------(structs)---*/
/*
struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };
*/
/*---------------------------------------------------------------(globals)---*/

#include <string.h>

int window = 12;
int downset = (window+1)/2 - 1;
int upset = window - downset;;
double locut = 2.2;
double hicut = 2.5;

int hilenmin = 0;
int overlaps = FALSE;
int hionly = FALSE;
int loonly = FALSE;
int entinfo = TRUE;
int singleseq = TRUE; //FALSE;
int prettyseq = FALSE;
int prettytree = FALSE; //TRUE;
int charline = 60;
int maxtrim = 100;

double getprob(), lnperm(), lnass();
void segseq(struct Sequence *seq, struct Segment **segs, int offset);
void usage();
void appendseg(struct Segment *segs, struct Segment *seg);
void freesegs(struct Segment *segs);
void seqout(struct Sequence *seq, int hilo, int begin, int end);
void pretreereport(struct Sequence *seq, struct Segment *segs);
void prettyreport(struct Sequence *seq, struct Segment *segs);
void singreport(struct Sequence *seq, struct Segment *segs);
void report(struct Sequence *seq, struct Segment *segs);
void mergesegs(struct Sequence *seq, struct Segment *segs);
void trim(struct Sequence *seq, int *leftend, int *rightend);
int hasdash(struct Sequence *win);
int findlo(int i, int limit, double* H);
int findhi(int i, int limit, double* H);


struct Sequence* test() {
    struct Sequence *seq = (struct Sequence *) malloc(sizeof(struct Sequence));
    
    seq->db = NULL;
    seq->parent = NULL;
    seq->root = NULL;
    seq->children = NULL;
    seq->rubberwin = FALSE;
    seq->floatwin = FALSE;

    seq->id = NULL;
    seq->name = NULL;
    seq->organism = NULL;

    seq->header = (char*) malloc(6);
    strcpy(seq->header, ">test");
    seq->seq = (char*) malloc(307);
    strcpy(seq->seq, "MELIFLGTNAGVPTKERNVTSIALDLHGIRPAFWLFDCGEGTQHQILRTPVKPGKLEKIFITHLHGDHLFGLPGLLCSRSMAGIETPLTLYGPAGLKTFVETTLALSGSWLTYPLEVVEVVPGTVCEDRQLRVTAHELSHTLYCVGYRIEERPKPGPLDVDKLAAEGIKPGAYFQQLKRGETVTLDDGRELNGWDYVGPGLPGKSLAIFGDTRPTPEALKLAAGVDVMVHEATLEGAMAERANERGHSTTLQTAAAARDAGAKRLIITHFSARYGQEDLERLRQECQTLFPATEVATDLATFRV");

    seq->length = strlen(seq->seq);
    seq->punctuation = FALSE;

    seq->entropy = -2.0;
    seq->state = NULL;
    seq->composition = NULL;
    seq->classvec = NULL;
    seq->scorevec = NULL;

    return seq;
}

/*------------------------------------------------------------------(main)---*/

int main(int argc, char** argv) {
   //struct Database *db;
   struct Sequence *seq;
   struct Segment *segs;
   //int ctime;

   genwininit();
/* readlencorr(); */                        /* #include lencorr file */
   downset = (window+1)/2 - 1;
   upset = window - downset;

   entropy_init(window);
/*
   if ((db=opendbase(argv[1]))==NULL)
     {
      fprintf(stderr, "Error opening file %s\n", argv[1]);
      exit(1);
     }

   for (seq=firstseq(db); seq!=NULL; seq=nextseq(db))
*/
     {
      seq = test();
      segs = (struct Segment *) NULL;
      segseq(seq, &segs, 0);
      mergesegs(seq, segs);

      //if (singleseq) 
        singreport(seq, segs);
      //else if (prettyseq) prettyreport(seq, segs);
      //else if (prettytree) pretreereport(seq, segs);
      //else report(seq, segs);

      freesegs(segs);
      closeseq(seq);
     }

//   closedbase(db);
   exit(0);
  }

/*---------------------------------------------------------------(segment)---*/

void segseq(struct Sequence *seq, struct Segment **segs, int offset)
  {struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int first, last, lowlim;
   int loi, hii, i;
   int leftend, rightend, lend, rend;
   double *H, *seqent();

   H = seqent(seq);
   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   for (i=first; i<=last; i++)
     {
      if (H[i]<=locut && H[i]!=-1)
        {
         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);
           }

         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = min(hii, rightend+downset);
         lowlim = i + 1;
        }
     }

   free(H);
   return;
  }

/*----------------------------------------------------------------(seqent)---*/

double *seqent(struct Sequence *seq)
  {struct Sequence *win;
   double *H;
   int i, first, last;

   if (window>seq->length)
     {
      return((double *) NULL);
     }

   H = (double *) malloc(seq->length*sizeof(double));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {H[i] = -1;
         shiftwin1(win);
         continue;}
      H[i] = win->entropy;
      shiftwin1(win);
     }

   closewin(win);
   return(H);
  }

/*---------------------------------------------------------------(hasdash)---*/

int hasdash(struct Sequence *win)
{
	register char	*seq, *seqmax;

	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		if (*seq++ == '-')
			return TRUE;
	}
	return FALSE;
}

/*----------------------------------------------------------------(findlo)---*/

int findlo(int i, int limit, double* H)
  {int j;

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*----------------------------------------------------------------(findhi)---*/

int findhi(int i, int limit, double *H)
  {int j;

   for (j=i; j<=limit; j++)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*------------------------------------------------------------------(trim)---*/

void trim(struct Sequence *seq, int *leftend, int *rightend)
  {struct Sequence *win;
   double prob, minprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
   for (len=seq->length; len>minlen; len--)
     {
      win = openwin(seq, 0, len);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len);
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin1(win);
         i++;
        }
      closewin(win);
     }

/* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

   closewin(seq);
   return;
  }

/*---------------------------------------------------------------(getprob)---*/

double getprob(int *sv, int total)
  {double ans, totseq;

#define LN20	2.9957322735539909
   totseq = ((double) total) * LN20;

   ans = lnass(sv) + lnperm(sv, total) - totseq;

   return(ans);
  }

/*----------------------------------------------------------------(lnperm)---*/

double lnperm(int *sv, int tot)
  {double ans;
   int i;

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfac[sv[i]];
     }

   return(ans);
  }

/*-----------------------------------------------------------------(lnass)---*/

double lnass(int *sv)
{
	double	ans;
	register int	svi, svim1;
	register int	class, total;
	register int    i;

	ans = lnfac[20];
	if (sv[0] == 0)
		return ans;

	total = 20;
	class = 1;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
	        if (++i==20) {
		        ans -= lnfac[class];
                        break;
		      }
		else if ((svi = *++sv) == svim1) {
			class++;
			continue;
		}
		else {
			total -= class;
			ans -= lnfac[class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
			}
			else {
				class = 1;
				continue;
			}
		}
	}

	return ans;
}

/*-------------------------------------------------------------(mergesegs)---*/

void mergesegs(struct Sequence *seq, struct Segment *segs)
  {struct Segment *seg, *nextseg;
   int len;

   if (overlaps) return;
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*------------------------------------------------------------(singreport)---*/

void singreport(struct Sequence *seq, struct Segment *segs)
{
	char	*proseq, *proseqmax;
	struct Segment	*seg;
	int	begin, end, i, ctr;

	proseq = seq->seq;
	proseqmax = proseq + seq->length;
	upper(proseq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 'x', end - begin +1);
	}

	fprintf(stdout, "%s\n", seq->header);

	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	putc('\n', stdout);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*----------------------------------------------------------------(seqout)---*/

#define HDRLEN 60
void seqout(struct Sequence *seq, int hilo, int begin, int end)
{
	char	*proseq, *proseqmax, *id, *header;
   char outbuf[HDRLEN+1];
   static int hi = 1;
   static int lo = 0;
   int i, ctr, iend;

   if (hionly && hilo==lo) return;
   if (loonly && hilo==hi) return;

   proseq = seq->seq;
   proseqmax = proseq + seq->length;
   id = seq->id;
   if (id==NULL) id = seq->parent->id;
   header = seq->header;
   if (header==NULL) header = seq->parent->header;

   iend = findchar(header, ' ');
   if (iend!=-1) header = header+iend;

   if (entinfo)
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
/*    if (iend!=-1 && strlen(header)<=HDRLEN) fprintf(stdout, "%s", header);
      else if (iend!=-1) for (i=0; i<HDRLEN; i++) putc(header[i], stdout); */
      fprintf(stdout, " complexity=%4.2f (%d/%4.2f/%4.2f)\n",
           seq->entropy, window, locut, hicut);
     }
   else
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
      if (iend!=-1)   /* fprintf(stdout, "%s\n", header); */
        {
		 i = MIN(HDRLEN, strlen(header));
		 fwrite(header, i, 1, stdout);
         putc('\n', stdout);
        }
      else putc('\n', stdout);
     }
   
   if (hilo==lo)
     {
      lower(proseq, seq->length);
     }
   else if (hilo==hi && seq->length>=hilenmin)
     {
      upper(proseq, seq->length);
     }
   else
     {
      lower(proseq, seq->length);
     }

   for (; proseq < proseqmax; proseq+=i) {
		i = MIN(charline, proseqmax - proseq);
		fwrite(proseq, i, 1, stdout);
		putc('\n', stdout);
	}

	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*-------------------------------------------------------------(appendseg)---*/

void appendseg(struct Segment *segs, struct Segment *seg)
  {struct Segment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*--------------------------------------------------------------(freesegs)---*/

void freesegs(struct Segment *segs)
  {struct Segment *temp;

   while (segs!=NULL)
     {
      temp = segs->next;
      free(segs);
      segs = temp;
     }
  }

/*-----------------------------------------------------------------(usage)---*/

void usage()

  {
   fprintf(stderr, "\
Usage: seg <file> <window> <locut> <hicut> <options>\n\
         <file>   - filename containing fasta-formatted sequence(s) \n\
         <window> - OPTIONAL window size (default 12) \n\
         <locut>  - OPTIONAL low (trigger) complexity (default 2.2) \n\
         <hicut>  - OPTIONAL high (extension) complexity (default locut + 0.3) \n\
	 <options> \n\
            -x  each input sequence is represented by a single output \n\
                sequence with low-complexity regions replaced by \n\
                strings of 'x' characters \n\
            -c <chars> number of sequence characters/line (default 60)\n\
            -m <size> minimum length for a high-complexity segment \n\
                (default 0).  Shorter segments are merged with adjacent \n\
                low-complexity segments \n\
            -l  show only low-complexity segments (fasta format) \n\
            -h  show only high-complexity segments (fasta format) \n\
            -a  show all segments (fasta format) \n\
            -n  do not add complexity information to the header line \n\
            -o  show overlapping low-complexity segments (default merge) \n\
            -t <maxtrim> maximum trimming of raw segment (default 100) \n\
            -p  prettyprint each segmented sequence (tree format) \n\
            -q  prettyprint each segmented sequence (block format) \n");
   exit(1);
  }

