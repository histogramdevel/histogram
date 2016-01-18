/*
Program to perform statistics and plot histograms
in text mode to ensure maximum portability
Written in ANSI C (I hope ;)
Written by L. Kral <http://kral.astronomy.cz/>, Feb. 2004
* 
* 
* 
Modified by M.Bohm, Dec.2015
Using additional graphical DISLIN library:
http://www.mps.mpg.de/dislin
For compiling and running on different platforms use clink -a c_source_code_name in command prompt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dislin.h"

char     version[7] = "0.8";

long     begindata;  /* position of first data in the file */
long     npoints = 0; /* number of points used to compute the histogram */
int      coln = 1; /* number of column of our interest in the file */
int      Wwidth = 79; /* histogram window width in characters */
int      Wheight = 45; /* histogram window height in characters */
int      nbars; /* number of histogram bars */
int      barwidth = 1; /* bar width in characters
                          -- must be 1 or a multiple of 2!!! */
int      lwidth = 15; /* width of legend to each column */
unsigned char
         histol = 220, histor = 223, histom = 219, lowc = '.',gaussc = '*';
         /* histogram and gauss style */
double   dmin,dmax,odmin,odmax; /* data min and max */
double   dx; /* x size of one histogram bar */
float    normfact,onormfact; /* norming factor for histogram */
long     nmax; /* maximum Y value in histogram */
/*Marek Bohm*/
float x_coordinates[200];
float zeroes[200];
float number_of_points[200];
FILE     *Fwork; /* working copy of input data file, can be cropped in range! */

struct Thist /* one bar of histogram */
{
  double x;     /* x coordinate of column middle */
  long   n;     /* counts */
  float  gauss; /* gaussian from statistics */
  float  zeroes;
} histoga[200];


int pause(void)
/* wait for pressing Enter */
{
 char s[10];
 printf("\n<Press Enter to continue>\n");
 gets(s);
 return 1;
}

int signum(double d)
{
 if (d>=0.0) return 1;
 else return -1;
}

char upcase(char c)
{
 if((c>='a') && (c<='z')) return c - 'a' + 'A';
 else return c;
}

int readchar(void)
/* reads first character from string entered to stdin, empties buffer
   Outputs 0 even if only Enter was pressed */

{
 char c;
 char s[100],ts[50];
 gets(s);
 if (sscanf(s,"%c%s",&c,ts)>=1) return c;
 else return 0;
}


void errhalt(char msg[100]) /* Stop program with error message */
{
 printf(msg);
 printf("\n");
 exit(1);
}

double round(double f)
/* round float to the NEAREST integer */
{
  float tmpi,tmpf;

  if (f>2147483647.49 || f<-2147483648.49)
  {
    printf("Overflow error rounding double to long: double out of range! Stopped.\n");
    exit(1);
  }

  tmpi = floor(fabs(f));
  tmpf = fabs(f) - tmpi;

  if (f>=0.0)
  {
    if (tmpf>=0.5) {return ceil(f);} else {return floor(f);}
  }
  else
  {
    if (tmpf>=0.5) {return floor(f);} else {return ceil(f);}
  }
}


double parserow(char rstr[200],int *err)
/* reads double from coln-th column from string rstr */
{
 double  row[12];

 if (sscanf(rstr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &row[0],&row[1],&row[2],&row[3],&row[4],&row[5],&row[6],&row[7],
        &row[8],&row[9],&row[10],&row[11])<(coln))
 {
   *err = 1;
   return 0;
 }
 else *err = 0;
 return row[coln-1];
}


int fxcopy(char src[40],char dest[40],double lmin,double lmax,int filter,
           long *npoints,double *mean,double *min,double *max,int verbose)
/* copies data file, but only lines containing values in selected range
   2 lines of header are expected!
   Perform basic statistics */
{
 FILE *fsrc,*fdest;
 char Stmp[1000];
 int  err;
 double Dtmp;

 if ((fsrc = fopen(src, "rt")) == NULL) return 0;
 if ((fdest = fopen(dest,"wt")) == NULL) return 0;

 if (verbose) printf("Copying file... ");

 /* copy 2 lines of header */
 if (fgets(Stmp,1000,fsrc)) fputs(Stmp,fdest);
 if (fgets(Stmp,1000,fsrc)) fputs(Stmp,fdest);

 *npoints = 0;
 *mean = 0.0;
 *max = 0.0;
 *min = 0.0;

 while (fgets(Stmp,1000,fsrc))
 {
   if (filter)  /* filtering */
   {
     if ((Dtmp = parserow(Stmp,&err))>=lmin && Dtmp<=lmax && !err)
     {
       fputs(Stmp,fdest);
       /* statistics */
       *npoints += 1;
       if (*npoints==1) *min = *max = Dtmp; /* first row -- initialize max and min */
       if (Dtmp>*max) *max = Dtmp;
       if (Dtmp<*min) *min = Dtmp;
       *mean += Dtmp;
     }
   }
   else /* copy everything */
   {
     fputs(Stmp,fdest); /*copy */
     Dtmp = parserow(Stmp,&err);
     /* statistics */
     *npoints += 1;
     if (*npoints==1) *min = *max = Dtmp; /* first row -- initialize max and min */
     if (Dtmp>*max) *max = Dtmp;
     if (Dtmp<*min) *min = Dtmp;
     *mean += Dtmp;
   }
 }

 *mean /= *npoints;

 fclose(fsrc);
 fclose(fdest);
 if (verbose) printf("Done.\n");
 return 1;
}


double readrow(FILE *Fin,int *eof)
/* outputs number in coln-column from file Fin */
{
 int     err;
 char    Stmp[200];
 double  Dtmp;

 if (fgets(Stmp,200,Fin) != NULL)
 {
   Dtmp = parserow(Stmp,&err);
   if(err)
   {
     printf("\nError reading file -- selected column not found or not complete! Stopped.\n");
     fclose(Fin);
     exit(1);
   }
   *eof = 0;
   return Dtmp;
 }
 else /* end of file */
 {
   *eof = 1;
   return 0;
 }
}

int findlimits(void)
/* finding data minimum and maximum */
{
 double Dtmp;
 int    eof;

 fseek(Fwork,begindata,0);
 Dtmp = readrow(Fwork,&eof);
 if (eof)
 {
  printf("Error -- file does not contain any data! Stopped.\n");
  exit(1);
 }

 dmax = dmin = Dtmp;

 printf("Scanning file for data range... ");
 do
 {
   Dtmp = readrow(Fwork,&eof);
   if (!eof)
   {
     if (Dtmp>dmax) dmax = Dtmp;
     if (Dtmp<dmin) dmin = Dtmp;
   }
 } while(!eof);
 if (dmax==dmin)
 {
   dmax += 1;
   dmin -= 1;
 }
 odmin = dmin; odmax = dmax; /* storing original values */
 printf("Done.\n");
 return 1;
}


int compnbars(void)
{
 /* number of histogram bars */
 nbars  = Wheight;
 nbars /= barwidth;
 return 1;
}

int compdx(void)
{
 /* x size of one bar */
 dx = (dmax - dmin)/(float)(nbars - 1);
 return 1;
}

int makehisto(void)
/* make histogram */
{
 int    i,eof;
 double Dtmp;
 /* computing x coordinates for histogram bars */
 for (i=0;i<nbars;i++)
 {
   histoga[i].x = dmin + (dx/(2.0))*(2*i+1);
 }

 for (i=0;i<nbars;i++)
 {
   x_coordinates[i] = histoga[i].x;
 }

 /* initializing histogram with zeros */
 for (i=0;i<200;i++)
 {
   histoga[i].n = 0;
   histoga[i].gauss = 0.0;
   histoga[i].zeroes = 0.0;
     number_of_points[i] = 0.0;
 }

 for (i=0;i<200;i++)
 {
  zeroes[i] = histoga[i].zeroes;
 }

 /* computing histogram values */
 printf("Computing histogram... ");
 fseek(Fwork,begindata,0);
 npoints = 0;
 do
 {
   Dtmp = readrow(Fwork,&eof);
   if (!eof)
   {
     npoints++;
     for (i=0;i<nbars;i++)
     {
       if (Dtmp>=(histoga[i].x-dx/2) && Dtmp<(histoga[i].x+dx/2))
       {
         histoga[i].n++;
         break;
       }
     }
   }
 } while(!eof);
 printf("Done.\n");

 /*Marek Bohm*/
 for (i=0;i<200;i++)
 {
  number_of_points[i] = histoga[i].n;
 }

 return 1;
}


int compnorm(void)
/* computing norming factor for histogram plotting */
{
 int i;

 nmax = 0;
 for (i=0; i<nbars; i++)
 {
   if (histoga[i].n>nmax) nmax = histoga[i].n;
 }
 normfact = Wwidth - 1 - lwidth;
 if (nmax==0)
 {
   printf("Error -- resulting histogram contains only zeros! Stopped.");
   exit(1);
 }
 onormfact = normfact /= nmax;   /* => #chars = normfact * n */
 return 1;
}

FILE *filetest(char fname[100])
/* tests if file exists, asks if overwrite, if successfull then returns
   pointer to file, if not returns NULL */
{
 FILE *ff;
 char c = '0';

 if ((ff = fopen(fname,"rt"))!=NULL) /* file already exists */
 {
   fclose(ff);
   printf("File %s already exists - overwrite (y/n)? ",fname);
   do {c = upcase(readchar());} while (c!='Y' && c!='N');
   if (c=='Y')
   {
     printf("\nFile will be overwritten.\n");
     if ((ff = fopen(fname,"wt"))!=NULL) return ff;
     else
     {
       printf("Error - cannot open file %s for writing!\n",fname);
       return NULL;
     }
   }
   else
   {
     printf("\nSaving cancelled!\n");
     return NULL;
   }
 }
 else
 {
   if ((ff = fopen(fname,"wt"))!=NULL) return ff;
   else
   {
     printf("Error - cannot open file %s for writing!\n",fname);
     return NULL;
   }
 }
}

/****************************************************************************/

int main(int argc,char *argv[])
{
 int      i;
 int      itern = 0; /* number of filtering interation */
 int      gflag;
 int      newnbars;
 float    Ns = -1; /* criterion for N*sigma filter */
 char     key = '0';
 char     Stmp[200];
 double   Dtmp,Ddum1,Ddum2,Ddum3;
 double   mean = 0,RMS = -1; /* statistical properties of data file */
 double   newdmin,newdmax; /* min and max after filtering */
 double   new_of_dmin; /* min after setting the ofset*/
 double   new_ext_dmax; /* new bax after setting the new number and x size of one histogram bar*/ 
 double   dlim1,dlim2;
 double   odx,newdx; /* original and new values of x size of one histogram bar*/
 long     Ldum;
 long     newnpoints1,newnpoints2; /* number of points in N*sigma filtered file */
 long     rempoints; /* number of points removed by N*sigma iteration */
 char     infn[200]; /* input data file name */
 char     workfn1[30] = "hwork1.tmp"; /* working data file name */
 char     workfn2[30] = "hwork2.tmp"; /* working data file name */
 char     filtfn1[30] = "hfilt1.tmp"; /* filtered data file name */
 char     filtfn2[30] = "hfilt2.tmp"; /* temporary file name for filtering */
 char     legend1[150] = "HISTOGRAM INFO: ";
 unsigned char plot[200][100]; /* ASCII 2D plot */
 unsigned char plot_s[200][100]; /* plot with statistics printed */
 FILE    *Fflt; /* n*sigma filtered copy of Fwork */
 FILE    *Ftmp;


 printf("\nHistogram plot v. %s    Author: Lukas Kral\n",version);
 printf("Two lines of header expected.\n");

 /* prompt for file name */
 if (argc<2)
 {
   printf("\nINPUT DATA FILENAME: ");
   fgets(Stmp,200,stdin);
   sscanf(Stmp,"%s\n",infn);
 }
 else strcpy(infn,argv[1]);

 /* checking for access to input file */
 if ((Fwork = fopen(infn, "rt")) == NULL)
 {
    fprintf(stderr, "Cannot open input file %s! Stopped.\n",infn);
    printf("q%sq, argc=%d",argv[1],argc);
    exit(1);
 }
 else fclose(Fwork);

 /* reading number of column */
 if (argc<3)
 {
   printf("\ncoordinate column?:  ");
   fgets(Stmp,200,stdin);
 }
 else strcpy(Stmp,argv[2]);

 if (sscanf(Stmp,"%d",&coln)!=1 || coln>12)
 {
   printf("Bad column number -- must be an integer number <12! Stopped.\n");
   exit(1);
 }

 /* copy input data file to working file */
 if (!fxcopy(infn,workfn1,0,0,0,&npoints,&mean,&dmin,&dmax,1)) errhalt("Error reading copying input file! Stopped.");
 odmin = dmin; odmax = dmax; /* storing original values */

 Fwork = fopen(workfn1,"rt");

 /* skipping 2 lines of expected header */
 fgets(Stmp,200,Fwork);
 fgets(Stmp,200,Fwork);
 begindata = ftell(Fwork);


 compnbars();
 compdx();
 makehisto();
 compnorm();
 
 
 /* printing legend */
   sprintf(legend1,"Bin size = %.3le   %ld points   Y-range 0...%.1lf counts   X-range %.2lf...%.2lf\n",dx,npoints,((Wwidth-1-lwidth)/normfact),dmin,dmax);
   printf(legend1);


   /* manually set histogram x offset */ 
 	 printf("HISTOGRAM OFFSET?:  ");
	 
     	 if ( scanf("%lf",&new_of_dmin) == 1 ) 
	 {
     odmin = dmin;
     odmax = dmax;
     dmin = new_of_dmin;
     compnbars();
     compdx();
     makehisto();
     compnorm();    
     } 
     
     else {
		   printf("Please enter a number. Exiting...  ");
           exit(1);
          }

     
    /* manually set number of bars */
	 printf("No. OF CELLS?:  ");
	 scanf("%d",&newnbars);
	 Wheight = newnbars + 1;
     compnbars();
     compdx();
     makehisto();
     compnorm();
     
     
    /* manual adjust of bin size */
     printf("CELL WIDTH?:  ");
	 scanf("%lf",&newdx);
	 odx = dx;
	 dx = newdx;
	 new_ext_dmax = dmin + newnbars*dx;
	 odmax = dmax;
	 dmax = new_ext_dmax;
     compnbars();
     compdx();
     makehisto();
     compnorm();
    
 do
 {

      if (key=='H') /* display histogram in separate window using DISLIN library */
   {
	   /* creating plot -> memory */
      //int nya = 2700, i;
      char   *ctit = "Histogram", cbuf[25];
      scrmod ("revers"); // background color set to white
      /* do not use setpag to set the height and width of the xwindow */
      page(4000,2100);     
      metafl ("CONS");
      disini ();
      nochek ();
      pagera (); //plot a border around the page
      complx ();
      ticks  (5, "x");
      axslen(3650,1700);
      axspos(200,1800);
      shdpat (6); // selects shading patterns
      graf   (dmin, dmax, dmin + (dx/(2.0)), 2*dx, 0.0, nmax + 100, 0.0, 100.0);
      labels ("second", "bars");
      labpos ("outside", "bars");
      labdig (1, "x");
      labdig (-2, "bars");
      color  ("red");
      bars(x_coordinates,zeroes,number_of_points, nbars);
      //barwth (-12);
      //barmod("fixed","width");
      //bargrp (nbars, 0);
      //barwth(1);
      color  ("fore");
      height (50);
      title  ();
      endgrf ();
      disfin ();
      //pause();
   }
   
   
         if (key=='C') /* display histogram in separate window using DISLIN library */
   {
	   /* creating plot -> memory */
      //int nya = 2700, i;
      char filename[256];
      printf("FILENAME?:  ");
      fgets(filename, sizeof(filename), stdin);
      printf(filename);
      char   *ctit = "Histogram", cbuf[25];
      scrmod ("revers"); // background color set to white
      /* do not use setpag to set the height and width of the xwindow */
      page(4000,2100);     
      metafl ("SVG");
      //setfil ("nostromo");
      disini ();
      nochek ();
      pagera (); //plot a border around the page
      complx ();
      ticks  (5, "x");
      axslen(3650,1700);
      axspos(200,1800);
      shdpat (6); // selects shading patterns
      graf   (dmin, dmax, dmin + (dx/(2.0)), 2*dx, 0.0, nmax + 100, 0.0, 100.0);
      labels ("second", "bars");
      labpos ("outside", "bars");
      labdig (1, "x");
      labdig (-2, "bars");
      color  ("red");
      bars(x_coordinates,zeroes,number_of_points, nbars);
      //barwth (-12);
      //barmod("fixed","width");
      //bargrp (nbars, 0);
      //barwth(1);
      color  ("fore");
      height (50);
      title  ();
      endgrf ();
      disfin ();
      //pause();
   }


   
  if (key=='R') /* redraw histogram */   
    {
     /* manually set histogram x offset */ 
 	 printf("HISTOGRAM OFFSET?:  ");
	 
     	 if ( scanf("%lf",&new_of_dmin) == 1 ) 
	 {
     odmin = dmin;
     odmax = dmax;
     dmin = new_of_dmin;
     compnbars();
     compdx();
     makehisto();
     compnorm();    
     } 
     
     else {
		   printf("Please enter a number. Exiting...  ");
           exit(1);
          }

     
     /* manually set number of bars */
	 printf("No. OF CELLS?:  ");
	 scanf("%d",&newnbars);
	 Wheight = newnbars + 1;
     compnbars();
     compdx();
     makehisto();
     compnorm();
     
     
     /* manual adjust of bin size */
     printf("CELL WIDTH?:  ");
	 scanf("%lf",&newdx);
	 odx = dx;
	 dx = newdx;
	 new_ext_dmax = dmin + newnbars*dx;
	 odmax = dmax;
	 dmax = new_ext_dmax;
     compnbars();
     compdx();
     makehisto();
     compnorm();
    }   
    
    /* printing legend */
    sprintf(legend1,"Bin size = %.3le   %ld points   Y-range 0...%.1lf counts   X-range %.2lf...%.2lf\n",dx,npoints,((Wwidth-1-lwidth)/normfact),dmin,dmax);
    printf(legend1);
    printf(" L - label | Q - quit | R - redraw | S - statistics | H - plot\n");
    printf(" C - save plot as SVG...\n");
    
    /* reading input from user - first character only */
    key = upcase(readchar());    

 } while(key!='Q');

 /* close input data file */
 fclose(Fwork);

 return 0;
}
