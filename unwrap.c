/*
 * Construct a CG trajectory from a AA trajectory.
 * Print the CG trajectory in reverse order.
 * By Zhiyong Zhang, 09/24/2013.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "pca.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

typedef struct 
{
  int num_bpi;
  int ndxlist[60];
} BP_NDX;

int main(int argc, char *argv[])
{

 FILE *ipar, *ocgj;
 XDRFILE *itrj;
 int nbp = 78;
 float cut = 1.0;
 char itrj_nm[30], ocgj_nm[30], *inp_nm;
 int i, j, k, m, n, result, unwrapn1,unwrapn2,bp_j;
 int ifram, nfram, nat, nDNA, nCore, step, *ndxDNA, *ndxCore, **unwrapn;
 rvec dx, *x, *xcg;
 float time, prec ;
 float dist2,dist1,mindist;
 char str[STRLEN];
 BP_NDX *bp_ndx;
 float max[3] = {-1e10, -1e10, -1e10};
 float min[3] = { 1e10,  1e10,  1e10};

 matrix box;

/* open the control file */
 if (argc != 2){ 
     fprintf(stderr,"Usage: unwrap [input parameters]\n");
     exit(1);
 }
 inp_nm = argv[1];
 fprintf(stderr,"Reading input parameter file: %s\n",inp_nm);
 ipar = fopen(inp_nm, "r");

/* control parameters */
 fgets2(str, STRLEN, ipar); // one-line statement of the file

/* Total Nr. of atoms, end DNAs & number of histone core CAs */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d%d%d",&nat, &nDNA, &nCore);

/* Nr. of frames, input AA trajectory and output data file*/
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d%s%s", &nfram, itrj_nm, ocgj_nm);
 fprintf(stderr,"%d frames traj %s to %s\n",nfram,itrj_nm,ocgj_nm);

/* get end-dna index */
 fgets2(str, STRLEN, ipar); // one-line comment
 snew(bp_ndx,  nbp);

 int current_bp,bpi,ndxi;
 current_bp = 0;
 j = 0;
 i = 0;
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d%d", &bpi,&ndxi);
 do {
    if (bpi==current_bp)
    {
        bp_ndx[current_bp].ndxlist[j] = ndxi;
        j = j+1;
        i = i+1;
        fgets2(str, STRLEN, ipar);
        sscanf(str, "%d%d", &bpi,&ndxi);
        /*fprintf(stderr,"%d %d %d\n",i,bpi,ndxi)*/;
    }
    else
    {
        bp_ndx[current_bp].num_bpi = j;
        current_bp += 1;
        j = 0;
    }
    //fprintf(stderr,"%d %d %d %d\n",i,j,ndxi,bp_ndx[current_bp].ndxlist[j]);
 } while (i<nDNA);

 bp_ndx[current_bp].num_bpi = j;

 /*for (i=0;i<nbp;i++) fprintf(stderr,"%d %d\n",i,bp_ndx[i].num_bpi);*/

/* get core heavy atom index */
 fgets2(str, STRLEN, ipar); // one-line comment
 snew(ndxCore,  nCore);

 for (i=0; i<nCore; i++)
 {
    fgets2(str, STRLEN, ipar);
    sscanf(str, "%d", &ndxCore[i]);
    ndxCore[i]--;
 }

/* close the control file */
 fclose(ipar);

 itrj = xdrfile_open(itrj_nm, "r");
 ocgj = fopen(ocgj_nm, "w");

/* allocate space */
 snew(x, nat);
 /*snew(unwrapn,nfram);*/

 float max_x,max_y,max_z,min_x,min_y,min_z;
/* start main loop of all the structures in trajectory */
 n = 0;
 do
 {
/* read frame */
   result = read_xtc(itrj,nat, &step, &time, box, x, &prec);

   if (exdrENDOFFILE != result)
   {

     i = 0;
     do {
       for (m=0; m<3; m++) max[m] = -1e10;
       for (m=0; m<3; m++) min[m] =  1e10;

       for (k=0; k<bp_ndx[i].num_bpi; k++)
       {
         bp_j = bp_ndx[i].ndxlist[k];
         for (m=0; m<3; m++)
         {
           if (x[bp_j][m] > max[m])
           {
             max[m] = x[bp_j][m];
           }
           if (x[bp_j][m] < min[m])
           {
             min[m] = x[bp_j][m];
           }
         }
       }

       for (m=0; m<3; m++) max[m] += cut;
       for (m=0; m<3; m++) min[m] -= cut;
       //fprintf(stderr,"%f %f\n",min[0],max[0]);

       mindist = 1e10;
       for (k=0; k<bp_ndx[i].num_bpi; k++)
       {
          bp_j = bp_ndx[i].ndxlist[k];
          for (j=0; j<nCore; j++) 
          {
             if (x[ndxCore[j]][0] >= min[0] && x[ndxCore[j]][0] <= max[0] &&
                 x[ndxCore[j]][1] >= min[1] && x[ndxCore[j]][1] <= max[1] &&
                 x[ndxCore[j]][2] >= min[2] && x[ndxCore[j]][2] <= max[2])
             {
                dist1 = sqrt(distance2(x[bp_j],x[ndxCore[j]]));
                if (dist1 < mindist) mindist = dist1;
             }
          }
          //fprintf(stderr,"%d %f\n",bp_ndx[i].ndxlist[k],dist1);
        
       }
       //fprintf(stderr,"%f %f\n",mindist,dist1);
       if (mindist <= cut) {
           break;
       }
        /*Unwrapping is determined if all the basepairs before this base pair unwraps*/
       //fprintf(stderr,"%d %f\n",i,mindist);
       /*fprintf(stderr,"%8d %8.3f %8.3f %8.3f\n",i,mindist,wrapratios[n][i],distref[i]);*/
       i++;
     } while(i<nbp/2);
     unwrapn1 = i;
     
     i = nbp-1;
     do {
        /*Unwrapping is determined if all the basepairs before this base pair unwraps*/
       for (m=0; m<3; m++) max[m] = -1e10;
       for (m=0; m<3; m++) min[m] =  1e10;

       //fprintf(stderr,"%d %d\n",i,bp_ndx[i].num_bpi,x[bp_ndx[i].ndxlist[k]][0]);
      
       for (k=0; k<bp_ndx[i].num_bpi; k++)
       {
         for (m=0; m<3; m++)
         {
           
       //    fprintf(stderr,"%f %f\n",min[m],x[bp_ndx[i].ndxlist[k]][m]);
           bp_j = bp_ndx[i].ndxlist[k];
           if (x[bp_j][m] > max[m])
           {
             max[m] = x[bp_j][m];
           }
           if (x[bp_j][m] < min[m])
           {
             min[m] = x[bp_j][m];
           }
         }
       }

       for (m=0; m<3; m++) max[m] += cut;
       for (m=0; m<3; m++) min[m] -= cut;

       mindist = 1e10;
       for (k=0; k<bp_ndx[i].num_bpi; k++)
       {
          bp_j = bp_ndx[i].ndxlist[k];
          for (j=0; j<nCore; j++) 
          {
             if (x[ndxCore[j]][0] >= min[0] && x[ndxCore[j]][0] <= max[0]  &&
                 x[ndxCore[j]][1] >= min[1] && x[ndxCore[j]][1] <= max[1]  &&
                 x[ndxCore[j]][2] >= min[2] && x[ndxCore[j]][2] <= max[2])
             {
                dist1 = sqrt(distance2(x[bp_j],x[ndxCore[j]]));
                if (dist1 < mindist) mindist = dist1;
             }
          }
        
        /*Unwrapping is determined if all the basepairs before this base pair unwraps*/
       }
       if (mindist <= cut) {
        break;
       }
       i--;
     } while(i >= nbp/2);
     unwrapn2 = nbp-1-i;


     /*total unwrapping*/
     /*snew(unwrapn[n],2);
     unwrapn[n][0] = unwrapn1;
     unwrapn[n][1] = unwrapn2;
     */

     printf("\rFrame current: %6d --> Total: %6d",n,nfram);
     fflush(stdout);

    /* output the CG trajectory */
     fprintf(ocgj, "%3d %3d\n", unwrapn1, unwrapn2);

     n++;
   }

 } while(result == exdrOK && n <= nfram);


/* close file */
 xdrfile_close(itrj);
 fclose(ocgj);

 printf("\n");
 fflush(stdout);
/* the end of main */
 return 0;
}
