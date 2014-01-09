/*
**      File:   nrutil.h
**      Purpose: Memory allocation routines borrowed from the
**              book "Numerical Recipes" by Press, Flannery, Teukolsky,
**              and Vetterling.
**              state sequence and probablity of observing a sequence
**              given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.h,v 1.2 1998/02/19 16:32:42 kanungo Exp kanungo $
*/
#ifndef NRUTIL_H_INCLUDED
#define NRUTIL_H_INCLUDED
float *Vector(int,int);
float **matrix(int,int,int,int);
float **convert_matrix(float*,int,int,int,int);
double *dvector(int,int);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
int *ivector(int,int);
int **imatrix(int,int,int,int);
float **submatrix(int,int,int,int);
void free_vector(float*,int,int);
void free_dvector(double*,int,int);
void free_ivector(int*,int,int);
void free_matrix(float*,int,int,int,int);
void free_dmatrix(double**,int,int,int,int);
void free_imatrix(int**,int,int,int,int);
void free_submatrix(float**,int,int,int,int);
void free_convert_matrix(float*,int,int,int,int);
void nrerror();
#endif
