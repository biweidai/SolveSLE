#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mkl_lapacke.h"
#define RowNumber 2999 
#define N RowNumber
#define LDA RowNumber
#define LDVL RowNumber
#define LDVR RowNumber
#define l 0
#define rmin 0.0
#define rmax 10000.0
#define step (double)(rmax-rmin)/(RowNumber+1.0)
#define a 1                                           //rmin=a*r[0] boundary condition
#define b 1                                           //rmax=b*r[RowNumber-1] boundary condition

extern void print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi );
extern void print_eigenvectors( char* desc, MKL_INT n, double* wi, double* v, MKL_INT ldv );

double r[RowNumber];
double Matrix[RowNumber*RowNumber]={0};

extern void ConstructMatrix();

int main()
{
	int i;
	for(i=0;i<RowNumber;i++)                                                                       //initialize r[]
	{
		r[i]=rmin+step*(i+1);
	}

	ConstructMatrix();

	printf("Matrix:\n%f %f\n",Matrix[0],Matrix[1]);                                                       //output Matrix
	for(i=1;i<(RowNumber-1);i++)
	{
		printf("%f %f %f\n",Matrix[i*RowNumber+i-1],Matrix[i*RowNumber+i],Matrix[i*RowNumber+i+1]);
	}
	printf("%f %f\n",Matrix[RowNumber*RowNumber-2],Matrix[RowNumber*RowNumber-1]);	
        
	/* Locals */
        MKL_INT n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
        /* Local arrays */

        double wr[N], wi[N], *vl=NULL, *vr=NULL; 
       
        /* Solve eigenproblem */
        info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'N', n, Matrix, lda, wr, wi,
                        vl, ldvl, vr, ldvr );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        print_eigenvalues( "Eigenvalues", n, wr, wi );
        //print_eigenvectors( "Right eigenvectors", n, wi, vr, ldvr );
        exit( 0 );
}

void SecondDerivative()                                                                    //construct the matrix of second derivative term
{
	Matrix[0]-=(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0])*(a-2.0);
	Matrix[1]-=(1.0+r[0])*(1.0+r[0])/(step*step)*log(1.0+r[0]);
	int i;
	for(i=1;i<RowNumber-1;i++)
	{
		Matrix[i-1+i*RowNumber]-=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i]);	
		Matrix[i+i*RowNumber]+=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i])*2.0;	
		Matrix[i+1+i*RowNumber]-=(1.0+r[i])*(1.0+r[i])/(step*step)*log(1.0+r[i]);	
	}
	Matrix[RowNumber*RowNumber-2]-=(1.0+r[RowNumber-1])*(1.0+r[RowNumber-1])/(step*step)*log(1.0+r[RowNumber-1]);
	Matrix[RowNumber*RowNumber-1]-=(1.0+r[RowNumber-1])*(1.0+r[RowNumber-1])/(step*step)*log(1.0+r[RowNumber-1])*(b-2.0);
}

void FirstDerivative()                                                                    // construct the matrix of first derivative term
{
	Matrix[0]+=(1.0+r[0])/(step)*a;
	Matrix[1]-=(1.0+r[0])/(step);
	int i;
	for(i=1;i<RowNumber-1;i++)
	{
		Matrix[i-1+i*RowNumber]+=(1.0+r[i])/(step);
		Matrix[i+1+i*RowNumber]-=(1.0+r[i])/(step);
	}
	Matrix[RowNumber*RowNumber-2]+=(1.0+r[RowNumber-1])/(step);
	Matrix[RowNumber*RowNumber-1]-=b*(1.0+r[RowNumber-1])/(step);
}

void ConstructMatrix()
{
	SecondDerivative();

	FirstDerivative();

	int i;
	for(i=0;i<RowNumber;i++) 
	{
		Matrix[i*RowNumber+i]+=1.0;//+l*(l+1)*(1+r[i])*(1+r[i])/(r[i]*r[i])*log(1+r[i]);
	}
}

/* Auxiliary routine: printing eigenvalues */
void print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi ) {
        MKL_INT j;
        printf( "\n %s\n", desc );
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (float)0.0 ) {
	if(wr[j]>15||wr[j]<0) continue;
         printf( " %6.2f", wr[j] );
      } else {
         printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
      }
   }
   printf( "\n" );
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, MKL_INT n, double* wi, double* v, MKL_INT ldv ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (float)0.0 ) {
            printf( " %6.2f", v[i*ldv+j] );
            j++;
         } else {
            printf( " (%6.2f,%6.2f)", v[i*ldv+j], v[i*ldv+(j+1)] );
            printf( " (%6.2f,%6.2f)", v[i*ldv+j], -v[i*ldv+(j+1)] );
            j += 2;
         }
      }
      printf( "\n" );
   }
}
