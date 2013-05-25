#include "boundary_val.h"
#include "NSDefinitions.h"

void boundaryvalues(int imax, int jmax, double **U, double **V, int wl, int wr, int wt, int wb)
{
	int i;
	int j;
	for(i=1; i<=imax; i++)
	{
		/* lower bounder */
		switch(wb)
		{
			case NO_SLIP:
				U[i][0] = -U[i][1];
				V[i][0] = 0.0;
				break;
			case FREE_SLIP:
				U[i][0] = U[i][1];
				V[i][0] = 0.0;
				break;
			case OUTFLOW:
				U[i][0] = U[i][1];
				V[i][0] = V[i][1];
				break;
		}

		/* upper bounder */
		switch(wt)
		{
			case NO_SLIP:
				U[i][jmax+1] = -U[i][jmax];
				V[i][jmax] = 0.0;
				break;
			case FREE_SLIP:
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = 0.0;
				break;
			case OUTFLOW:
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = V[i][jmax-1];
				break;
		}
	}
	for(j=1; j<=jmax; j++)
	{
		/* left border */
		switch(wl)
		{
			case NO_SLIP:
				U[0][j] = 0.0;
				V[0][j] = -V[1][j];
				break;
			case FREE_SLIP:
				U[0][j] = 0.0;
				V[0][j] = V[1][j];
				break;
			case OUTFLOW:
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
				break;
		}

		/* right border */
		switch(wb)
		{
			case NO_SLIP:
				U[imax][j] = 0.0;
				V[imax+1][j] = -V[imax][j];
				break;
			case FREE_SLIP:
				U[imax][j] = 0.0;
				V[imax+1][j] = V[imax][j];
				break;
			case OUTFLOW:
				U[imax][j] = U[imax-1][j];
				V[imax+1][j] = V[imax][j];
				break;
		}
	}
}


void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V)
{
	int i,j;
	double Re, dp;
	/* TODO: wo kommen Re und P her? */
	Re = 10;
	dp = 4;
	if(strcmp(problem, "Driven cavity"))
	{
		for(i=0; i<imax+1; i++)
		{
			U[i][jmax] = -U[i][jmax-1] + 2;
		}
	}
	else if(strcmp(problem, "Karman vortex street"))
	{
		for(j=0; j<jmax+2; ++j)
		{
			U[0][j] = 1;
			V[0][j] = 0;
		}
	}
	else if(strcmp(problem, "Plane shear flow"))
	{
		for(j=0; j<jmax+2; ++j)
		{
			U[0][j] = -0.5 * Re * dp * j * (j-jmax);
			V[0][j] = 0;
		}
	}
	else if(strcmp(problem, "Flow over a step"))
	{
		for(j=0; j<(jmax+2)/2; ++j)
		{
			U[0][j] = 1;
			V[0][j] = 0;
		}
		for(j=(jmax+2)/2+1; j<jmax+2; ++j)
		{
			U[0][j] = 0;
			V[0][j] = 0;
		}
	}
}
