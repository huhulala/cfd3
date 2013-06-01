#include "boundary_val.h"
#include "NSDefinitions.h"
#include <string.h>

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax,int jmax, int wl, int wr, int wt, int wb, double **U, double **V,
	 double **F, double **G,double **P, int** Flag)
{
    int i = 0;
    int j = 0;

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
		switch(wr)
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


    /* set the values of the inner obstacles */
	/* treat obstacle boundaries, see 1.4 to 1.6 TODO Unperformant!? */
	for(j = 1; j <= jmax; j++)
	for(i = 1; i <= imax; i++)
	{
		switch(Flag[i][j])
		{
		case B_N:
			V[i][j] = 0;
			U[i-1][j] = -U[i-1][j+1];
			U[i][j] = -U[i][j+1];
			break;
		case B_S:
			V[i][j-1] = 0;
			U[i-1][j] = -U[i-1][j-1];
			U[i][j] = -U[i][j-1];
			break;
		case B_W:
			U[i-1][j] = 0;
			V[i][j-1] = -V[i-1][j-1];
			V[i][j] = -V[i-1][j];
			break;
		case B_O:
			U[i][j] = 0;
			V[i][j-1] = -V[i+1][j-1];
			V[i][j] = -V[i+1][j];
			break;
		case B_NO:
			U[i][j] = 0;
			U[i-1][j] = -U[i-1][j+1];
			V[i][j] = 0;
			V[i][j-1] = -V[i+1][j-1];
			break;
		case B_NW:
			U[i-1][j] = 0;
			U[i][j] = -U[i][j+1];
			V[i][j] = 0;
			V[i][j-1] = -V[i-1][j-1];
			break;
		case B_SO:
			U[i][j] = 0;
			U[i-1][j] = -U[i-1][j-1];
			V[i][j-1] = 0;
			V[i][j] = -V[i+1][j];
			break;
		case B_SW:
			U[i-1][j] = 0;
			U[i][j] = -U[i][j-1];
			V[i][j-1] = 0;
			V[i][j] = -V[i-1][j];
			break;
		}
	}
}

void spec_boundary_val(char *problem,int imax, int jmax, double dx, double dy,
    double Re, double deltaP, double **U, double **V, double **P )
{
    int j = 0;
    if (strcmp(problem, "karman") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
        	/* set karman inflow */
            U[0][j] = 1.0;
            V[0][j] = 0.0;
        }
    }
    else if(strcmp(problem, "plane")==0)
    {
    	for(j=1; j<jmax+1; ++j)
    	{
    		U[0][j] = -0.5 * Re * deltaP * j * (j-jmax);
    		V[0][j] = 0;
    	}
    }
	else if(strcmp(problem, "step")==0)
	{
		for(j=1; j<(jmax+2)/2; ++j)
		{
			U[0][j] = 0.0;
			V[0][j] = 0.0;
		}
		for(j=(jmax+2)/2+1; j<jmax+1; ++j)
		{
			U[0][j] = 1.0;
			V[0][j] = 0.0;
		}
	}
}

