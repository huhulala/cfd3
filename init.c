#include "helper.h"
#include "init.h"
#include "NSDefinitions.h"

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                                        /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    int    *wl,
                    int    *wr,
                    int    *wt,
                    int    *wb,
        		    double *dt_value,            /* time for output */
                    double *deltaP
)
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );
   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );
   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *dt_value );
   READ_DOUBLE( szFileName, *deltaP );

   READ_INT( szFileName, *wl );
   READ_INT( szFileName, *wr );
   READ_INT( szFileName, *wt );
   READ_INT( szFileName, *wb );
   READ_INT   ( szFileName, *itermax );
   return 1;
}



/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvp(double UI, double VI, double PI, int imax, int jmax,
		char* problem, double **U, double **V, double **P)
{
	int i;
	int j;
	for(j=1; j<=jmax; j++)
	for(i=1; i<=imax; i++)
	{
		U[i][j] = UI;
		V[i][j] = VI;
		P[i][j] = PI;
	}

    if (strcmp(problem, "step") == 0)
    {
        for (i = 0; i <=imax + 1; ++i)
        for (j = 0; j < jmax/2; ++j)
        	U[i][j] = 0;
    }
}


/**  Init the flag field **/
int init_flag(int **Problem,int imax,int jmax, int **Flag)
{
    int i, j;
    /* normalize problem - due to read_pgm returns not always 1 and zero,
     * the problem needs to be normalized */
    for (i = 1; i < imax+1; i++)
    for (j = 1; j < jmax+1; j++)
    {
    	if(Problem[i][j] != 0 )
    		Problem[i][j] = 1;
    }


    for(i = 1; i < imax+1; i++)
    for(j = 1; j < jmax+1; j++)
    {
        if(Problem[i][j] == 1) /* all fluid cells to C_F */
        {
        	Flag[i][j] = C_F;
        }
        else   /* 8*east + 4*west + 2*south + 1*north */
        {
         	Flag[i][j] = 8 * Problem[i][j-1] + 4 * Problem[i][j+1] + 2 * Problem[i+1][j] + 1 * Problem[i-1][j];
            /* is flag valid */
            if(Flag[i][j] == 3 || Flag[i][j] == 7 || Flag[i][j] == 11 || Flag[i][j] == 12 || Flag[i][j] == 13 || Flag[i][j] == 14 || Flag[i][j] == 15)
            	return -1;
        }
    }

    for(i = 1; i < imax+1; i++)
    {
        Flag[i][0] = 8 * Problem[i][1];
        Flag[i][jmax+1] = 4 * Problem[i][jmax];
    }

    for(j = 1; j < jmax+1; j++)
    {
        Flag[0][j] = 2 * Problem[1][j];
        Flag[imax+1][j] = 1 * Problem[imax][j];
    }
    return 1;
}
