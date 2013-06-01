#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include "fd.h"
#include <math.h>
#include "NSDefinitions.h"


void calculate_fg(double Re, double GX, double GY,double alpha, double dt,double dx,double dy, int imax, int jmax,
		double **U,double **V,double **F,double **G, int **Flag)
{
	/* see formulas 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		/********** calculate F **********/
		/* calculate f and g only between two fluid cells */
		if ((i <= imax - 1) && (Flag[i+1][j] == C_F && Flag[i][j] == C_F))
		{
				F[i][j] = U[i][j] + dt * (
				/* 1/Re * (d²u/dx² + d²u/dy²) */
				1 / Re * (d2udx2(i, j, U, dx) + d2udy2(i, j, U, dy))
				/* - du²/dx */
				- du2dx(i, j, U, dx, alpha)
				/* - duv/dy */
				- duvdy(i, j, U, V, dy, alpha) + GX);
		}
		else
			/* F aus 1.4 bis 1.6  */
			F[i][j] = U[i][j];
		/********** calculate G **********/
		if ((j <= jmax - 1) && (Flag[i][j+1] == C_F && Flag[i][j] == C_F))
		{
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY);
		}
		/* G aus 1.4 bis 1.6  */
		else
			G[i][j] = V[i][j];
	}

	/* calculate boundary values -  see formula 17 */
	for (j = 1; j <= jmax; j++)
	{
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (i = 1; i <= imax; i++)
	{
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V) {
	/* See formula 13 */
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
    double dtcon, dxcon, dycon;
    double min;
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		if (fabs(U[i][j]) > umax)umax = fabs(U[i][j]);
		if (fabs(V[i][j]) > vmax)vmax = fabs(V[i][j]);
	}

	/* conditions */
	dtcon = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
	dxcon = dx/fabs(umax);
	dycon = dy/fabs(vmax);

	/* determine smalles condition */
    min = dtcon;
	if(min > dxcon)min = dxcon;
	if(min > dycon) min = dycon;

	/* calculate dt */
	*dt = tau * min;
}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int **Flag
)
{
    int i, j;
    double dtodx, dtody;

    /******** Calculate dt/dx and dt/dy (it is the same for each element) ********/
    dtodx = dt/dx;
    dtody = dt/dy;

    /******** Calculate u in step next step ********/
    for(i = 1; i < imax; i++)
    {
        for(j = 1; j < jmax+1; j++)
        {
            if (Flag[i][j] >= C_F && Flag[i+1][j] >= C_F){
                 U[i][j] = F[i][j] - dtodx*(P[i+1][j] - P[i][j]);
            }
        }
    }

    /******** Calculate v in step next step ********/
    for(i = 1; i < imax + 1; i++)
    {
        for(j = 1; j < jmax; j++)
        {
            if (Flag[i][j] >= C_F && Flag[i][j+1] >= C_F){
                V[i][j] = G[i][j] - dtody*(P[i][j+1] - P[i][j]);
            }
        }
    }
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
  double **F, double **G,double **RS, int **Flag)
{
	/* right hand side of formula 11 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		if (Flag[i][j] == C_F && Flag[i][j+1] == C_F)
			RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j]
					   - G[i][j - 1]) / dy);
	}
}
