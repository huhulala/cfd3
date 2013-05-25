#include "uvp.h"
#include <math.h>
#include "fd.h"
#include "NSDefinitions.h"

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G, int **Flag) {
	/* see formula 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		/********** calculate F **********/
		/* calculate f and g only between two fluid cells */
		if (i <= imax - 1 && ((Flag[i][j] & 24) == 24))
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
			/* F aus 1.4 bis 1.6 TODO: ? */
			F[i][j] = U[i][j];
		/********** calculate G **********/
		if (j <= jmax - 1 && ((Flag[i][j] & 17) == 17))
		{
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, dy, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY);
		}
		/* G aus 1.4 bis 1.6 TODO: ? */
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

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS) {

	/* right hand side of formula 11 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j]
					- G[i][j - 1]) / dy);
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

void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P, int **Flag) {
	/* see formula 7 and 8 */
	int i,j;
	double dtdx = dt / dx;
	double dtdy = dt / dy;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		/* claculate for neighbouring fluid cells */
		if (i <= imax - 1 && ((Flag[i][j] & 24) == 24))
		{
				U[i][j] = F[i][j] - dtdx * (P[i + 1][j] - P[i][j]);
		}
		if (j <= jmax - 1 && ((Flag[i][j] & 17) == 17))
		{
			V[i][j] = G[i][j] - dtdy * (P[i][j + 1] - P[i][j]);
		}
	}
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
