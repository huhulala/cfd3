#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int wl,
  int wr,
  int wt,
  int wb
);

/**
 * Set additional boundary condidtions according to problem
 * Re and dp as extra parameters, because they are needed for the
 * Plane shear flow borders
 */
void spec_boundary_val(
		char *problem,
		int imax,
		int jmax,
		double **U,
		double **V,
		double Re,
		double dp
		);

#endif
