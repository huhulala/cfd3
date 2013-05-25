/*
 * NSDefinitions.h
 *
 */

#ifndef NSDEFINITIONS_H_
#define NSDEFINITIONS_H_

/* boundary values */
#define NO_SLIP		1
#define FREE_SLIP	2
#define OUTFLOW 	3
/* cell values fluid/obstacle (Bits: center|east|west|south|north),
 * 1 if corresponding cell is a fluid cell, 0 if it is an obstacle */
#define C_F			16
#define C_B			0
#define B_O			8
#define B_W			4
#define B_S			2
#define B_N			1
#define B_SO		10
#define B_NO		9
#define B_SW		6
#define B_NW		5

#endif
