#ifndef CONSTANTS
#define CONSTANTS

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** Infinity as used in minimization routines */
#define INF 10000000 /* (INT_MAX/10) */

#define EMAX (INF/10)
/** forbidden */
#define FORBIDDEN 9999
/** bonus contribution */
#define BONUS 10000
/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30

/** The maximum interior loop length */
#define MAX_INTLOOP  (static_cast<unsigned int>(std::min(30, MAXLOOP)))

#define UNIT 100

#define MINPSCORE -2 * UNIT

#endif