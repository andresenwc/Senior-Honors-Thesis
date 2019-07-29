/* config.h */

// Rectangular dimensions of the mic array
#define X_DIM 32.0
#define Y_DIM 32.0

// mic0 coordinates
#define M0X 0.0
#define M0Y 0.0
// mic1 coordinates
#define M1X 0.0
#define M1Y 32.0
// mic2 coordinates
#define M2X 32.0
#define M2Y 0.0
// mic3 coordinates
#define M3X 32.0
#define M3Y 32.0

// file to be parsed
#define FILENAME "2D-01-3-samples-per-point.log"

// Must define one (and only one) of the three signal processing
// algorithms. Documentation for each is within ImpactFinder.C.
//#define FIRST_SUSTAINED_HIT
//#define FIRST_HIT
#define CROSS_CORRELATION

// Must define one (and only one) of the two final choice
// algorithms. Outliers can be removed by DBSCAN (CLUSTER)
// or standard deviation (SD).
#define CLUSTER
//#define STAT

// Must define one (and only one) speed of sound.
// 700.0 for CROSS_CORRELATION or FIRST_HIT
//#define SPEEDOFSOUND 700.0
// 350.0 for FIRST_SUSTAINED_HIT
#define SPEEDOFSOUND 700.0

// Epsilon dictates at what point to stop bisection
#define EPSILON .00001