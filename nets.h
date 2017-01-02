// Sun-il Kim and Steve Lumetta
// Modified: 12/12/2006

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_LAMBDA 750
#define HOP_LIMIT 30
#define TOTAL_LPS (N_NODES*(N_NODES-1))/2
#define MAX_STREAMS TOTAL_LPS

#define VERBOSE 0
#define DEBUG 0
#define PRINT_NUM_PATHS 0
#define SHOW_TIME 1

#define ARPANET  0
#define NJ_LATA  1
#define LATA_X   2
#define NATIONAL 3
#define COST     4
#define COST_M   5
#define NSFNET   6
#define TEMP     7
#define SPARS    9

#define NETWORK NJ_LATA

#define LINK_ONLY 1

#define NOT_FOUND    0
#define FOUND        1


/*************************************************************************/
#if (NETWORK==SPARS)
#define N_NODES 15
#define N_EDGES 21
#define MAX_DEG 5
#define MAX_PATHS 20000

static int base_edges[][2] = {
  {0, 2},
  {0, 6},
  {0, 5},
  {0, 12},
  {0, 10},
  {1, 5},
  {1, 6},
  {2, 3},
  {3, 12},
  {4, 5},  /* 10th*/
  {4, 12},
  {5, 6},
  {5, 7},
  {5, 14},
  {6, 7},
  {6, 8},
  {8, 9},
  {9, 10},
  {11, 12},
  {11, 13},
  {13, 14}, /* 20 */
  {-1,-1}
};

#elif (NETWORK==NSFNET)
#define N_NODES 14
#define N_EDGES 20
#define MAX_DEG 4
#define MAX_PATHS 15000

static int base_edges[][2] = {
    { 0, 1},
    { 1, 2},
    { 0, 2},
    { 1, 5},
    { 0, 3},
    { 2, 8},
    { 3, 6},
    { 5, 11},
    { 5, 6},
    { 6, 7},
    { 7, 8},
    { 8, 9},
    { 9, 10},
    { 11, 12},
    { 9, 12},
    { 10, 13},
    { 3, 13},
    { 3, 4},
    { 4, 9},
    { 13, 12},
    {-1, -1}
};


static short node_deg[N_NODES] = {
  3,
  3,
  3,
  4,
  2,
  3,
  3,
  2,
  3,
  4,
  3,
  3,
  3,
  3,
};


#elif (NETWORK == TEMP)

#define N_NODES 3
#define N_EDGES 3
#define MAX_DEG 2
#define MAX_PATHS 2

static int base_edges[][2] = {
    { 0, 1},
    { 1, 2},
    { 2, 0},
    {-1, -1}
};

#elif (NETWORK==ARPANET)

#define N_NODES 20
#define N_EDGES 32
#define MAX_DEG 4
#define MAX_PATHS 2000

static int base_edges[][2] = {
    { 0, 17},
    {17, 13},
    {13, 11},
    {11,  9},
    { 9, 10},
    {10,  7},
    { 7,  6},
    { 6,  5},
    { 5,  4},
    { 4,  2},
    { 2,  1},
    { 1,  3},
    { 3,  8},
    { 8, 12},
    {12, 18},
    {18, 14},
    {14, 15},
    {15, 16},
    {16, 19},
    {19,  0},
    { 8, 10},
    { 9,  8},
    {10, 19},
    {14, 12},
    { 0, 14},
    {19, 17},
    {18, 16},
    { 6,  1},
    { 2,  5},
    { 4, 11},
    {13, 15},
    { 3,  7},

    {-1, -1}
};


#elif (NETWORK==NJ_LATA)

#define N_NODES 11
#define N_EDGES 23
#define MAX_DEG 8
#define MAX_PATHS 1300

static int base_edges[][2] = {
    {0, 1},
    {0, 2}, //{0,7} -> {0,2}
    {0, 4},//{4, 0}, -> {0,4}
    {1, 2},
    {2, 4},
    {3, 7},
    {3 ,4},//{4, 3}, -> {3, 4}
    {3, 0},
    {4, 6},
    {4, 10},
    {4, 8},//{8, 4}, -> {4,8}
    {5, 0},
    {5, 4},//{4, 5}, -> {5, 4}
    {6, 7},
    {7, 5},
    {7, 4},
    {7, 8},
    {7, 0},//{2, 0}, -> {7, 0}
    {7, 10},
    {8, 3},
    {9, 8},
    {10, 9},
    {10, 8},

    {-1, -1}
};

#elif (NETWORK==NJ_LATA_2)

#define N_NODES 11
#define N_EDGES 23
#define MAX_DEG 8
#define MAX_PATHS 1000

static int base_edges[][2] = {
    {0, 1},
    {1, 2},
    {2, 4},
    {3, 7},
    {4, 6},
    {4, 10},
    {5, 0},
    {6, 7},
    {7, 5},
    {8, 3},
    {9, 8},
    {10, 9},
    {4, 7},
    {0, 7},
    {4, 0},
    {4, 3},
    {3, 0},
    {7, 8},
    {4, 8},
    {2, 0},
    {4, 5},
    {7, 10},
    {10, 8},

    {-1, -1}
};


#elif (NETWORK==LATA_X)

#define N_NODES 28
#define N_EDGES 47
#define MAX_DEG 8
#define MAX_PATHS 100000

static int base_edges[][2] = {
    { 0,  1},
    { 1,  2},
    { 2,  3},
    { 3, 12},
    {12, 13},
    {13,  7},
    { 7, 11},
    {11, 14},
    {14, 20},
    {20, 15},
    {15, 16},
    {16, 17},
    {17, 18},
    {18, 19},
    {19, 23},
    {23, 21},
    {21, 22},
    {22, 26},
    {26, 27},
    {27, 25},
    {25, 24},
    {24, 10},
    {10,  9},
    { 9,  6},
    { 6,  5},
    { 5,  8},
    { 8,  4},
    { 4,  0},
    { 1,  6},
    {10, 21},
    {11, 10},
    {21, 20},
    {20, 13},
    { 7,  3},
    { 6,  7},
    {13, 14},
    {25, 26},
    {21, 25},
    {22, 23},
    {23, 15},
    {15, 13},
    {20, 22},
    {18, 15},
    {15, 19},
    {19, 17},
    { 5,  4},
    { 8,  9},

    {-1, -1}
};


#elif (NETWORK==NATIONAL)

#define N_NODES 24
#define N_EDGES 44
#define MAX_DEG 8
#define MAX_PATHS 30000

static int base_edges[][2] = {
    { 0,  9},
    { 9,  1},
    { 1,  5},
    { 5, 10},
    {10, 11},
    {11, 12},
    {12,  4},
    { 4,  2},
    { 2,  3},
    { 3,  8},
    { 8,  7},
    { 7,  6},
    { 6, 13},
    {13, 14},
    {14, 15},
    {15, 16},
    {16, 17},
    {17, 18},
    {18, 19},
    {19, 20},
    {20, 21},
    {21, 22},
    {22, 23},
    {23,  0},
    { 0,  8},
    { 1,  2},
    { 3,  4},
    { 3,  7},
    { 3,  9},
    { 4,  5},
    { 4,  6},
    { 4,  7},
    { 6, 12},
    { 6, 14},
    { 6, 15},
    { 7, 15},
    { 7, 16},
    { 7, 17},
    { 8,  9},
    { 8, 17},
    { 8, 22},
    {10, 12},
    {12, 13},
    {17, 22},
    {-1, -1}
};

#elif (NETWORK==COST)

#define N_NODES 19
#define N_EDGES 37
#define MAX_DEG 7
#define MAX_PATHS 25000
#define MIN_LINKS 20

static int base_edges[][2] = {
  {1,0},
  {0,2},
  {2,8},
  {8,18},
  {18,17},
  {17,16},
  {16,15},
  {15,7},
  {7,5},
  {5,6},
  {6,11},
  {11,10},
  {10,9},
  {9,12},
  {12,14},
  {14,13},
  {13,1},
  {1,3},
  {3,4},
  {4,2},
  {1,2},
  {1,5},
  {2,5},
  {1,11},
  {5,11},
  {2,7},
  {7,8},
  {8,13},
  {8,15},
  {8,17},
  {11,15},
  {11,13},
  {13,15},
  {13,16},
  {1,9},
  {10,12},
  {12,13},

  {-1,-1}
};

#elif (NETWORK==COST_M)

#define N_NODES 11
#define N_EDGES 24
#define MAX_DEG 5
#define MAX_PATHS 15000
#define MIN_LINKS 20

static int base_edges[][2] = {
 {0, 1},
 {0, 3},
 {0, 9},
 {0, 10},
 {1, 2},
 {1, 3},
 {1, 4},
 {2, 3},
 {2, 4},
 {2, 7},
 {2, 9},
 {3, 5},
 {3, 10},
 {4, 5},
 {4, 6},
 {4, 7},
 {5, 6},
 {5, 10},
 {6, 7},
 {6, 8},
 {6, 9},
 {7, 8},
 {8, 9},
 {8, 10},

  {-1,-1}
};


#elif (NETWORK == COST2)
/* this fixes the directions so that
   there is a loop 0-1-3-4-2 */
#define N_NODES 19
#define N_EDGES 37
#define MAX_DEG 7
#define MAX_PATHS 25000
#define MIN_LINKS 20

static int base_edges[][2] = {
  {1,0},
  {0,2},
  {2,4},
  {4,3},
  {3,1},
  {2,1},
  {8,2},
  {8,18},
  {18,17},
  {17,16},
  {16,15},
  {15,7},
  {7,5},
  {5,6},
  {6,11},
  {11,10},
  {10,9},
  {9,12},
  {12,14},
  {14,13},
  {13,1},
  {1,5},
  {2,5},
  {1,11},
  {5,11},
  {2,7},
  {7,8},
  {8,13},
  {8,15},
  {8,17},
  {11,15},
  {11,13},
  {13,15},
  {13,16},
  {1,9},
  {10,12},
  {12,13},

  {-1,-1}
};

#else
#error "unknown network"
#endif

/***********************************************************************************/
typedef int graph_t[N_NODES][MAX_DEG + 1][2];
typedef int edges_t[N_EDGES][2];

#if (N_EDGES <= 32)
typedef long unsigned path_t;
#define PATH_EMPTY(p)       ((*(p)) = 0)
#define PATH_SET(p,n)       ((*(p)) |= (1UL << (n)))
#define PATH_CONTAINS(p,n)  (((*(p)) & (1UL << (n))) != 0)
#define PATHS_OVERLAP(p,q)  (((*(p)) & (*(q))) != 0)
#define PATH_INTERSECT(p,q) ((*(p)) &= (*(q)))
#define PATH_KEEPALL(p,q) ((*(p)) |= (*(q)))
#elif (N_EDGES <= 64)
typedef long long unsigned path_t;
#define PATH_EMPTY(p)      ((*(p)) = 0)
#define PATH_SET(p,n)      ((*(p)) |= (1ULL << (n)))
#define PATH_CONTAINS(p,n) (((*(p)) & (1ULL << (n))) != 0)
#define PATHS_OVERLAP(p,q) (((*(p)) & (*(q))) != 0)
#define PATH_INTERSECT(p,q) ((*(p)) &= (*(q)))
#define PATH_KEEPALL(p,q) ((*(p)) |= (*(q)))
#endif

#if (N_NODES <= 32)
typedef long unsigned nodeset_t;
#define NSET_EMPTY(p)       ((*(p)) = 0)
#define NSET_FLIP(p,n)      ((*(p)) ^= (1UL << (n)))
#define NSET_SET(p,n)       ((*(p)) |= (1UL << (n)))
#define NSET_CLEAR(p,n)     ((*(p)) &= ~(1UL << (n)))
#define NSET_CONTAINS(p,n)  (((*(p)) & (1UL << (n))) != 0)
#define NSETS_OVERLAP(p,q)  (((*(p)) & (*(q))) != 0)
#define NSET_INTERSECT(p,q) ((*(p)) &= (*(q)))
#define NSET_EQUAL(p,q)     ((*(p)) == (*(q)))
#define NSET_KEEPALL(p,q) ((*(p)) |= (*(q)))
#else
typedef long long unsigned nodeset_t;
#define NSET_EMPTY(p)       ((*(p)) = 0)
#define NSET_FLIP(p,n)      ((*(p)) ^= (1ULL << (n)))
#define NSET_SET(p,n)       ((*(p)) |= (1ULL << (n)))
#define NSET_CLEAR(p,n)     ((*(p)) &= ~(1ULL << (n)))
#define NSET_CONTAINS(p,n)  (((*(p)) & (1ULL << (n))) != 0)
#define NSETS_OVERLAP(p,q)  (((*(p)) & (*(q))) != 0)
#define NSET_INTERSECT(p,q) ((*(p)) &= (*(q)))
#define NSET_EQUAL(p,q)     ((*(p)) == (*(q)))
#define NSET_KEEPALL(p,q) ((*(p)) |= (*(q)))
#endif



typedef struct {
  short hops;
  short links[N_EDGES];   /* links# */
  short dir[N_EDGES];     /* -1 unset, 0 & 1-same as edges[] */
} route_t;

typedef struct {
  long unsigned nodes_seen;
  short n_num;
  short depth;
  path_t path;
  route_t route;
} queue_t;

typedef struct {
  short s, d; /* source and destination nodes */
  short color_p; /* wavelength used */
  short color_b; /* wavelength used */
  path_t path_p, path_b;
  short dir;
} conn_t;
