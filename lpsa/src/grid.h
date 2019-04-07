#include head.h

// Parameters of environment grid

#define GRID_SIZE 1.0	// size of a grid, unit: A
// GRID_SIZE should be related to MAX_INCLUDE, the maximum number of atoms
// that could be mapped onto a single grid, roughtly, 1.0A-40, 2.0A-60, etc.

#define MAX_INCLUDE 40	// maximum number of atoms mapped onto a single grid
#define GRIDMAX_X 47	// maximum number of grids in X axis direction
#define GRIDMAX_Y 43	// maximum number of grids in Y axis direction
#define GRIDMAX_Z 42	// maximum number of grids in Z axis direction
#define Rvdw 1.8	// maximum van der Waal radius
#define Dvdw 3.6	// maximum van der Waal overlaping distance between 
			// any atomic pairs
				
/* grid
   ----
   define a grid to map all protein atoms onto it
*/

class grid
{
public:
  float x,y,z;
  unsigned int include_num;
  unsigned int include[MAX_INCLUDE];

  grid();
  void set_coor(float a, float b, float c);
  void add(int i);
};

