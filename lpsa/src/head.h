/*============================================================================

                    #       #####    ####     ##
                    #       #    #  #        #  #
                    #       #    #   ####   #    #
                    #       #####        #  ######
                    #       #       #    #  #    #
                    ######  #        ####   #    #


      modeling Loops of Proteins by Simulated Annealing method

                     WRITTEN BY DR. HONGYU ZHANG

                     CARB, University of Maryland

                           August, 1998

                           VERSION 1.0b

============================================================================*/



#include "string.h"	// for string manipulate
#include "stdlib.h"	// for standard library function
#include "math.h"	// for math
#include "fstream.h"	// for file i/o
#include "iostream.h"	// for i/o
#include "iomanip.h"	// for i/o format
#include "time.h"	// for random seed initialization

#define MAX_ATOM 1500	// maximum protein atom number in system
#define K_OVERLAP 10.0	// force weights for soft-sphere potential
#define K_ANCHOR 100.0	// force weights for anchoring potential

#define GRID_SIZE 1.0	// size of a grid, unit: A
// GRID_SIZE should be related to MAX_INCLUDE, the maximum number of atoms
// that could be mapped onto a single grid, roughtly, 1.0A-40, 2.0A-60, etc.

#define MAX_INCLUDE 40	// maximum number of atoms mapped onto a single grid
#define GRIDMAX_X 45	// maximum number of grids in X axis direction
#define GRIDMAX_Y 55	// maximum number of grids in Y axis direction
#define GRIDMAX_Z 48	// maximum number of grids in Z axis direction
#define Rvdw 1.8	// maximum van der Waal radius
#define Dvdw 3.6	// maximum van der Waal overlaping distance between 
			// any atomic pairs
				
#define MAX_RES int(MAX_ATOM/5)	// maximum residue number of a protein, 
//	roughtly, it will be less than MAX_ATOM divided by 5
#define MAX_SIDE_TOR 5	// maximum mumber of side-chain 
			// torsion angles of any residue
#define INIT_TEMP 50.	// intial temperture of annealing
#define MAX_LP_LENGTH int(MAX_ATOM/2)	// maximum number of loop atoms
#define MAX_TORSION 100	// maximum torsion angles in loop
#define MAX_NONBOND_PAIRS 100000	// maximum nonbonded atomic pair number
#define ROT_NUM 500	// maximum number of loop atoms to be rotated
#define RADIUS 9.*9.	// not used anymore
#define QUEUE_SIZE 150	// maximum queue size for class que
#define MAX_NEIGHBOR 100	// not used anymore
#define MAX(a,b) a>b? a:b
#define MIN(a,b) a>b? b:a

#define CONV (180./acos(-1.))
#define GROWRATE 10		//growing rate of queue
#define STEP 180.		//step of torsion angle in annealing


/* ---------- CLASS DEFINITION --------------- */ 

/* PDBdata
   ------- 
   include PDB data and functions(I/O, torsion angle calculating,
   linkage, operator etc., where the  rotating function is provided 
   by Luo Zhaowen
*/

class PDBdata
{
  char    atom[7];
  char    res[5];
  char	atomNum[6];
  char	resNum[5];
  char	chain[2];
  float   xyz[3];
  short int     bond_num;
  int	bond[4]; 

public:
  PDBdata();				//constructor

  void link(int i);		//link two atoms into a bond

  void cut(int j);		//cut off a bond

  PDBdata operator=(PDBdata pp);	//project equal

  float operator-(PDBdata pp);	//distance between two atoms

  friend float bondlength(const PDBdata aa, const PDBdata bb);//standard bondlen

  friend ofstream &operator<<(ofstream &ostrm, PDBdata pp);  //overlode <<, >>

  friend ifstream &operator>>(ifstream &istrm, PDBdata &pp);

  friend short int getlink(PDBdata pp, int bd[4]); //link relations of pp

  // calculate angle, aa is bond_start, bb is end, cc is pre_bond, dd is next.
  friend float torsionAngle(const PDBdata aa, const PDBdata bb,
	const PDBdata cc, const PDBdata dd);

  // rotating function, aa is bond_start, bb is end, _cc is input, cc is output
  friend void Rotate(float theta, const PDBdata aa, const PDBdata bb,
	const PDBdata _cc, PDBdata &cc);

  friend float energy(PDBdata aa, PDBdata bb);	//overlap energy

  friend float Aenergy(PDBdata aa, PDBdata bb);	//anchor energy

  void getInfo(char atomName[7],char resName[5],float coor[3]);

  void storeInfo(char atomName[7],char resName[5],float coor[3]);

  PDBdata operator<=(PDBdata bb);//vectors form atom bb to ...

  PDBdata operator+(PDBdata bb);//vectors adding

  //angle between two vectors
  friend float angle(PDBdata aa, PDBdata bb, PDBdata cc);

  PDBdata operator*(PDBdata b);	//vector multiply

  friend PDBdata center(PDBdata aa, PDBdata bb);

  //generating next loop atom coordinates
  friend PDBdata nextAtom(PDBdata p1,PDBdata p2,PDBdata p3);

};

/* queue
   -----
Define a queue.
*/

class queue
{
  short int spos,rpos;
  int que[QUEUE_SIZE];
public:
  queue();
  void init();
  void store(int i);
  int retrieve();
};

/*	RANDOM NUMBER GENERATOR    GGL
	original reference for algorithm
	P. Lewis, A. Goodman, J. Miller,IBM Sys. J., 2, 136 (1969)
*/

typedef double TP;
class GGL
{
  TP Seed;
public:
  GGL(TP initSeed=31415926);
  void SetSeed(TP newSeed=(TP)time(NULL));
  TP operator()();
};

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


/* loop
   ----
Monte-Carlo simulated annealing to generate a candidate conformation for a loop.
*/

class loop
{
  PDBdata pdb[MAX_ATOM];
  PDBdata anchor[2];
  float vertex[3][2];	//vertex of grid cube
//  double Etotal[3];
  grid 	gr[GRIDMAX_X][GRIDMAX_Y][GRIDMAX_Z];	// grid matrix
  int	torsion[MAX_TORSION][5], torsion_num,	  // torsion angle
  	lp_start, lp_end, csnec,	// loop terminal position
	random0,random1,random2,	//random seed
	nbPair[MAX_NONBOND_PAIRS][2], nbPair_num, //all nonbond pair atoms
        rot[MAX_TORSION][ROT_NUM],      //to be rotated atoms for each rotating
        fix[MAX_TORSION][ROT_NUM],      //to be fixed   atoms for each rotating
        rot_num[MAX_TORSION], fix_num[MAX_TORSION],
	np[MAX_RES],	//sequence position of N
	total_atom_num;	//total atom number in PDB

public:
  loop();
  void input(); 	// input data 
  void NPos();		// determine N position of each residue
  void environment(); 	// find out environment atoms
  void torsionGen(); 	// generation of torsion angle array
  void construct(); 	// construct initial conformation
  void connect(); 	// identifing bond connections
  void nonbond(); 	// find out nonbond pairs 
//  void TotalEnergy(); 	// evaluate total overlap energy
  void anneal(); 	// Monte Carlo Annealing 
};

/* ------------- public function definition ---------------------- */

void MultipleMatrix(float *p, double mat[3][3], double t[3]);
void RotateCoor(float theta, float _rot[3], float _coor[3],
		float coor[3], float _tran[3]);
int In_neighbor(float coor[3],float gr00[3],float gr11[3]);
