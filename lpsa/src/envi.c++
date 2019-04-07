#include "head.h"

int In_neighbor(float coor[3],float gr00[3],float gr11[3])
{

  float R2=Dvdw*Dvdw;
  float z2=(coor[2]-gr00[2])*(coor[2]-gr00[2]);
  if(z2>R2) return 0;

  float x00_2=(coor[0]-gr00[0])*(coor[0]-gr00[0]);
  float y00_2=(coor[1]-gr00[1])*(coor[1]-gr00[1]);
  float x11_2=(coor[0]-gr11[0])*(coor[0]-gr11[0]);
  float y11_2=(coor[1]-gr11[1])*(coor[1]-gr11[1]);

  if(
   (gr00[0]<coor[0])&&(gr11[0]>coor[0])&&(gr00[1]<coor[1])&&(gr11[1]>coor[1]) 
	|| ((z2+x00_2)<R2&&gr00[1]<coor[1]&&gr11[1]>coor[1])
	|| ((z2+y00_2)<R2&&gr00[0]<coor[0]&&gr11[0]>coor[0])
	|| ((z2+x11_2)<R2&&gr00[1]<coor[1]&&gr11[1]>coor[1])
	|| ((z2+y11_2)<R2&&gr00[0]<coor[0]&&gr11[0]>coor[0])
	|| (x00_2+y00_2+z2)<R2
	|| (x00_2+y11_2+z2)<R2
	|| (x11_2+y00_2+z2)<R2
	|| (x11_2+y11_2+z2)<R2)
   return 1;

  return 0;

}
