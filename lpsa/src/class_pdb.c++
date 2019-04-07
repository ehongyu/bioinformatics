/* ==== The function of class PDBdata was defined in this file. ==== */

#include "head.h"

PDBdata::PDBdata()
{
  bond_num=0;
  for(register int i=0;i<4;i++) bond[i]=0;	//for loop maybe wrong !!!!
} 

void PDBdata::getInfo(char atomName[7],char resName[5],float coor[3])
{
  strcpy(resName, res);
  strcpy(atomName, atom);
  coor[0]=xyz[0];
  coor[1]=xyz[1];
  coor[2]=xyz[2];
}

void PDBdata::storeInfo(char atomName[7],char resName[5],float coor[3])
{
  strcpy(res, resName);
  strcpy(atom,atomName);
  xyz[0]=coor[0];
  xyz[1]=coor[1];
  xyz[2]=coor[2];
}

void PDBdata::link(int i)
{
  bond[bond_num++]=i;	//x++ not ++x
}

void PDBdata::cut(int j)
{
  for(register int i=0;i<4;i++){
     if(bond[i]==j){
	bond[i]=0;
  	bond_num--;
	for(register int k=i;k<3;k++) bond[k]=bond[k+1];
     }
  }
}
	
//inline PDBdata PDBdata::operator=(PDBdata pp)
PDBdata PDBdata::operator=(PDBdata pp)
{
  strcpy(atom,pp.atom);
  strcpy(res,pp.res);
  xyz[0]=pp.xyz[0]; 
  xyz[1]=pp.xyz[1]; 
  xyz[2]=pp.xyz[2]; 
  return *this;
}

//inline float PDBdata::operator-(PDBdata pp)	//square of distance
float PDBdata::operator-(PDBdata pp)	//square of distance
{
  return        (xyz[0]-pp.xyz[0])*(xyz[0]-pp.xyz[0])+ 
		(xyz[1]-pp.xyz[1])*(xyz[1]-pp.xyz[1])+ 
		(xyz[2]-pp.xyz[2])*(xyz[2]-pp.xyz[2]); 
}

//inline float bondlength(const PDBdata aa, const PDBdata bb) //up standard
float bondlength(const PDBdata aa, const PDBdata bb) //up standard
{
  return aa.atom[2]=='S'||bb.atom[2]=='S'? 3.61:2.89;
}

ofstream &operator<<(ofstream &ostrm, PDBdata pp) //output of PDB data
{
  static int i=0;       // atom num
  static short int j=0; // residue num
  i++;
  if(!strcmp(pp.atom,"  N   ")||!strcmp(pp.atom,"  C1  ")) j++;
  ostrm<<"ATOM";
  ostrm<<setiosflags(ios::right);
  ostrm<<setw(7)<<i;
  ostrm<<pp.atom<<pp.res;
  ostrm<<setw(5)<<j<<"    ";
  ostrm.setf(ios::showpoint|ios::fixed);
  ostrm<<setprecision(3)<<setw(8)<<pp.xyz[0];
  ostrm<<setprecision(3)<<setw(8)<<pp.xyz[1];
  ostrm<<setprecision(3)<<setw(8)<<pp.xyz[2];
  ostrm<<"  1.00  0.00     PDB1"<<endl;
  return ostrm;

}

ifstream &operator>>(ifstream &istrm, PDBdata &pp) //input of PDB
{
  istrm.ignore(2);
  istrm.getline(pp.atomNum,6);
  istrm.getline(pp.atom,7);
  istrm.getline(pp.res,5);
  istrm.getline(pp.chain,2);
  istrm.getline(pp.resNum,5);
//  istrm.ignore(10);
  istrm>>pp.xyz[0]>>pp.xyz[1]>>pp.xyz[2];
  istrm.ignore(80,'\n');
  return istrm;
}

short int getlink(PDBdata pp, int bd[4])
{
  if(pp.bond_num){
    for(register int i=0;i<4;i++) bd[i]=pp.bond[i];
  }
  return pp.bond_num;
}
  
void Rotate(float theta, const PDBdata aa, const PDBdata bb,
			 const PDBdata _cc, PDBdata &cc)
{
  float _rot[3];
  float _coor[3];
  float coor[3];
  float _tran[3];

  for(register int i=0;i<3;i++) _tran[i]=aa.xyz[i];
  for(i=0;i<3;i++) _rot[i]=bb.xyz[i]-aa.xyz[i];
  float rr=_rot[0]*_rot[0]+_rot[1]*_rot[1]+_rot[2]*_rot[2];
  rr=fsqrt(rr);
  for(i=0;i<3;i++) _rot[i]=_rot[i]/rr;
  for(i=0;i<3;i++) _coor[i]=_cc.xyz[i];
  RotateCoor(theta, _rot, _coor, coor, _tran);
  for(i=0;i<3;i++) cc.xyz[i]=coor[i];
}

void MultipleMatrix(float *p, double mat[3][3], double t[3])
{
        int i,j;

        for(i=0;i<=2;i++){
                t[i]=0.0;
                for(j=0;j<=2;j++)
                        t[i]=t[i]+*(p+j)*mat[j][i];
        }
}

void RotateCoor(float theta, float _rot[3], float _coor[3],
		float coor[3], float _tran[3])
{
        double mat[3][3],p[3];
        double a,b;
        int j;
	
        a=sin(theta/CONV);
        b=cos(theta/CONV);
        /* translate the coordinate */
        for(j=0;j<=2;j++)
                coor[j]=_coor[j]-_tran[j];

        mat[0][0]=_rot[0]*_rot[0]+(1-_rot[0]*_rot[0])*b;
        mat[1][0]=_rot[0]*_rot[1]*(1-b)-_rot[2]*a;
        mat[2][0]=_rot[0]*_rot[2]*(1-b)+_rot[1]*a;
        mat[0][1]=_rot[0]*_rot[1]*(1-b)+_rot[2]*a;
        mat[1][1]=_rot[1]*_rot[1]+(1-_rot[1]*_rot[1])*b;
        mat[2][1]=_rot[1]*_rot[2]*(1-b)-_rot[0]*a;
        mat[0][2]=_rot[0]*_rot[2]*(1-b)-_rot[1]*a;
        mat[1][2]=_rot[1]*_rot[2]*(1-b)+_rot[0]*a;
        mat[2][2]=_rot[2]*_rot[2]+(1-_rot[2]*_rot[2])*b;

	MultipleMatrix(&coor[0],mat,p);

        /* reverse the translation */
        for(j=0;j<=2;j++)
                coor[j]=p[j]+_tran[j];
}

float torsionAngle(const PDBdata aa, const PDBdata bb,
	const PDBdata cc, const PDBdata dd)
//this subroutine is translated from FORTRAN code of Deng QiaoLing
{
  float    xij=cc.xyz[0]-aa.xyz[0];
  float    yij=cc.xyz[1]-aa.xyz[1];
  float    zij=cc.xyz[2]-aa.xyz[2];
  float    xkj=bb.xyz[0]-aa.xyz[0];
  float    ykj=bb.xyz[1]-aa.xyz[1];
  float    zkj=bb.xyz[2]-aa.xyz[2];
  float    xkl=bb.xyz[0]-dd.xyz[0];
  float    ykl=bb.xyz[1]-dd.xyz[1];
  float    zkl=bb.xyz[2]-dd.xyz[2];
  float    dx=yij*zkj-zij*ykj;
  float    dy=zij*xkj-xij*zkj;
  float    dz=xij*ykj-yij*xkj;
  float    gx=ykj*zkl-zkj*ykl;
  float    gy=zkj*xkl-xkj*zkl;
  float    gz=xkj*ykl-ykj*xkl;
  float    bi=dx*dx+dy*dy+dz*dz;
  float    bk=gx*gx+gy*gy+gz*gz;
  float    ct=dx*gx+dy*gy+dz*gz;
  if(bi<0.0001) bi=0.0001;
  if(bk<0.0001) bk=0.0001;
  bi=sqrtf(bi);
  bk=sqrtf(bk);
  float    z1=1./bi;
  float    z2=1./bk;

// ct is the cosine of the angle between the normals to the two atom planes
  ct=ct*z1*z2;
  if(ct>1.) ct=1.;
  if(ct<(-1.)) ct=-1.;
// ap is the angle, in radians, between the two planes
  float   ap=acosf(ct);
/*  --------------------------------------------------------------------  C
  the vector perpendicular to the normals to the two planes is compared
  with the direction of the central bond vector to determine the sign of 
  the torsion
C  -------------------------------------------------------------------- */ 
  float si=xkj*(dy*gz-dz*gy)+ykj*(dz*gx-dx*gz)+zkj*(dx*gy-dy*gx);
  if(si<0.) ap=-ap;

  ap=ap*CONV;
  return ap;
}

//inline float energy(PDBdata aa, PDBdata bb)
float energy(PDBdata aa, PDBdata bb)
{
     float ra, rb;
     switch(aa.atom[2]) {
           case 'C':
                ra=1.70;
                break;
           case 'N':
                ra=1.55;
                break;
           case 'O':
                ra=1.52;
                break;
           default:
                ra=1.80;
                break;
     }

     switch(bb.atom[2]) {
           case 'C':
                rb=1.70;
                break;
           case 'N':
                rb=1.55;
                break;
           case 'O':
                rb=1.52;
                break;
           default:
                rb=1.80;
                break;
     }

     float rr0=(ra+rb)*(ra+rb);
     float rr=aa-bb;
     return rr>rr0? 0:K_OVERLAP*(rr0-rr);
}

float Aenergy(PDBdata aa, PDBdata bb)
{
     float rr=aa-bb;
     return K_ANCHOR*rr;
}

PDBdata PDBdata::operator<=(PDBdata bb)	//vector from bb to ...
{
  PDBdata temp=*this;
  for(int i=0;i<3;i++) temp.xyz[i]=xyz[i]-bb.xyz[i];
  return temp;
}

PDBdata PDBdata::operator+(PDBdata bb)	//vectors adding
{
  PDBdata temp;
  for(int i=0;i<3;i++) temp.xyz[i]=this->xyz[i]+bb.xyz[i];
  return temp;
}

float angle(PDBdata aa, PDBdata bb, PDBdata cc)	//angle <abc
{
  PDBdata a,b;
  a=aa<=bb;
  b=cc<=bb;
  return (acosf((a.xyz[0]*b.xyz[0]+a.xyz[1]*b.xyz[1]+a.xyz[2]*b.xyz[2])/
	 	sqrtf(a.xyz[0]*a.xyz[0]+a.xyz[1]*a.xyz[1]+a.xyz[2]*a.xyz[2])/
		sqrtf(b.xyz[0]*b.xyz[0]+b.xyz[1]*b.xyz[1]+b.xyz[2]*b.xyz[2]))); 
} 

PDBdata PDBdata::operator*(PDBdata b)	//vector multiply
{
  PDBdata temp;

  temp.xyz[0]=this->xyz[1]*b.xyz[2]-this->xyz[2]*b.xyz[1];
  temp.xyz[1]=-this->xyz[0]*b.xyz[2]+this->xyz[2]*b.xyz[0];
  temp.xyz[2]=this->xyz[0]*b.xyz[1]-this->xyz[1]*b.xyz[0];
  return temp;
}

PDBdata center(PDBdata aa, PDBdata bb)
{
  PDBdata temp=aa+bb;
  for(int i=0;i<3;i++) temp.xyz[i]=temp.xyz[i]/2;
  return temp;
}

PDBdata nextAtom(PDBdata p1,PDBdata p2,PDBdata p3)
{
PDBdata temp;
float afa, bet, r01, r12, r13;

if(strcmp(p1.atom,"  N   "))
{
     afa=(180-123.5)/CONV;
     bet=(180-113.9)/CONV;
     r01=1.345;
     r12=1.225;
     r13=1.515;
     strcpy(temp.atom,"  N   ");
}
else
{
     afa=(180-117.6)/CONV;
     bet=(180-152.9)/CONV;
     r01=1.46;
     r12=1.345;
     r13=2.4;
     strcpy(temp.atom,"  CA  ");
}

float x12=p1.xyz[0]-p2.xyz[0];
float y12=p1.xyz[1]-p2.xyz[1];
float z12=p1.xyz[2]-p2.xyz[2];
float x13=p1.xyz[0]-p3.xyz[0];
float y13=p1.xyz[1]-p3.xyz[1];
float z13=p1.xyz[2]-p3.xyz[2];

float A=(r01*r12*y13*cos(afa)-r01*r13*y12*cos(bet))/(x12*y13-x13*y12);
float B=(y12*z13-y13*z12)/(x12*y13-x13*y12);
float C=(r01*r12*x13*cos(afa)-r01*r13*x12*cos(bet))/(y12*x13-y13*x12);
float D=(z12*x13-z13*x12)/(x12*y13-x13*y12);

float A1=y12*z13-y13*z12;
float B1=x13*z12-x12*z13;
float C1=x12*y13-x13*y12;

float z01=(-A*A1-B1*C)/(A1*B+B1*D+C1);
float x01=A+B*z01;
float y01=C+D*z01;

temp.xyz[0]=p1.xyz[0]+x01;
temp.xyz[1]=p1.xyz[1]+y01;
temp.xyz[2]=p1.xyz[2]+z01;

return temp;

}
