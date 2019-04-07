/* ================ Function of class loop ==================== */
#include "head.h"


/* ----------- input data ----------------- */

loop::loop()
  {}
void loop::input()
{
  ifstream in1("INNAME");
  ifstream in2("LOOPSE");
  ifstream in3("RANDOM");

  if(!in1)
  {
    cout<<"Cannot open file INNAME\n";
    exit(1);
  }

  if(!in2)
  {
    cout<<"Cannot open file LOOPSE\n";
    exit(1);
  }

  if(!in3)
  {
    cout<<"Cannot open file RANDOM\n";
    exit(1);
  }

  char begin[5], buf[200];
  for(register int i=0;i<MAX_ATOM&&!in1.eof();i++) {
	in1.getline(begin, 5);
	if(!strcmp(begin,"ATOM")) {
//		in1.seekg(-4, ios::cur);
		in1>>pdb[i];
	}
	else {
		in1.getline(buf, 200);
		if(!in1.eof()) i--;
	}
  }
  total_atom_num=--i;

  in2>>lp_start>>lp_end>>csnec;

  in3>>random0>>random1>>random2;
  in1.close();
  in2.close();
  in3.close();

}

/* ----------- determine environmental atoms by a grid -------------- */

void loop::environment()
{
float coor[3],coor1[3];
char Name[7];
int grid_num[3];

ifstream gp("GRID_PARA");
// generate grid parameter file if not existed
if(!gp)
 {

 ofstream gp("GRID_PARA");
  pdb[0].getInfo(Name, Name, coor);
  vertex[0][0]=coor[0];
  vertex[0][1]=coor[0];
  vertex[1][0]=coor[1];
  vertex[1][1]=coor[1];
  vertex[2][0]=coor[2];
  vertex[2][1]=coor[2];

// Firstly, select a box that could contain all the environmental atoms,
// i.e., those of protein atoms except loop atoms and the two or three anchor
// atoms( here is two atoms N and CA).

  for(int i=0;i<total_atom_num;i++){

    if(i>=np[lp_start-1]&&i<np[lp_end]+2)  continue;	//If i in loop, ignore!
    pdb[i].getInfo(Name, Name, coor);
    for(int j=0;j<3;j++){
	if( coor[j]-Dvdw<vertex[j][0] ) vertex[j][0]=coor[j]-Dvdw;
	if( coor[j]+Dvdw>vertex[j][1] ) vertex[j][1]=coor[j]+Dvdw;
    }
  }
  for(i=0;i<3;i++) grid_num[i]=int((vertex[i][1]-vertex[i][0])/GRID_SIZE)+1;
    // to plus 2 means considering the one more shell for the grid 
    // --- changed in Jan 8, 1996
  if(grid_num[0]>GRIDMAX_X||grid_num[1]>GRIDMAX_Y||grid_num[2]>GRIDMAX_Z)
  {
    cout<<" grid number overflow! GRIDMAX in head.h need to be enlarged "<<endl;
    cout<<" Given GRID_SIZE = "<<GRID_SIZE <<" anstrom, it's necessary that:\n";
    cout<<"	GRIDMAX_X > "<< grid_num[0] <<endl;
    cout<<" 	GRIDMAX_Y > "<< grid_num[1] <<endl;
    cout<<" 	GRIDMAX_Z > "<< grid_num[2] <<endl;
    exit(1);
  }
  for(i=0;i<grid_num[0];i++) {
    for(int j=0;j<grid_num[1];j++) {
      for(int k=0;k<grid_num[2];k++) {
	gr[i][j][k].set_coor(GRID_SIZE*i+vertex[0][0], 
	GRID_SIZE*j+vertex[1][0], GRID_SIZE*k+vertex[2][0]); 
      }
    }
  }

// Secondly, map all evironmental atoms into grid

  int grid_x,grid_y,grid_z;
  for(i=0;i<total_atom_num;i++){

    if(i>=np[lp_start-1]&&i<np[lp_end]+2)  continue;	//whether i is in loop
    pdb[i].getInfo(Name, Name, coor);

    for(grid_x=0;grid_x<grid_num[0]-1;grid_x++){
     for(grid_y=0;grid_y<grid_num[1]-1;grid_y++){
      for(grid_z=0;grid_z<grid_num[2]-1;grid_z++){

     //for the six cube plane of the grid, do
	float gr00[3],gr11[3]; //two diagonal vertex of a rectangular plane

     //001
	gr00[0]=gr[grid_x][grid_y][grid_z+1].x;
	gr00[1]=gr[grid_x][grid_y][grid_z+1].y;
	gr00[2]=gr[grid_x][grid_y][grid_z+1].z;
	gr11[0]=gr[grid_x+1][grid_y+1][grid_z+1].x;
	gr11[1]=gr[grid_x+1][grid_y+1][grid_z+1].y;
	gr11[2]=gr[grid_x+1][grid_y+1][grid_z+1].z;
	if(In_neighbor(coor,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

     //00-1
	gr00[0]=gr[grid_x][grid_y][grid_z].x;
	gr00[1]=gr[grid_x][grid_y][grid_z].y;
	gr00[2]=gr[grid_x][grid_y][grid_z].z;
	gr11[0]=gr[grid_x+1][grid_y+1][grid_z].x;
	gr11[1]=gr[grid_x+1][grid_y+1][grid_z].y;
	gr11[2]=gr[grid_x+1][grid_y+1][grid_z].z;
	if(In_neighbor(coor,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

     //100
	gr00[0]=gr[grid_x+1][grid_y][grid_z].x;
	gr00[1]=gr[grid_x+1][grid_y][grid_z].y;
	gr00[2]=gr[grid_x+1][grid_y][grid_z].z;
	gr11[0]=gr[grid_x+1][grid_y+1][grid_z+1].x;
	gr11[1]=gr[grid_x+1][grid_y+1][grid_z+1].y;
	gr11[2]=gr[grid_x+1][grid_y+1][grid_z+1].z;
	//swap x<->z
	coor1[0]=coor[2];
	coor1[1]=coor[1];
	coor1[2]=coor[0];
	float tmp;
	tmp=gr00[2];
	gr00[2]=gr00[0];
	gr00[0]=tmp;
	tmp=gr11[2];
	gr11[2]=gr11[0];
	gr11[0]=tmp;
	if(In_neighbor(coor1,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

     //-100
	gr00[0]=gr[grid_x][grid_y][grid_z].x;
	gr00[1]=gr[grid_x][grid_y][grid_z].y;
	gr00[2]=gr[grid_x][grid_y][grid_z].z;
	gr11[0]=gr[grid_x][grid_y+1][grid_z+1].x;
	gr11[1]=gr[grid_x][grid_y+1][grid_z+1].y;
	gr11[2]=gr[grid_x][grid_y+1][grid_z+1].z;
	//swap x<->z
	coor1[0]=coor[2];
	coor1[1]=coor[1];
	coor1[2]=coor[0];
	tmp=gr00[2];
	gr00[2]=gr00[0];
	gr00[0]=tmp;
	tmp=gr11[2];
	gr11[2]=gr11[0];
	gr11[0]=tmp;
	if(In_neighbor(coor1,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

     //010
	gr00[0]=gr[grid_x][grid_y+1][grid_z].x;
	gr00[1]=gr[grid_x][grid_y+1][grid_z].y;
	gr00[2]=gr[grid_x][grid_y+1][grid_z].z;
	gr11[0]=gr[grid_x+1][grid_y+1][grid_z+1].x;
	gr11[1]=gr[grid_x+1][grid_y+1][grid_z+1].y;
	gr11[2]=gr[grid_x+1][grid_y+1][grid_z+1].z;
	//swap y<->z
	coor1[0]=coor[0];
	coor1[1]=coor[2];
	coor1[2]=coor[1];
	tmp=gr00[2];
	gr00[2]=gr00[1];
	gr00[1]=tmp;
	tmp=gr11[2];
	gr11[2]=gr11[1];
	gr11[1]=tmp;
	if(In_neighbor(coor1,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

     //-010
	gr00[0]=gr[grid_x][grid_y][grid_z].x;
	gr00[1]=gr[grid_x][grid_y][grid_z].y;
	gr00[2]=gr[grid_x][grid_y][grid_z].z;
	gr11[0]=gr[grid_x+1][grid_y][grid_z+1].x;
	gr11[1]=gr[grid_x+1][grid_y][grid_z+1].y;
	gr11[2]=gr[grid_x+1][grid_y][grid_z+1].z;
	//swap y<->z
	coor1[0]=coor[0];
	coor1[1]=coor[2];
	coor1[2]=coor[1];
	tmp=gr00[2];
	gr00[2]=gr00[1];
	gr00[1]=tmp;
	tmp=gr11[2];
	gr11[2]=gr11[1];
	gr11[1]=tmp;
	if(In_neighbor(coor1,gr00,gr11)){
	  gr[grid_x][grid_y][grid_z].add(i);
	  continue;
	}

      }
     }
    }
  }

  //write to file GRID_PARA
  gp<<grid_num[0]<<" "<<grid_num[1]<<" "<<grid_num[2]<<endl;
  gp<<vertex[0][0]<<" "<<vertex[0][1]<<" "<<vertex[1][0]<<endl;
  gp<<vertex[1][1]<<" "<<vertex[2][0]<<" "<<vertex[2][1]<<endl;

  for(i=0;i<grid_num[0]-1;i++){
    for(int j=0;j<grid_num[1]-1;j++){
      for(int k=0;k<grid_num[2]-1;k++){
	gp<<gr[i][j][k].x<<" "<<gr[i][j][k].y<<" "<<gr[i][j][k].z<<endl;
	gp<<gr[i][j][k].include_num<<endl;
	for(int kk=0;kk<gr[i][j][k].include_num;kk++) 
		gp<<" "<<gr[i][j][k].include[kk];  
	gp<<endl;
      }
    }
  }
 }

// if existed, input grid parameters

  gp>>grid_num[0]>>grid_num[1]>>grid_num[2];
  gp>>vertex[0][0]>>vertex[0][1]>>vertex[1][0];
  gp>>vertex[1][1]>>vertex[2][0]>>vertex[2][1];
  for(int i=0;i<grid_num[0]-1;i++){
    for(int j=0;j<grid_num[1]-1;j++){
      for(int k=0;k<grid_num[2]-1;k++){
	gp>>gr[i][j][k].x>>gr[i][j][k].y>>gr[i][j][k].z;
	gp>>gr[i][j][k].include_num;
	for(int kk=0;kk<gr[i][j][k].include_num;kk++) 
		gp>>gr[i][j][k].include[kk];  
      }
    }
  }

}

/* ----------- identifing bond connections ------------ */

void loop::connect()
{
  for(register int i=np[lp_start-1];i<np[lp_end]-1;i++){
     for(register int j=i+1;j<np[lp_end];j++){
	if((pdb[i]-pdb[j])<bondlength(pdb[i],pdb[j])){
	  pdb[i].link(j);
	  pdb[j].link(i);
	}
     }
  }
}


/* ----------- generation of torsion angle array ------------- */

void loop::torsionGen()
{
  char resType[20][5];
  short int tor[20][MAX_SIDE_TOR][4], torNum[20];
  char atomName[7],resName[5];
  float coor[3];
  register int i,j,k,ii,jj;

// --------- open standard topology file
  ifstream in2("/hosts/tabasco/usr/people/hyzhang/work/lpsa/example/TORTEMP");
  if(!in2)
  {
    cout<<"Cannot open file TORTEMP in directory /hosts/tabasco/usr/people/hyzhang/work/lpsa/example/\n";
    exit(1);
  }

  for(i=0;i<20;i++){
    in2>>resType[i]>>torNum[i];
    in2.ignore(80,'\n');
    for(j=0;j<torNum[i];j++){
	 in2>>tor[i][j][0]>>tor[i][j][1]>>tor[i][j][2]>>tor[i][j][3];
    }
  }

// ---------- determine torsion angle atoms of both main and side chain
  j=0;
  for(i=lp_start-1;i<lp_end;i++){
    pdb[np[i]].getInfo(atomName,resName,coor);
    if(strcmp("PRO ", resName))
    {
    	torsion[j][0]=np[i];
    	torsion[j][1]=np[i]+1;
    	torsion[j][2]=np[i-1]+2;
    	torsion[j][3]=np[i]+2;
    	torsion[j][4]=0;	//main chain
    	j++;
    	torsion[j][0]=np[i]+1;
    	torsion[j][1]=np[i]+2;
    	torsion[j][2]=np[i];
    	torsion[j][3]=np[i+1];
    	torsion[j][4]=0;	//main chain
    	j++;
    }
    else
    {
    	torsion[j][0]=np[i]+4;
    	torsion[j][1]=np[i]+2;
    	torsion[j][2]=np[i];
    	torsion[j][3]=np[i+1];
    	torsion[j][4]=0;	//main chain
    	j++;
    }

    for(jj=0;jj<20;jj++){
	char tmp[4];
	strncpy(tmp,resName,3);
	if(!strcmp(resType[jj],tmp)){
    	  for(ii=0;ii<torNum[jj];ii++){
	     for(k=0;k<4;k++) torsion[j][k]=np[i]+tor[jj][ii][k];
    	     torsion[j][4]=1;	//side chain
    	     j++;
	  }
	}
    }
    torsion_num=j;
    in2.close();	
  }
}

/* ------- find out nonbond pairs to each atom of loop ,
	   determine the total nonbond pairs ----------------- */

void loop::nonbond()
{
  register int i,j,k,ii,jj;

  // ---------  neighboring atoms related with bond, angle and dihedral  
  short int neighbor_num=0;
  int  neighbor[64]; //4*4*4
  short int visited[MAX_ATOM];	// being searched or not
  int bd[4];		// linked atoms

  nbPair_num=0;

  for(i=np[lp_start-1];i<np[lp_end];i++){
    neighbor_num=0;
    for(j=0;j<total_atom_num;j++) visited[j]=0;
    neighbor[0]=i;
    visited[i]=1;
    short int tail=0;
    for(jj=0;jj<2;jj++){  // jj<2 means bond,angle, excluding dihedrals
      for(ii=0;ii<tail+1;ii++){
	short int aa=getlink(pdb[neighbor[ii]], bd);
	for(j=0; j<aa; j++){
	   if(!visited[bd[j]]){
		neighbor[++neighbor_num]=bd[j]; 	//++x not x++
		visited[bd[j]]=1;
	   }
	}
      }
      tail=neighbor_num;
    }

    for(j=np[lp_start-1];j<np[lp_end];j++){	// pairs in loop
	if(!visited[j]&&i<j){
	   nbPair[nbPair_num][0]=i; 
	   nbPair[nbPair_num][1]=j;
	   nbPair_num++;
	}
    }
  }

  // --------  warning informations of array overflow
  if(nbPair_num>MAX_NONBOND_PAIRS){
    cout<<"'nbPair' overflow ! 'MAX_NONBOND_PAIRS' should be larger than "
	<<nbPair_num+1<<endl;
    exit(1);
  }

  // ---- determine the rotating and fixed parts around each torsion bond
  queue	search;			//queue for searching rotating atoms

  for(k=0;k<torsion_num;k++){

    // -------- initialize
    pdb[torsion[k][1]].cut(torsion[k][0]);
    for(j=0;j<total_atom_num;j++) visited[j]=0;
    i=torsion[k][1];	//end atom of bond
    search.init();
    search.store(i);
    visited[i]=1;
    rot_num[k]=fix_num[k]=0;

    // -------- search rotating part, excluding bond end terminal atom
    int seed;
    while((seed=search.retrieve())>=0){
   	short int aa=getlink(pdb[seed], bd);
   	for(j=0; j<aa;j++){
	   if(!visited[bd[j]]){
		rot[k][rot_num[k]++]=bd[j]; 	//x++ not ++x
		visited[bd[j]]=1;
    		search.store(bd[j]);
	   }
	}
    }

    // --------  detach fixed part, excludingg bond terminal atoms
    for(j=np[lp_start-1];j<np[lp_end];j++){
	if(!visited[j]&&(j!=torsion[k][0])&&(j!=i)) fix[k][fix_num[k]++]=j;
    }

    // --------  warning informations of array overflow
    if(rot_num[k]>ROT_NUM||fix_num[k]>ROT_NUM){
	cout<<"'rot or fix' overflow ! 'ROT_NUM' should be larger than "
	<<(MAX(rot_num[k]+1,fix_num[k]+1))<<endl;
	exit(1);
    }

    pdb[torsion[k][1]].link(torsion[k][0]);
  }

/* ------------- output rotating atom ------------------- */
/*
  ofstream ou1("tmp");

  for(k=0;k<torsion_num;k++){
	ou1<<"torsion "<<k+1<<endl;
	ou1<<torsion[k][0]+1<<" "<<torsion[k][1]+1<<" "
	<<torsion[k][2]+1<<" "<<torsion[k][3]+1<<endl;
	ou1<<rot_num[k]<<" "<<fix_num[k]<<endl;
	for(j=0;j<rot_num[k];j++){
	    ou1<<rot[k][j]+1<<" ";
	    if(j%16==15) ou1<<endl;
	}
	ou1<<endl;
	for(j=0;j<fix_num[k];j++){
	    ou1<<fix[k][j]+1<<" ";
	    if(j%16==15) ou1<<endl;
	}
	ou1<<endl<<endl;
  } 
*/
}

/* ------------- determine N postion ----------------------- */

void loop::NPos()
{
  char atomName[7],resName[5];
  float coor[3];

  int j=0;
  for(int i=0;i<total_atom_num;i++){
    pdb[i].getInfo(atomName,resName,coor);
    if(!strcmp(atomName,"  N   ")||!strcmp(atomName,"  C1  ")){
	np[j++]=i;
    }
  }
  for(i=0;i<2;i++) anchor[i]=pdb[np[lp_end]+i];
}

/* ------------- construct initial conformation --------------------- */

void loop::construct()
{
  register int i;
  int lastNum=np[lp_start-2], lpStartNum=np[lp_start-1];
  PDBdata last_CA,last_C,last_O,dummy_N,dummy_CA,lpStart_N,lpStart_CA,lpStart_C;
  PDBdata lk,lk1,lk2,mut_vec;
  float CA_N_TCA;
  float coor[3];
  char atomName[7],resName[5];
  
  last_C=pdb[lastNum+2];
  last_O=pdb[lastNum+3];
  pdb[lastNum].getInfo(atomName,resName,coor);
  if(strcmp("PRO ", resName)) 
    { last_CA=pdb[lastNum+1]; }
  else
    { last_CA=pdb[lastNum+4]; }
  dummy_N=nextAtom(last_C,last_O,last_CA);
  dummy_CA=nextAtom(dummy_N,last_C,last_CA);

  lpStart_N=pdb[lpStartNum];
  lk=lpStart_N<=dummy_N;
  for(i=lpStartNum;i<np[lp_end];i++) pdb[i]=pdb[i]<=lk;
  lpStart_N=pdb[lpStartNum];
  lpStart_N.getInfo(atomName,resName,coor);
  if(strcmp("PRO ", resName)) 
    { lpStart_CA=pdb[lpStartNum+1]; }
  else
    { lpStart_CA=pdb[lpStartNum+4]; }
  lk1=lpStart_CA<=dummy_N;
  lk2=dummy_CA<=dummy_N;
 /* 
	X=nextAtom(A,B,C)

       B
	\
	 \
	  A ----- X
	 /	   \ lk
	/	    v
	C	    X'

	    mut_vec
		^
		|   lk2
	dummy_N	|-------->dummy_CA
	   (N)  \
	 	 \lk1
		  v
		  CA
*/
  mut_vec=lk1*lk2+dummy_N;
  CA_N_TCA=CONV*angle(lpStart_CA,dummy_N,dummy_CA);
  for(i=lpStartNum+1;i<np[lp_end];i++)
  { Rotate(CA_N_TCA, dummy_N, mut_vec, pdb[i], pdb[i]); }

  if(!strcmp("PRO ", resName)) 
  {
    lpStart_CA=pdb[lpStartNum+4];
    lpStart_C=pdb[lpStartNum+2];
    float Fi_of_Pro=torsionAngle(lpStart_N,lpStart_CA,last_C,lpStart_C);
    for(i=lpStartNum+1;i<np[lp_end];i++)
    { Rotate(-Fi_of_Pro-75.0, lpStart_N, lpStart_CA, pdb[i], pdb[i]); }
//    lpStart_C=pdb[lpStartNum+2];
//    lpStart_CA=pdb[lpStartNum+4];
//    float Fi=torsionAngle(lpStart_N,lpStart_CA,last_C,lpStart_C);
//    cout<<Fi;
  }

  int lpEndNum=np[lp_end-1];
  anchor[0]=nextAtom(pdb[lpEndNum+2],pdb[lpEndNum+3],pdb[lpEndNum+1]);
  anchor[1]=nextAtom(anchor[0],pdb[lpEndNum+2],pdb[lpEndNum+1]);
  pdb[lpEndNum].getInfo(atomName,resName,coor);
  if(!strcmp("PRO ", resName))
  {
    anchor[1]=nextAtom(anchor[0],pdb[lpEndNum+2],pdb[lpEndNum+4]);
  }
}

/* ------------- Monte Carlo Annealing ------------------- */

void loop::anneal()
{
  GGL 		random[3]={random0,random1,random2};
  PDBdata 	loop[MAX_LP_LENGTH], // Minimum conformation
	  	anchorRot[2], pdbRot[MAX_ATOM];  // Rotated conformation
  int 		lp_length=np[lp_end]-np[lp_start-1];
  char 		Name[7];
  float 	Etotal[3],Emin[3], coor[3];
  double 	de[3];
  register int 	i,j,k,ii,jj,kk,iii,jjj;
  int 		grid_x,grid_y,grid_z;
  double 	uab[MAX_ATOM][MAX_ATOM]; 
			//energy matrix, upright is before rotating
			//lowerleft is post rotating.
  double 	enviE[MAX_ATOM][2]; 
			//total interatomic interaction to one atom, 
			//[0] is before rotating, [1] is post rotating.

  ofstream ou("OUNAME");
  if(!ou)
  {
    cout<<"Cannot open file OUNAME\n";
    exit(1);
  }

//* ------------ initialize energy

//----- first, total overlap energy

  Etotal[1]=Etotal[2]=0.;
  for(i=0;i<total_atom_num;i++){
    enviE[i][0]=0;
    enviE[i][1]=0;
  }

  //----- intra-action in loop
  for(i=0;i<nbPair_num;i++){
     ii=nbPair[i][0];
     jj=nbPair[i][1];
     iii=MIN(ii,jj);
     jjj=MAX(ii,jj);	
     uab[iii][jjj]=energy(pdb[ii],pdb[jj]);
     Etotal[1]+=uab[iii][jjj];
  }

  //----- inter-action between loop and environment
  for(i=np[lp_start-1];i<np[lp_end];i++){
    pdb[i].getInfo(Name,Name,coor);
    if( coor[0]>vertex[0][0] && coor[0]<vertex[0][1] &&
	coor[1]>vertex[1][0] && coor[1]<vertex[1][1] &&
	coor[2]>vertex[2][0] && coor[2]<vertex[2][1] )
    {
    	grid_x=(coor[0]-vertex[0][0])/GRID_SIZE;
    	grid_y=(coor[1]-vertex[1][0])/GRID_SIZE;
    	grid_z=(coor[2]-vertex[2][0])/GRID_SIZE;
    	for(kk=0;kk<gr[grid_x][grid_y][grid_z].include_num;kk++)
	{
	enviE[i][0]+=energy(pdb[i],pdb[gr[grid_x][grid_y][grid_z].include[kk]]);
	}
    	Etotal[1]+=enviE[i][0];
    }
  }
     
// --- second, anchoring energy
  PDBdata anchor_ref[2];
  char atomName[7],resName[5];

  anchor_ref[0]=pdb[np[lp_end]];
  anchor_ref[1]=pdb[np[lp_end]+1];
  pdb[np[lp_end]].getInfo(atomName,resName,coor);
  if(!strcmp("PRO ", resName)) anchor_ref[1]=pdb[np[lp_end]+4];

  for(i=0;i<2;i++) Etotal[2]+=Aenergy(anchor_ref[i],anchor[i]);
  Etotal[0]=Etotal[1]+Etotal[2];

//* ------- initializing the annealing parameters
  short int iterate=0;
  float t=INIT_TEMP;
  float accept_ratio=0.;
  for(i=0;i<3;i++) Emin[i]=Etotal[i];
  for(i=0;i<total_atom_num;i++) pdbRot[i]=pdb[i];
  for(i=0;i<lp_length;i++) loop[i]=pdb[np[lp_start-1]+i];
  ou<<"CSNEC= "<<csnec<<endl;
  ou<<"RANDOM SEED= "<<random0<<" "<<random1<<" "<<random2<<endl;
  ou<<"INITIAL ENERGY: E_total= "<<Etotal[0]<<", E_softsphere= "<<Etotal[1]<<", E_anchoring= "<<Etotal[2]<<endl;

  // ------------ iteration
  while(iterate<100&&t>.1){
     int count=0;
     for(i=0;i<csnec;i++){
	k=int(torsion_num*random[0]());
        float degree=STEP*(2*random[1]()-1);
	for(j=0;j<3;j++) de[j]=0.0;

	for(j=0;j<rot_num[k];j++){
	    iii=rot[k][j];

	    // ------------ rotate post-bond atoms
	    Rotate(degree, pdb[torsion[k][0]], pdb[torsion[k][1]],
		pdb[iii], pdbRot[iii]);

	    // ------------ calculate energy variation
	    // --- first, overlap energy
	    // intra-loop terms
	    for(kk=0;kk<fix_num[k];kk++){
	 	jjj=fix[k][kk];
		ii=MIN(iii,jjj);
		jj=MAX(iii,jjj);	
		uab[jj][ii]=energy(pdbRot[iii],pdb[jjj]);
		de[1]+=uab[jj][ii]-uab[ii][jj];
	    }
	    // loop-enviroment terms
	    enviE[iii][1]=0;
	    pdbRot[iii].getInfo(Name,Name,coor);
	    if( coor[0]>vertex[0][0] && coor[0]<vertex[0][1] &&
		coor[1]>vertex[1][0] && coor[1]<vertex[1][1] &&
		coor[2]>vertex[2][0] && coor[2]<vertex[2][1] )
	    {
		grid_x=(coor[0]-vertex[0][0])/GRID_SIZE;
		grid_y=(coor[1]-vertex[1][0])/GRID_SIZE;
		grid_z=(coor[2]-vertex[2][0])/GRID_SIZE;
		for(kk=0;kk<gr[grid_x][grid_y][grid_z].include_num;kk++)
		{
		  enviE[iii][1]+=energy(pdbRot[iii],
		  pdb[gr[grid_x][grid_y][grid_z].include[kk]]);
		}
            }
	    de[1]+=enviE[iii][1]-enviE[iii][0];
	}
	// --- second, anchoring energy
	if(!torsion[k][4]){
            // ------------ rotate anchor atoms if it's main chain dihedral
	    for(kk=0;kk<2;kk++){
		Rotate(degree, pdb[torsion[k][0]], pdb[torsion[k][1]],
		anchor[kk], anchorRot[kk]);
		de[2]+=Aenergy(pdb[np[lp_end]+kk],anchorRot[kk]);
	    }
	    de[2]-=Etotal[2];
	}
	de[0]=de[1]+de[2];

	//---------------------------------- if accepted
	if(exp(-de[0]/t)>random[2]()){
	    if(!torsion[k][4]){
	    	for(kk=0;kk<2;kk++) anchor[kk]=anchorRot[kk];
	    }
	    for(j=0;j<rot_num[k];j++){
		iii=rot[k][j];
	        pdb[iii]=pdbRot[iii];
		// intra-action
		for(kk=0;kk<fix_num[k];kk++){
	 	   jjj=fix[k][kk];
		   ii=MIN(iii,jjj);
		   jj=MAX(iii,jjj);	
		   uab[ii][jj]=uab[jj][ii];
		}
		// inter-action
	        enviE[iii][0]=enviE[iii][1];
	    }
	    for(j=0;j<3;j++) Etotal[j]+=de[j];
	    if(Etotal[0]<Emin[0]){ // then reserve minimum conformation
		for(j=0;j<lp_length;j++) loop[j]=pdb[np[lp_start-1]+j];
	        for(j=0;j<3;j++) Emin[j]=Etotal[j];
	    }
	    count++;
	}
     }
     accept_ratio=float(count)/float(csnec);
     iterate++;
     if(iterate==1&&accept_ratio<.8)   //if it is initial iteration
     {
	   ou<<"ACCEPTION RATE="<<accept_ratio;
	   ou<<", THE TEMPERATURE "<<t<<" IS NOT HIGH ENOUGH\n";
	   t=t*2;
	   iterate=0;
     }
     else
     {
	ou<<"ITERATION:"<<iterate<<"  TEMP:"<<t<<"  ACCEPTION RATIO:";
	ou<<accept_ratio<<endl;
	ou<<"	     ENERGY:"<<Etotal[0]<<"  "<<Etotal[1]<<"  ";
	ou<<Etotal[2]<<endl;
	t=t/1.2;
     }
  }
  ou<<"THE MINIMAL ENERGY BEING FOUND IS "<<Emin[0]<<" "<<Emin[1]<<" ";
  ou<<Emin[2]<<endl;
  ou<<"  SAVED\n";
  ou<<total_atom_num<<endl;
  for(j=0;j<lp_length;j++) pdb[np[lp_start-1]+j]=loop[j];
  for(i=0;i<total_atom_num;i++) ou<<pdb[i];

  ou.close();
}

