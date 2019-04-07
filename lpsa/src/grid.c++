#include "head.h"

  grid::grid()
  {
    include_num=0;
  }

  void grid::set_coor(float a, float b, float c) 
  {
    x=a;
    y=b;
    z=c;
  }

  void grid::add(int i)
  {
    include[include_num++]=i;
    if(include_num==MAX_INCLUDE){
	cout<<" MAX_INCLUDE overflow ! "<<endl;
	exit(1);
    }
  }
