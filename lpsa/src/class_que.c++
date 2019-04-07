/* =============	Function of class queue 	============ */
#include "head.h"

queue::queue()
{
  spos=rpos=0;
}

void queue::init()
{
  spos=rpos=0;
}

void queue::store(int i)
{
  if(spos==QUEUE_SIZE){
    cout<<"queue space overflow ! QUEUE_SIZE at least > "<<spos<< endl;
    exit(1);
  }
  que[spos]=i;
  spos++;
}

int queue::retrieve()
{
  if(rpos==spos) return -1;
  rpos ++;
  return que[rpos-1];
}
