/* =================  Main Program of Program LPSA =================== */

#include "head.h"

main()
{
  loop lp;
  lp.input();
  lp.NPos();
  lp.torsionGen();
  lp.construct();
  lp.environment();
  lp.connect();
  lp.nonbond();
  lp.anneal();

  return 0;
}
