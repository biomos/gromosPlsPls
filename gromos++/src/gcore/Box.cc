// gcore_Box.cc

#include <cassert>
#include "../gromos/Exception.h"
#include "Box.h"


void gcore::Box::update_triclinic()
{
  d_K_L_M = K().cross(L()).dot(M());
  d_cross_K_L_M[0] = L().cross(M())/ -d_K_L_M;
  d_cross_K_L_M[1] = K().cross(M())/  d_K_L_M;
  d_cross_K_L_M[2] = K().cross(L())/ -d_K_L_M;
}

void gcore::Box::setNtb(boxshape_enum b)
{
  d_ntb = b;
  switch(b){
    case vacuum:
    case rectangular:
    case truncoct:
      d_boxformat = box96;
      break;
    case triclinic:
      d_boxformat = triclinicbox;
      break;
    default:
      throw gromos::Exception("Box", "Boxshape not handled");
  }
}
