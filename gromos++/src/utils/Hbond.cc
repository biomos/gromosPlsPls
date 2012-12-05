#include <cassert>
#include "../args/Arguments.h"

#include "Hbond.h"

using utils::HB;
using utils::Pairlist;

HB::HB(gcore::System &sys, args::Arguments &args, HBPara2c hbparas2c, HBPara3c hbparas3c) :
hb2c_calc(hbparas2c), hb3c_calc(hbparas3c) {
  hbpara2c = hbparas2c;
  hbpara3c = hbparas3c;
  do_native = false;
  hb2c_calc.setval(sys, args);
  if (args.count("threecenter") >= 0) {
    do3c = true;
    hb3c_calc.setval(sys, args);
  } else
    do3c = false;
}//end HB::HB()

void HB::init() {
  hb2c_calc.init();
  if (do3c)
    hb3c_calc.init();
}//end HB::init()

void HB::clear() {
  do_native = true;
  hb2c_calc.clear();
  hb3c_calc.clear();
}//end HB::clear()

void HB::calc() {
  if (do_native) {
    hb2c_calc.calc_native();
    if (do3c) {
      hb3c_calc.calc_native();
    }
  } else {
    hb2c_calc.calc2();
    if (do3c) {
      hb3c_calc.calc3();
    }
  }
}//end HB::calc()

void HB::settime(double times) {
  hb2c_calc.settime(times);
  if (do3c)
    hb3c_calc.settime(times);
}//end HB::settime()

void HB::printstatistics() {
  std::cout << "#\n"
          << "# 2-Centered hydrogen bond D-H..A counted if:\n"
          << "#     Distance H..A is at most " << hbpara2c.maxdist << "\n"
          << "#     Angle    D-H..A is at least " << hbpara2c.minangle << "\n";
  std::cout << "#\n#\n";
  if (do3c) {
    std::cout << "#\n"
            << "# 3-Centered hydrogen bond D-H..A1\n"
            << "#                              \\A2 counted if:\n"
            << "#     Distance H..A1 is at most " << hbpara3c.maxdist << "\n"
            << "#     Distance H..A2 is at most " << hbpara3c.maxdist << "\n"
            << "#     Angle    D-H..A1 is at least " << hbpara3c.minangle << "\n"
            << "#     Angle    D-H..A2 is at least " << hbpara3c.minangle << "\n"
            << "#     Sum of angles D-H..A1, D-H..A2, A1..H..A2 is at least "
            << hbpara3c.minanglesum << "\n"
            << "#     Dihedral angle D..A1..A2..H is at most " << hbpara3c.maxdihedral << "\n";
    std::cout << "#\n#\n";
  }
  hb2c_calc.printstatistics();
  if (do3c) {
    hb3c_calc.printstatistics(hb2c_calc);
  }
}//end HB::printstatistics()

utils::HBPara2c HB::mk_hb2c_paras(const vector<double> &hbparas) {
  utils::HBPara2c hbparas2c;
  hbparas2c.maxdist = hbparas[DISTANCE];
  hbparas2c.maxdist2 = hbparas[DISTANCE] * hbparas[DISTANCE];
  hbparas2c.minangle = hbparas[ANGLE];
  return hbparas2c;
}// end HB::mk_hb2c_paras()

utils::HBPara3c HB::mk_hb3c_paras(const vector<double> &hbparas) {
  utils::HBPara3c hbparas3c;
  hbparas3c.maxdist = hbparas[DISTANCE];
  hbparas3c.maxdist2 = hbparas[DISTANCE] * hbparas[DISTANCE];
  hbparas3c.minangle = hbparas[ANGLE];
  hbparas3c.minanglesum = hbparas[SUM];
  hbparas3c.maxdihedral = abs(hbparas[DIHEDRAL]);
  return hbparas3c;
}// end HB::mk_hb3c_paras()