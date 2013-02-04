#include <cassert>
#include <iomanip>
#include <algorithm>
#include "../args/Arguments.h"
#include "../bound/Boundary.h"

#include "Hbond_calc_2c.h"
#include "Hbond.h"

using utils::HB2c_calc;
using utils::HB2c;

HB2c_calc::HB2c_calc(HBPara2c para) :
HB_calc(para.maxdist2, para.minangle), hbpara(para) {

}

void HB2c_calc::setval(gcore::System& sys, args::Arguments& args) {
  HB_calc::setval(sys, args);
  //open timeseries file
  opents("Hbts.out", "Hbnumts.out");
}//end HB2c_calc::setval()

void HB2c_calc::init() {
  HB_calc::init();
  numHb = 0;
  maxindex = acceptors.size() * acceptors.size() + acceptors.size();
}//end HB2c_calc::init()

void HB2c_calc::clear() {
  frames--;
  tstime.resize(0);
  tsnum.resize(0);
  tsnumHB.resize(0);
  hb2cc.clear();
}//end HB2c_calc::clear()

void HB2c_calc::opents(string fi1, string fi2) {
  timeseriesHB.open(fi1.c_str());
  timeseriesHBtot.open(fi2.c_str());
}//end HB2c_calc::opents()

void HB2c_calc::calc_native() {
  frames++;
  numHb = 0;
  // loop over hydrogen bonds that we already have only
  //Hb2c_container::const_iterator iter = hb2cc.begin(), to = hb2cc.end();
  for (int k = 0; k < maxindex; k++) {
    int num = hb2cc[k].getnum();
    if (num > 0) {
      int i = k / acceptors.size();
      int j = k % acceptors.size();
      gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
              sys -> box());
      calc2c(i, j, bound_i);
    }
  }
  tsnumHB.push_back(numHb);
  if (numHb == 0)
    tstime.push_back(time);

}//end HB2c_calc::calc_native()

void HB2c_calc::calc2() {
  frames++;
  numHb = 0;
#ifdef OMP
#pragma omp parallel for
#endif
  // loop over possible hydrogen bonds
  // first A -> B
  for (int i = 0; i < num_A_donors; i++) {
    gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
            sys -> box());
    for (int j = num_A_acceptors; j < acceptors.size(); j++)
      calc2c(i, j, bound_i);
  }
#ifdef OMP
#pragma omp parallel for
#endif
  // then A -> B
  for (int i = num_A_donors; i < donors.size(); i++) {
    gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
            sys -> box());
    for (int j = 0; j < num_A_acceptors; j++)
      calc2c(i, j, bound_i);
  }
  tsnumHB.push_back(numHb);
  tstimeHB.push_back(time);
}//end HB2c_calc::calc2()

void HB2c_calc::calc2c(int i, int j, gmath::Vec &bound_i) {
  double dist, angles;
  gmath::Vec acceptor_j, vec_ij;
  acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j),
          sys->box());
  vec_ij = acceptor_j - donors.pos(i);
  if (neighbour(i, j) && distances(dist, vec_ij) && angle(i, angles, bound_i, vec_ij)) {
    int index = i * acceptors.size() + j;
#ifdef OMP
#pragma omp critical
#endif
    {
      numHb++;
      tstime.push_back(time);
      tsnum.push_back(index);
      hb2cc[index].incr_num();
      hb2cc[index].setdist(dist);
      hb2cc[index].setangle(angles);
    }
  }
}//end HB2c_calc::calc2c()

void HB2c_calc::printstatistics() {
  std::cout << "# Statistics of the run:" << endl;
  std::cout << endl << "# Two-centered hydrogen bonds:\n";
  std::cout << "#"
          << setw(8) << "HB"
          << setw(18) << "Donor"
          << setw(15) << "Acceptor"
          << setw(19) << "D -"
          << setw(14) << "H ..."
          << setw(10) << "A"
          << setw(15) << "DIST"
          << setw(8) << "ANGLE"
          << setw(8) << "OCCUR"
          << setw(8) << "%" << endl;
  int count = 0;
  int i_d, i_a;
  vector<int> totnum, realnum;
  std::cout.setf(ios::floatfield, ios::fixed);
  Hb2c_container::iterator it = hb2cc.begin(), to = hb2cc.end();
  for (; it != to; it++) {
    HB2c &hb2cprint = it->second;
    int occur = hb2cprint.getnum();
    if (occur > 0) {
      ++count;
      totnum.push_back(it->first);
      realnum.push_back(count);
      hb2cprint.calcmean();
      i_d = getdon_ind(it->first);
      i_a = getacc_ind(it->first);
      std::cout << setw(8) << count;
      if (donors.mol(i_d) < 0)
        std::cout << setw(8) << " ";
      else
        std::cout << setw(8) << donors.mol(i_d) + 1;
      std::cout << setw(4) << donors.resnum(i_d) + 1
              << setw(6) << donors.resname(i_d)
              << setw(2) << "-";
      if (acceptors.mol(i_a) < 0)
        std::cout << setw(4) << " ";
      else
        std::cout << setw(4) << acceptors.mol(i_a) + 1;
      std::cout << setw(4) << acceptors.resnum(i_a) + 1
              << setw(6) << acceptors.resname(i_a)
              << setw(11) << bound.atom(i_d) + 1
              << setw(6) << bound.name(i_d)
              << setw(2) << "-"
              << setw(6) << donors.atom(i_d) + 1
              << setw(6) << donors.name(i_d)
              << setw(2) << "-"
              << setw(6) << acceptors.atom(i_a) + 1
              << setw(6) << acceptors.name(i_a);
      std::cout.precision(3);
      std::cout << setw(13) << hb2cprint.getmeandist();
      std::cout.precision(3);
      std::cout << setw(8) << hb2cprint.getmeanangle();
      std::cout.precision(0);
      std::cout << setw(8) << occur;
      std::cout.precision(2);
      std::cout << setw(8) << ((occur / (double) frames)*100)
              << endl;
    }
  }
  std::vector<int>::const_iterator iter;
  std::vector<int> tmp;
  for (unsigned int ia = 0; ia < tstimeHB.size(); ia++) {
    timeseriesHBtot.precision(9);
    timeseriesHBtot << setw(15) << tstimeHB[ia];
    timeseriesHBtot << setw(10) << tsnumHB[ia] << endl;
  }
  for (unsigned int ib = 0; ib < tstime.size(); ib++) {
    int find = tsnum[ib];
    iter = std::find(totnum.begin(), totnum.end(), find);
    tmp.push_back(realnum[iter - totnum.begin()]);
    if (tstime[ib] != tstime[ib + 1]) {
      sort(tmp.begin(), tmp.end());
      for (int j = 0; j < tmp.size(); j++) {
        timeseriesHB.precision(9);
        timeseriesHB << setw(15) << tstime[ib] << " ";
        timeseriesHB << setw(10) << tmp[j] << endl;
      }
      tmp.clear();
    }
  }
}//end HB2c_calc::printstatistics()

int HB2c_calc::getnum(int index) {
  return hb2cc[index].getnum();
}//end HB2c_calc::getnum()

void HB2c::clear() {
  dist = 0;
  angle = 0;
  num = 0;
}//end HB2c::clear()

void HB2c::incr_num() {
  num++;
}//end HB2c::incr_num()

void HB2c::setdist(double &distance) {
  dist += sqrt(distance);
}//end HB2c::setdist()

void HB2c::setangle(double &ang) {
  angle += ang;
}//end HB2c::setangle()

void HB2c::calcmean() {
  mean_dist = dist / num;
  mean_angle = angle / num;
}//end HB2c::calcmean()

int HB2c::getnum() {
  return num;
}//end HB2c::getnum()

double HB2c::getmeandist() {
  return mean_dist;
}//end HB2c::getmeandist()

double HB2c::getmeanangle() {
  return mean_angle;
}//end HB2c::getmeanangle()
