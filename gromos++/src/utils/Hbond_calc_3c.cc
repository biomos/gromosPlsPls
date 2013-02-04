#include <cassert>
#include <iomanip>
#include <algorithm>
#include "../args/Arguments.h"
#include "../bound/Boundary.h"

#include "Hbond_calc_3c.h"
#include "Hbond.h"

using utils::HB3c_calc;
using utils::HB3c;

HB3c_calc::HB3c_calc(HBPara3c para) :
HB_calc(para.maxdist2, para.minangle), min_angle_sum(para.minanglesum),
max_dihedral(para.maxdihedral), hbpara(para) {
}

void HB3c_calc::setval(gcore::System &sys, args::Arguments &args) {
  HB_calc::setval(sys, args);
  opents3c("Hb3cts.out", "Hb3cnumts.out");
}//end HB3c_calc::setval()

void HB3c_calc::init() {
  HB_calc::init();
  numHb3c = 0;
}//end HB3c_calc::init()

void HB3c_calc::clear() {
  frames--;
  tstime3c.resize(0);
  tsnum3c.resize(0);
  tsnumHB3c.resize(0);
  Hb3c_container::iterator iter = hb3cc.begin(), to = hb3cc.end();
  for (; iter != to; ++iter)
    iter->second.clear();
}//end HB3c_calc::clear()

void HB3c_calc::opents3c(string fi1, string fi2) {
  timeseriesHB3c.open(fi1.c_str());
  timeseriesHB3ctot.open(fi2.c_str());
}//end HB3c_calc::opents3c()

void HB3c_calc::calc_native() {
  frames++;
  numHb3c = 0;
  // loop over hydrogen bonds that we already have only
  //Hb3c_container::const_iterator iter = hb3cc.begin(), to = hb3cc.end();
  Hb3c_container::iterator iter = hb3cc.begin(), to = hb3cc.end();
  for (; iter != to; ++iter) {
    // reset(/*0*/);
    int i = getdon_ind(iter->first);
    int j = getacc_indj(iter->first);
    int k = getacc_indk(iter->first);
    gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
            sys -> box());
    double angles1, dist1;
    gmath::Vec acceptor_j, vec_ij;
    acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j),
            sys->box());
    vec_ij = acceptor_j - donors.pos(i);
    if (neighbour(i, j) && distances(dist1, vec_ij) && angle(i, angles1, bound_i, vec_ij)) {
      double angles2, dist2;
      gmath::Vec acceptor_k, vec_ik;
      acceptor_k = pbc->nearestImage(donors.pos(i), acceptors.pos(k),
              sys->box());
      vec_ik = acceptor_k - donors.pos(i);
      if (neighbour(i, k) && distances(dist2, vec_ik) && angle(i, angles2, bound_i, vec_ik)) {
        double angle_sum, dihedral;
        angle_sum = angles1 + angles2;
        if (anglesum(i, angle_sum, acceptor_j, acceptor_k) &&
                dihedrals(i, dihedral, bound_i, acceptor_j, acceptor_k)) {
          tstime3c.push_back(time);
          tsnum3c.push_back(iter->first);
          numHb3c++;
          hb3cc[iter->first].incr_num();
          hb3cc[iter->first].setdist(dist1, dist2);
          hb3cc[iter->first].setangle(angles1, angles2, angle_sum);
          hb3cc[iter->first].setdihed(dihedral);
        }
      }
    }
  }
  tsnumHB3c.push_back(numHb3c);
}//end HB3c_calc::calc_native()

void HB3c_calc::calc3() {
  frames++, numHb3c = 0;
  // loop over possible hydrogen bonds
  // first A -> B
  int k_lim = acceptors.size();
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < num_A_donors; i++) {
    gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
            sys -> box());
    for (int j = num_A_acceptors; j < k_lim; j++) {
      calc3c(i, j, k_lim, bound_i);
    }
  }
  k_lim = num_A_acceptors;
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = num_A_donors; i < donors.size(); i++) {
    gmath::Vec bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i),
            sys -> box());
    for (int j = 0; j < k_lim; j++) {
      calc3c(i, j, k_lim, bound_i);
    }
  }
  tsnumHB3c.push_back(numHb3c);
  tstimeHB3c.push_back(time);
}//end HB3c_calc::calc3()

void HB3c_calc::calc3c(int i, int j, int k_lim, gmath::Vec &bound_i) {
  double angles1, dist1;
  gmath::Vec acceptor_j, vec_ij;
  acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j),
          sys->box());
  vec_ij = acceptor_j - donors.pos(i);
  if (neighbour(i, j) && distances(dist1, vec_ij) && angle(i, angles1, bound_i, vec_ij)) {
    for (int k = j + 1; k < k_lim; k++) {
      double angles2, dist2;
      gmath::Vec acceptor_k, vec_ik;
      acceptor_k = pbc->nearestImage(donors.pos(i), acceptors.pos(k),
              sys->box());
      vec_ik = acceptor_k - donors.pos(i);
      if (neighbour(i, k) && distances(dist2, vec_ik) && angle(i, angles2, bound_i, vec_ik)) {
        double angle_sum, dihedral;
        angle_sum = angles1 + angles2;
        if (anglesum(i, angle_sum, acceptor_j, acceptor_k) &&
                dihedrals(i, dihedral, bound_i, acceptor_j, acceptor_k)) {
          int index = i * acceptors.size() * acceptors.size() + j * acceptors.size() + k;
#ifdef OMP
#pragma omp critical
#endif
          {
            tstime3c.push_back(time);
            tsnum3c.push_back(index);
            numHb3c++;

            hb3cc[index].incr_num();
            hb3cc[index].setdist(dist1, dist2);
            hb3cc[index].setangle(angles1, angles2, angle_sum);
            hb3cc[index].setdihed(dihedral);
          }
        }
      }
    }
  }
}//end HB3c_calc::calc3c()

void HB3c_calc::printstatistics(HB2c_calc &hb2c) {
  std::cout << endl << endl
          << "# Three-centered hydrogen bonds:" << endl
          << "#"
          << setw(82) << " "
          << setw(29) << "2-CENTER"
          << setw(39) << "3-CENTER" << endl
          << "#"
          << setw(8) << "HB"
          << setw(18) << "Donor"
          << setw(15) << "Acceptor"
          << setw(17) << "D -"
          << setw(14) << "H ..."
          << setw(10) << "A"
          << setw(15) << "DIST"
          << setw(8) << "ANGLE"
          << setw(8) << "OCCUR"
          << setw(8) << "%"
          << setw(15) << "SUM"
          << setw(8) << "DIHED."
          << setw(8) << "OCCUR"
          << setw(8) << "%" << endl;
  int count = 0;
  int i_a1, i_a2, i_d;
  vector<int> totnum, realnum;
  //Hb3c_container::const_iterator it = hb3cc.begin(), to = hb3cc.end();
  Hb3c_container::iterator it = hb3cc.begin(), to = hb3cc.end();
  std::cout.setf(ios::floatfield, ios::fixed);
  for (; it != to; ++it) {
    HB3c &hb3cprint = it->second;
    if (hb3cprint.getnum() > 0) {
      ++count;
      totnum.push_back(it->first);
      realnum.push_back(count);
      hb3cprint.calcmean();
      i_d = getdon_ind(it->first);
      i_a1 = getacc_indj(it->first);
      i_a2 = getacc_indk(it->first);
      std::cout << setw(8) << count;
      if (donors.mol(i_d) < 0) std::cout << setw(8) << " ";
      else std::cout << setw(8) << donors.mol(i_d) + 1;
      std::cout << setw(4) << donors.resnum(i_d) + 1
              << setw(6) << donors.resname(i_d)
              << setw(2) << "-";
      if (acceptors.mol(i_a1) < 0) std::cout << setw(4) << " ";
      else std::cout << setw(4) << acceptors.mol(i_a1) + 1;
      std::cout << setw(4) << acceptors.resnum(i_a1) + 1
              << setw(6) << acceptors.resname(i_a1);
      std::cout << setw(11) << bound.atom(i_d) + 1
              << setw(4) << bound.name(i_d)
              << setw(2) << "-"
              << setw(6) << donors.atom(i_d) + 1
              << setw(6) << donors.name(i_d)
              << setw(2) << "-"
              << setw(6) << acceptors.atom(i_a1) + 1
              << setw(6) << acceptors.name(i_a1);
      std::cout.precision(3);
      std::cout << setw(13) << hb3cprint.getmeandist(0);
      std::cout.precision(3);
      std::cout << setw(8) << hb3cprint.getmeanangle(0);

      // get the occurrence of this HB as two centered HB
      index = acceptors.size() * i_d + i_a1;
      int occur1 = hb2c.getnum(index);

      std::cout.precision(0);
      std::cout << setw(8) << occur1;
      std::cout.setf(ios::floatfield, ios::fixed);
      std::cout.precision(2);
      std::cout << setw(8) << ((occur1 / (double) frames)*100);
      std::cout.precision(3);
      std::cout << setw(15) << hb3cprint.getmeanangle_sum();
      std::cout.precision(3);
      std::cout << setw(8) << hb3cprint.getmeandihedral();
      std::cout.precision(0);
      std::cout << setw(8) << hb3cprint.getnum();
      std::cout.setf(ios::floatfield, ios::fixed);
      std::cout.precision(2);
      std::cout << setw(8) << ((hb3cprint.getnum() / (double) frames)*100)
              << endl;

      // and the second line
      std::cout << setw(27) << " "
              << setw(2) << "\\";
      if (acceptors.mol(i_a2) < 0) std::cout << setw(4) << " ";
      else std::cout << setw(4) << acceptors.mol(i_a2) + 1;
      std::cout << setw(4) << acceptors.resnum(i_a2) + 1
              << setw(6) << acceptors.resname(i_a2)
              << setw(27) << " "
              << setw(4) << "\\"
              << setw(6) << acceptors.atom(i_a2) + 1
              << setw(6) << acceptors.name(i_a2);

      std::cout.precision(3);
      std::cout << setw(13) << hb3cprint.getmeandist(1);
      std::cout.precision(3);
      std::cout << setw(8) << hb3cprint.getmeanangle(1);
      // get the occurrence of this HB as two centered HB
      index = acceptors.size() * i_d + i_a2;
      int occur2 = hb2c.getnum(index);

      std::cout.precision(0);
      std::cout << setw(8) << occur2;
      std::cout.setf(ios::floatfield, ios::fixed);
      std::cout.precision(2);
      std::cout << setw(8) << ((occur2 / (double) frames)*100);
      std::cout << endl;
    }
  }
  std::vector<int>::const_iterator iter;
  std::vector<int> tmp;
  for (unsigned int ia = 0; ia < tstimeHB3c.size(); ia++) {
    timeseriesHB3ctot.precision(9);
    timeseriesHB3ctot << setw(15) << tstimeHB3c[ia];
    timeseriesHB3ctot << setw(10) << tsnumHB3c[ia] << endl;
  }
  for (unsigned int ib = 0; ib < tstime3c.size(); ib++) {
    int find = tsnum3c[ib];
    iter = std::find(totnum.begin(), totnum.end(), find);
    tmp.push_back(realnum[iter - totnum.begin()]);
    if (tstime3c[ib] != tstime3c[ib + 1]) {
      sort(tmp.begin(), tmp.end());
      for (int j = 0; j < tmp.size(); j++) {
        timeseriesHB3c.precision(9);
        timeseriesHB3c << setw(15) << tstime3c[ib] << " ";
        timeseriesHB3c << setw(10) << tmp[j] << endl;
      }
      tmp.clear();
    }
  }
}//end HB3c_calc::printstatistics()

void HB3c::clear() {
  distance.clear();
  distance.resize(2);
  angle.clear();
  angle.resize(2);
  dihed = 0.0;
  angletot = 0.0;
  num = 0;
}//end HB3c::clear()

void HB3c::incr_num() {
  num++;
}//end HB3c::incr_num()

void HB3c::setdist(double &dist1, double &dist2) {
  distance[0] += sqrt(dist1);
  distance[1] += sqrt(dist2);
}//end HB3c::setdist()

void HB3c::setangle(double &angles1, double &angles2, double &angle_sum) {
  angle[0] += angles1;
  angle[1] += angles2;
  angletot += angle_sum;
}//end HB3c::setangle()

void HB3c::setdihed(double dihedral) {
  dihed += dihedral;
}//end HB3c::setdihed()

void HB3c::calcmean() {
  mean_distance.clear();
  mean_angle.clear();
  mean_distance.push_back(distance[0] / num);
  mean_distance.push_back(distance[1] / num);
  mean_angle.push_back(angle[0] / num);
  mean_angle.push_back(angle[1] / num);
  angletot_mean = angletot / num;
  dihed_mean = dihed / num;
}//end HB3c::calcmean()

int HB3c::getnum() {
  return num;
}//end HB3c::getnum()

double HB3c::getmeandist(int i) {
  return mean_distance[i];
}//end HB3c::getmeandist()

double HB3c::getmeanangle(int i) {
  return mean_angle[i];
}//end HB3c::getmeanangle()

double HB3c::getmeanangle_sum() {
  return angletot_mean;
}//end HB3c::getmeanangle_sum()

double HB3c::getmeandihedral() {
  return dihed_mean;
}//end HB3c::getmeandihedral()
