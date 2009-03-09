// gcore_Box.cc

#include <cassert>
#include "../gromos/Exception.h"
#include "Box.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"

void gcore::Box::update_triclinic() {
  d_K_L_M = K().cross(L()).dot(M());
  if (d_K_L_M == 0.0) return; // weird box, vacuum or similar
  d_cross_K_L_M[0] = L().cross(M()) / -d_K_L_M;
  d_cross_K_L_M[1] = K().cross(M()) / d_K_L_M;
  d_cross_K_L_M[2] = K().cross(L()) / -d_K_L_M;
}

void gcore::Box::setNtb(boxshape_enum b) {
  d_ntb = b;
}

gcore::Box::Box(gcore::Box::boxshape_enum bound, double a, double b, double c, double alpha, double beta, double gamma,
        double phi, double theta, double psi) {
  d_boxformat = gcore::Box::genbox;
  d_dim.resize(3);
  d_cross_K_L_M.resize(3);
  d_ntb = bound;

  if (d_ntb == gcore::Box::vacuum) {

    K() = gmath::Vec(0.0, 0.0, 0.0);
    L() = gmath::Vec(0.0, 0.0, 0.0);
    M() = gmath::Vec(0.0, 0.0, 0.0);
    update_triclinic();

    return;
  }

  if (d_ntb == gcore::Box::rectangular || d_ntb == gcore::Box::truncoct) {
    if (alpha != 90.0 || beta != 90.0 || gamma != 90.0)
      throw gromos::Exception("GENBOX", "For rectangular and truncated octahedral boxes, alpha, beta"
            " and gamma should be 90 degrees");
    if (phi != 0.0 || theta != 0.0 || psi != 0.0)
      throw gromos::Exception("GENBOX","For rectangular and truncated octahedral boxes, phi, theta"
            " and phi should be 0 degrees");
    K() = gmath::Vec(a, 0.0, 0.0);
    L() = gmath::Vec(0.0, b, 0.0);
    M() = gmath::Vec(0.0, 0.0, c);
    update_triclinic();
    return;
  }
  alpha *= M_PI / 180.0;
  beta *= M_PI / 180.0;
  gamma *= M_PI / 180.0;
  phi *= -M_PI / 180.0;
  theta *= -M_PI / 180.0;
  psi *= -M_PI / 180.0;

  const double cosphi = cos(phi);
  const double sinphi = sin(phi);
  const double costheta = cos(theta);
  const double sintheta = sin(theta);
  const double cospsi = cos(psi);
  const double sinpsi = sin(psi);

  K() = gmath::Vec(1, 0, 0);
  L() = gmath::Vec(cos(gamma), sin(gamma), 0);
  const double cosbeta = cos(beta);
  const double cosalpha = cos(alpha);

  M() = gmath::Vec(cosbeta, 1,
          sqrt((cosbeta * L()[0] + L()[1]) *
          (cosbeta * L()[0] + L()[1]) / (cosalpha * cosalpha)
          - cosbeta * cosbeta - 1));
  M() = M().normalize();
  K() *= a;
  L() *= b;
  M() *= c;

  gmath::Vec x(cospsi * cosphi - costheta * sinphi*sinpsi,
          -sinpsi * cosphi - costheta * sinphi*cospsi,
          sintheta * sinphi);
  gmath::Vec y(cospsi * sinphi + costheta * cosphi*sinpsi,
          -sinpsi * sinphi + costheta * cosphi*cospsi,
          -sintheta * cosphi);
  gmath::Vec z(sinpsi*sintheta,
          cospsi*sintheta,
          costheta);
  gmath::Matrix mat(x, y, z);
  K() = mat * K();
  L() = mat * L();
  M() = mat * M();
  update_triclinic();
}

