#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>

#include "../gcore/Box.h"
#include "Mesh.h"

using namespace gmath;

template<typename T>
gmath::Mesh<T>::Mesh(const Mesh<T> & mesh) {
  resize(mesh.size()[0], mesh.size()[1], mesh.size()[2]);
  copy(mesh.data.begin(), mesh.data.end(), data.begin());
  setBox(mesh.m_box);
}

template<typename T>
void gmath::Mesh<T>::write(std::ostream& os) const {
  os << title() << std::endl
          << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << std::endl
          << "    1    0.000000    0.000000    0.000000" << std::endl;
  os.precision(6);
  os << std::setw(5) << m_size[0]
          << std::setw(12) << m_box.K()[0]
          << std::setw(12) << m_box.K()[1]
          << std::setw(12) << m_box.K()[2] << std::endl;
  os << std::setw(5) << m_size[1] 
          << std::setw(12) << m_box.L()[0]
          << std::setw(12) << m_box.L()[1]
          << std::setw(12) << m_box.L()[2] << std::endl;
  os << std::setw(5) << m_size[2] 
          << std::setw(12) << m_box.M()[0]
          << std::setw(12) << m_box.M()[1]
          << std::setw(12) << m_box.M()[2] << std::endl;
  os << "    1    0.000000    0.000000    0.000000" << std::endl;

  for (int i = 0, x = 0; x < m_size[0]; ++x) {
    for (int y = 0; y < m_size[1]; ++y) {
      for (int z = 0; z < m_size[2]; ++z, ++i) {
        os << std::setw(12) << (*this)(x, y, z) << " ";
        if (i % 6 == 5)
          os << std::endl;

      }
      os << std::endl;
    }
  }
}

template class Mesh<double>;
template class Mesh<float>;
template class Mesh<int>;
template class Mesh<short>;
template class Mesh<char>;
