/**
 * @file Mesh.h
 * a 3D mesh class
 */

#ifndef INCLUDED_GMATH_MESH_H
#define	INCLUDED_GMATH_MESH_H

#include "Vec.h"


namespace gmath
{
  /**
   * A generic 3D mesh class which can be used for spatial distributions
   * and can be written in formats readable by VMD
   *
   * @class Mesh
   * @author N. Schmid
   * @ingroup gmath
   */
  template<typename T>
  class Mesh {
    std::vector<T> data;
    int m_size[3];
    std::string m_title;
    gcore::Box m_box;
    T dummy;
  public:
    /**
     * default constructor
     */
    Mesh() : m_title("GROMOS Mesh") { clear(); }

    /**
     * copy constructor
     * @param mesh the mesh to copy
     */
    Mesh(const Mesh<T> & mesh);
    /**
     * assignment copy
     * @param mesh the mesh to copy
     * @return reference to the mesh
     */
    Mesh<T> & operator=(const Mesh<T> & mesh) {
      if (this != &mesh) {
        new(this) Mesh(mesh);
      }
      return *this;
    }

    Mesh<T> & operator=(const T & val) {
      for (typename std::vector<T>::iterator it = data.begin(), to = data.end();
              it != to; ++it) {
        *it = val;
      }
      return *this;
    }

    /**
     * clear the mesh, free the data, set the size to zero
     */
    void clear() {
      data.clear();
      m_size[0] = m_size[1] = m_size[2] = 0;
    }

    /**
     * accessor to size
     * @return a three membered array containing the dimensions along x, y and z.
     */
    const int* size() const {
      return m_size;
    }

    /**
     * accessor to the title
     * @return the title of the mesh
     */
    const std::string & title() const {
      return m_title;
    }

    /**
     * setter of the title
     * @param title the new title
     */
    void setTitle(const std::string & title) {
      m_title = title;
    }

    /**
     * the the dimensionality of the grid
     * @param x dimension along x
     * @param y dimension along y
     * @param z dimension along z
     */
    void resize(int x, int y, int z) {
      data.resize(x * y * z);
      m_size[0] = x;
      m_size[1] = y;
      m_size[2] = z;
    }

    /**
     * set the box
     * @param box the box
     */
    void setBox(const gcore::Box & box) {
      m_box = box;
    }

    /**
     * access a data element (const version)
     * @param x x index
     * @param y y index
     * @param z z index
     * @return the data element at position (x,y,z)
     */
    inline const T & operator()(int x, int y, int z) const {
      assert(x < m_size[0] && y < m_size[1] && z < m_size[2]);
      return data[z + m_size[2]*(y + m_size[1] * x)];
    }

    /**
     * access a data element
     * @param x x index
     * @param y y index
     * @param z z index
     * @return the data element at position (x,y,z)
     */
    inline T & operator()(int x, int y, int z) {
      return const_cast<T&>((*this)(x, y, z));
    }

    /**
     * access a data element (const version) at a position in space.
     * A box is needed for this function. Set that first.
     * @param v position
     * @return the data element at position v
     */
    inline const T & operator()(const gmath::Vec & v) const {
      int x = int(m_box.K().dot(v) * m_size[0]);
      int y = int(m_box.L().dot(v) * m_size[1]);
      int z = int(m_box.M().dot(v) * m_size[2]);

      if (x < 0 || x >= m_size[0] ||
              y < 0 || y >= m_size[1] ||
              z < 0 || z >= m_size[2])
        return dummy;

      return (*this)(x, y, z);
    }

    /**
     * access a data element at a position in space.
     * A box is needed for this function. Set that first.
     * @param v position
     * @return the data element at position v
     */
    inline T & operator()(const gmath::Vec & v) {
      return const_cast<T&>((*this)(v));
    }

    /**
     * write the mesh to a file in Gaussians CUBE format.
     * A box is required for this.
     * In addition one dummy hydrogen atom at the origin will be added.
     * @param os the stream to write the stuff to.
     */
    void write(std::ostream & os) const;

  };
}

#endif	/* INCLUDED_GMATH_MESH_H */

