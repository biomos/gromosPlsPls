/**
 * @file AtomicRadii.h
 * AtomicRadii methods
 */

#ifndef INCLUDED_UTILS_ATOMICRADII
#define INCLUDED_UTILS_ATOMICRADII

namespace gcore {
  class System;
  class GromosForceField;
}

namespace utils
{
  /**
   * Compute the atomic radii as the minimal Lennard Jones energy distance of the
   * atoms and the probe particle minus the radius of that particle.
   *
   * @param probe_iac the IAC of the probe
   * @param probe_radius the radius of the probe
   * @param sys the system
   * @param gff the GROMOS force field parameters
   */
  void compute_atomic_radii_vdw(int probe_iac, double probe_radius, gcore::System & sys, const gcore::GromosForceField & gff);
}

#endif

