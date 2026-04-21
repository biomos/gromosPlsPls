/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// gio_OutG96S.h

#ifndef INCLUDED_GIO_OUTG96S
#define INCLUDED_GIO_OUTG96S

#include <string>

#include "OutCoordinates.h"

namespace gcore{
  class System;
}

namespace utils{
  class AtomSpecifier;
}

namespace gio{
  class OutG96S_i;
  /**
   * Class OutG96S
   * is of type OutCoordinates and defines how a single coordinate file 
   * is printed out (POSITION block)
   *
   * @class OutG96S
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @ingroup gio
   */
  class OutG96S: public OutCoordinates{
    OutG96S_i *d_this;
    bool posres;
    // prevent copying and assignment
    OutG96S(const OutG96S &);
    OutG96S &operator=(const OutG96S&);
  public:
    OutG96S(bool posres = false);
    OutG96S(std::ostream &os, bool posres = false);
    ~OutG96S();
    void select(const std::string &thing) override;
    void open(std::ostream &os) override;
    void close() override;
    void writeTitle(const std::string &title) override;
    void writeTimeFrame(const gcore::System &sys, int timeframe) override;
    void writeTimeFrame(const utils::AtomSpecifier & atoms, int timeframe) override;
    void writeTimestep(const int step, const double time) override;
    OutG96S &operator<<(const gcore::System & sys) override;
    OutG96S &operator<<(const utils::AtomSpecifier & atoms) override;
  };
}

#endif
