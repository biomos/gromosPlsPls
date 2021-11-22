/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AmberTopology.h
 * Author: bschroed
 *
 * Created on March 7, 2018, 3:19 PM
 */



#ifndef AMBERTOPOLOGY_H
#define AMBERTOPOLOGY_H
//for .h - here
#include <map>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <cstdio>
#include <cmath>
#include <c++/4.8.2/iosfwd>

#include "Ginstream.h"

#include "../gcore/GromosForceField.h"
#include "../gcore/LinearTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/LJException.h"
#include "../gcore/LJType.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MassType.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/Improper.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/DihedralType.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/LinearTopology.h"

#include "../gmath/Physics.h"
#include "../utils/StringOps.h"
#include "../gromos/Exception.h"


using namespace std;
using namespace gcore;

namespace gio{
    
    class AmberTopology  : public Ginstream {
            AmberTopology &operator=(const AmberTopology&);
            
      public:
        GromosForceField d_gff;
        map<string, vector<string>> d_blocks;

          /**
         * Give Constructor the path<string> to amber top file!
         */
        AmberTopology();
        AmberTopology(const AmberTopology&);
        AmberTopology(string s);  // : d_blocks()
        virtual ~AmberTopology();

        /**
         * the init function reads in the whole file into the map of blocks and
         * reads in the topology version
         */
        void init();
        /**
         * parseForceField takes all blocks that end up in the forcefield
         * and stores the information in... d_gff
         */
        void parseFile(LinearTopology &lt, double ljscaling);

        const GromosForceField & forceField() const {
          return d_gff;
        }

        private:
            static string readAmber(string inputfile);

    };
}

#endif /* AMBERTOPOLOGY_H */

