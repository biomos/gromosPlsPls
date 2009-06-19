/**
 * @file InBFactorOccupancy.cc
 * implemenation of the b-factor and occupancy reader
 */

#include "InBFactorOccupancy.h"
#include <sstream>
#include <iostream>
#include "../gromos/Exception.h"

gio::InBFactorOccupancy::InBFactorOccupancy(std::string file) {
  open(file);
}

gio::InBFactorOccupancy::~InBFactorOccupancy() {
  close();
}

void gio::InBFactorOccupancy::open(std::string file) {
  file_stream.open(file);
}

void gio::InBFactorOccupancy::close() {
  file_stream.close();
}

std::vector<gio::BFactorOccupancyData> gio::InBFactorOccupancy::getData() {
  std::vector<std::string> buffer;
  std::vector<gio::BFactorOccupancyData> bfoccu;
  gio::BFactorOccupancyData sdata;
  // Get Block
  file_stream.getblock(buffer);
  if (buffer[0] != "BFACTOROCCUPANCY")
    throw gromos::Exception("InBFactorOccupancy",
          "library file does not contain a BFACTOROCCUPANCY block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("InBFactorOccupancy", "No END in BFACTOROCCUPANCY"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    // Data-Vars
    std::istringstream ss(buffer[i]);
    ss >> sdata.b_factor >> sdata.occupancy;
    if (ss.fail())
      throw gromos::Exception("InBFactorOccupancy", "bad line in BFACTOROCCUPANCY"
            " block.");
    // Push back to vetor
    bfoccu.push_back(sdata);
  }
  return bfoccu;
}




