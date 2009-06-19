/**
 * @file InIACElementNameMapping.cc
 * implemenation of the element mapping reader
 */

#include "InIACElementNameMapping.h"
#include <sstream>
#include <iostream>
#include <vector>
#include "../gromos/Exception.h"

gio::InIACElementNameMapping::InIACElementNameMapping(std::string file) {
  open(file);
}

gio::InIACElementNameMapping::~InIACElementNameMapping() {
  close();
}

void gio::InIACElementNameMapping::open(std::string file) {
  file_stream.open(file);
}

void gio::InIACElementNameMapping::close() {
  file_stream.close();
}

std::map<int,std::string> gio::InIACElementNameMapping::getData() {
  std::vector<std::string> buffer;
  std::map<int, std::string> gactoele;

  file_stream.getblock(buffer);
  if (buffer[0] != "ELEMENTMAPPING")
    throw gromos::Exception("getElementMapping",
          "library file does not contain a ELEMENTMAPPING block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("InIACElementNameMapping", "No END in ELEMENTMAPPING"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    int gac;
    std::string element;
    std::istringstream ss(buffer[i]);
    ss >> gac >> element;
    if (ss.fail())
      throw gromos::Exception("InIACElementNameMapping", "bad line in ELEMENTMAPPING"
            " block.");
    gactoele[gac - 1] = element;
  }
  return gactoele;
}



