/**
 * @file InCIF.cc
 * implemenation of the CIF reader
 */

#include "InCIF.h"
#include "StringTokenizer.h"
#include <sstream>
#include <iostream>

gio::InCIF::InCIF(std::string file) {
  open(file);
}

gio::InCIF::~InCIF() {
  if (isOpen()) close();
}

void gio::InCIF::open(std::string file) {
  file_stream.open(file.c_str());
  if (!file_stream.is_open())
    throw gromos::Exception("InCIF", "Cannot open CIF file " + file + ".");
}

void gio::InCIF::close() {
  if (isOpen())
    file_stream.close();
}

bool gio::InCIF::isOpen() {
  return file_stream.is_open();
}

std::vector<gio::CIFData> gio::InCIF::getData() {
  std::vector<gio::CIFData> cifdata;
  gio::CIFData tcif;
  std::string line;
  int refln_nr = 0, hnr = 0, knr = 0, lnr = 0, Fnr = 0, Fsnr = 0;
  while (true) {
    if (file_stream.eof()) {
      break;
    }
    std::getline(file_stream, line);
    if (line.substr(0, 6) == "_refln") {
      refln_nr++;
      // check index of desired data
      if (line == "_refln.index_h") {
        hnr = refln_nr;
      } else if (line == "_refln.index_k") {
        knr = refln_nr;
      } else if (line == "_refln.index_l") {
        lnr = refln_nr;
      } else if (line == "_refln.F_meas_au") {
        Fnr = refln_nr;
      } else if (line == "_refln.F_meas_sigma_au") {
        Fsnr = refln_nr;
      }
    } else {
      // Check if before or after _refln-Block
      if (refln_nr == 0) {
        // do nothing / del lines
      } else {
        // check if data or not
        gio::StringTokenizer tok(line);
        std::vector<std::string> split = tok.tokenize();
        if (split.size() < 4) {
          // no valid line...skip
        } else {
          std::istringstream is1(split[hnr - 1]); is1 >> tcif.h;
          std::istringstream is2(split[knr - 1]); is2 >> tcif.k;
          std::istringstream is3(split[lnr - 1]); is3 >> tcif.l;
          std::istringstream is4(split[Fnr - 1]); is4 >> tcif.f_obs;
          if (Fsnr != 0) {
            std::istringstream is(split[Fsnr - 1]); is >> tcif.stddev_f_obs;
          } else {
            tcif.stddev_f_obs = 0.0;
          }
          if (is1.fail() || is2.fail() || is3.fail() || is4.fail()) 
            std::cerr << "Skipping bad line." << std::endl;
          else
            cifdata.push_back(tcif);
        }
      }
    }
  }
  return cifdata;
}


