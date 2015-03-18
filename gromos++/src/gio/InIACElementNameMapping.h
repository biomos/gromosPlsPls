/**
 * @file InIACElementNameMapping.h
 * a file to read element mapping files
 */

#ifndef INCLUDED_INIACELEMENTNAMEMAPPING_H
#define	INCLUDED_INIACELEMENTNAMEMAPPING_H

#include <string>
#include <map>

#include "Ginstream.h"

namespace gio {

  /**
   * @class InIACElementNameMapping
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief reads a IAC to element mapping file
   * A class to handle and read file that map IACs to element names
   *
   * Format of the IAC to element-name mapping file:
   * @verbatim
TITLE
map IAC to element name for xxxx force-field
END
ELEMENTMAPPING
# IAC ELEMENTNAME
1 O
11 N
12 C
26 Fe
END
@endverbatim
   */
  class InIACElementNameMapping {
  public:
    /**
     * constructor
     */
    InIACElementNameMapping() {}
    /**
     * construct from file name
     * @param file file name
     */
    InIACElementNameMapping(std::string file);
    /**
     * destructor
     */
    ~InIACElementNameMapping();
    /**
     * open a file
     * @param file file name
     */
    void open(std::string file);
    /**
     * close the file
     */
    void close();
    /**
     * get the mapping data
     * @return a map containing the IAC to element mapping
     */
    std::map<int, std::string> getData();
  protected:
    Ginstream file_stream;
  };
}


#endif	/* INCLUDED_INIACELEMENTNAMEMAPPING_H */

