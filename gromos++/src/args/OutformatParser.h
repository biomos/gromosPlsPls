#ifndef INCLUDED_ARGS_OUTFORMATPARSER
#define	INCLUDED_ARGS_OUTFORMATPARSER

#include <string>

namespace gio {
  class OutCoordinates;
}

namespace args {
  class Arguments;
  /**
   * @class OutformatParser
   * @ingroup args
   * @author N. Schmid
   * @date 04.11.2010
   *
   * Used to parse the coordinate output format argument (usually named \@outformat).
   * The following formats are supported:
   *
   * <table>
   * <tr><th>Argument</th><th>Description</th></tr>
   * <tr><td>cnf</td><td>Configuration format containing the POSITION block.</td></tr>
   * <tr><td>trc</td><td>Trajectory format containing the POSITIONRED block.</td></tr>
   * <tr><td>por</td><td>Position restraints specification format.</td></tr>
   * <tr><td>pdb [&lt;factor to convert length unit to Angstrom, 10.0&gt; and/or &lt;"renumber" keyword to start numbering at 1 at each molecule&gt;] </td><td>Protein Data Bank (PDB) format.</td></tr>
   * <tr><td>vmdam [&lt;factor to convert length unit to Angstrom, 10.0&gt;]</td><td> VMD's Amber Coordinates format.</td></tr>
   * <\table>
   */
  class OutformatParser {
  private:
    OutformatParser(const OutformatParser & op);
    OutformatParser & operator=(const OutformatParser & op);
  public:
    /**
     * parse the arguments and reaturn the OutCoordinates for the correct format
     * @param args the arguments
     * @param ext the extension as a string (e.g. '.cnf')
     * @param argname the argument name used
     * @return pointer to the OutCoordinates
     */
    static gio::OutCoordinates * parse(Arguments & args,
            std::string & ext,
            std::string argname = "outformat");
  };
}

#endif	/* INCLUDED_ARGS_OUTFORMATPARSER */

