// utils_Cluster.h


#ifndef INCLUDED_UTILS_CLUSTER
#define INCLUDED_UTILS_CLUSTER

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace utils{
  /**
   * Class Cluster
   * I have no clue whatsoever
   *
   * @class Cluster
   * @author V. Kraeutler
   * @ingroup utils
   * @todo Find out what this actually does and write the documentation
   */
  class Cluster{
    public:
      int center;
      bool is_taken;
      std::vector<int> neighbors;
      Cluster::Cluster(): is_taken(false) {};
  };
}
#endif
