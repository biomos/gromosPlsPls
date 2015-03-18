// gcore_Weight.h

#ifndef INCLUDED_GCORE_WEIGHT
#define INCLUDED_GCORE_WEIGHT

#include <queue>

namespace args {
  class Arguments;
}

namespace gcore {
  /**
   * Class Weight
   * Purpose: representing a weight of a configuration
   *
   * Description:
   * Represents the weight of a configuration. By default the weight is 1 but
   * when re-weighting is applied this can change.
   *
   * @class Weight
   * @version $Date: Jul 07 2011
   * @author N. Schmid
   * @ingroup gcore
   * @sa gcore::System
   */
  class Weight {
  public:
    /**
     * Construct the weights and initialise summed weights to zero
     */
    Weight() : has_weights(false), first(true), sum_weights(0.0) {}
    /**
     * copy constructor
     */
    Weight(const Weight & w) : has_weights(w.has_weights), first(w.first),
    weights(w.weights), sum_weights(w.sum_weights) {}
    /**
     * Detects whether weighting is applied and and reads the weights from the
     * argument name specified.
     * @param args the arguments
     * @param argname the name of the argument used
     */
    void read(const args::Arguments & args, const std::string argname = std::string("weights"));
    /**
     * Pops a weight from the queue. The current weight is adjusted and added to
     * the summed weights. 
     * Call pop everytime a new frame is read.
     */
    void pop();
    /**
     * get the current weight
     * @return the current weight
     */
    double getWeight() const;
    /**
     * @sa getWeight
     */
    operator double() const { return getWeight(); }
    /**
     * get the summed weights
     * @return the sum of the weights
     */
    double getSummedWeights() const { return sum_weights; }
    
  private:
    bool has_weights, first;
    std::queue<double> weights;
    double sum_weights;
  };
}

#endif
