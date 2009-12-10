/* 
 * File:   RdcFuncs.h
 * Author: jallison
 *
 * Created on November 24, 2009, 4:09 PM
 */

#ifndef INCLUDED_UTILS_RDCFUNCS
#define INCLUDED_UTILS_RDCFUNCS
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif
#ifndef INCLUDED_ARGS_ARGUMENTS
#include "../args/Arguments.h"
#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace gcore {
    class System;
}
namespace args {
    class Arguments;
}

// define class for storing RDC data and the corresponding atoms and parameters
namespace utils {

    /**
     * Class RdcFuncs
     * Defines a data structure for storing RDCs and related information
     * and some functions for SVD-fitting to and back-calculating RDCs.
     *
     * @class RdcFuncs
     * @author J. Allison
     */

    /**
     * Class RDCData
     * A class to store RDC data read in from a RDC specification file
     *
     * @class RDCData
     * @author J. Allison
     * @ingroup utils
     */
    class RDCData {
    public:

        // struct for storing RDC data
        struct rdcparam {
            // atoms defining inter-nuclear vector
            int mol;
            int i;
            int j;
            // parameters
            double w; // weight factor
            double exp; // experimental rdc
            double gi; // gyromagnetic ratio of atom i
            double gj; // gyromagnetic ratio of atom j
            double dmax; // maximum possible rdc

            rdcparam() : mol(0), i(0), j(0) {
            }

            rdcparam(const rdcparam & rdcp) : mol(rdcp.mol), i(rdcp.i), j(rdcp.j),
            w(rdcp.w), exp(rdcp.exp), gi(rdcp.gi), gj(rdcp.gj), dmax(rdcp.dmax) {
            }

            rdcparam & operator=(const rdcparam & rdcp) {
                mol = rdcp.mol;
                i = rdcp.i;
                j = rdcp.j;
                w = rdcp.w;
                exp = rdcp.exp;
                gi = rdcp.gi;
                gj = rdcp.gj;
                dmax = rdcp.dmax;
                return *this;
            }
        };

        /**
         * Constructor
         */
        RDCData() {
        }

        /**
         *  RDCData copy constructor
         */
        RDCData(const RDCData &rdcdata) {
            m_data = rdcdata.m_data;
        }

        /**
         *  RDCData deconstructor
         */
        ~RDCData() {
        }

        /**
         * const accessor to rdcparam data
         */
        const std::vector<rdcparam> & data() const {
            return m_data;
        }

        /**
         * accessor to rdcparam data
         */
        std::vector<rdcparam> & data() {
            return m_data;
        }
        
    private:
        std::vector<rdcparam> m_data;

    };

    /**
     * Class RDCWeights
     * A class to store weights for individual frames of a trajectory
     *
     * @class RDCWeights
     * @author J. Allison
     * @ingroup utils
     */
    class RDCWeights {
    public:


        // struct for storing weights (for frames of a trajectory)

        struct weights {
            int frame;
            double weight;

            weights() : frame(0), weight(0) {
            }

            weights(const weights & w) : frame(w.frame), weight(w.weight) {
            }

            weights & operator=(const weights & w) {
                frame = w.frame;
                weight = w.weight;
                return *this;
            }

        };

        /**
         * Constructor
         */
        RDCWeights() {
        }

        /**
         *  RDCWeights copy constructor
         */
        RDCWeights(const RDCWeights &weights) {
            m_weights = weights.m_weights;
        }

        /**
         *  RDCWeights deconstructor
         */
        ~RDCWeights() {
        }

        /**
         * const accessor to weight data
         */
        const std::vector<weights> & wdata() const {
            return m_weights;
        }

        /**
         * accessor to weights data
         */
        std::vector<weights> & wdata() {
            return m_weights;
        }

    private:

        std::vector<weights> m_weights;

    };

    /**
     * Class RdcFuncs
     * A class of functions for fitting to and calculating RDCs
     *
     * @class RdcFuncs
     * @author J. Allison
     * @ingroup utils
     */
    class RdcFuncs {
        /**
         * copy constructor
         * not implemented
         */
        RdcFuncs(const RdcFuncs&);

        /**
         * operator =
         * not implemented
         */
        RdcFuncs & operator=(const RdcFuncs&);

        public:

        /**
         * RdcFuncs Constructor
         */
        RdcFuncs(gcore::System &sys, args::Arguments &args);

        /**
         * RdcFuncs Destructor
         */
        ~RdcFuncs();

        // Methods

        // function to read in rdc data
        void read_rdc(std::vector<std::string> buffer, const gcore::System &sys,
                std::vector<utils::RDCData::rdcparam> &rdcp);

        // function to compute dmax from first set of coordinates
        void calc_dmax(const gcore::System &sys,
                std::vector<utils::RDCData::rdcparam> &rdcp, int nrdc);

        // function to read weights for individual frames from file
        void read_weights(std::vector<std::string> buffer, std::vector<RDCWeights::weights> &weight_data);

        // compute the coefficients of the matrix describing bond vector fluctuations
        void calc_coef(const gcore::System &sys, std::vector<utils::RDCData::rdcparam> &fit_data,
                gsl_matrix *coef_mat, int nrdc, double w);

        // fill a vector with the normalised RDCs (scaled by Dmax)
        void fill_rdcvec(const std::vector<utils::RDCData::rdcparam> &R, gsl_vector *v,
        int nrdc);

        // calculate Q value (goodness of fit)
        double calc_Q(gsl_vector *calc, gsl_vector *expt);

        // function to calculate unnormalised RDCs
        void unnorm_rdcvec(const std::vector<utils::RDCData::rdcparam> &R,
                gsl_vector *v, gsl_vector *w);

        /**
         * @struct Exception
         * Throws an exception if something is wrong
         */
        struct Exception : public gromos::Exception {

            /**
             * @exception If called says RdcFuncs, followed by the argument
             * @param what The string that is thrown
             */
            Exception(const std::string & what) :
            gromos::Exception("RdcFuncs", what) {
            }
        };

    }; // end class RdcFuncs

}

#endif	/* _RDCFUNCS_H */

