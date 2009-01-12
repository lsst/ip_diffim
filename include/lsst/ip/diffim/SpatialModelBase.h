// -*- lsst-c++ -*-
/**
 * @file SpatialModelBase.h
 *
 * @brief Base class for members of SpatialModelCell _modelPtrList
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup afw
 */

#ifndef LSST_IP_DIFFIM_SPATIALMODELBASE_H
#define LSST_IP_DIFFIM_SPATIALMODELBASE_H

#include <boost/shared_ptr.hpp>

namespace lsst {
namespace ip {
namespace diffim {

    /** @class
     * 
     * @brief Base class for members of SpatialModelCell _modelPtrList
     * 
     * @todo Include Sdqa classes as class members; make notion of being "good"
     * or usable a function of Sdqa classes.
     *
     * An instance of SpatialModelCell will contain a list of instances of
     * SpatialModelBase in _modelPtrList.  The methods that are needed by
     * SpatialModelCell are defined here and should be overriden in each class
     * derived from SpatialModelBase (e.g. SpatialModelKernel).
     *
     * Derived members of this class need to know how to build themselves from
     * their own member variables, implmented in the buildModel() method.
     */
    template <typename ImageT, typename MaskT>
    class SpatialModelBase {
    public: 
        typedef boost::shared_ptr<SpatialModelBase<ImageT, MaskT> > Ptr;

        /** Empty constructor
         */
        SpatialModelBase(){;};

        /** Destructor
         */
        virtual ~SpatialModelBase() {;};

        /** Set col centroid of Model; range -1 to 1
         *
         * @param colc  Column center
         */
        void setColcNorm(double colc) {_colcNorm = colc;};

        /** Get col centroid of Model; range -1 to 1
         */
        double getColcNorm() {return _colcNorm;};

        /** Set row centroid of Model; range -1 to 1
         *
         * @param rowc  Row center
         */
        void setRowcNorm(double rowc) {_rowcNorm = rowc;};

        /** Get row centroid of Model; range -1 to 1
         */
        double getRowcNorm() {return _rowcNorm;};

        /** Execute the time-consuming process of building the local model
         */
        bool buildModel();

        /** Return Sdqa rating
         */
        double returnSdqaRating();

        /** Get its build status
         */
        bool getBuildStatus() {return _isBuilt;};

        /** Set its Sdqa status
         * 
         * @param status  Boolean status of model
         */
        void setSdqaStatus(bool status) {_isGood = status;};

        /** Get its Sdqa status
         */
        bool getSdqaStatus() {return _isGood;};

        /** Get its Sdqa status
         */
        bool isGood() {return _isGood;};

    private: 
        /** Set its build status
         *
         * @param status  Boolean status of build
         */
        void _setBuildStatus(bool built) {_isBuilt = built;};

        double _colcNorm; ///< Effective col position of model in overall image
        double _rowcNorm; ///< Effective col position of model in overall image

        bool _isBuilt;    ///< Model has been built
        bool _isGood;     ///< Passes local and/or Sdqa requirments

    }; // end of class

}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELBASE_H
    
