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

#include <lsst/afw/image/Mask.h>
#include <lsst/pex/policy/Policy.h>

namespace lsst {
namespace ip {
namespace diffim {

    /** 
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
     *
     * An abstract class.  Derived classes gain access to private variables
     * through the protected methods.
     */
    template <typename ImageT>
    class SpatialModelBase {
    public: 
        typedef boost::shared_ptr<SpatialModelBase<ImageT> > Ptr;

        /** Empty constructor
         */
        SpatialModelBase();

        /** Destructor
         */
        virtual ~SpatialModelBase() {;};

        /** Set col centroid of Model; range -1 to 1
         *
         * @param colc  Column center
         */
        void setColc(double colc) {_colc = colc;};

        /** Get col centroid of Model; range -1 to 1
         */
        double getColc() {return _colc;};

        /** Set row centroid of Model; range -1 to 1
         *
         * @param rowc  Row center
         */
        void setRowc(double rowc) {_rowc = rowc;};

        /** Get row centroid of Model; range -1 to 1
         */
        double getRowc() {return _rowc;};

        /** Execute the time-consuming process of building the local model
         * 
         * Virtual function that must be overridden by derived class
         */
        virtual bool buildModel() = 0;

        /** Return Sdqa rating
         * 
         * Virtual function that must be overridden by derived class
         */
        virtual double returnSdqaRating(lsst::pex::policy::Policy &policy2) = 0;

        /** Set its build status
         *
         * @param built  Boolean status of build
         */
        void setBuildStatus(bool built) {_isBuilt = built;};

        /** Get its build status
         */
        bool getBuildStatus() {return _isBuilt;};

        /** Get its build status
         */
        bool isBuilt() {return _isBuilt;};

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

        /** Set running ID
         *
         * @param id  integer that represents the id
         */
        void setID(int id) {_id = id;};

        /** Get running ID
         */
        int getID() {return _id;};


    private: 
        int _id;          ///< Running ID
        bool _isBuilt;    ///< Model has been built
        bool _isGood;     ///< Passes local and/or Sdqa requirments
        double _colc;     ///< Effective col position of model in overall image
        double _rowc;     ///< Effective col position of model in overall image


    }; // end of class

}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELBASE_H
    
