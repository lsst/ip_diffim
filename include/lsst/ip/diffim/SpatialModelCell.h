// -*- lsst-c++ -*-
/**
 * @file SpatialModelCell.h
 *
 * @brief Class to ensure constraints for spatial modeling
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_SPATIALMODELCELL_H
#define LSST_IP_DIFFIM_SPATIALMODELCELL_H

#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image/Mask.h>
#include <lsst/detection/Footprint.h>

#include <lsst/ip/diffim/SpatialModelBase.h>
#include <lsst/ip/diffim/SpatialModelKernel.h>

namespace lsst {
namespace ip {
namespace diffim {

    /** 
     * 
     * @brief Class to ensure constraints for spatial modeling
     * 
     * A given MaskedImage will be divided up into cells, with each cell
     * represented by an instance of this class.  Each cell itself contains a
     * list of instances of classes derived from SpatialModelBase.  One member
     * from each cell will be fit for a spatial model.  In case of a poor fit,
     * the next SpatialModelBase instance in the list will be fit for.  If all
     * instances in a list are rejected from the spatial model, the best one
     * will be used.
     */
    template <typename ImageT, typename MaskT = lsst::afw::image::maskPixelType>
    class SpatialModelCell {
        
    public:
        typedef boost::shared_ptr<SpatialModelCell<ImageT, MaskT> > Ptr;
        typedef std::vector<typename SpatialModelCell<ImageT, MaskT>::Ptr> SpatialModelCellList;

        typedef typename SpatialModelBase<ImageT, MaskT>::Ptr SpatialModel;
        typedef std::vector<typename SpatialModelBase<ImageT, MaskT>::Ptr> ModelPtrList;
        typedef std::vector<lsst::detection::Footprint::PtrType> FpPtrList;
        
        /** Constructor
         *
         * @param label  string representing "name" of cell
         * @param fpPtrList  vector of pointers to footprints within the cell
         * @param modelPtrList  vector of pointers to models of the function you are fitting for
         */
        SpatialModelCell(std::string label,
                         FpPtrList fpPtrList, 
                         ModelPtrList modelPtrList);
        
        /** Constructor
         *
         * @param label  string representing "name" of cell
         * @param colC  effective location of column center of cell within overall MaskedImage
         * @param rowC  effective location of row center of cell within overall MaskedImage
         * @param fpPtrList  vector of pointers to footprints within the cell
         * @param modelPtrList  vector of pointers to models of the function you are fitting for
         */
        SpatialModelCell(std::string label, int colC, int rowC, 
                         FpPtrList fpPtrList,
                         ModelPtrList modelPtrList);

        SpatialModelCell(ModelPtrList modelPtrList);
        SpatialModelCell(std::vector<typename SpatialModelKernel<ImageT, MaskT>::Ptr> modelPtrList);

        /** Destructor
         */
        virtual ~SpatialModelCell() {;};

        /** Get current footprint
         */
        lsst::detection::Footprint::PtrType getCurrentFootprint() {return _fpPtrList[_currentID];};
        
        /** Get current model
         */
        SpatialModel getCurrentModel() {return _modelPtrList[_currentID];};

        /** Get footprint in list
         * 
         * @param i  index of footprint you want to retrieve
         */
        lsst::detection::Footprint::PtrType getFootprint(int i);

        /** Get model in list
         * 
         * @param i  index of model you want to retrieve
         */
        SpatialModel getModel(int i);

        /** Get vector of all footprints
         */
        FpPtrList getFootprints() {return _fpPtrList;};

        /** Get vector of all models
         */
        ModelPtrList getModels() {return _modelPtrList;};

        /** Get number of models
         */
        int  getNModels() {return _nModels;};

        /** Select best model based upon QA assesment
         *
         * @param fix  optionally fix model as the one to use for this cell
         */
        void selectBestModel(bool fix);

        /** Get index of current model
         */
        int  getCurrentID() {return _currentID;};

        /** Set index of current model
         *
         * @param id  index to use; will build model if not built
         */
        void setCurrentID(int id);

        /** Set label
         *
         * @param label  string that represents the label
         */
        void setLabel(std::string label) {_label = label;};

        /** Get label
         */
        std::string getLabel() {return _label;};

        /** Go to next model in list; call its buildModel() method
         */
        bool incrementModel();

        /** Is cell usable for spatial fit; false if no members or all are bad
         */
        bool isUsable();

        /** Is cell usable but the model is fixed?
         */
        bool isFixed() {return _modelIsFixed;};

    private:
        /** @todo Implement method _orderFootprints
         */
        void _orderFootprints();

        std::string _label;         ///< Name of cell for logging/trace
        int _colC;                  ///< Effective col position of cell in overall image
        int _rowC;                  ///< Effective row position of cell in overall image

        FpPtrList _fpPtrList;       ///< List of footprints in cell
        ModelPtrList _modelPtrList; ///< List of models associated with the footprints

        int _nModels;               ///< Number of entries; len(_fpPtrList)
        int _currentID;             ///< Which entry is being used; 0 <= _currentID < _nModels
        bool _modelIsFixed;         ///< Use model _currentID no matter what

    }; // end of class
 
}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELCELL_H
    
