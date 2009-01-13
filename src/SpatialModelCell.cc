// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of SpatialModelCell class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include <lsst/ip/diffim/SpatialModelCell.h>
#include <lsst/ip/diffim/SpatialModelKernel.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>

namespace lsst {
namespace ip {
namespace diffim {

template <typename ImageT, typename MaskT, class ModelT>
SpatialModelCell<ImageT, MaskT, ModelT>::SpatialModelCell(
    std::string                                      label,
    std::vector<lsst::detection::Footprint::PtrType> fpPtrList,
    std::vector<ModelT>                              modelPtrList) :
    lsst::daf::data::LsstBase(typeid(this)),
    _label(label),
    _colC(0),
    _rowC(0),
    _fpPtrList(fpPtrList),
    _modelPtrList(modelPtrList),
    _nModels(fpPtrList.size()),
    _currentID(-1),
    _modelIsFixed(false)
{
    SpatialModelCell(label, 0, 0, fpPtrList, modelPtrList);
}

template <typename ImageT, typename MaskT, class ModelT>
SpatialModelCell<ImageT, MaskT, ModelT>::SpatialModelCell(
    std::string                                      label,
    int                                              colC, 
    int                                              rowC, 
    std::vector<lsst::detection::Footprint::PtrType> fpPtrList,
    std::vector<ModelT>                              modelPtrList) :
    lsst::daf::data::LsstBase(typeid(this)),
    _label(label),
    _colC(colC),
    _rowC(rowC),
    _fpPtrList(fpPtrList),
    _modelPtrList(modelPtrList),
    _nModels(fpPtrList.size()),
    _currentID(-1),
    _modelIsFixed(false)
{
    if (!(fpPtrList.size() == modelPtrList.size())) {
        throw lsst::pex::exceptions::DomainError("SpatialModelCell : footprint list and model list are not the same size");
    }
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.SpatialModelCell.SpatialModelCell", 
                                  "Cell %s : created with %d models", this->_label.c_str(), this->_nModels);
    this->_orderFootprints();
}

/** Order footprints; currently this just sets their IDs 
 *
 * @note Synchronously modify fpPtrList and modelPtrList
 */
template <typename ImageT, typename MaskT, class ModelT>
void SpatialModelCell<ImageT, MaskT, ModelT>::_orderFootprints() {
    for (int i = 0; i < this->_nModels; i++) {
        this->_modelPtrList[i]->setID(i);
    }
    // Initialize first model
    this->incrementModel();
}

/** Select best model; if no good models, set Cell as fixed with ID=-1.
 *
 * @note Currently this does *not* use Sdqa objects, but it will.  It
 * only selects the first "good" model.
 * 
 * @note Optionally, if all models are *really* bad (this needs to
 * defined) set Cell as fixed with ID=-1
 */
template <typename ImageT, typename MaskT, class ModelT>
void SpatialModelCell<ImageT, MaskT, ModelT>::selectBestModel(bool fix) {
    bool found = false;
    for (int i = 0; i < this->_nModels; i++) {
        if (this->_modelPtrList[i]->isGood()) {
            this->setCurrentID(i);
            this->_modelIsFixed = fix;
            found = true;
            break;
        }
    }

    if (found == false) {
        this->_currentID    = -1;
        this->_modelIsFixed = true;
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelCell.selectBestModel", 
                                      "Cell %s : Locking with no good models", this->_label.c_str());
    }
}

/** Determine if cell has a usable model
 *
 * @note On the condition the cell is "fixed" and has "currentID =
 * -1", the cell is not usable.
 */
template <typename ImageT, typename MaskT, class ModelT>
bool SpatialModelCell<ImageT, MaskT, ModelT>::isUsable() {
    if ( (this->_currentID == -1) && (this->_modelIsFixed) ) {
        return false;
    }
    return true;
}

template <typename ImageT, typename MaskT, class ModelT>
lsst::detection::Footprint::PtrType SpatialModelCell<ImageT, MaskT, ModelT>::getFootprint(int i) {
    if ( (i < 0) || (i >= this->_nModels) ) {
        throw lsst::pex::exceptions::DomainError("Index out of range");
    }        
    return this->_fpPtrList[i];
}

template <typename ImageT, typename MaskT, class ModelT>
ModelT SpatialModelCell<ImageT, MaskT, ModelT>::getModel(int i) {
    if ( (i < 0) || (i >= this->_nModels) ) {
        throw lsst::pex::exceptions::DomainError("Index out of range");
    }        
    return this->_modelPtrList[i];
}

/** Move to the next model in the list
 *
 * @note  Method setCurrentID() actually builds the model
 */
template <typename ImageT, typename MaskT, class ModelT>
bool SpatialModelCell<ImageT, MaskT, ModelT>::incrementModel() {
    /* If the Cell has a fixed Kernel */
    if (this->_modelIsFixed) {
        return false;
    }

    if (this->_currentID == -1) {
        /* Its the first time through */

        if ((this->_nModels) == 0) {
            /* There are 0 Footprints */
            this->_modelIsFixed = true;
            return false;
        }
        else {
            /* There are at least 1 */
            this->setCurrentID(0);
            return true;
        }            
    }
    else {
        if ( (this->_currentID) == ((this->_nModels) - 1) ) {
            /* You are at the last one */
            this->selectBestModel(true);
            return false;
        }
        else {
            /* Standard increment */
            this->setCurrentID(this->_currentID + 1);
            return true;
        }
    }
}

template <typename ImageT, typename MaskT, class ModelT>
void SpatialModelCell<ImageT, MaskT, ModelT>::setCurrentID(int id) {
    if ( (id < 0) || (id >= this->_nModels) ) {
        throw lsst::pex::exceptions::DomainError("Index out of range");
    }        

    this->_currentID = id;
    lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelCell.setCurrentID", 
                                  "Cell %s : Building Footprint %d / %d", this->_label.c_str(), this->_currentID+1, this->_nModels);

    // If the model does not build for some reason, move on to the next model
    if (! (this->_modelPtrList[this->_currentID]->isBuilt()) )
        if (! (this->_modelPtrList[this->_currentID]->buildModel()) ) 
            this->incrementModel();
}

// Explicit instantiations
template class SpatialModelCell<float, lsst::afw::image::maskPixelType, 
                                lsst::ip::diffim::SpatialModelKernel<float, lsst::afw::image::maskPixelType>::Ptr >;
template class SpatialModelCell<double, lsst::afw::image::maskPixelType, 
                                lsst::ip::diffim::SpatialModelKernel<double, lsst::afw::image::maskPixelType>::Ptr >;

}}} // end of namespace lsst::ip::diffim

