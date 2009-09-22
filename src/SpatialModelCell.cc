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

#include <algorithm>

#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/Image.h>

#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>

#include <lsst/ip/diffim/SpatialModelCell.h>
#include <lsst/ip/diffim/SpatialModelKernel.h>

namespace lsst {
namespace ip {
namespace diffim {

#if 0
template <typename ClassT>
SpatialModelCell<ClassT>::SpatialModelCell(
    std::string  label,
    ModelPtrList modelPtrList) :
    _label(label),
    _colC(0),
    _rowC(0),
    _modelPtrList(modelPtrList),
    _nModels(modelPtrList.size()),
    _currentId(-1),
    _modelIsFixed(false)
{
    SpatialModelCell(label, 0, 0, modelPtrList);
}

template <typename ClassT>
SpatialModelCell<ClassT>::SpatialModelCell(
    std::string  label,
    int          colC, 
    int          rowC, 
    ModelPtrList modelPtrList) :
    _label(label),
    _colC(colC),
    _rowC(rowC),
    _modelPtrList(modelPtrList),
    _nModels(modelPtrList.size()),
    _currentId(-1),
    _modelIsFixed(false)
{
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.SpatialModelCell.SpatialModelCell", 
                                  "Cell %s : created with %d models", this->_label.c_str(), this->_nModels);

    this->_orderModels();
    this->incrementModel();
}

/** Order models
 *
 */
template <typename ClassT>
void SpatialModelCell<ClassT>::_orderModels() {
    sort (this->_modelPtrList.begin(), this->_modelPtrList.end(), cmpSpatialModels<ClassT>());
    reverse(this->_modelPtrList.begin(), this->_modelPtrList.end());
}

/** Select best model; if no good models, set Cell as fixed with Id=-1.
 *
 * @note Currently this does *not* use Sdqa objects, but it will.  It
 * only selects the first "good" model.
 * 
 * @note Optionally, if all models are *really* bad (this needs to
 * defined) set Cell as fixed with Id=-1
 *
 * @note For now lets just make it bad
 */
template <typename ClassT>
void SpatialModelCell<ClassT>::selectBestModel(bool fix) {
    bool found = false;
    /*
    for (int i = 0; i < this->_nModels; i++) {
        if (this->_modelPtrList[i]->isGood()) {
            this->setCurrentId(i);
            this->_modelIsFixed = fix;
            found = true;
            break;
        }
    }
    */

    if (found == false) {
        this->_currentId    = -1;
        this->_modelIsFixed = true;
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelCell.selectBestModel", 
                                      "Cell %s : Locking with no good models", this->_label.c_str());
    }
}

/** Determine if cell has a usable model
 *
 * @note On the condition the cell is "fixed" and has "currentId =
 * -1", the cell is not usable.
 */
template <typename ClassT>
bool SpatialModelCell<ClassT>::isUsable() {
    if ( (this->_currentId == -1) && (this->_modelIsFixed) ) {
        return false;
    }
    return true;
}

template <typename ClassT>
typename SpatialModelCell<ClassT>::SpatialModel SpatialModelCell<ClassT>::getModel(int i) {
    if ( (i < 0) || (i >= this->_nModels) ) {
        throw LSST_EXCEPT(lsst::pex::exceptions::Exception, 
                          "Index out of range");
    }        
    return this->_modelPtrList[i];
}

/** Move to the next model in the list
 *
 * @note  Method setCurrentId() actually builds the model
 */
template <typename ClassT>
bool SpatialModelCell<ClassT>::incrementModel() {
    /* If the Cell has a fixed Model */
    if (this->_modelIsFixed) {
        return false;
    }

    if (this->_currentId == -1) {
        /* Its the first time through */

        if ((this->_nModels) == 0) {
            /* There are 0 Models */
            this->_modelIsFixed = true;
            return false;
        }
        else {
            /* There are at least 1 */
            this->setCurrentId(0);
            return true;
        }            
    }
    else {
        if ( (this->_currentId) == ((this->_nModels) - 1) ) {
            /* You are at the last one */
            this->selectBestModel(true);
            return false;
        }
        else {
            /* Standard increment */
            this->setCurrentId(this->_currentId + 1);
            return true;
        }
    }
}

template <typename ClassT>
void SpatialModelCell<ClassT>::setCurrentId(int id) {
    if ( (id < 0) || (id >= this->_nModels) ) {
        throw LSST_EXCEPT(lsst::pex::exceptions::Exception, 
                          "Index out of range");
    }        

    this->_currentId = id;
    lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelCell.setCurrentId", 
                                  "Cell %s : Building Model %d / %d", this->_label.c_str(), this->_currentId+1, this->_nModels);

    // If the model does not build for some reason, move on to the next model
    if (! (this->_modelPtrList[this->_currentId]->isBuilt()) )
        if (! (this->_modelPtrList[this->_currentId]->buildModel()) ) 
            this->incrementModel();
}

// Explicit instantiations
template class SpatialModelCell<SpatialModelKernel<float> >;
template class SpatialModelCell<SpatialModelKernel<double> >;

template class cmpSpatialModels<SpatialModelKernel<float> >;
template class cmpSpatialModels<SpatialModelKernel<double> >;

#endif

}}} // end of namespace lsst::ip::diffim

