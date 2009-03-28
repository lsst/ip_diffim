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

namespace lsst {
namespace ip {
namespace diffim {

    /** 
     * 
     * @brief Class to ensure constraints for spatial modeling
     * 
     * A given MaskedImage will be divided up into cells, with each cell
     * represented by an instance of this class.  Each cell itself contains a
     * list of instances of classes on which this is templated.  One class
     * member from each cell will be fit for a spatial model.  In case of a poor
     * fit, the next class instance in the list will be fit for.  If all
     * instances in a list are rejected from the spatial model, the best one
     * will be used.
     *
     * Requirements on the class needed by SpatialModelCell are :
     * 
     * Method : double returnCellRating()
     *     A way to rank the object w.r.t. other instances of the class
     *       **before** construction of the model.  E.g. total flux.
     * 
     * Method : bool isBuilt()
     *     Has the model been build yet?
     *     
     * Method : bool buildModel()
     *     Build the model.
     */
    template <typename ClassT>
    class SpatialModelCell {
        
    public:
        typedef boost::shared_ptr<SpatialModelCell<ClassT> > Ptr;
        typedef std::vector<Ptr> SpatialModelCellList;

        typedef boost::shared_ptr<ClassT> SpatialModel;
        typedef std::vector<SpatialModel> ModelPtrList;

        /** Constructor
         *
         * @param label  string representing "name" of cell
         * @param modelPtrList  vector of pointers to models of the function you are fitting for
         */
        SpatialModelCell(std::string label,
                         ModelPtrList modelPtrList);

        /** Constructor
         *
         * @param label  string representing "name" of cell
         * @param colC  effective location of column center of cell within overall MaskedImage
         * @param rowC  effective location of row center of cell within overall MaskedImage
         * @param modelPtrList  vector of pointers to models of the function you are fitting for
         */
        SpatialModelCell(std::string label, int colC, int rowC, 
                         ModelPtrList modelPtrList);

        /** Destructor
         */
        virtual ~SpatialModelCell() {};

        /** Get Col position in overall image
         */
        int getCol() const { return _colC; }

        /** Get row position in overall image
         */
        int getRow() const { return _rowC; }
        
        /** Get current model
         */
        SpatialModel getCurrentModel() {return _modelPtrList[_currentId];};

        /** Get model in list
         * 
         * @param i  index of model you want to retrieve
         */
        SpatialModel getModel(int i);

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
        int  getCurrentId() {return _currentId;};

        /** Set index of current model
         *
         * @param id  index to use; will build model if not built
         */
        void setCurrentId(int id);

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
        std::string _label;         ///< Name of cell for logging/trace
        int _colC;                  ///< Effective col position of cell in overall image
        int _rowC;                  ///< Effective row position of cell in overall image

        ModelPtrList _modelPtrList; ///< List of models associated with the cell

        int _nModels;               ///< Number of entries; len(_modelPtrList)
        int _currentId;             ///< Which entry is being used; 0 <= _currentId < _nModels
        bool _modelIsFixed;         ///< Use model _currentId no matter what

        void _orderModels();        ///< Orders the models

    }; // end of class

    /**
     * @brief Class for sorting spatial models based on their cell rating
     */
    template <typename ClassT>
    class cmpSpatialModels {
    public:
        bool operator() (boost::shared_ptr<ClassT> sm1,
                         boost::shared_ptr<ClassT> sm2) const
            {
                return (sm1->returnCellRating() < sm2->returnCellRating());
            }
    }; 


}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELCELL_H
    
