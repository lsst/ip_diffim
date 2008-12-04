// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of Image Subtraction
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */

#ifndef LSST_IMAGEPROC_IMAGESUBTRACT_H
#define LSST_IMAGEPROC_IMAGESUBTRACT_H

#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math/Function.h>
#include <lsst/detection/Footprint.h>

namespace lsst {
namespace ip {
namespace diffim {
    
    /**
     * @brief Class to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup diffim
     */
    template <typename ImageT, typename MaskT>
    class DifferenceImageStatistics : public lsst::daf::base::Persistable,
                                      public lsst::daf::data::LsstBase {
    public:
        DifferenceImageStatistics();
        DifferenceImageStatistics(const lsst::afw::image::MaskedImage<ImageT, MaskT> differenceMaskedImage);

        virtual ~DifferenceImageStatistics() {};
        void setResidualMean(double mean) {_residualMean = mean;}
        void setResidualStd(double std) {_residualStd = std;}
        double getResidualMean() {return _residualMean;}
        double getResidualStd() {return _residualStd;}

        bool evaluateQuality(lsst::pex::policy::Policy &policy);
    private:
        double _residualMean;
        double _residualStd;
    };

    /**
     * @brief Class containing Image Subtraction information for a given
     * Footprint
     * 
     * For each set of input template,science MaskedImages, the Detection
     * pipeline or a database query will return a list of acceptible coordinates
     * around which to build a difference imaging Kernel.  Each of these
     * positions will be assigned a DifferenceImageFootprintInformation instance
     * containing pointers to the Footprint generated at this position, the
     * row,col position itself, and associated subimages in the MaskedImages.
     * 
     * @ingroup diffim
     */
    template <typename ImageT, typename MaskT>
    class DifferenceImageFootprintInformation : public lsst::daf::base::Persistable,
                                                public lsst::daf::data::LsstBase {
    public:
        typedef boost::shared_ptr<DifferenceImageFootprintInformation<ImageT, MaskT> > Ptr;
        typedef std::vector<typename DifferenceImageFootprintInformation<ImageT, MaskT>::Ptr> DifiList;

        /* SWIG DOES NOT LIKE BELOW */
        //typedef lsst::afw::image::MaskedImage<ImageT, MaskT> MaskedImage;
        //typedef boost::shared_ptr<MaskedImage> MaskedImagePtr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > MaskedImagePtr; 

        DifferenceImageFootprintInformation(lsst::detection::Footprint::PtrType footprintPtr,
                                            MaskedImagePtr imageToConvolvePtr,
                                            MaskedImagePtr imageToNotConvolvePtr);
        virtual ~DifferenceImageFootprintInformation() {};

        void setID(int id) {_id = id;};
        int getID() {return _id;};
        void setColcNorm(double colc) {_colcNorm = colc;};
        void setRowcNorm(double rowc) {_rowcNorm = rowc;};
        double getColcNorm() {return _colcNorm;};
        double getRowcNorm() {return _rowcNorm;};

        lsst::detection::Footprint::PtrType getFootprintPtr() {return _footprintPtr;};
        MaskedImagePtr getImageToNotConvolvePtr() {return _imageToNotConvolvePtr;};
        MaskedImagePtr getImageToConvolvePtr() {return _imageToConvolvePtr;};

        void setSingleKernelPtr(boost::shared_ptr<lsst::afw::math::Kernel> ptr) {_singleKernelPtr = ptr;};
        boost::shared_ptr<lsst::afw::math::Kernel> getSingleKernelPtr() {return _singleKernelPtr;};
        void setSingleKernelErrorPtr(boost::shared_ptr<lsst::afw::math::Kernel> ptr) {_singleKernelErrorPtr = ptr;};
        boost::shared_ptr<lsst::afw::math::Kernel> getSingleKernelErrorPtr() {return _singleKernelErrorPtr;};

        void setSingleBackground(double background) {_singleBackground = background;};
        double getSingleBackground() {return _singleBackground;};
        void setSingleBackgroundError(double backgroundError) {_singleBackgroundError = backgroundError;};
        double getSingleBackgroundError() {return _singleBackgroundError;};

        void setSingleStats(DifferenceImageStatistics<ImageT, MaskT> stats) {_singleKernelStats = stats;};
        DifferenceImageStatistics<ImageT, MaskT> getSingleStats() {return _singleKernelStats;};

        DifferenceImageStatistics<ImageT, MaskT> computeImageStatistics(boost::shared_ptr<lsst::afw::math::Kernel> kernelPtr,
                                                                        double background);

        void setStatus(bool status) {_isGood = status;};
        bool getStatus() {return _isGood;};

    private:
        /* running ID */
        int _id;

        /* position of the Footprint in the image, -1 to 1 */
        double _colcNorm;
        double _rowcNorm;

        /* footprint assocated with the object we're building the kernels around */
        lsst::detection::Footprint::PtrType _footprintPtr;

        /* subimages associated with the Footprint */
        MaskedImagePtr _imageToConvolvePtr;    /* Typically the template image */
        MaskedImagePtr _imageToNotConvolvePtr; /* Typically the science image */

        /* results from individual kernel fit */
        boost::shared_ptr<lsst::afw::math::Kernel> _singleKernelPtr;
        boost::shared_ptr<lsst::afw::math::Kernel> _singleKernelErrorPtr;
        double _singleBackground;
        double _singleBackgroundError;
        DifferenceImageStatistics<ImageT, MaskT> _singleKernelStats;

        bool _isGood;
    };
    
    template <typename ImageT, typename MaskT>
    typename DifferenceImageFootprintInformation<ImageT, MaskT>::DifiList
    getGoodFootprints(std::vector<boost::shared_ptr<DifferenceImageFootprintInformation<ImageT, MaskT> > > & difiList );

    template <typename ImageT, typename MaskT>
    bool evaluateDiffimQuality(boost::shared_ptr<DifferenceImageFootprintInformation<ImageT, MaskT> > & difiPtr, 
                               lsst::pex::policy::Policy &policy
        );
    
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateDeltaFunctionKernelSet(
        unsigned int nCols,
        unsigned int nRows
        );

    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateAlardLuptonKernelSet(
        unsigned int nCols,
        unsigned int nRows,
        std::vector<double> const &sigGauss,
        std::vector<double> const &degGauss
        );

    template <typename ImageT, typename MaskT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::Kernel> const &convolutionKernelPtr,
        double background
        );

    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<double> computePsfMatchingKernelForFootprint(
        double &background,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList, ///< Input kernel basis set
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<std::pair<double,double> > computePsfMatchingKernelForFootprint2(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &varianceImage,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<std::pair<double,double> > computePsfMatchingKernelForFootprintGSL(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &varianceImage,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
        lsst::pex::policy::Policy &policy
        );

    template <typename MaskT>
    bool maskOk(
        lsst::afw::image::Mask<MaskT> const &inputMask,
        MaskT const badPixelMask
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage,
        MaskT const badPixelMask
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage
        );

    template <typename ImageT>
    void calculateImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::Image<ImageT> const &inputImage
        );

    template <typename VectorT>
    void calculateVectorStatistics(
        vw::Vector<VectorT> const &inputVector,
        double &mean,
        double &variance
        );

    template <typename PixelT, typename FunctionT>
    void addFunctionToImage(
        lsst::afw::image::Image<PixelT> &image,
        lsst::afw::math::Function2<FunctionT> const &function
        );


    /* This will replace difi */
    template <typename ImageT, typename MaskT>
    class KernelModelQa {
    public: 
        typedef boost::shared_ptr<KernelModelQa<ImageT, MaskT> > Ptr;
        typedef std::vector<typename KernelModelQa<ImageT, MaskT>::Ptr> KernelModelQaList;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > MaskedImagePtr; 

        //KernelModelQa(){;};
        virtual ~KernelModelQa() {};

        KernelModelQa(lsst::detection::Footprint::PtrType fpPtr,
                      MaskedImagePtr miToConvolveParentPtr,
                      MaskedImagePtr miToNotConvolveParentPtr,
                      lsst::afw::math::KernelList<lsst::afw::math::Kernel> kBasisList,
                      lsst::pex::policy::Policy &policy,
                      bool build=false
            );

        void setID(int id) {_id = id;};
        int getID() {return _id;};

        void setColcNorm(double colc) {_colcNorm = colc;};
        double getColcNorm() {return _colcNorm;};

        void setRowcNorm(double rowc) {_rowcNorm = rowc;};
        double getRowcNorm() {return _rowcNorm;};

        void setFootprintPtr(lsst::detection::Footprint::PtrType fpPtr) {_fpPtr = fpPtr;};
        lsst::detection::Footprint::PtrType getFootprintPtr() {return _fpPtr;};

        void setMiToConvolvePtr(MaskedImagePtr miPtr) {_miToConvolvePtr = miPtr;};
        MaskedImagePtr getMiToConvolvePtr() {return _miToConvolvePtr;};

        void setMiToNotConvolvePtr(MaskedImagePtr miPtr) {_miToNotConvolvePtr = miPtr;};
        MaskedImagePtr getMiToNotConvolvePtr() {return _miToNotConvolvePtr;};

        void setKernelPtr(boost::shared_ptr<lsst::afw::math::Kernel> kPtr) {_kPtr = kPtr;};
        boost::shared_ptr<lsst::afw::math::Kernel> getKernelPtr() {return _kPtr;};

        void setKernelErrPtr(boost::shared_ptr<lsst::afw::math::Kernel> kPtr) {_kErrPtr = kPtr;};
        boost::shared_ptr<lsst::afw::math::Kernel> getKernelErrPtr() {return _kErrPtr;};

        void setKernelSum(double kSum) {_kSum = kSum;};
        double getKernelSum() {return _kSum;};

        void setBg(double bg) {_bg = bg;};
        double getBg() {return _bg;};

        void setBgErr(double bgErr) {_bgErr = bgErr;};
        double getBgErr() {return _bgErr;};

        /* Ideally this stuff below will be replaced by SDQA metrics */
        void setStats(DifferenceImageStatistics<ImageT, MaskT> kStats) {_kStats = kStats;};
        DifferenceImageStatistics<ImageT, MaskT> getStats() {return _kStats;};

        /* Using a different kernel and background, compute stats for this Footprint */
        //DifferenceImageStatistics<ImageT, MaskT> computeImageStats(boost::shared_ptr<lsst::afw::math::Kernel> kPtr,
        //double bg);

        void setQaStatus(bool status) {_isGood = status;};
        bool getQaStatus() {return _isGood;};
        bool isGood() {return _isGood;};

        void setBuildStatus(bool built) {_isBuilt = built;};
        bool getBuildStatus() {return _isBuilt;};

        /* actually execute building the Kernel; method needed for SpatialModelCell */
        bool buildModel();

        /* Requirement to be used with SpatialModelCell */
        double returnRating();

    private: 
        /* objects needed to build itself; only initializable from constructor */
        MaskedImagePtr _miToConvolveParentPtr;    /* Typically the template image */
        MaskedImagePtr _miToNotConvolveParentPtr; /* Typically the template image */
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> _kBasisList;
        lsst::pex::policy::Policy _policy;

        /* running ID */
        int _id;

        /* position of the Footprint in the image, -1 to 1 */
        double _colcNorm;
        double _rowcNorm;

        /* footprint assocated with the object we're building the kernels around */
        lsst::detection::Footprint::PtrType _fpPtr;

        /* subimages associated with the Footprint */
        MaskedImagePtr _miToConvolvePtr;    /* Typically the template image */
        MaskedImagePtr _miToNotConvolvePtr; /* Typically the science image */

        /* results from individual kernel fit */
        double _kSum;
        boost::shared_ptr<lsst::afw::math::Kernel> _kPtr;
        boost::shared_ptr<lsst::afw::math::Kernel> _kErrPtr;
        double _bg;
        double _bgErr;

        /* SdqaRating _rating; */
        /* Until we are merged with SDQA, use DifferenceImageStatistics class */
        DifferenceImageStatistics<ImageT, MaskT> _kStats;

        /* Is built */
        bool _isBuilt;

        /* Is usable */
        bool _isGood;
    };


    template <typename ImageT, typename MaskT, class ModelT>
    class SpatialModelCell: public lsst::daf::base::Persistable,
                            public lsst::daf::data::LsstBase {
        
    public:
        typedef boost::shared_ptr<SpatialModelCell<ImageT, MaskT, ModelT> > Ptr;
        typedef std::vector<typename SpatialModelCell<ImageT, MaskT, ModelT>::Ptr> SpatialModelCellList;

        SpatialModelCell(std::string label,
                         std::vector<lsst::detection::Footprint::PtrType> fpPtrList, 
                         std::vector<ModelT> modelPtrList);
        SpatialModelCell(std::string label, int colC, int rowC, 
                         std::vector<lsst::detection::Footprint::PtrType> fpPtrList,
                         std::vector<ModelT> modelPtrList);
        virtual ~SpatialModelCell() {};

        lsst::detection::Footprint::PtrType getCurrentFootprint() {return _fpPtrList[_currentID];};
        ModelT                              getCurrentModel() {return _modelPtrList[_currentID];};

        lsst::detection::Footprint::PtrType getFootprint(int i);
        ModelT                              getModel(int i);

        std::vector<lsst::detection::Footprint::PtrType> getFootprints() {return _fpPtrList;};
        std::vector<ModelT>                              getModels() {return _modelPtrList;};

        void selectBestModel(bool fix); /* select best model based upon QA assesment; optionally fix this for the Cell */

        int  getNModels() {return _nModels;};

        int  getCurrentID() {return _currentID;};
        void setCurrentID(int id);      /* choose a particular entry */

        void setLabel(std::string label) {_label = label;};
        std::string getLabel() {return _label;};

        bool increment();               /* go to the next one in the list; call ModelT's buildModel() method */
        bool isUsable();
        bool isFixed() {return _modelIsFixed;};

    private:
        void _orderFootprints();        /* based upon brightness/proximity */

        std::string _label;

        /* optional : position of the grid point */
        int _colC;
        int _rowC;

        /* synchronized lists of footprints and models */
        std::vector<lsst::detection::Footprint::PtrType> _fpPtrList;
        std::vector<ModelT>                           _modelPtrList;

        /* number of entries; len(_fpPtrList) */
        int _nModels;

        /* which entry we are using; 0 <= _currentID < _nModels */
        /* -1 at initilization */
        int _currentID;

        /* we are using _currentID no matter what */
        bool _modelIsFixed;

    };

}}}

#endif



