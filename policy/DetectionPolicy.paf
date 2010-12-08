definitions: {
    ######
    #
    # pure detection 
    #
    detThreshold: {
        type: double
        description: "Value of footprint detection threshold"
        maxOccurs: 1
        default: 10.
    }

    detThresholdType: {
        type: string
        description: "Type of detection threshold"
        maxOccurs: 1
        default: "stdev"
        allowed: {
           value:        "value"
           description:  "Use counts as the detection threshold type"
        }
        allowed: {
           value:        "stdev"
           description:  "Use standard deviation as the detection threshold type"
        }
        allowed: {
           value:        "variance"
           description:  "Use variance as the detection threshold type"
        }
    }


    detOnTemplate: {
        type: bool
        description: "If true run detection on the template (imageToConvolve);
            if false run detection on the science image (imageToNotConvolve)"
        maxOccurs: 1
        default: true
    }

    detBadMaskPlanes: {
        type: string
        description: "Mask planes that lead to an invalid detection.
            Options: EDGE SAT BAD CR INTRP
            E.g. : EDGE SAT BAD allows CR-masked and interpolated pixels"
        default: "EDGE" "SAT" "BAD"
    }

    ######
    #
    # specializations for psf matching selection
    #
    fpNpixMin: {
        type: int
        description: "Minimum number of pixels in an acceptable Footprint"
        maxOccurs: 1
        default: 5
    }

    fpNpixMax: {
        type: int
        description: "Maximum number of pixels in an acceptable Footprint;
            too big and the subsequent convolutions become unwieldy"
        maxOccurs: 1
        default: 500
    }

    fpGrowFwhmScaling: {
        type: double
        description: "Grow the footprint based on the Psf Fwhm;
            should be larger than kernelRadiusFwhmScaling"
        maxOccurs: 1
        default: 10.
    }

    fpGrowPix: {
        type: int
        description: "Grow each raw detection footprint by
            this many pixels.  The smaller the faster; however the kernel sum does
            not converge if the stamp is too small; and the kernel is not
            constrained at all if the stamp is the size of the kernel.  Rule of
            thumb is at least 1.5 times the kernel size.  The grown stamp is
            ~2*fpGrowPix pixels larger in each dimension."
        maxOccurs: 1
        default: 30
    }

    fpGrowMin: {
        type: int
        description: "Minimum amount to grow the footprint"
        maxOccurs: 1
        default: 20
    }

    fpGrowMax: {
        type: int
        description: "Maximum amount to grow the footprint"
        maxOccurs: 1
        default: 40
    }
}
