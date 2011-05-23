definitions: {
    algorithm: {
        type: string
        description: "How to interpolate the background values"
        maxOccurs: 1
        default: "NATURAL_SPLINE" 
        allowed: {
            value: "CONSTANT"
        }
        allowed: {
            value: "LINEAR"
        }
        allowed: {
            value: "NATURAL_SPLINE"
        }
        allowed: {
            value: "CUBIC_SPLINE"
        }
        allowed: {
            value: "CUBIC_SPLINE_PERIODIC"
        }
        allowed: {
            value: "AKIMA_SPLINE"
        }
        allowed: {
            value: "AKIMA_SPLINE_PERIODIC"
        }
    }

    binsize: {
        type: int
        description: "How large of regions should be used for each background point"
        maxOccurs: 1
        default: 2048
    }

    undersample: {
        type: string
        description: "What to do if there are not enough regions for the interpolation"
        maxOccurs: 1
        default: "REDUCE_INTERP_ORDER"
        allowed: {
            value: "THROW_EXCEPTION"
        }
        allowed: {
            value: "REDUCE_INTERP_ORDER"
        }
        allowed: {
            value: "INCREASE_NXNYSAMPLE"
        }
    }
}

