from __future__ import print_function
import gdb
import re
import sys

try:
    import gdb.printing

    # example from meas_alg
    class CRPixelPrinter(object):
        "Print a CRPixel"

        def __init__(self, val):
            self.val = val

        def to_string(self):
            return "{id=%d (%d, %d)}" % (self.val["id"], self.val["col"], self.val["row"])

    printers = []

    def register(obj):
        "Register my pretty-printers with objfile Obj."

        if obj is None:
            obj = gdb

        for p in printers:
            gdb.printing.register_pretty_printer(obj, p)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def build_ip_diffim_dictionary():
        printer = gdb.printing.RegexpCollectionPrettyPrinter("ip_diffim")

        # example from meas_alg
        #printer.add_printer('lsst::meas::algorithms::CRPixel',
        #                    '^lsst::meas::algorithms::CRPixel', CRPixelPrinter)

        return printer

    printers.append(build_ip_diffim_dictionary())

except ImportError:
    def register(*args, **kwargs):
        print("Your version of gdb is too old to load the ip.diffim python pretty printers", file=sys.stderr)
        pass

    pass
