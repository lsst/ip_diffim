# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv("imageproc",
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/tickets/7/SConstruct $",
                    [["boost", "boost/version.hpp", "boost_filesystem:C++"],
                     ["vw", "vw/Core.h", "vw:C++"],
                     ["fw", "lsst/fw/MaskedImage.h lsst/fw/Trace.h", "fw"],
                     ["python", "Python.h"],
                     ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"],
                     ["wcstools", "wcs.h", "wcs", "wcscat"],
                     ["xpa", "xpa.h", "xpa", "XPAPuts"],
                     ["minuit", "Minuit/FCNBase.h"]
                     ])

#
# Libraries that I need to link things.  This should be handled better
#
env.libs = dict([
    ("boost",	Split("boost_filesystem")),
    ("fits",	Split("fitsio")),
    ("vw",	Split("vw vwCore vwFileIO")),
    ("fw",	Split("fw")),
    ("minuit", Split("liblcg_Minuit")),
    ])

#
# Build/install things
#
for d in Split("lib src tests examples"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", env.Install(env['prefix'], "python"))
Alias("install", env.Install(env['prefix'], "include"))
Alias("install", env.Install(env['prefix'], "lib"))
Alias("install", env.Install(env['prefix'] + "/bin", glob.glob("bin/*.py")))
Alias("install", env.InstallEups(env['prefix'] + "/ups", glob.glob("ups/*.table")))

scons.CleanTree(r"*~ core *.so *.os *.o")

env.Declare()
env.Help("""
LSST Image Processing Pipeline packages
""")

