# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv("imageproc",
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/tickets/7/SConstruct $",
                    [["boost", "boost/version.hpp", "boost_filesystem:C++"],
                     ["boost", "boost/regex.hpp", "boost_regex:C++"],
                     ["vw", "vw/Core.h", "vw:C++"],
                     ["python", "Python.h"],
                     ["m", "math.h", "m", "sqrt"],
                     ["cfitsio", "fitsio.h", "cfitsio", "ffopen"],
                     ["wcslib", "wcslib/wcs.h", "wcs"],
                     ["xpa", "xpa.h", "xpa", "XPAPuts"],
                     ["mwi", "lsst/mwi/data.h", "mwi:C++"],
                     ["fw", "lsst/fw/MaskedImage.h", "fw"],
                     ["minuit", "Minuit/MnMigrad.h", "lcg_Minuit"],
                     ]
                    )

env.libs["imageproc"] +=  env.getlibs("mwi fw boost wcslib")

#
# Build/install things
#
#for d in Split("doc include/lsst/%s lib src tests" % env['eups_product']):
for d in Split("examples include/lsst/fw lib src tests"):
    SConscript(os.path.join(d, "SConscript"))

for d in map(lambda str: os.path.join("python/lsst/%s" % env['eups_product'], str), Split(".")):
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

