# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv("ip_diffim",
                    r"$HeadURL$",
                    [["boost", "boost/version.hpp", "boost_filesystem:C++"],
                     ["boost", "boost/regex.hpp", "boost_regex:C++"],
                     ["vw", "vw/Core.h", "vw:C++"],
                     ["vw", "vw/Core.h", "vwCore:C++"],
                     ["vw", "vw/Math.h", "vwMath:C++"],
                     ["lapack", None, "lapack", "dgesdd_"],
                     ["daf_base", "lsst/daf/base/Citizen.h", "daf_base:C++"],
                     ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
                     #["pex_harness", "lsst/pex/harness/Stage.h", "pex_harness:C++"],
                     ["afw", "lsst/afw/MaskedImage.h", "fw"],
                     ["detection", "lsst/detection/Footprint.h", "detection"],
                     ["python", "Python.h"],
                     ["m", "math.h", "m", "sqrt"],
                     ["cfitsio", "fitsio.h", "cfitsio", "ffopen"],
                     ["wcslib", "wcslib/wcs.h", "wcs"],
                     ["xpa", "xpa.h", "xpa", "XPAPuts"],
                     ["minuit", "Minuit/FCNBase.h", "lcg_Minuit:C++"],
                     ])

env.libs["ip_diffim"] = env.getlibs("boost vw wcslib fw cfitsio daf_base daf_data pex_logging minuit detection") + env.libs["ip_diffim"]
env.libs["ip_diffim"] += ["lapack"]     # bug in scons 1.16; getlibs("lapack") fails as lapack isn't in eups

#
# Build/install things
#
for d in Split("doc include/lsst/ip/diffim lib src tests examples"):
    SConscript(os.path.join(d, "SConscript"))

# python
for d in map(lambda str: os.path.join("python/lsst/%s" % env['eups_product'], str), Split(".")):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", env.Install(env['prefix'], "python"))
Alias("install", env.Install(env['prefix'], "include"))
Alias("install", env.Install(env['prefix'], "lib"))
Alias("install", env.Install(env['prefix'], "pipeline"))
Alias("install", env.Install(env['prefix'] + "/bin", glob.glob("bin/*.py")))
Alias("install", env.InstallEups(env['prefix'] + "/ups", glob.glob("ups/*.table")))

scons.CleanTree(r"*~ core *.so *.os *.o")

files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST Image Processing Pipeline packages
""")

