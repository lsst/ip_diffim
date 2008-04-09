# -*- python -*-
#
# Setup our environment
#
import glob, os.path
import lsst.SConsUtils as scons

env = scons.makeEnv(
    "ip_diffim",
    r"$HeadURL$",
    [
        ["boost", "boost/version.hpp", "boost_filesystem:C++"],
        ["boost", "boost/regex.hpp", "boost_regex:C++"],
        ["boost", "boost/serialization/base_object.hpp", "boost_serialization:C++"],
        ["vw", "vw/Core.h", "vw:C++"],
        ["vw", "vw/Core.h", "vwCore:C++"],
        ["vw", "vw/Math.h", "vwMath:C++"],
        ["vw", "vw/FileIO.h", "vwFileIO:C++"],
        ["vw", "vw/Image.h", "vwImage:C++"],
        ["python", "Python.h"],
        ["m", "math.h", "m", "sqrt"],
        ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"], # remove m once SConsUtils bug fixed
        ["wcslib", "wcslib/wcs.h", "m wcs"], # remove m once SConsUtils bug fixed
        ["xpa", "xpa.h", "xpa", "XPAPuts"],
        ["minuit", "Minuit/FCNBase.h", "lcg_Minuit:C++"],
        ["lapack", None, "lapack", "dgesdd_"],
        ["utils", "lsst/utils/Utils.h", "utils:C++"],
        ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
        ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
        ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
        ["security", "lsst/security/Security.h", "security:C++"],
        ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
        ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
        ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
        ["afw", "lsst/afw.h", "afw:C++"],
        ["mpich2", "mpi.h", "mpich:C++"],
        ["pex_harness", "lsst/pex/harness/Stage.h", "pex_harness:C++"],
        ["detection", "lsst/detection/Footprint.h", "detection"],
    ],
)
env.libs["ip_diffim"] = env.getlibs("boost vw lapack wcslib fw cfitsio daf_base daf_data daf_persistence pex_logging pex_exceptions pex_logging minuit detection")
env.libs["ip_diffim"] += ["lapack"]     # bug in scons 1.16; getlibs("lapack") fails as lapack isn't in eups

#
# Build/install things
#
for d in Split("doc include/lsst/ip/diffim lib python/lsst/ip/diffim tests examples"):
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

