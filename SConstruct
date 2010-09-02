# -*- python -*-
#
# Setup our environment
#
import glob, os.path, sys
import lsst.SConsUtils as scons

dependencies = [
    ["boost", "boost/version.hpp", "boost_system:C++"],
    ["boost", "boost/version.hpp", "boost_filesystem:C++"],
    ["boost", "boost/regex.hpp", "boost_regex:C++"],
    ["boost", "boost/filesystem.hpp", "boost_system:C++"],
    ["boost", "boost/serialization/base_object.hpp", "boost_serialization:C++"],
    ["boost", "boost/test/unit_test.hpp", "boost_unit_test_framework:C++"],
    ["python", "Python.h"],
    ["cfitsio", "fitsio.h", "m cfitsio", "ffopen"],
    ["wcslib", "wcslib/wcs.h", "m wcs"], # remove m once SConsUtils bug fixed
    ["xpa", "xpa.h", "xpa", "XPAPuts"],
    ["minuit2", "Minuit2/FCNBase.h", "Minuit2:C++"],
    ["gsl", "gsl/gsl_matrix.h", "gslcblas gsl"],
    ["eigen", "Eigen/Core.h"],
    ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
    ["utils", "lsst/utils/Utils.h", "utils:C++"],
    ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
    ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
    ["security", "lsst/security/Security.h", "security:C++"],
    ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
    ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
    ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
    ["base", "lsst/base.h"],
    ["afw", "lsst/afw.h", "afw:C++"],
    ["sdqa", "lsst/sdqa/SdqaMetric.h", "sdqa:C++"],
    ]

env = scons.makeEnv(
    "ip_diffim",
    r"$HeadURL$",
    dependencies)
#
# Libraries needed to link libraries/executables
#
env.libs["ip_diffim"] += env.getlibs("boost wcslib cfitsio minuit2 utils daf_base daf_data daf_persistence pex_exceptions pex_logging pex_policy security afw gsl eigen sdqa")

#
# Build/install things
#
for d in (
    ".",
    "doc",
    "examples",
    "lib",
    "python/lsst/ip/diffim",
    "tests",
):
    if d != ".":
        try:
            SConscript(os.path.join(d, "SConscript"))
        except Exception, e:
            print >> sys.stderr, "%s: %s" % (os.path.join(d, "SConscript"), e)
    Clean(d, Glob(os.path.join(d, "*~")))
    Clean(d, Glob(os.path.join(d, "*.pyc")))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.Install(env['prefix'], "examples"),
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.Install(env['prefix'], "pipeline"),
    env.Install(env['prefix'], "python"),
    env.Install(env['prefix'], "src"),
    env.Install(env['prefix'], "tests"),
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"), os.path.join("doc", "htmlDir")),
    env.InstallEups(env['prefix'] + "/ups"),
])

scons.CleanTree(r"*~ core *.so *.os *.o")
#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST Image Processing Pipeline Difference Imaging packages
""")
