# -*- python -*-
#
# Setup our environment
#
import glob, os.path, sys
import lsst.SConsUtils as scons

try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.SconsUtils.ConfigureDependentProducts

env = scons.makeEnv("ip_diffim",
                    r"$HeadURL$",
                    scons.ConfigureDependentProducts("ip_diffim"))

env.libs["ip_diffim"] +=  env.getlibs("boost eigen afw sdqa")

#env.libs["ip_diffim"] += env.getlibs("boost wcslib cfitsio minuit2 utils daf_base daf_data daf_persistence pex_exceptions pex_logging pex_policy security afw gsl eigen sdqa")

if True:
    #
    # Workaround SConsUtils failure to find numpy .h files. Fixed in sconsUtils >= 3.3.2
    #
    import numpy
    env.Append(CCFLAGS = ["-I", numpy.get_include()])

#
# Build/install things
#
for d in Split("doc examples lib src tests python/lsst/ip/diffim"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [
    env.Install(env['prefix'], "doc"),
    env.Install(env['prefix'], "etc"),
    env.Install(env['prefix'], "examples"),
    env.Install(env['prefix'], "include"),
    env.Install(env['prefix'], "lib"),
    env.Install(env['prefix'], "policy"),
    env.Install(env['prefix'], "python"),
    env.Install(env['prefix'], "src"),
    env.Install(env['prefix'], "tests"),
    env.InstallEups(os.path.join(env['prefix'], "ups")),
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
