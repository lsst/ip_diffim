# -*- python -*-
#
# Setup our environment
#
import glob, sys
import LSST.SConsUtils as scons

env = scons.makeEnv("imageproc",
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/branches/scons/SConstruct $",
                    [["support"]])
#
# Build/install things
#
env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$)"

Alias("install", env.Install(env['prefix'], "python"))
Alias("install", env.InstallEups(env['prefix'] + "/ups", glob.glob("ups/*.table")))

env.Declare()
env.Help("""
Help, help, me do.
""")

