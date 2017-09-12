cd ~/lsstsw
source setup.csh
setup -t b2983 lsst_apps
setup -t b2983 pipe_tasks
setup -t b2983 obs_decam

cd ~/GIT_REPOS/TICKETS/DM-3704/temp

git clone https://github.com/lsst/pipe_tasks
cd pipe_tasks
git checkout u/djreiss/DM-3704
setup -r . -t b2983
scons
cd ..

git clone https://github.com/lsst/ip_diffim
cd ip_diffim
git checkout u/djreiss/DM-3704
setup -r . -t b2983
scons
cd ..
