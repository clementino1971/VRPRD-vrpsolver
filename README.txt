# The VRPRD application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPRD_APP:/VRPRD bapdock /VRPRD/src/run.jl /VRPRD/data/100/C203.txt -u 588.8 --cfg /VRPRD/config/VRPRD_set_2.cfg -o /VRPRD/sol/C203.sol 

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_VRPRD_APP:/VRPRD bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_VRPRD_APP:/VRPRD bapdock /VRPRD/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPRD_APP:/VRPRD bapdock /VRPRD/src/run.jl -b /VRPRD/solomon.batch

# The application directory (/ABSOLUTE_PATH_TO_VRPRD_APP) was mounted with -v as /VRPRD inside the container. Also, it is possible to mount a different directory to read/write solutions:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPRD_APP:/VRPRD -v /ABSOLUTE_PATH_TO_OUTPUT:/OUT bapdock /VRPRD/src/run.jl /VRPRD/data/100/C203.txt -u 588.8 --cfg /VRPRD/config/VRPRD_set_2.cfg -o /OUT/C203.sol

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/100/C203.txt -u 588.8 --cfg config/VRPRD_set_2.cfg -o sol/C203.sol

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

