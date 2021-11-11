# The VRPTW application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPTW_APP:/VRPTW bapdock /VRPTW/src/run.jl /VRPTW/data/100/C203.txt -u 588.8 --cfg /VRPTW/config/VRPTW_set_2.cfg -o /VRPTW/sol/C203.sol 

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_VRPTW_APP:/VRPTW bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_VRPTW_APP:/VRPTW bapdock /VRPTW/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPTW_APP:/VRPTW bapdock /VRPTW/src/run.jl -b /VRPTW/solomon.batch

# The application directory (/ABSOLUTE_PATH_TO_VRPTW_APP) was mounted with -v as /VRPTW inside the container. Also, it is possible to mount a different directory to read/write solutions:
docker run --rm -v /ABSOLUTE_PATH_TO_VRPTW_APP:/VRPTW -v /ABSOLUTE_PATH_TO_OUTPUT:/OUT bapdock /VRPTW/src/run.jl /VRPTW/data/100/C203.txt -u 588.8 --cfg /VRPTW/config/VRPTW_set_2.cfg -o /OUT/C203.sol

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/100/C203.txt -u 588.8 --cfg config/VRPTW_set_2.cfg -o sol/C203.sol

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see solomon.batch before for adjustments)
./VRPSolver -b solomon.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
