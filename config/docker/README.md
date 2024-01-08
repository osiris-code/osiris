
## Using Docker + Osiris
This is the Docker configuration for GCC that the build server uses... (currently, GCC8.4). It contains all the needed HDF5, OpenMPI, and FFTW stuffs Osiris needs. It is a turn-key, easy-to-use all-in-one way to build and run Osiris. 

See the _Super Easy I Don't Care About Details Method_ section below to quick-start getting Osiris up and running in a matter of minutes. 

As a next step, if interested, it is highly recommended that you checkout the _Details About Osiris Development Using Docker_. It will only take a few minutes to understand the details, but will really help you to customize the process - especially if you wish to frequently compile and use Osiris.

Note: This commit contains the actual Dockerfile used to make the Docker image - but it's not needed. The image has been uploaded to _DockerHub_ and all will be automatically downloaded upon first use. See the _Docker Image Details_ section below for more details.

## Super Easy I Don't Care About Details Method
So you wanna build Osiris and run it _right now_ and don't care about learning anything about Linux or Docker? Or maybe you want to compile and run Osiris with the exact same settings as the build server uses to see + correct any errors reported? No problem. Just run the following commands from the root directory of the Osiris repository. (These instruction show building a 2D production build. Change thing every _2_ to your desired dimension as needed.)

First, make sure Docker is installed and works. The easiest way to test that Docker is properly installed and ready to use is to run the following:
```
docker run hello-world
```
Now, let's compile Osiris! Run the following 2 commands:
```
./config/docker/osiris -s docker.gnu -d 2 -t production
./config/docker/osiris make
```
Osiris should now be compiled and ready to use. Now run a test simulation! (We will run one of the provided test simulations provided in /decks/test/base2d)
```
./config/docker/osiris bin/osiris-2D.e decks/test/base-2d
```
Finally, to run Osiris on multiple CPU cores via MPI, so something like this (this example runs on 4 cores, change to suit your needs):
```
./config/docker/osiris mpirun -n 4 bin/osiris-2D.e path/to/your/input/deck
```

_One final note: If you are using Docker because you wish to compile/run Osiris using the exact same configuration as the Buildserver, then use the configuration file ```docker.gnu.buildserver``` in debug mode. That is, run the commands_
```
./config/docker/osiris -s docker.gnu.buildserver -d 2 -t debug
./config/docker/osiris make
```
_note: the buildserver compiles using 'debug' mode for speed reasons. This is likely change in the future._

**_Done! Fim! Finis! Das ist alles!_**


## Details About Osiris Developement Using Docker
Let's get into the details about what is actually going on.

I am first going to show the full Docker command sued to run things... it is long and clunky. We will simplify things below.

Assuming Docker is installed, you use the container by invoking:
```
docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic ANY_NORMAL_COMMAND
```
```ANY_NORMAL_COMMAND``` is literally any normal Unix-y command.. and that command will be run from inside the Docker container. Here's a quotidian one:
```
# this will run the humble 'ls' command to list all the files in the current directory.. except
#  that 'ls' will be executed from inside the container.. which makes no difference and will do exactly
#  the same thing 'ls' has always done.
docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic ls
```
(For a more full explanation of the purpose of each part of this Docker command, see the _Docker Command Details_ section below)
### Osiris
Again, I will show the full clunky commands.  We will simplify things below.

Here's how to use this Docker container to build and run Osiris. All these command assume that you are in the root directory of the Osiris repository and that we are making a 2D, production build.
```
# configure the build
./configure -s docker.gnu -d 2 -t production
# run make to make, but from inside the Docker container....
docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic make
# run the Osiris you just compiled.
docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic bin/osiris-2D.e decks/test/base-2d
# if you want to use mutiple CPU cores (via MPI)  - say 4 cores - then:
docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic mpirun -n 4 bin/osiris-2D.e path/to/your/input/deck
```

### More Pleasing Osiris Building
Rather then that whole long command each time, we can make a script to help out or define a Bash function to simplify things.

Choose ONE of the following (or make a better way!)

1. The file _config/docker/osiris_ is a small script.. put this into you path. It runs the long Docker command above and passes in any commands you specify.
2. Run the script _config/docker/install_helper.sh_ ONCE. This install_helper adds the following 3 lines to your _.bashrc_
3. Any other Linux-y way to run the above commands in an easy way.
```bash
osiris() {
    docker run -it -v$(pwd):/osiris:delegated --user $(id -u):$(id -g) -w="/osiris" atableman/osirispic $1 $2 $3 $4 $5 $6
}
```
Which ever route you choose (I like choice 2 personally), you will now have a command _osiris_ at your disposal.. Then the build process looks like this:
```bash 
# run 'ls' again from inside the container!!!!
osiris ls
# configure the build
./configure -s docker.gnu -d 2 -t production
osiris make
osiris bin/osiris-2D.e decks/test/base-2d
osiris mpirun -n 4 bin/osiris-2D.e path/to/your/deck
```
## Docker Command Details
Why do we use the specific options above in the Docker command?
#### -it
Tells docker to run programs in _interactive_ mode (basically so that we can see the output of the programs as they run)
#### -v$(pwd):/osiris
Tells Docker to map the folder _/osiris_ inside the container to the current directory.  
This is part of a larger design decision. Writing files to the docker container itself is  _slow and wasteful_ (this is due to Docker using a layered files system)_._  If we wrote to the container, then every file generated when compiling.. every output HDF5 file.. every file of any kind.. would be added to the container and never removed.. making a container that would grow in size forever until taking the entire hard drive (_important note:  this container grows without bound even if you delete files!! 'Deleting' files in a layered file doesn't actually delete them. A layered file-system tracks the difference between each change written. So, in fact, deleting files actually makes the container grow in size!_ ) 

No problem. To get around this, the container never writes files to itself. It is always run so that any files created while compiling or running Osiris are written to folders  _outside_  the container (i.e. to normal folders on the hard drive, not inside any docker containers). Thus, from the container's perspective, it writes all output files to a folder ```/osiris```. But the  ```-v$(pwd):/osiris``` option ensures that all writes to ```/osiris``` get redirected to the current working directory (i.e. to the normal current directory on the hard drive, not inside any docker containers).
#### :delegated
This is an important (new) option that can dramatically increase performance, especially on Mac OSX. It sets the caching strategy used in ensuring that the state of the file system is consistent between the Docker container and the host Operating System. See https://docs.docker.com/docker-for-mac/osxfs-caching](https://docs.docker.com/docker-for-mac/osxfs-caching for details.

#### --user $(id -u): $(id -g)
This option force Docker to run all programs inside the container as the _current user_ .

Without this option, Docker defaults to running programs inside the container as the _root_ user.  This is not only dangerous but also really annoying. Why annoying? Because all output files (e.g. any Osiris executable, any HDF5 files, etc) will also belong to _root._ And then you will have to run any post-processing Python scripts as  _root_. Then you would have to run any movie/image viewers as _root_. Yuck.

#### -w="/osiris"
Tells Docker to use the directory (from the container's perspective) _/osiris_ as the current working directory when starting up.

### atableman/osirispic
This tells Docker which container to run! I already built the image (with all the GCC, MPI, HDF5, FFTW etc) and uploaded to the public _Dockerhub_ website. The first time this image is run, it will be automatically downloaded from the _Dockerhub_ website. If you would like to run your own version of the Docker image, see the _Building the Docker Image_ section below.

## Docker Image Details
Here, we'll discuss details about the Docker image. The specific versions used are:
```
GCC:      8.3.0
FFTW:     3.3O.7
OpenMPI:  2.2.1
HDF5:     1.10.5
```

#### Locations

Note, all these locations are from the perspective of a program running within a Docker container.

This Docker container has GCC and 3 libraries: HDF5 (with parallel support), OpenMPI and FFTW (Fastest FFT in the West). For reference, these three libraries are installed in the following locations within the Docker container:
```
HDF5:    /osiris_libs/hdf5
OpenMPI: /osiris_libs/mpi
FFTW:    /osiris_libs/fftw
```
The container assumes that all source code files (when compiling) and all output files (when running Osiris) are located with the ```/osiris``` directory. This directory is then mapped to a real directory in the host operating system (usually the current working directory).
#### Building the Docker Image
Most times, the Docker image will not need to be built. It will be  automatically downloaded from the _Dockerhub_ website. But occasionally, you might want to alter and build the image yourself. 

 - The Dockerfile is located in this directory and is named ```Dockerfile.osiris_gcc83```
 - Make any addition/changes you wish
 - Build the image using the command  
 ```
 docker build -t pick_a_name -f Dockerfile.osiris_gcc83 .
 ```
 (then replace ```atableman/osirispic``` in the above commands with ```pick_a_name```)
 - Note: for advanced users: add the ```--squash``` option to produce a smaller image (shrinks the image from about 1.9GB to 1.19GB)


