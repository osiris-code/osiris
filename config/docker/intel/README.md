## Introduction
Running this Intel Docker image is very similar to running the GNU one in ```/config/docker```. The only main difference is that the Intel ```osiris``` file in this directory 
(```/config/docker/intel/osiris```) runs a different docker image from the GNU version (```/config/docker/osiris```):

```
GNU:   /config/docker/osiris       --> Docker Image: atableman/osirispic
INTEL: /config/docker/intel/osiris --> Docker Image: atableman/osirispic-intel
```

### Quickstart

The following is very similar to that seen in /config/docker/README.md. Read that file for more information.

Run the following two commands. These two commands only need to run _one_ time in order to download the Intel (v2020 Update 1) image to your local machine. 
```
docker login -u atableman -p 876ee8bc-02bb-4403-8671-c8d7093ff213
docker pull atableman/osirispic-intel
```

Now here's an example of how to build and run Osiris (a 2D 'debug' build.. change the various settings as needed)
```
# configure the build. Notice that the system configuration
#   to use is 'docker.intel'
./config/docker/intel/osiris ./configure -s docker.intel -d 2 -t debug
# build
./config/docker/intel/osiris make
# run, using an example input deck in ./decks
./config/docker/intel/osiris ./bin/osiris-2D.e ./decks/test/base-2d
```


## Advanced
### Building new image from scratch
Building a new image from the Dockerfile is a annoying. Ask me about it.. TOSO: document it here one day.

