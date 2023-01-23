#FROM debian:latest
FROM ubuntu:18.04


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

COPY config.cfg /tmp/icc-config.cfg
COPY license.lic /tmp/icc-license.lic

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git cpio build-essential python3 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Copy from local OR
#COPY parallel_studio_xe_2020_update1_cluster_edition.tgz /tmp/intel.tgz 
#COPY parallel_studio_xe_2020_update1_cluster_edition_online.tgz /tmp/intel.tgz
COPY my_parallel_studio_xe_2020_update1_cluster_edition.tgz /tmp/intel.tgz
# OR download:
# RUN wget -O intel.tgz http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16526/parallel_studio_xe_2020_update1_cluster_edition.tgz

#RUN tar -xvzf intel.tgz 

# && \

RUN tar -xvzf intel.tgz && \
  cd /tmp/parallel_studio_xe_* && \
  bash ./install.sh --cli-mode --silent=/tmp/icc-config.cfg && \ 
  cd /tmp && \
  rm -rf parallel_studio_xe_* intel.tgz && \
  rm /tmp/icc-config.cfg

#RUN echo "source /opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020/linux/pkg_bin/compilervars.sh intel64 \&\& source /opt/intel/parallel_studio_xe_2020/compilers_and_libraries_2020/linux/mkl/bin/mklvars.sh intel64" >> ~/.bashrc

RUN echo "echo Heloo Adam!!!" >> ~/.bashrc
#RUN echo "source /opt/intel/bin/compilervars.sh intel64" >> ~/.bashrc
RUN echo "env" >> ~/.bashrc
RUN echo "echo ---------------------" >> ~/.bashrc

#RUN source ~/.bashrc
ENV CC=icc CXX=icpc FC=ifort

##ENV TINI_VERSION v0.16.1
##ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
##RUN chmod +x /usr/bin/tini

#ENTRYPOINT [ "/usr/bin/tini", "--" ]
#CMD [ "/bin/bash" ]
# Make the final directories where all our compiled libraries will reside
RUN mkdir /osiris_libs
RUN mkdir /osiris_libs/hdf5
#RUN mkdir /osiris_libs/mpi
RUN mkdir /osiris_libs/fftw

#
# Configure + Build FFTW
#
# Note: We build FFTW twice. The first build is for Double precision floating point.
#                            The second build is for Single precision floating point.
#
WORKDIR /root
RUN wget http://www.fftw.org/fftw-3.3.7.tar.gz
RUN tar xzvf fftw-3.3.7.tar.gz
WORKDIR /root/fftw-3.3.7
#RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && pwd && env && icc --version"
#RUN /bin/bash -c "./configure --prefix=/osiris_libs/fftw && make && make install"
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && ./configure CFLAGS=\"-gcc-sys\" --prefix=/osiris_libs/fftw && make && make install"
RUN make clean
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && ./configure CFLAGS=\"-gcc-sys\" --enable-float --prefix=/osiris_libs/fftw && make && make install"

#
# Configure + Build OpenMPI
#
#WORKDIR /root
#RUN wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz
#RUN tar xvzf openmpi-2.1.1.tar.gz
#WORKDIR /root/openmpi-2.1.1
#RUN ./configure --enable-mpi-fortran=all --enable-static --prefix=/osiris_libs/mpi
#RUN make
#RUN make install
#RUN ldconfig

#
# Setup path for MPI
#   (needed so that we can use 'mpicc', 'mpif90' so that HDF5 can be compiled with MPI support)
#
ENV H5_ROOT=/osiris_libs/hdf5
#ENV OPENMPI_ROOT=/osiris_libs/mpi
#ENV OPENMPI_ROOT=/opt/intel/compilers_and_libraries_2019.4.243/linux/mpi/intel64
#ENV PATH="/opt/intel/compilers_and_libraries_2019.4.243/linux/mpi/intel64/bin:${PATH}"

#
# Configure + Build HDF5 (with MPI support)
#
WORKDIR /root
RUN wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar
RUN tar xvf hdf5-1.10.5.tar
WORKDIR /root/hdf5-1.10.5
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && ./configure --enable-fortran --enable-parallel --enable-shared --prefix=/osiris_libs/hdf5 CC=mpicc FC=mpiifort F90=mpiifort CXX=icpc"
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && make"
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh intel64 && make install"

#
# Final setup and cleanup
#
WORKDIR /root
#RUN ln -s /usr/local/bin/gcc /usr/local/bin/gcc-4.8
ENV LD_LIBRARY_PATH=/osiris_libs/hdf5/lib
run chmod -R 757 /osiris_libs
#ENV PATH="/opt/intel/compilers_and_libraries_2019.4.243/linux/mpi/intel64/bin:${PATH}"
ENV TERM xterm-256color

# Remove build files.. when the image is built with the --squash option,
#  removing these files can reduce the image size by 40 percent.
#RUN rm -rf hdf5* && rm -rf openmpi* && rm -rf fftw*
RUN rm -rf hdf5* && rm -rf fftw*
RUN ln -d /usr/bin/python3 /usr/bin/python
RUN echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/osiris_libs/hdf5/lib:/osiris_libs/fftw/lib" >> ~/.bashrc 

ENTRYPOINT ["/bin/bash", "-c", "source /opt/intel/bin/compilervars.sh intel64 && \"$@\"", "-s"]

#
#
