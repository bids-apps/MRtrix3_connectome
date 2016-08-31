FROM ubuntu:14.04
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>

# Core system capabilities required
RUN apt-get update && apt-get install -y git python software-properties-common wget

# Now that we have software-properties-common, can use add-apt-repository to get to g++ version 5, which is required by JSON for Modern C++
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y g++-5

# Additional dependencies for MRtrix3 compilation
RUN apt-get install -y libeigen3-dev zlib1g-dev

# Additional neuroimaging software dependencies required
RUN wget -O- http://neuro.debian.net/lists/trusty.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && apt-get update
RUN apt-get install -y ants
RUN apt-get install -y fsl-5.0-core
RUN apt-get install -y fsl-5.0-eddy-nonfree
RUN apt-get install -y fsl-first-data
RUN wget -qO- ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.3.0-HCP/freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0-HCP.tar.gz | tar zxv -C /opt
# TODO Acquire own FreeSurfer license file
RUN echo $'fsa@brain.org.au\n13477\n*CCll67nWRT9.\n' > /opt/freesurfer/license

# MRtrix3 setup
# NOTE: After the Stanford Coding Sprint, the command "git checkout stanford" will likely be removed, since prerequisite changes should be merged into master
ENV CXX=/usr/bin/g++-5
RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 && cd mrtrix3 && git checkout stanford && python configure -nogui -verbose && python build

# Acquire script to be executed
RUN cd ../ && wget https://raw.githubusercontent.com/BIDS-Apps/MRtrix3_connectome/master/run.py && chmod 775 run.py

# Setup environment variables
ENV FREESURFER_HOME=/opt/freesurfer
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLMULTIFILEQUIT=TRUE
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV PATH=/opt/freesurfer/bin:/usr/lib/fsl/5.0:/usr/lib/ants:/mrtrix3/release/bin:/mrtrix3/scripts:$PATH
ENV PYTHONPATH=/mrtrix3/scripts

ENTRYPOINT ["/run.py"]
