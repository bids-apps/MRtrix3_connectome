FROM ubuntu:18.04
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>

# Core system capabilities required
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    nano \
    nodejs \
    npm \
    perl-modules \
    python \
    tar \
    unzip \
    wget

# ANTs installs tzdata as a dependency; however its installer is interactive
# Therefore we need to do some shenanigans here to force it though
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y tzdata

# NeuroDebian setup
RUN wget -qO- http://neuro.debian.net/lists/artful.au.full | \
    tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9 && apt-get update

# Additional dependencies for MRtrix3 compilation
RUN apt-get install -y \
    libeigen3-dev \
    libfftw3-dev \
    libtiff5-dev \
    zlib1g-dev

# Attempt to install CUDA 8.0 for eddy_cuda
#RUN wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1404/x86_64/cuda-repo-ubuntu1404_8.0.61-1_amd64.deb && \
#    dpkg -i cuda-repo-ubuntu1404_8.0.61-1_amd64.deb && \
#    apt-get update && apt-get install -y cuda && \
#    rm -f cuda-repo-ubuntu1404_8.0.61-1_amd64.deb

# Neuroimaging software / data dependencies
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | \
    tar zx -C /opt \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage5' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt'
RUN apt-get install -y ants
RUN apt-get install -y fsl-5.0-core
RUN apt-get install -y fsl-first-data
RUN apt-get install -y fsl-mni152-templates
# FSL installer appears to not yet be ready for use; encountered a lot of dict key errors
#RUN wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
#    chmod 775 fslinstaller.py
#RUN apt-get install -y fsl-5.0-eddy-nonfree # Use direct download instead (below) - more up-to-date version
RUN rm -f `which eddy`
RUN mkdir /opt/eddy/
RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.11/centos6/eddy_openmp -O /opt/eddy/eddy_openmp
#RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.11/centos6/eddy_cuda8.0 -O /opt/eddy/eddy_cuda8.0
RUN wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | \
    tar zx -C /opt
RUN npm install -g bids-validator

# Download additional data for neuroimaging software, e.g. templates / atlases
RUN wget -qO- http://www.gin.cnrs.fr/AAL_files/aal_for_SPM12.tar.gz | \
    tar zx -C /opt
RUN wget -qO- http://www.gin.cnrs.fr/AAL2_files/aal2_for_SPM12.tar.gz | \
    tar zx -C /opt
#RUN wget -q http://www.nitrc.org/frs/download.php/4499/sri24_anatomy_nifti.zip -O sri24_anatomy_nifti.zip && \
#    unzip -qq -o sri24_anatomy_nifti.zip -d /opt/ && \
#    rm -f sri24_anatomy_nifti.zip
#RUN wget -q http://www.nitrc.org/frs/download.php/4502/sri24_anatomy_unstripped_nifti.zip -O sri24_anatomy_unstripped_nifti.zip && \
#    unzip -qq -o sri24_anatomy_unstripped_nifti.zip -d /opt/ && \
#    rm -f sri24_anatomy_unstripped_nifti.zip
#RUN wget -q http://www.nitrc.org/frs/download.php/4508/sri24_labels_nifti.zip -O sri24_labels_nifti.zip && \
#    unzip -qq -o sri24_labels_nifti.zip -d /opt/ && \
#    rm -f sri24_labels_nifti.zip
RUN wget -q https://github.com/AlistairPerry/CCA/raw/master/parcellations/512inMNI.nii -O /opt/512inMNI.nii
#RUN wget -q https://ndownloader.figshare.com/files/3133832 -O oasis.zip && \
#    unzip -qq oasis.zip -d /opt/ && \
#    rm -f oasis.zip
RUN wget -qO- http://www.nitrc.org/frs/download.php/5906/ADHD200_parcellations.tar.gz | \
    tar zx -C /opt
RUN wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528816/lh.HCPMMP1.annot" -O /opt/freesurfer/subjects/fsaverage/label/lh.HCPMMP1.annot
RUN wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528819/rh.HCPMMP1.annot" -O /opt/freesurfer/subjects/fsaverage/label/rh.HCPMMP1.annot

# Make ANTS happy
ENV ANTSPATH=/usr/lib/ants
ENV PATH=/usr/lib/ants:$PATH

# Make FreeSurfer happy
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/mni/bin:$PATH
ENV OS Linux
ENV SUBJECTS_DIR /opt/freesurfer/subjects
ENV FSF_OUTPUT_FORMAT nii.gz
ENV MNI_DIR /opt/freesurfer/mni
ENV LOCAL_DIR /opt/freesurfer/local
ENV FREESURFER_HOME /opt/freesurfer
ENV FSFAST_HOME /opt/freesurfer/fsfast
ENV MINC_BIN_DIR /opt/freesurfer/mni/bin
ENV MINC_LIB_DIR /opt/freesurfer/mni/lib
ENV MNI_DATAPATH /opt/freesurfer/mni/data
ENV FMRI_ANALYSIS_DIR /opt/freesurfer/fsfast
ENV PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
ENV MNI_PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
RUN echo "cHJpbnRmICJyb2JlcnQuc21pdGhAZmxvcmV5LmVkdS5hdVxuMjg1NjdcbiAqQ3FLLjFwTXY4ZE5rXG4gRlNvbGRZRXRDUFZqNlxuIiA+IC9vcHQvZnJlZXN1cmZlci9saWNlbnNlLnR4dAo=" | base64 -d | sh

# Make FSL happy
ENV PATH=/usr/lib/fsl/5.0:/opt/eddy:$PATH
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLMULTIFILEQUIT=TRUE
ENV FSLOUTPUTTYPE=NIFTI
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
RUN chmod 775 /opt/eddy/eddy_openmp
#RUN chmod 775 /opt/eddy/eddy_cuda8.0

# Make ROBEX happy
ENV PATH=/opt/ROBEX:$PATH

# MRtrix3 setup
RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 && \
    cd mrtrix3 && \
    git checkout 3.0_RC3 && \
    python configure -nogui && \
    python build -persistent -nopaginate && \
    git describe --tags > /mrtrix3_version
#RUN echo $'FailOnWarn: 1\n' > /etc/mrtrix.conf

# Setup environment variables for MRtrix3
ENV PATH=/mrtrix3/bin:$PATH
ENV PYTHONPATH=/mrtrix3/lib:$PYTHONPATH

# Acquire script to be executed
COPY run.py /run.py
RUN chmod 775 /run.py

COPY version /version

ENTRYPOINT ["/run.py"]
