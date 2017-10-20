FROM ubuntu:14.04
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>

# Core system capabilities required
RUN apt-get update && apt-get install -y curl git nano perl-modules python software-properties-common tar unzip wget
RUN curl -sL https://deb.nodesource.com/setup_4.x | bash -
RUN apt-get install -y nodejs

# Now that we have software-properties-common, can use add-apt-repository to get to g++ version 5, which is required by JSON for Modern C++
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y g++-5

# Additional dependencies for MRtrix3 compilation
RUN apt-get install -y libeigen3-dev libfftw3-dev libtiff5-dev zlib1g-dev

# Neuroimaging software / data dependencies
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | tar zxv -C /opt \
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
RUN wget -qO- http://neuro.debian.net/lists/trusty.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && apt-get update
RUN apt-get install -y ants
RUN apt-get install -y fsl-5.0-core
RUN apt-get install -y fsl-first-data
RUN apt-get install -y fsl-mni152-templates
# FSL installer appears to not yet be ready for use; encountered a lot of dict key errors
#RUN wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && chmod 775 fslinstaller.py
#RUN apt-get install -y fsl-5.0-eddy-nonfree # Use direct download instead (below) - more up-to-date version
RUN rm -f `which eddy`
RUN mkdir /opt/eddy/
RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/eddy_openmp -O /opt/eddy/eddy_openmp
RUN chmod 775 /opt/eddy/eddy_openmp
RUN wget -qO- http://www.gin.cnrs.fr/AAL_files/aal_for_SPM12.tar.gz | tar zxv -C /opt
RUN wget -qO- http://www.gin.cnrs.fr/AAL2_files/aal2_for_SPM12.tar.gz | tar zxv -C /opt
RUN wget -qO- http://www.nitrc.org/frs/download.php/4499/sri24_anatomy_nifti.zip -O sri24_anatomy_nifti.zip && unzip -qq -o sri24_anatomy_nifti.zip -d /opt/ && rm -f sri24_anatomy_nifti.zip
RUN wget -qO- http://www.nitrc.org/frs/download.php/4508/sri24_labels_nifti.zip -O sri24_labels_nifti.zip && unzip -qq -o sri24_labels_nifti.zip -d /opt/ && rm -f sri24_labels_nifti.zip
RUN wget -q https://github.com/AlistairPerry/CCA/raw/master/parcellations/512inMNI.nii -O /opt/512inMNI.nii
RUN wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | tar zxv -C /opt

RUN npm install -g bids-validator

# Make ANTS happy
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
# Note: Would prefer NIFTI, but need to stick to compressed for now due to FSL Ubuntu not honoring this variable. May be able to revert once fsl.checkFirst() is merged in.
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0

# Make ROBEX happy
ENV PATH=/opt/ROBEX:$PATH

# MRtrix3 setup
ENV CXX=/usr/bin/g++-5
# Note: Current commit being checked out includes various fixes that have been necessary to get test data working; eventually it will instead point to a release tag that includes these updates
RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 && cd mrtrix3 && git checkout fc3fc21 && python configure -nogui && NUMBER_OF_PROCESSORS=1 python build && git describe --tags > /mrtrix3_version
#RUN echo $'FailOnWarn: 1\n' > /etc/mrtrix.conf

# Setup environment variables for MRtrix3
ENV PATH=/mrtrix3/bin:$PATH
ENV PYTHONPATH=/mrtrix3/lib

# Acquire script to be executed
COPY run.py /run.py
RUN chmod 775 /run.py

COPY version /version

ENTRYPOINT ["/run.py"]
