FROM ubuntu:18.04
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>

# Core system capabilities required
RUN apt-get update && apt-get install -y \
    bc=1.07.1-2 \
    build-essential=12.4ubuntu1 \
    curl=7.58.0-2ubuntu3.9 \
    dc=1.07.1-2 \
    git=1:2.17.1-1ubuntu0.7 \
    libegl1-mesa-dev=19.2.8-0ubuntu0~18.04.3 \
    libopenblas-dev=0.2.20+ds-4 \
    nano=2.9.3-2 \
    perl-modules-5.26=5.26.1-6ubuntu0.3 \
    python2.7=2.7.17-1~18.04ubuntu1 \
    python3=3.6.7-1~18.04 \
    python3-pip=9.0.1-2.3~ubuntu1.18.04.1 \
    tar=1.29b-2ubuntu0.1 \
    tcsh=6.20.00-7 \
    unzip=6.0-21ubuntu1 \
    wget=1.19.4-1ubuntu2.2


RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y tzdata=2020a-0ubuntu0.18.04

# PPA for newer version of nodejs, which is required for bids-validator
RUN curl -sL https://deb.nodesource.com/setup_12.x -o nodesource_setup.sh && \
    bash nodesource_setup.sh && \
    rm -f nodesource_setup.sh
RUN apt-get install -y nodejs=12.18.2-1nodesource1

# NeuroDebian setup
RUN wget -qO- http://neuro.debian.net/lists/bionic.au.full | \
    tee /etc/apt/sources.list.d/neurodebian.sources.list
COPY neurodebian.gpg /neurodebian.gpg
RUN apt-key add /neurodebian.gpg && \
    apt-get update

# Additional dependencies for MRtrix3 compilation
RUN apt-get update && apt-get install -y \
    libeigen3-dev=3.3.4-4 \
    libfftw3-dev=3.3.7-1 \
    libpng-dev=1.6.34-1ubuntu0.18.04.2 \
    libtiff5-dev=4.0.9-5ubuntu0.3 \
    zlib1g-dev=1:1.2.11.dfsg-0ubuntu2

# Python dependencies for eddyqc
RUN pip3 install --upgrade pip==20.1.1
RUN pip3 install \
    cycler==0.10.0 \
    kiwisolver==1.2.0 \
    matplotlib==3.3.0rc1 \
    nibabel==2.5.2 \
    numpy==1.16 \
    pandas==1.0.5 \
    Pillow==7.2.0 \
    pyparsing==2.4.7 \
    PyPDF2==1.26.0 \
    python-dateutil==2.8.1 \
    pytz==2020.1 \
    scipy==1.5.0 \
    seaborn==0.10.1 \
    setuptools==39.0.1 \
    six==1.11.0

# Neuroimaging software / data dependencies
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.1.0/freesurfer-linux-centos8_x86_64-7.1.0.tar.gz | \
    tar zx -C /opt \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt'
RUN echo "cHJpbnRmICJyb2JlcnQuc21pdGhAZmxvcmV5LmVkdS5hdVxuMjg1NjdcbiAqQ3FLLjFwTXY4ZE5rXG4gRlNvbGRZRXRDUFZqNlxuIiA+IC9vcHQvZnJlZXN1cmZlci9saWNlbnNlLnR4dAo=" | base64 -d | sh
RUN FREESURFER_HOME=/opt/freesurfer /bin/bash -c 'source /opt/freesurfer/SetUpFreeSurfer.sh'
RUN apt-get install -y ants=2.2.0-1ubuntu1
# FSL installer appears to now be ready for use with version 6.0.0
# eddy is also now included in FSL6
RUN wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
    chmod 775 fslinstaller.py
RUN python2 /fslinstaller.py -d /opt/fsl -V 6.0.3 -q
RUN rm -f /fslinstaller.py
RUN which immv || ( rm -rf /opt/fsl/fslpython && /opt/fsl/etc/fslconf/fslpython_install.sh -f /opt/fsl )
RUN git clone https://git.fmrib.ox.ac.uk/matteob/eddy_qc_release.git /opt/eddyqc && \
    cd /opt/eddyqc && git checkout v1.0.2 && python3 ./setup.py install && cd /
RUN wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | \
    tar zx -C /opt
RUN npm install -gq bids-validator@1.5.3

# apt cleanup to recover as much space as possible
RUN apt-get remove -y libegl1-mesa-dev && apt-get autoremove -y
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

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
RUN wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528816/lh.HCPMMP1.annot" \
    -O /opt/freesurfer/subjects/fsaverage/label/lh.HCPMMP1.annot
RUN wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528819/rh.HCPMMP1.annot" \
    -O /opt/freesurfer/subjects/fsaverage/label/rh.HCPMMP1.annot
RUN wget -q "http://ddl.escience.cn/f/IiyU?func=download&rid=8135438" -O /opt/freesurfer/average/rh.BN_Atlas.gcs
RUN wget -q "http://ddl.escience.cn/f/IiyP?func=download&rid=8135433" -O /opt/freesurfer/average/lh.BN_Atlas.gcs
RUN wget -q "http://ddl.escience.cn/f/PC7Q?func=download&rid=9882718" -O /opt/freesurfer/average/BN_Atlas_subcortex.gca
RUN mkdir /opt/brainnetome && \
    wget -q "http://ddl.escience.cn/f/PC7O?func=download&rid=9882716" -O /opt/brainnetome/BN_Atlas_246_LUT.txt
RUN wget -q "http://ddl.escience.cn/f/Bvhg?func=download&rid=6516020" -O /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz
RUN cp /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/freesurfer/
RUN wget -qO- "https://github.com/ThomasYeoLab/CBIG/archive/v0.11.1-Wu2017_RegistrationFusion.tar.gz" | \
    tar zx -C /opt && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/*h.Yeo2011_*Networks_N1000.split_components.annot /opt/freesurfer/subjects/fsaverage5/label/ && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/project_to_individual/Yeo2011_*networks_Split_Components_LUT.txt /opt/freesurfer/ && \
    mkdir /opt/Yeo2011 && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_*Networks_N1000.split_components.FSL_MNI152_*mm.nii.gz /opt/Yeo2011/ && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/*Networks_ColorLUT_freeview.txt /opt/Yeo2011/ && \
    rm -rf /opt/CBIG-0.11.1-Wu2017_RegistrationFusion


# Make ANTS happy
ENV ANTSPATH /usr/lib/ants
ENV PATH /usr/lib/ants:$PATH

# Make FreeSurfer happy
ENV PATH /opt/freesurfer/bin:/opt/freesurfer/mni/bin:$PATH
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

# Make FSL happy
ENV FSLDIR /opt/fsl
ENV PATH $FSLDIR/bin:$PATH
RUN /bin/bash -c 'source /opt/fsl/etc/fslconf/fsl.sh'
ENV FSLMULTIFILEQUIT TRUE
ENV FSLOUTPUTTYPE NIFTI

# Make ROBEX happy
ENV PATH /opt/ROBEX:$PATH

# MRtrix3 setup
RUN git clone -b 3.0.1 --depth 1 https://github.com/MRtrix3/mrtrix3.git mrtrix3 && \
    cd mrtrix3 && \
    python3 configure -nogui && \
    python3 build -persistent -nopaginate && \
    git describe --tags > /mrtrix3_version && \
    rm -rf cmd/ core/ src/ testing/ tmp/ && \
    cd /

# Setup environment variables for MRtrix3
ENV PATH /mrtrix3/bin:$PATH
ENV PYTHONPATH /mrtrix3/lib:$PYTHONPATH

# Acquire script to be executed
COPY mrtrix3_connectome.py /mrtrix3_connectome.py
RUN chmod 775 /mrtrix3_connectome.py

COPY version /version

ENTRYPOINT ["/mrtrix3_connectome.py"]
