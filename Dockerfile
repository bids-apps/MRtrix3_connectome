FROM ubuntu:18.04
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>

# Core system capabilities required
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    bc \
    build-essential \
    curl \
    dc \
    git \
    libegl1-mesa-dev \
    libopenblas-dev \
    nano \
    perl-modules-5.26 \
    python2.7 \
    python3 \
    tar \
    tcsh \
    tzdata \
    unzip \
    wget

# PPA for newer version of nodejs, which is required for bids-validator
RUN curl -sL https://deb.nodesource.com/setup_12.x -o nodesource_setup.sh && \
    bash nodesource_setup.sh && \
    rm -f nodesource_setup.sh && \
    apt-get install -y nodejs

# NeuroDebian setup
COPY neurodebian.gpg /neurodebian.gpg
RUN wget -qO- http://neuro.debian.net/lists/bionic.au.full | \
    tee /etc/apt/sources.list.d/neurodebian.sources.list && \
    apt-key add /neurodebian.gpg && \
    apt-get update

# Additional dependencies for MRtrix3 compilation
RUN apt-get update && apt-get install -y \
    libeigen3-dev \
    libfftw3-dev \
    libpng-dev \
    libtiff5-dev \
    zlib1g-dev

# Neuroimaging software / data dependencies
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.1.1/freesurfer-linux-centos8_x86_64-7.1.1.tar.gz | \
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
RUN apt-get install -y ants=2.2.0-1ubuntu1
# FSL installer appears to now be ready for use with version 6
# eddy is also now included in FSL6
RUN wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
    chmod 775 fslinstaller.py && \
    python2 /fslinstaller.py -d /opt/fsl -V 6.0.4 -q && \
    rm -f /fslinstaller.py
RUN which immv || ( echo "FSLPython not properly configured; re-running" && rm -rf /opt/fsl/fslpython && /opt/fsl/etc/fslconf/fslpython_install.sh -f /opt/fsl || ( cat /tmp/fslpython*/fslpython_miniconda_installer.log && exit 1 ) )
RUN wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | \
    tar zx -C /opt
RUN npm install -gq bids-validator@1.5.3

# apt cleanup to recover as much space as possible
RUN apt-get remove -y libegl1-mesa-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Download additional data for neuroimaging software, e.g. templates / atlases
RUN wget -q https://object.cscs.ch/v1/AUTH_4791e0a3b3de43e2840fe46d9dc2b334/ext-d000035_AAL1Atlas_pub/Release2018_SPM12/aal_for_SPM12.zip && \
    unzip aal_for_SPM12.zip -d /opt && \
    rm -f aal_for_SPM12.zip && \
    wget -qO- http://www.gin.cnrs.fr/wp-content/uploads/aal2_for_SPM12.tar.gz | \
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
    -O /opt/freesurfer/subjects/fsaverage/label/lh.HCPMMP1.annot && \
    wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528819/rh.HCPMMP1.annot" \
    -O /opt/freesurfer/subjects/fsaverage/label/rh.HCPMMP1.annot
RUN mkdir /opt/brainnetome && \
    ( wget -q "http://ddl.escience.cn/f/IiyU?func=download&rid=8135438" -O /opt/freesurfer/average/rh.BN_Atlas.gcs || \
    wget -q "https://osf.io/e6zkg/download" -O /opt/freesurfer/average/rh.BN_Atlas.gcs ) && \
    ( wget -q "http://ddl.escience.cn/f/IiyP?func=download&rid=8135433" -O /opt/freesurfer/average/lh.BN_Atlas.gcs || \
    wget -q "https://osf.io/af9ut/download" -O /opt/freesurfer/average/lh.BN_Atlas.gcs ) && \
    ( wget -q "http://ddl.escience.cn/f/PC7Q?func=download&rid=9882718" -O /opt/freesurfer/average/BN_Atlas_subcortex.gca || \
    wget -q "https://osf.io/k2cd8/download" -O /opt/freesurfer/average/BN_Atlas_subcortex.gca ) && \
    ( wget -q "http://ddl.escience.cn/f/PC7O?func=download&rid=9882716" -O /opt/brainnetome/BN_Atlas_246_LUT.txt || \
    wget -q "https://osf.io/eb7pm/download" -O /opt/brainnetome/BN_Atlas_246_LUT.txt ) && \
    ( wget -q "http://ddl.escience.cn/f/Bvhg?func=download&rid=6516020" -O /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz || \
    wget -q "https://osf.io/dbqep/download" -O /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz ) && \
    cp /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/freesurfer/
RUN wget -qO- "https://github.com/ThomasYeoLab/CBIG/archive/v0.11.1-Wu2017_RegistrationFusion.tar.gz" | \
    tar zx -C /opt && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/*h.Yeo2011_*Networks_N1000.split_components.annot /opt/freesurfer/subjects/fsaverage5/label/ && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/project_to_individual/Yeo2011_*networks_Split_Components_LUT.txt /opt/freesurfer/ && \
    mkdir /opt/Yeo2011 && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_*Networks_N1000.split_components.FSL_MNI152_*mm.nii.gz /opt/Yeo2011/ && \
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/*Networks_ColorLUT_freeview.txt /opt/Yeo2011/ && \
    rm -rf /opt/CBIG-0.11.1-Wu2017_RegistrationFusion

# Setup envvars
ENV ANTSPATH=/usr/lib/ants \
    FREESURFER_HOME=/opt/freesurfer \
    FMRI_ANALYSIS_DIR=/opt/freesurfer/fsfast \
    FSF_OUTPUT_FORMAT=nii.gz \
    FSFAST_HOME=/opt/freesurfer/fsfast \
    LOCAL_DIR=/opt/freesurfer/local \
    MINC_BIN_DIR=/opt/freesurfer/mni/bin \
    MINC_LIB_DIR=/opt/freesurfer/mni/lib \
    MNI_DATAPATH=/opt/freesurfer/mni/data \
    MNI_DIR=/opt/freesurfer/mni \
    MNI_PERL5LIB=/opt/freesurfer/mni/lib/perl5/5.8.5 \
    OS=Linux \
    PERL5LIB=/opt/freesurfer/mni/lib/perl5/5.8.5 \
    SUBJECTS_DIR=/opt/freesurfer/subjects \
    FSLDIR=/opt/fsl \
    FSLOUTPUTTYPE=NIFTI \
    FSLMULTIFILEQUIT=TRUE \
    FSLTCLSH=/opt/fsl/bin/fsltclsh \
    FSLWISH=/opt/fsl/bin/fslwish \
    LD_LIBRARY_PATH=/opt/fsl/lib:$LD_LIBRARY_PATH \
    PATH=/opt/mrtrix3/bin:/usr/lib/ants:/opt/freesurfer/bin:/opt/freesurfer/mni/bin:/opt/fsl/bin:/opt/ROBEX:$PATH \
    PYTHONPATH=/opt/mrtrix3/lib:$PYTHONPATH

# MRtrix3 setup
# Commitish is 3.0.2 plus relevant hotfix
RUN git clone https://github.com/MRtrix3/mrtrix3.git /opt/mrtrix3 && \
    cd /opt/mrtrix3 && \
    git checkout 4ab54489f40997f7da1e1915c2adde3373cf6039 && \
    python3 configure -nogui && \
    python3 build -persistent -nopaginate && \
    git describe --tags > /mrtrix3_version && \
    rm -rf .git/ cmd/ core/ src/ testing/ tmp/ && \
    cd /

# Acquire extra MRtrix3 data
RUN wget -q "https://osf.io/v8n5g/download" -O /opt/mrtrix3/share/mrtrix3/labelconvert/Yeo2011_7N_split.txt && \
    wget -q "https://osf.io/ug2ef/download" -O /opt/mrtrix3/share/mrtrix3/labelconvert/Yeo2011_17N_split.txt

# Acquire script to be executed
COPY mrtrix3_connectome.py /mrtrix3_connectome.py
RUN chmod 775 /mrtrix3_connectome.py

COPY version /version

ENTRYPOINT ["/mrtrix3_connectome.py"]
