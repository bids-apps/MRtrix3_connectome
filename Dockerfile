ARG MAKE_JOBS="1"
ARG DEBIAN_FRONTEND="noninteractive"

FROM ubuntu:24.04 AS base

# Core system capabilities that should be present in all images
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    python3 \
    tar \
    unzip \
    wget

FROM base AS aal-downloader
WORKDIR /opt/aal
RUN wget https://data.kg.ebrains.eu/zip?container=https://data-proxy.ebrains.eu/api/v1/buckets/p4791e-ext-d000035_AAL1Atlas_pub?prefix=Release2018_SPM12 \
       -O aal1_for_SPM12.zip && \
   unzip aal1_for_SPM12.zip && \
   rm -f aal1_for_SPM12.zip && \
   unzip aal_for_SPM12.zip && \
   rm -f aal_for_SPM12.zip && \
   mv aal_for_SPM12/* . && \
   rmdir aal_for_SPM12 && \
   rm -rf __MACOSX && \
   wget --no-check-certificate -qO- http://www.gin.cnrs.fr/wp-content/uploads/aal2_for_SPM12.tar.gz | \
   tar zx --strip-components=1

FROM base AS adhd200-downloader
WORKDIR /opt
RUN wget -qO- http://www.nitrc.org/frs/download.php/5906/ADHD200_parcellations.tar.gz | \
    tar zx && \
    rm -f ADHD200_parcellations.tar.gz

FROM base AS ants-installer
WORKDIR /opt/ants
RUN wget -q https://github.com/ANTsX/ANTs/releases/download/v2.6.2/ants-2.6.2-ubuntu18.04-X64-gcc.zip && \
   unzip ants-2.6.2-ubuntu18.04-X64-gcc.zip

FROM base AS brainnetome-downloader
WORKDIR /opt/brainnetome
RUN \
   # freesurfer/average/rh.BN_Atlas.gcs
   #( wget -q "http://ddl.escience.cn/f/IiyU?func=download&rid=8135438" -O rh.BN_Atlas.gcs || \
   ( curl https://pan.cstcloud.cn/unode/stor/downloadByUrl?downloadId=1.eyJidWNrZXQiOiJkZWZhdWx0IiwibGVuIjozODIzMDY1Nywic2l6ZSI6MzgyMzA2NTcsInBvcyI6MCwibmFtZSI6InJoLkJOX0F0bGFzLmdjcyIsImN0aW1lIjoxNzYxNjE3OTkyLCJrZXkiOiJzTEhJaHcxYkZYZ0RTeGRld3I0ZWtxRjJuQ0VBQUFBQ1IxcUIiLCJhZ2UiOjg2NDAwfQ.3670708009 \
       -o rh.BN_Atlas.gcs || \
       wget -q "https://osf.io/e6zkg/download" -O rh.BN_Atlas.gcs) && \
   # freesurfer/average/lh.BN_Atlas.gcs
   #( wget -q "http://ddl.escience.cn/f/IiyP?func=download&rid=8135433" -O lh.BN_Atlas.gcs || \
   ( curl https://pan.cstcloud.cn/unode/stor/downloadByUrl?downloadId=1.eyJidWNrZXQiOiJkZWZhdWx0IiwibGVuIjozNjQ5NzM5OSwic2l6ZSI6MzY0OTczOTksInBvcyI6MCwibmFtZSI6ImxoLkJOX0F0bGFzLmdjcyIsImN0aW1lIjoxNzYxNjE3OTEwLCJrZXkiOiJsMzFISVJNZ3FScXZlREg4S2RPVjRqUXNVRnNBQUFBQ0xPZjMiLCJhZ2UiOjg2NDAwfQ.4292748252 \
       -o lh.BN_Atlas.gcs || \
       wget -q "https://osf.io/af9ut/download" -O lh.BN_Atlas.gcs ) && \
   # freesurfer/average/BN_Atlas_subcortex.gca
   #( wget -q "http://ddl.escience.cn/f/PC7Q?func=download&rid=9882718" -O BN_Atlas_subcortex.gca || \
   ( curl https://pan.cstcloud.cn/unode/stor/downloadByUrl?downloadId=1.eyJidWNrZXQiOiJkZWZhdWx0IiwibGVuIjozNzAxODQwOCwic2l6ZSI6MzcwMTg0MDgsInBvcyI6MCwibmFtZSI6IkJOX0F0bGFzX3N1YmNvcnRleC5nY2EiLCJjdGltZSI6MTc2MTYxODA2Niwia2V5IjoicEhiTXdwbjBKRUd4cmJ2Z0NScDFTdzNsS2RrQUFBQUNOTnNvIiwiYWdlIjo4NjQwMH0.1891919902 \
       -o BN_Atlas_subcortex.gca || \
       wget -q "https://osf.io/k2cd8/download" -O BN_Atlas_subcortex.gca ) && \
   # brainnetome/BN_Atlas_246_LUT.txt
   #( wget -q "http://ddl.escience.cn/f/PC7O?func=download&rid=9882716" -O BN_Atlas_246_LUT.txt || \
   ( curl https://pan.cstcloud.cn/unode/stor/downloadByUrl?downloadId=1.eyJidWNrZXQiOiJkZWZhdWx0IiwibGVuIjo1ODAyLCJzaXplIjo1ODAyLCJwb3MiOjAsIm5hbWUiOiJCTl9BdGxhc18yNDZfTFVULnR4dCIsImN0aW1lIjoxNzYxNjE4MDkwLCJrZXkiOiJ1bGV4dG1Kd0hJeGNjdVFDZWwzdTZnUVJXUThBQUJhcSIsImFnZSI6ODY0MDAsInBhcnRPbmUiOnsic2l6ZSI6NTgwMiwiZm4iOiJoa194MmRaLVFYZy0wLTU4MDIiLCJjcmMzMiI6ODk4MzQwOTA0LCJiaWQiOjEsImNpZCI6MX19.3302844208 \
       -o BN_Atlas_246_LUT.txt || \
       wget -q "https://osf.io/eb7pm/download" -O BN_Atlas_246_LUT.txt ) && \
   # brainnetome/BNA_MPM_thr25_1.25mm.nii.gz
   #( wget -q "http://ddl.escience.cn/f/Bvhg?func=download&rid=6516020" -O BNA_MPM_thr25_1.25mm.nii.gz || \
   ( curl https://pan.cstcloud.cn/unode/stor/downloadByUrl?downloadId=1.eyJidWNrZXQiOiJkZWZhdWx0IiwibGVuIjoxNzEyMzMsInNpemUiOjE3MTIzMywicG9zIjowLCJuYW1lIjoiQk5BX01QTV90aHIyNV8xLjI1bW0ubmlpLmd6IiwiY3RpbWUiOjE3NjE2MTk3NzcsImtleSI6IlFrbzBFajRIRTEyS1U4N0tBOXdVZlVMX1RIY0FBcHpoIiwiYWdlIjo4NjQwMCwicGFydE9uZSI6eyJzaXplIjoxNzEyMzMsImZuIjoibDlsX3ZVVkRRMDAtMC0xNzEyMzMiLCJjcmMzMiI6MjY1NjI1NTM0NywiYmlkIjoxLCJjaWQiOjF9fQ.2692886444 \
       -o BNA_MPM_thr25_1.25mm.nii.gz || \
       wget -q "https://osf.io/dbqep/download" -O BNA_MPM_thr25_1.25mm.nii.gz )
   # cp /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/freesurfer/

FROM base AS freesurfer-installer
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-centos8_x86_64-7.4.1.tar.gz | \
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
RUN wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528816/lh.HCPMMP1.annot" \
   -O /opt/freesurfer/subjects/fsaverage/label/lh.HCPMMP1.annot && \
   wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528819/rh.HCPMMP1.annot" \
   -O /opt/freesurfer/subjects/fsaverage/label/rh.HCPMMP1.annot
RUN echo "cHJpbnRmICJyb2JlcnQuc21pdGhAZmxvcmV5LmVkdS5hdVxuMjg1NjdcbiAqQ3FLLjFwTXY4ZE5rXG4gRlNvbGRZRXRDUFZqNlxuIiA+IC9vcHQvZnJlZXN1cmZlci9saWNlbnNlLnR4dAo=" | base64 -d | sh

FROM base AS fsl-installer
# Installer script needs to not reside in destination installation location
WORKDIR /
RUN wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
   chmod 775 fslinstaller.py && \
   python3 /fslinstaller.py -d /opt/fsl -V 6.0.7.18

FROM base AS mni512-downloader
WORKDIR /opt
RUN wget -q https://github.com/AlistairPerry/CCA/raw/master/parcellations/512inMNI.nii

#FROM base AS mrtrix-30x-builder
#WORKDIR /opt/mrtrix3
#RUN apt-get install -y \
#    build-essential \
#    git \
#    libeigen3-dev \
#    libfftw3-dev \
#    zlib1g-dev
## Commitish is 3.0.5 plus relevant changes for dwicat and -export_grad_fsl hotfix
#RUN git clone https://github.com/MRtrix3/mrtrix3.git . && \
#    git checkout 906730011b5e21f1449cc7d60ec145375de07479 && \
#    python3 configure -nogui && \
#    python3 build -persistent -nopaginate && \
#    git describe --tags > /mrtrix3_version && \
#    rm -rf .git/ cmd/ core/ src/ testing/ tmp/ && \
#    # Some extra lookup tables for parcellations in use
#    wget -q "https://osf.io/v8n5g/download" -O share/mrtrix3/labelconvert/Yeo2011_7N_split.txt && \
#    wget -q "https://osf.io/ug2ef/download" -O share/mrtrix3/labelconvert/Yeo2011_17N_split.txt
#
#FROM base AS mrtrix-31x-builder
#WORKDIR /opt/dev
#RUN apt-get install -y \
#    build-essential \
#    cmake \
#    git \
#    libfftw3-dev \
#    ninja-build \
#    pkg-config \
#    zlib1g-dev
## Committish is tip of MRtrix3 Issue#3029 as at 2025-06-19
## TODO Is it possible to do a limited build followed by a cmake --install?
## TODO This branch doesn't yet have targets per Python command;
##   just remove all unwanted Python commands
#RUN git clone https://github.com/MRtrix3/mrtrix3.git . && \
#    git checkout 3bb025a0f8b9edbed187510fa81c9cce422311d3 && \
#    cmake -Bbuild -GNinja --preset=release -DMRTRIX_BUILD_GUI=OFF && \
#    cmake --build build --target dwidenoise2 && \
#    #cmake --build build --target dwibiasnormmask && \
#    #cmake --build build --target dwicat
#    rm -f build/bin/5ttgen build/bin/dwi2response build/bin/dwibiascorrect build/bin/dwinormalise build/bin/population_template build/bin/dwifslpreproc build/bin/dwigradcheck build/bin/dwishellmath build/bin/for_each build/bin/labelsgmfirst build/bin/mask2glass build/bin/mrtrix_cleanup build/bin/responsemean

FROM base AS mrtrix-builder
WORKDIR /opt/mrtrix3
RUN apt-get install -y \
    build-essential \
    cmake \
    git \
    libfftw3-dev \
    ninja-build \
    pkg-config \
    zlib1g-dev
# Committish is tip of MRtrix3 #3029 as at 2025-10-21
RUN git clone https://github.com/MRtrix3/mrtrix3.git . && \
    git checkout 26965d57b374a733ac0c583d3b92bad17923128a
# Main tip as at 2025-10-21
RUN git clone https://github.com/Lestropie/dwidenoise2.git dwidenoise2 && \
    cd dwidenoise2 && \
    git checkout 37c9b70cf7e1c67ac846311ddd0df925424797f5 && \
    cd ../ && \
    cp -r dwidenoise2/cpp .
# Since external project compilation may not yet be working on 3.1.0,
#   just dump the code contents of this repository into the appropriate location,
#   and the build process of MRtrix3 itself should deal with the issue
# TODO Could also make use of cmake --install
COPY mrtrix3_connectome/ /opt/mrtrix3/python/mrtrix3/commands/mrtrix3_connectome
RUN cmake -Bbuild -GNinja --preset=release -DMRTRIX_BUILD_GUI=OFF && \
    cmake --build build

FROM base AS robex-installer
WORKDIR /opt/robex
RUN wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | \
    tar zx

FROM base AS yeo2011-downloader
WORKDIR /opt/Yeo2011
RUN wget -qO- "https://github.com/ThomasYeoLab/CBIG/archive/v0.11.1-Wu2017_RegistrationFusion.tar.gz" | \
    tar zx && \
    mkdir -p freesurfer/subjects/fsaverage5/label && \
    cp CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/*h.Yeo2011_*Networks_N1000.split_components.annot freesurfer/subjects/fsaverage5/label/ && \
    cp CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/project_to_individual/Yeo2011_*networks_Split_Components_LUT.txt freesurfer/ && \
    cp CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_*Networks_N1000.split_components.FSL_MNI152_*mm.nii.gz . && \
    cp CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/*Networks_ColorLUT_freeview.txt . && \
    rm -rf CBIG-0.11.1-Wu2017_RegistrationFusion

FROM base AS final
# Install runtime system dependencies
RUN apt-get -qq update && \
    apt-get install -yq --no-install-recommends \
    bc \
    dc \
    libfftw3-single3 \
    libfftw3-double3 \
    nano \
    nodejs \
    npm \
    python3 \
    tcsh && \
    # apt cleanup to recover as much space as possible
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY --from=aal-downloader /opt/aal /opt/aal
COPY --from=adhd200-downloader /opt/ADHD200_parcellate_200.nii.gz /opt/ADHD200_parcellate_200.nii.gz
COPY --from=adhd200-downloader /opt/ADHD200_parcellate_400.nii.gz /opt/ADHD200_parcellate_400.nii.gz
COPY --from=ants-installer /opt/ants/ants-2.6.2 /opt/ants
COPY --from=brainnetome-downloader /opt/brainnetome/rh.BN_Atlas.gcs /opt/freesurfer/average/rh.BN_Atlas.gcs
COPY --from=brainnetome-downloader /opt/brainnetome/lh.BN_Atlas.gcs /opt/freesurfer/average/lh.BN_Atlas.gcs
COPY --from=brainnetome-downloader /opt/brainnetome/BN_Atlas_subcortex.gca /opt/freesurfer/average/BN_Atlas_subcortex.gca
COPY --from=brainnetome-downloader /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/brainnetome/BN_Atlas_246_LUT.txt
COPY --from=brainnetome-downloader /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/freesurfer/BN_Atlas_246_LUT.txt
COPY --from=brainnetome-downloader /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz
COPY --from=freesurfer-installer /opt/freesurfer /opt/freesurfer
COPY --from=fsl-installer /opt/fsl /opt/fsl
COPY --from=mni512-downloader /opt/512inMNI.nii /opt/512inMNI.nii
#COPY --from=mrtrix-30x-builder /opt/mrtrix3 /opt/mrtrix3
#COPY --from=mrtrix-31x-builder /opt/dev/build /opt/dev
COPY --from=mrtrix-builder /opt/mrtrix3 /opt/mrtrix3
COPY --from=robex-installer /opt/robex /opt/robex
COPY --from=yeo2011-downloader /opt/Yeo2011 /opt/Yeo2011

RUN mv /opt/Yeo2011/freesurfer/subjects/fsaverage5/label/* /opt/freesurfer/subjects/fsaverage5/label && \
   mv /opt/Yeo2011/freesurfer/*.* /opt/freesurfer/ && \
   rm -rf /opt/Yeo2011/freesurfer

# PPA for newer version of nodejs, which is required for bids-validator
#RUN curl -sL https://deb.nodesource.com/setup_12.x -o nodesource_setup.sh && \
#    bash nodesource_setup.sh && \
#    rm -f nodesource_setup.sh && \
#    apt-get install -y nodejs && \
#    npm install -gq bids-validator@1.5.3
RUN npm install -gq bids-validator@1.15.0

# Setup envvars
ENV ANTSPATH=/opt/ants \
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
    PATH=/opt/mrtrix3/build/bin:/opt/ants/bin:/opt/freesurfer/bin:/opt/freesurfer/mni/bin:/opt/fsl/bin:/opt/ROBEX:$PATH \
    PYTHONPATH=/opt/mrtrix3/lib

## Acquire script to be executed
#COPY mrtrix3_connectome.py /mrtrix3_connectome.py
#RUN chmod 775 /mrtrix3_connectome.py

COPY version /version

ENTRYPOINT ["/opt/mrtrix3/build/bin/mrtrix3_connectome"]
