Bootstrap: debootstrap
MirrorURL: http://us.archive.ubuntu.com/ubuntu/
OSVersion: bionic
Include: apt bc build-essential dc file git gnupg libfftw3-dev libpng-dev libtiff5-dev nano python python-numpy tar tzdata unzip wget zlib1g-dev

%labels
MAINTAINER Robert E. Smith <robert.smith@florey.edu.au>
HARDWARE gpu

%files
    mrtrix3_connectome.py /mrtrix3_connectome.py
    neurodebian.gpg /neurodebian.gpg
    version /version

%environment

# ANTs
    ANTSPATH=/usr/lib/ants
    export ANTSPATH

# FreeSurfer
    OS=Linux
    SUBJECTS_DIR=/opt/freesurfer/subjects
    FSF_OUTPUT_FORMAT=nii.gz
    MNI_DIR=/opt/freesurfer/mni
    LOCAL_DIR=/opt/freesurfer/local
    FREESURFER_HOME=/opt/freesurfer
    FSFAST_HOME=/opt/freesurfer/fsfast
    MINC_BIN_DIR=/opt/freesurfer/mni/bin
    MINC_LIB_DIR=/opt/freesurfer/mni/lib
    MNI_DATAPATH=/opt/freesurfer/mni/data
    FMRI_ANALYSIS_DIR=/opt/freesurfer/fsfast
    PERL5LIB=/opt/freesurfer/mni/lib/perl5/5.8.5
    MNI_PERL5LIB=/opt/freesurfer/mni/lib/perl5/5.8.5
    export OS SUBJECTS_DIR FSF_OUTPUT_FORMAT MNI_DIR LOCAL_DIR FREESURFER_HOME FSFAST_HOME MINC_BIN_DIR MINC_LIB_DIR MNI_DATAPATH FMRI_ANALYSIS_DIR PERL5LIB MNI_PERL5LIB
# FSL
    FSLDIR=/opt/fsl
    FSLMULTIFILEQUIT=TRUE
    FSLOUTPUTTYPE=NIFTI
    export FSLDIR FSLMULTIFILEQUIT FSLOUTPUTTYPE
# MRtrix3
    PYTHONPATH=/mrtrix3/lib:$PYTHONPATH
    export PYTHONPATH
# All
    LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda/bin:/.singularity.d/libs:/usr/lib/$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    PATH=/usr/local/cuda/bin:/usr/lib/ants:/opt/freesurfer/bin:/opt/freesurfer/mni/bin:$FSLDIR/bin:/opt/ROBEX:/mrtrix3/bin:$PATH
    export PATH

%post
# Non-interactive installation of packages
    export DEBIAN_FRONTEND=noninteractive

# Grab additional repositories
    sed -i 's/main/main restricted universe multiverse/g' /etc/apt/sources.list
    apt update
    apt upgrade -y

# NeuroDebian setup
    wget -qO- http://neuro.debian.net/lists/bionic.au.full | tee -a /etc/apt/sources.list
    apt-key add /neurodebian.gpg
    apt update
    apt upgrade -y

# Packages that coulnd't be installed upfront
    apt install -y clang libeigen3-dev libopenblas-dev nodejs npm perl-modules tcsh

# CUDA setup
# 18.04 only has CUDA 10, for which there is no eddy_cuda compilation;
# couldn't get CUDA 8 or 9 to install on 16.04
#    echo "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64 /" > /etc/apt/sources.list.d/cuda.list
#    apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
# Try installing 9.1 via runfile instead of repo
#    wget -q http://developer.nvidia.com/compute/cuda/9.1/Prod/local_installers/cuda_9.1.85_387.26_linux --no-check-certificate -O cuda_9.1.85_387.26_linux.run
#    chmod +x cuda_9.1.85_387.26_linux.run
#    sh cuda_9.1.85_387.26_linux.run --silent --override --driver --toolkit --toolkitpath=/usr/local/cuda-9.1
#    ln -s /usr/local/cuda-9.1 /usr/local/cuda
#    rm -f cuda_9.1.85_387.26_linux.run
#
# Attempt to install CUDA 9.1 for eddy_cuda
#    apt install -y cuda-9-1 cuda-runtime-9-1 cuda-demo-suite-9-1 cuda-drivers
#
# Trying duplication of apparently working formula from DKP
    cd /tmp
    wget -q http://developer.nvidia.com/compute/cuda/9.1/Prod/local_installers/cuda_9.1.85_387.26_linux --no-check-certificate -O cuda_9.1.85_387.26_linux.run
    mkdir -p nvidia_installers
    chmod +x cuda_9.1.85_387.26_linux.run
    ./cuda_9.1.85_387.26_linux.run -extract=/tmp/nvidia_installers
    rm -f cuda_9.1.85_387.26_linux.run
    cd nvidia_installers
    rm cuda-samples.9.1.85-23083092-linux.run
    rm NVIDIA-Linux-x86_64-387.26.run
    ls
    ./cuda-linux.9.1.85-23083092.run -noprompt
    cd ..
    rm -rf nvidia_installers
    cd ..
    ln -s /usr/local/cuda-9.1 /usr/local/cuda


# Neuroimaging software / data dependencies
    wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.1/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.1.tar.gz | tar zx -C /opt --exclude='freesurfer/trctrain' --exclude='freesurfer/subjects/fsaverage_sym' --exclude='freesurfer/subjects/fsaverage3' --exclude='freesurfer/subjects/fsaverage4' --exclude='freesurfer/subjects/fsaverage6' --exclude='freesurfer/subjects/cvs_avg35' --exclude='freesurfer/subjects/cvs_avg35_inMNI152' --exclude='freesurfer/subjects/bert' --exclude='freesurfer/subjects/V1_average' --exclude='freesurfer/average/mult-comp-cor' --exclude='freesurfer/lib/qt'
    echo "cHJpbnRmICJyb2JlcnQuc21pdGhAZmxvcmV5LmVkdS5hdVxuMjg1NjdcbiAqQ3FLLjFwTXY4ZE5rXG4gRlNvbGRZRXRDUFZqNlxuIiA+IC9vcHQvZnJlZXN1cmZlci9saWNlbnNlLnR4dAo=" | base64 -d | sh
    apt install -y ants
    wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py
    chmod 775 fslinstaller.py
    /fslinstaller.py -d /opt/fsl -V 6.0.1 -q
    FSLDIR=/opt/fsl /bin/bash -c 'source /opt/fsl/etc/fslconf/fsl.sh'
    wget -qO- "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" | tar zx -C /opt
    npm install -gq bids-validator

# apt cleanup to recover as much space as possible
    apt clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Download additional data for neuroimaging software, e.g. templates / atlases
    wget -qO- http://www.gin.cnrs.fr/AAL_files/aal_for_SPM12.tar.gz | tar zx -C /opt
    wget -qO- http://www.gin.cnrs.fr/AAL2_files/aal2_for_SPM12.tar.gz | tar zx -C /opt
    wget -q https://github.com/AlistairPerry/CCA/raw/master/parcellations/512inMNI.nii -O /opt/512inMNI.nii
    wget -qO- http://www.nitrc.org/frs/download.php/5906/ADHD200_parcellations.tar.gz | tar zx -C /opt
    wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528816/lh.HCPMMP1.annot" -O /opt/freesurfer/subjects/fsaverage/label/lh.HCPMMP1.annot
    wget -q "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/5528819/rh.HCPMMP1.annot" -O /opt/freesurfer/subjects/fsaverage/label/rh.HCPMMP1.annot
    wget -q "http://ddl.escience.cn/f/IiyU?func=download&rid=8135438" -O /opt/freesurfer/average/rh.BN_Atlas.gcs
    wget -q "http://ddl.escience.cn/f/IiyP?func=download&rid=8135433" -O /opt/freesurfer/average/lh.BN_Atlas.gcs
    wget -q "http://ddl.escience.cn/f/PC7Q?func=download&rid=9882718" -O /opt/freesurfer/average/BN_Atlas_subcortex.gca
    mkdir /opt/brainnetome && wget -q "http://ddl.escience.cn/f/PC7O?func=download&rid=9882716" -O /opt/brainnetome/BN_Atlas_246_LUT.txt
    wget -q "http://ddl.escience.cn/f/Bvhg?func=download&rid=6516020" -O /opt/brainnetome/BNA_MPM_thr25_1.25mm.nii.gz
    cp /opt/brainnetome/BN_Atlas_246_LUT.txt /opt/freesurfer/
    wget -qO- "https://github.com/ThomasYeoLab/CBIG/archive/v0.11.1-Wu2017_RegistrationFusion.tar.gz" | tar zx -C /opt
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/*h.Yeo2011_*Networks_N1000.split_components.annot /opt/freesurfer/subjects/fsaverage5/label/
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/project_to_individual/Yeo2011_*networks_Split_Components_LUT.txt /opt/freesurfer/
    mkdir /opt/Yeo2011
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/Yeo2011_*Networks_N1000.split_components.FSL_MNI152_*mm.nii.gz /opt/Yeo2011/
    cp /opt/CBIG-0.11.1-Wu2017_RegistrationFusion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/MNI152/*Networks_ColorLUT_freeview.txt /opt/Yeo2011/
    rm -rf /opt/CBIG-0.11.1-Wu2017_RegistrationFusion

# MRtrix3 setup
    git clone https://github.com/MRtrix3/mrtrix3.git && cd mrtrix3 && git checkout 81036fcc6dc11222515fc6cc1b2403585560bfcb && python configure -nogui && python build -persistent -nopaginate && git describe --tags > /mrtrix3_version

# MRtrix3_connectome script
    chmod 775 /mrtrix3_connectome.py

# Mount points
    mkdir /bids_dataset
    mkdir /output

%runscript
    exec /usr/bin/python /mrtrix3_connectome.py "$@"