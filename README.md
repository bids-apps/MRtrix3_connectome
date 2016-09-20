### Description
This BIDS App enables generation and subsequent group analysis of structural connectomes generated from diffusion MRI data. The analysis pipeline relies primarily on the *MRtrix3* software package, and includes a number of state-of-the-art methods for image processing, tractography reconstruction, connectome generation and inter-subject connection density normalisation.

**NOTE**: App is still under development; script is not guaranteed to be operational for all use cases.

### Documentation
Please use the official [*MRtrix3* documentation](http://mrtrix.readthedocs.org) for reference. Additional information may be found in the [online *MRtrix3* community forum](http://community.mrtrix.org).

### Error Reporting
Experiencing problems? You can post a message on the [*MRtrix3* community forum](http://community.mrtrix.org); or if you are confident that what you are experiencing is a genuine issue, you can report it directly to the [GitHub issues list](https://github.com/MRtrix3/mrtrix3/issues). In both cases, please include as much information as possible.

### Acknowledgement
When using this pipeline, please use the following snippet to acknowledge the relevant work (amend as appropriate depending on options used):

Structural connectomes were generated using tools provided in the MRtrix3 software package (http://mrtrix.org). This included: DWI denoising (Veraart et al., 2016), pre-processing (Andersson et al., 2003; Andersson and Sotiropoulos, 2015) and bias field correction (Tustison et al., 2010); inter-modal registration (Bhushan et al., 2015); T1 tissue segmentation (Zhang et al., 2001; Smith, 2002; Patenaude et al., 2011; Smith et al., 2012); spherical deconvolution (Tournier et al., 2004; Jeurissen et al., 2014); probabilistic tractography (Tournier et al., 2010) utilizing ACT (Smith et al., 2012) and dynamic seeding (Smith et al., 2015); SIFT2 (Smith et al., 2015); T1 parcellation (Tzourio-Mazoyer et al., 2002 OR (Dale et al., 1999 AND (Desikan et al., 2006 OR Destrieux et al., 2010) ) ); robust structural connectome construction (Yeh et al., 2016).

```
Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888
Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078
Bhushan, C.; Haldar, J. P.; Choi, S.; Joshi, A. A.; Shattuck, D. W. & Leahy, R. M. Co-registration and distortion correction of diffusion and anatomical images based on inverse contrast normalization. NeuroImage, 2015, 115, 269-280
Dale, A. M.; Fischl, B. & Sereno, M. I. Cortical Surface-Based Analysis: I. Segmentation and Surface Reconstruction NeuroImage, 1999, 9, 179-194
Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data NeuroImage, 2014, 103, 411-426
Desikan, R. S.; Ségonne, F.; Fischl, B.; Quinn, B. T.; Dickerson, B. C.; Blacker, D.; Buckner, R. L.; Dale, A. M.; Maguire, R. P.; Hyman, B. T.; Albert, M. S. & Killiany, R. J. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest NeuroImage, 2006, 31, 968-980
Destrieux, C.; Fischl, B.; Dale, A. & Halgren, E. Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature NeuroImage, 2010, 53, 1-15
Patenaude, B.; Smith, S. M.; Kennedy, D. N. & Jenkinson, M. A Bayesian model of shape and appearance for subcortical brain segmentation. NeuroImage, 2011, 56, 907-922
Smith, S. M. Fast robust automated brain extraction. Human Brain Mapping, 2002, 17, 143-155
Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. Anatomically-constrained tractography: Improved diffusion MRI streamlines tractography through effective use of anatomical information. NeuroImage, 2012, 62, 1924-1938
Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography. NeuroImage, 2015, 119, 338-351
Tournier, J.-D.; Calamante, F., Gadian, D.G. & Connelly, A. Direct estimation of the fiber orientation density function from diffusion-weighted MRI data using spherical deconvolution. NeuroImage,     2004, 23, 1176-1185
Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
Tustison, N.; Avants, B.; Cook, P.; Zheng, Y.; Egan, A.; Yushkevich, P. & Gee, J. N4ITK: Improved N3 Bias Correction. IEEE Transactions on Medical Imaging, 2010, 29, 1310-1320
Tzourio-Mazoyer, N.; Landeau, B.; Papathanassiou, D.; Crivello, F.; Etard, O.; Delcroix, N.; Mazoyer, B. & Joliot, M. Automated Anatomical Labeling of activations in SPM using a Macroscopic Anatomical Parcellation of the MNI MRI single-subject brain. NeuroImage, 15(1), 273–289
Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping using random matrix theory Magn. Res. Med., 2016, early view, doi:10.1002/mrm.26059
Yeh, C.H.; Smith, R.E.; Liang, X.; Calamante, F.; Connelly, A. Correction for diffusion MRI fibre tracking biases: The consequences for structural connectomic metrics. Neuroimage, 2016, doi: 10.1016/j.neuroimage.2016.05.047
Zhang, Y.; Brady, M. & Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Transactions on Medical Imaging, 2001, 20, 45-57
```

### Usage

Command-line usage of the processing script `run.py` is as follows (also accessible by running the script without any command-line options):

#### Synopsis

    run.py [ options ] bids_dir output_dir analysis_level

-  *bids_dir*: The directory with the input dataset formatted according to the BIDS standard.
-  *output_dir*: The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.
-  *analysis_level*: Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: participant

#### Description

Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in MRtrix3

#### Options

#### Options that are relevant to participant-level analysis

+ **--parc**<br>The choice of connectome parcellation scheme. Options are: aal, aal2, fs_2005, fs_2009

+ **--streamlines**<br>The number of streamlines to generate for each subject

#### Options specific to the batch processing of subject data

+ **---participant_label**<br>The label(s) of the participant(s) that should be analyzed. The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a space-separated list.

#### Standard options

+ **--continue <TempDir> <LastFile>**<br>Continue the script from a previous execution; must provide the temporary directory path, and the name of the last successfully-generated file

+ **--force**<br>Force overwrite of output files if pre-existing

+ **--help**<br>Display help information for the script

+ **--nocleanup**<br>Do not delete temporary files during script, or temporary directory at script completion

+ **--nthreads/-n_cpus number**<br>Use this number of threads in MRtrix multi-threaded applications (0 disables multi-threading)

+ **--tempdir /path/to/tmp/**<br>Manually specify the path in which to generate the temporary directory

+ **--quiet**<br>Suppress all console output during script execution

+ **--verbose**<br>Display additional information and progress for every command invoked

+ **--debug**<br>Display additional debugging information over and above the verbose output

#### optional arguments

+ **--v/--version**<br>show program's version number and exit

---

**Author:** Robert E. Smith (robert.smith@florey.edu.au)

**Copyright:** Copyright (c) 2008-2016 the MRtrix3 contributors

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/

MRtrix is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

For more details, see www.mrtrix.org

### Instructions

*** WARNING: Incomplete ***

The [bids/MRtrix3_connectome](https://hub.docker.com/r/bids/mrtrix3_connectome/) Docker container enables users to generate structural connectomes from diffusion MRI data using state-of-the-art techniques. The pipeline requires that data be organized in accordance with the [BIDS specification](http://bids.neuroimaging.io).

In your terminal, type:
```{bash}
$ docker pull bids/mrtrix3_connectome
```

To run the script in participant level mode (for processing one subject only), use e.g.:

```
docker run -i --rm \
    -v /Users/yourname/data/ds005:/bids_dataset \
    -v /Users/yourname/outputs:/outputs \
    bids/example \
    /bids_dataset /outputs participant --participant_label 01 -parc fs_2005
```

Watch this space for the addition of group-level analysis capability.

