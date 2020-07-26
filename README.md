## Description

This BIDS App enables generation and subsequent group analysis of structural connectomes
generated from diffusion MRI data. The analysis pipeline relies primarily on the *MRtrix3*
software package, and includes a number of state-of-the-art methods for image processing,
tractography reconstruction, connectome generation and inter-subject connection density
normalisation.

**NOTE**: App is still under development; script is not guaranteed to be operational
for all use cases.

## Requirements

Due to use of the Anatomically-Constrained Tractography (ACT) framework, correction of
EPI susceptibility distortions is a prerequisite for this pipeline. Currently, this is
only possible within this pipeline through use of the FSL tool `topup`, which relies
on the presence of spin-echo EPI images with differences in phase encoding to estimate
the causative inhomogeneity field. In the absence of such data, this pipeline is not
currently applicable; though recommendations for alternative mechanisms for such
correction in the Issues page are welcome, and development of novel techniques for
performing this correction are additionally underway.

While many common DICOM conversion software are capable of providing data characterising
the phase and slice encoding performed in the acquisition protocol, which are subsequently
used by this pipeline to automate DWI data pre-processing, for some softwares and/or
some data (particularly those not acquired on a Siemens platform), such data may not be
present in the sidecar JSON files for files in the BIDS `dwi/` and `fmap/` directories.
In this circumstance, it will be necessary for users to manually enter the relevant
information into these files in order for this script to be capable of processing the
data. Every JSON file in these two directories should contain the BIDS fields
`PhaseEncodingDirection` and `TotalReadoutTime`. For DWI data, it is also preferable to
provide the `SliceEncodingDirection` and `SliceTiming` fields. More information on these
data can be found in the [BIDS documentation](http://bids.neuroimaging.io/).

## Instructions

This script can be utilised in one of three ways:

### 1. As a stand-alone *MRtrix3* script

The script ``mrtrix3_connectome.py`` can additionally be used *outside* of this Docker
container, as a stand-alone Python script build against the *MRtrix3* Python libraries.
Using the script in this way requires setting the ``PYTHONPATH`` environment variable to
include the path to the *MRtrix3* ``lib/`` directory where it is installed on your local
system, as described [here](https://community.mrtrix.org/t/the-mrtrix3-python-script-library/2243).
When used in this way, the command-line interface of the script will be more consistent
with the rest of *MRtrix3*. Note that this usage will require version `3.0.0` of *MRtrix3*
to be installed and configured appropriately on your local system.

### 2. As a Docker container

The [bids/MRtrix3_connectome](https://hub.docker.com/r/bids/mrtrix3_connectome/) Docker
container enables users to generate structural connectomes from diffusion MRI data using
state-of-the-art techniques. The pipeline requires that data be organized in accordance
with the [BIDS specification](http://bids.neuroimaging.io).

In your terminal, type:
```{bash}
$ docker pull bids/mrtrix3_connectome
```

To query the help page of the tool:

```{bash}
$ docker run -i --rm bids/mrtrix3_connectome
```

To run the script in participant1 level mode (for processing one subject only), use e.g.:

```{bash}
$ docker run -i --rm \
      -v /Users/yourname/data:/bids_dataset \
      -v /Users/yourname/output:/output \
      bids/mrtrix3_connectome \
      /bids_dataset /output participant1 --participant_label 01 --parcellation desikan
```

Following processing of all participants, the script can be run in group analysis mode
using e.g.:

```{bash}
$ docker run -i --rm \
      -v /Users/yourname/data:/bids_dataset \
      -v /Users/yourname/output:/output \
      bids/mrtrix3_connectome \
      /bids_dataset /output group
```

### 3. As a Singularity container

The *MRtrix3_connectome* BIDS App can also be built locally as a Singularity container.
This is particularly useful for subsequent utilisation on high-performance computing
hardware, as unlike Docker there are no super-user privileges or user group
memberships required for execution. I have also personally been able to utilise the
CUDA version of FSL's `eddy` command within this tool when running on a computing
cluster node with GPU capability (though this can require explicit configuration;
speak to your system administrator).

Within the location in which the *MRtrix3_connectome* source code has been cloned,
type:
```{bash}
$ sudo singularity build MRtrix3_connectome.sif Singularity
```

The resulting container file "`MRtrix3_connectome.sif`" can be run as a stand-alone
executable, as long as the system on which the file is executed has a version of
Singularity installed that is compatible with that of the system used to build the
container.

## Documentation

The help page of the tool itself can be generated by executing the script without
providing any command-line options. The help page is additionally presented at the
bottom of this README page for reference. Documentation regarding the underlying
*MRtrix3* tools can be found in the official
[*MRtrix3* documentation](http://mrtrix.readthedocs.org). Additional information
may be found in the [online *MRtrix3* community forum](http://community.mrtrix.org).

## Error Reporting

Experiencing problems? You can either post a private message to me on the
[*MRtrix3* community forum](http://community.mrtrix.org/u/rsmith), or you can report it
directly to the [GitHub issues list](https://github.com/BIDS-Apps/MRtrix3_connectome/issues).
In both cases, please include as much information as possible; this may include re-running
the script using the ``--debug`` option, which will provide additional information at the
terminal, and preserve temporary files generated by the script within your target output
directory, which can be forwarded to the developer.

## Acknowledgements

Development of this tool was made possible through funding from the National Health
and Medical Research Council (NHMRC) of Australia.

The developer acknowledges the facilities and scientific and technical assistance of
the National Imaging Facility, a National Collaborative Research Infrastructure Strategy
(NCRIS) capability, at the Florey Institute of Neuroscience and Mental Health.

The Florey Institute of Neuroscience and Mental Health acknowledges support from the
Victorian Government and in particular the funding from the Operational Infrastructure
Support Grant.

Robert Smith is supported by fellowship funding from the National Imaging Facility (NIF),
an Australian Government National Collaborative Research Infrastructure Strategy (NCRIS)
capability.

## Citation

When using this pipeline, please use the following snippet to acknowledge the relevant
work (amend as appropriate depending on options used):

Structural connectomes were generated using the *MRtrix3_connectome* BIDS App (Smith
et al., 2019), which operates principally using tools provided in the *MRtrix3*
software package (Tournier et al., 2019; http://mrtrix.org). This included: DWI
denoising (Veraart et al., 2016), Gibbs ringing removal (Kellner et al., 2016),
pre-processing (Andersson et al., 2003; Andersson and Sotiropoulos, 2016; Andersson
et al., 2016; (IF USING EDDY_CUDA: Andersson et al., 2017)); and bias field correction
(Tustison et al., 2010 OR Zhang et al., 2001); inter-modal registration (Bhushan et al.,
2015); brain extraction (Smith, 2002 OR Iglesias et al., 2011), T1 tissue segmentation
(Zhang et al., 2001; Smith, 2002; Patenaude et al., 2011; Smith et al., 2012); spherical
deconvolution (Tournier et al., 2004; Jeurissen et al., 2014); probabilistic tractography
(Tournier et al., 2010) utilizing Anatomically-Constrained Tractography (Smith et al.,
2012) and dynamic seeding (Smith et al., 2015b); SIFT2 (Smith et al., 2015b); T1
parcellation ((((Avants et al., 2008 AND Tustison et al., 2013) OR Andersson et al.,
2010) AND (Tzourio-Mazoyer et al., 2002 OR Yeo et al., 2011 OR Craddock et al., 2012
2011) OR Fan et al., 2016 OR (Zalesky et al., 2010 AND Perry et al., 2017))) OR
2012) (Dale et al., 1999 AND (Desikan et al., 2006 OR Destrieux et al., 2010 OR
2013) Glasser et al., 2016))); robust structural connectome construction (Smith et al.,
2014) 2015a; Yeh et al., 2016); inter-subject connection density normalisation
2015) (Smith et al., 2020).

```
Smith, R. E.; Connelly, A. MRtrix3_connectome: A BIDS Application for quantitative structural connectome construction. In Proc OHBM, 2019, W610
Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888
Andersson, J. L. R.; Jenkinson, M. & Smith, S. Non-linear registration, aka spatial normalisation. FMRIB technical report, 2010, TR07JA2
Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2016, 125, 1063-1078
Andersson, J. L. R. & Graham, M. S. & Zsoldos, E. & Sotiropoulos, S. N. Incorporating outlier detection and replacement into a non-parametric framework for movement and distortion correction of diffusion MR images. NeuroImage, 2016, 141, 556-572
Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H.; Filippini, N. & Bastiani, M. Towards a comprehensive framework for movement and distortion correction of diffusion MR images: Within volume movement. NeuroImage, 2017, 152, 450-466
Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H. & Campbell, J. Susceptibility-induced distortion that varies due to motion: Correction in diffusion MR without acquiring additional data. NeuroImage, 2018, 171, 277-295
Avants, B. B.; Epstein, C. L.; Grossman, M. & Gee, J. C. Symmetric diffeomorphic image registration with cross-correlation: Evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 2008, 12, 26-41
Bhushan, C.; Haldar, J. P.; Choi, S.; Joshi, A. A.; Shattuck, D. W. & Leahy, R. M. Co-registration and distortion correction of diffusion and anatomical images based on inverse contrast normalization. NeuroImage, 2015, 115, 269-280
Craddock, R. C.; James, G. A.; Holtzheimer, P. E.; Hu, X. P.; Mayberg, H. S. A whole brain fMRI atlas generated via spatially constrained spectral clustering. Human Brain Mapping, 2012, 33(8), 1914-1928
Dale, A. M.; Fischl, B. & Sereno, M. I. Cortical Surface-Based Analysis: I. Segmentation and Surface Reconstruction. NeuroImage, 1999, 9, 179-194
Desikan, R. S.; Ségonne, F.; Fischl, B.; Quinn, B. T.; Dickerson, B. C.; Blacker, D.; Buckner, R. L.; Dale, A. M.; Maguire, R. P.; Hyman, B. T.; Albert, M. S. & Killiany, R. J. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. NeuroImage, 2006, 31, 968-980
Destrieux, C.; Fischl, B.; Dale, A. & Halgren, E. Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature. NeuroImage, 2010, 53, 1-15
Fan, L.; Li, H.; Zhuo, J.; Zhang, Y.; Wang, J.; Chen, L.; Yang, Z.; Chu, C.; Xie, S.; Laird, A.R.; Fox, P.T.; Eickhoff, S.B.; Yu, C.; Jiang, T. The Human Brainnetome Atlas: A New Brain Atlas Based on Connectional Architecture. Cerebral Cortex, 2016, 26 (8), 3508-3526
Glasser, M. F.; Coalson, T. S.; Robinson, E. C.; Hacker, C. D.; Harwell, J.; Yacoub, E.; Ugurbil, K.; Andersson, J.; Beckmann, C. F.; Jenkinson, M.; Smith, S. M. & Van Essen, D. C. A multi-modal parcellation of human cerebral cortex. Nature, 2016, 536, 171-178
Iglesias, J. E.; Liu, C. Y.; Thompson, P. M. & Tu, Z. Robust Brain Extraction Across Datasets and Comparison With Publicly Available Methods. IEEE Transactions on Medical Imaging, 2011, 30, 1617-1634
Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426
Kellner, E.; Dhital, B.; Kiselev, V. G.; Reisert, M. Gibbs-ringing artifact removal based on local subvoxel-shifts. Magnetic Resonance in Medicine, 2006, 76(5), 1574-1581
Patenaude, B.; Smith, S. M.; Kennedy, D. N. & Jenkinson, M. A Bayesian model of shape and appearance for subcortical brain segmentation. NeuroImage, 2011, 56, 907-922
Perry, A.; Wen, W.; Kochan, N. A.; Thalamuthu, A.; Sachdev, P. S.; Breakspear, M. The independent influences of age and education on functional brain networks and cognition in healthy older adults. Human Brain Mapping, 2017, 38(10), 5094-5114
Smith, S. M. Fast robust automated brain extraction. Human Brain Mapping, 2002, 17, 143-155
Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. Anatomically-constrained tractography: Improved diffusion MRI streamlines tractography through effective use of anatomical information. NeuroImage, 2012, 62, 1924-1938
Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015a, 104, 253-265
Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography. NeuroImage, 2015b, 119, 338-351
Smith, R. E.; Calamante, F.; Connelly, A. Mapping connectomes with diffusion MRI: Deterministic or probabilistic tractography? Magnetic Resonance in Medicine, 2020, 83(3), 787-790
Tournier, J.-D.; Calamante, F., Gadian, D.G. & Connelly, A. Direct estimation of the fiber orientation density function from diffusion-weighted MRI data using spherical deconvolution. NeuroImage, 2004, 23, 1176-1185
Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
Tournier, J.-D.; Smith, R. E.; Raffelt, D. A.; Tabbara, R.; Dhollander, T.; Pietsch, M; Christiaens, D.; Jeurissen, B.; Y, C.-H.; Connelly, A. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 2019, 202, 116137
Tustison, N.; Avants, B.; Cook, P.; Zheng, Y.; Egan, A.; Yushkevich, P. & Gee, J. N4ITK: Improved N3 Bias Correction. IEEE Transactions on Medical Imaging, 2010, 29, 1310-1320
Tustison, N.; Avants, B. Explicit B-spline regularization in diffeomorphic image registration. Frontiers in Neuroinformatics, 2013, 7, 39
Tzourio-Mazoyer, N.; Landeau, B.; Papathanassiou, D.; Crivello, F.; Etard, O.; Delcroix, N.; Mazoyer, B. & Joliot, M. Automated Anatomical Labeling of activations in SPM using a Macroscopic Anatomical Parcellation of the MNI MRI single-subject brain. NeuroImage, 15(1), 273–289
Veraart, J.; Novikov, D. S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406
Yeh, C.-H.; Smith, R. E.; Liang, X.; Calamante, F. & Connelly, A. Correction for diffusion MRI fibre tracking biases: The consequences for structural connectomic metrics. NeuroImage, 2016, 142, 150-162
Yeo, B.T.; Krienen, F.M.; Sepulcre, J.; Sabuncu, M.R.; Lashkari, D.; Hollinshead, M.; Roffman, J.L.; Smoller, J.W.; Zollei, L.; Polimeni, J.R.; Fischl, B.; Liu, H. & Buckner, R.L. The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol, 2011, 106(3), 1125-1165
Zalesky, A.; Fornito, A.; Harding, I. H.; Cocchi, L.; Yücel, M.; Pantelis, C. & Bullmore, E. T. Whole-brain anatomical networks: Does the choice of nodes matter? NeuroImage, 2010, 50, 970-983
Zhang, Y.; Brady, M. & Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Transactions on Medical Imaging, 2001, 20, 45-57
```




## Help page

The following help page can equivalently be generated by executing the tool without
providing any command-line arguments (within a container environment; note interface
is slightly different if run natively).

---

## Synopsis

Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in *MRtrix3*

### Usage

    mrtrix3_connectome.py bids_dir output_dir analysis_level [ options ]

-  *bids_dir*: The directory with the input dataset formatted according to the BIDS standard.

-  *output_dir*: The directory where the output files should be stored.

-  *analysis_level*: Level of analysis that will be performed; options are: preproc, participant, group.

### Description

While preproc-level analysis only requires data within the BIDS directory, participant-level analysis requires that the output directory be pre-populated with the results from preproc-level processing; similarly, group-level analysis requires that the output directory be pre-populated with the results from participant-level analysis.

The operations performed by each of the three levels of analysis are as follows:

"preproc": DWI: Denoising; Gibbs ringing removal; motion, eddy current and EPI distortion correction and outlier detection & replacement; brain masking, bias field correction and intensity normalisation; rigid-body registration & transformation to T1-weighted image. T1-weighted image: bias field correction; brain masking.

"participant": DWI: Response function estimation; FOD estimation. T1-weighted image (if -parcellation is not none): Tissue segmentation; grey matter parcellation. Combined (if -parcellation is not none, or -streamlines is provided): Whole-brain streamlines tractography; SIFT2; connectome construction.

"group": Generation of FA-based population template; warping of template-based white matter mask to subject spaces; calculation of group mean white matter response function; scaling of connectomes based on white matter b=0 intensity, response function used during participant-level analysis, and SIFT model proportioinality coefficient; generation of group mean connectome.

The label(s) provided to the -participant_label and -session_label options correspond(s) to sub-<participant_label> and ses-<session_label> from the BIDS spec (so they do _not_ include "sub-" or "ses-"). Multiple participants / sessions can be specified with a space-separated list.

For both preproc-level and participant-level analyses, if no specific participants or sessions are nominated by the user (or the user explicitly specifies multiple participants / sessions), the script will process each of these in series. It is additionally possible for the user to invoke multiple instances of this script in order to process multiple subjects at once in parallel, ensuring that no single participant / session is being processed in parallel, and that preproc-level output data are written fully before commencing participant-level analysis.

The -output_verbosity option principally affects the participant-level analysis, modulating how many derivative files are written to the output directory. Permitted values are from 1 to 4: 1 writes only those files requisite for group-level analysis; 2 additionally writes files typically useful for post-hoc analysis (the default); 3 additionally generates files for enhanced connectome visualisation and copies the entire whole-brain tractogram; 4 additionally generates a full copy of the script scratch directory (with all intermediate files retained) to the output directory (and this applies to all analysis levels)

If running participant-level analysis using the script as a standalone tool rather than inside the provided container, data pertaining to atlas parcellations can no longer be guaranteed to be stored at a specific location on the filesystem. In this case, the user will most likely need to manually specify the location where the corresponding parcellation is stored using the -atlas_path option.

### Options

+ **--output_verbosity**<br>The verbosity of script output (number from 1 to 4).

#### Options that are relevant to participant-level analysis

+ **--parcellation**<br>The choice of connectome parcellation scheme (compulsory for participant-level analysis); options are: aal, aal2, brainnetome246fs, brainnetome246mni, craddock200, craddock400, desikan, destrieux, hcpmmp1, none, perry512, yeo7fs, yeo7mni, yeo17fs, yeo17mni.

+ **--streamlines**<br>The number of streamlines to generate for each subject (will be determined heuristically if not explicitly set).

+ **--template_reg software**<br>The choice of registration software for mapping subject to template space; options are: ants, fsl.

#### Options that are relevant to both preproc-level and participant-level analyses

+ **--t1w_preproc path**<br>Provide a path by which pre-processed T1-weighted image data may be found for the processed participant(s) / session(s)

#### Options specific to the batch processing of participant data

+ **--participant_label**<br>The label(s) of the participant(s) that should be analyzed.

+ **--session_label**<br>The session(s) within each participant that should be analyzed.

#### Standard options

+ **-d/--debug**<br>display debugging messages.

+ **-h/--help**<br>display this information page and exit.

+ **-n/--n_cpus number**<br>use this number of threads in multi-threaded applications (set to 0 to disable multi-threading).

+ **--scratch /path/to/scratch/**<br>manually specify the path in which to generate the scratch directory.

+ **--skip-bids-validator**<br>Skip BIDS validation

+ **-v/--version**<br>display version information and exit.

---

**Author:** Robert E. Smith (robert.smith@florey.edu.au)

**Copyright:** Copyright (c) 2016-2020 The Florey Institute of Neuroscience
and Mental Health.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Covered Software is provided under this License on an "as is"
basis, without warranty of any kind, either expressed, implied, or
statutory, including, without limitation, warranties that the
Covered Software is free of defects, merchantable, fit for a
particular purpose or non-infringing.
See the Mozilla Public License v. 2.0 for more details.

---
