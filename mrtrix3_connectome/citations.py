def add_citations(cmdline, option_prefix='-'):
    cmdline.add_citation(
        'Smith, R. E.; Connelly, A. '
        'MRtrix3_connectome: A BIDS Application for quantitative structural '
        'connectome construction. '
        'In Proc OHBM, 2019, W610',
        is_external=False)
    cmdline.add_citation(
        'Andersson, J. L.; Skare, S. & Ashburner, J. '
        'How to correct susceptibility distortions in spin-echo echo-planar '
        'images: application to diffusion tensor imaging. '
        'NeuroImage, 2003, 20, 870-888',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Andersson, J. L. R.; Jenkinson, M. & Smith, S. '
        'Non-linear registration, aka spatial normalisation. '
        'FMRIB technical report, 2010, TR07JA2',
        condition='If performing participant-level analysis '
                  + f'using {option_prefix}template_reg fsl',
        is_external=True)
    cmdline.add_citation(
        'Andersson, J. L. & Sotiropoulos, S. N. '
        'An integrated approach to correction for off-resonance effects and '
        'subject movement in diffusion MR imaging. '
        'NeuroImage, 2015, 125, 1063-1078',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Andersson, J. L. R.; Graham, M. S.; Zsoldos, E. '
        '& Sotiropoulos, S. N. '
        'Incorporating outlier detection and replacement into a '
        'non-parametric framework for movement and distortion correction of '
        'diffusion MR images. '
        'NeuroImage, 2016, 141, 556-572',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H.; '
        'Filippini, N. & Bastiani, M. '
        'Towards a comprehensive framework for movement and distortion '
        'correction of diffusion MR images: Within volume movement. '
        'NeuroImage, 2017, 152, 450-466',
        condition='If performing preproc-level analysis '
                  'utilising slice-to-volume motion correction',
        is_external=True)
    cmdline.add_citation(
        'Andersson, J. L. R.; Graham,, M. S.; Drobnjak, I.; Zhang, H. & '
        'Campbell, J. '
        'Susceptibility-induced distortion that varies due to motion: '
        'Correction in diffusion MR without acquiring additional data. '
        'NeuroImage, 2018, 171, 277-295',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Avants, B. B.; Epstein, C. L.; Grossman, M.; Gee, J. C. '
        'Symmetric diffeomorphic image registration with cross-correlation: '
        'Evaluating automated labeling of elderly and neurodegenerative brain. '
        'Medical Image Analysis, 2008, 12, 26-41',
        condition='If performing participant-level analysis, '
        f'using {option_prefix}parcellation [ aal, aal2, '
        'brainnetome246mni, craddock200, craddock400, perry512, '
        'yeo7mni or yeo17mni ], '
        f'and not using {option_prefix}template_reg fsl',
        is_external=True)
    cmdline.add_citation(
        'Bhushan, C.; Haldar, J. P.; Choi, S.; Joshi, A. A.; Shattuck, D. W. & '
        'Leahy, R. M. '
        'Co-registration and distortion correction of diffusion and anatomical '
        'images based on inverse contrast normalization. '
        'NeuroImage, 2015, 115, 269-280',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Craddock, R. C.; James, G. A.; Holtzheimer, P. E.; Hu, X. P.; '
        'Mayberg, H. S. '
        'A whole brain fMRI atlas generated via spatially constrained '
        'spectral clustering. '
        'Human Brain Mapping, 2012, 33(8), 1914-1928',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation '
        '[ craddock200 or craddock400 ]',
        is_external=True)
    cmdline.add_citation(
        'Dale, A. M.; Fischl, B. & Sereno, M. I. '
        'Cortical Surface-Based Analysis: '
        'I. Segmentation and Surface Reconstruction. '
        'NeuroImage, 1999, 9, 179-194',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation '
        '[ brainnetome246fs, desikan, destrieux, hcpmmp1, '
        'yeo7fs or yeo17fs ]',
        is_external=True)
    cmdline.add_citation(
        'Desikan, R. S.; Segonne, F.; Fischl, B.; Quinn, B. T.; '
        'Dickerson, B. C.; Blacker, D.; Buckner, R. L.; Dale, A. M.; '
        'Maguire, R. P.; Hyman, B. T.; Albert, M. S. & Killiany, R. J. '
        'An automated labeling system for subdividing the human cerebral '
        'cortex on MRI scans into gyral based regions of interest. '
        'NeuroImage, 2006, 31, 968-980',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation desikan',
        is_external=True)
    cmdline.add_citation(
        'Destrieux, C.; Fischl, B.; Dale, A. & Halgren, E. '
        'Automatic parcellation of human cortical gyri and sulci using '
        'standard anatomical nomenclature. '
        'NeuroImage, 2010, 53, 1-15',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation destrieux',
        is_external=True)
    cmdline.add_citation(
        'Fan, L.; Li, H.; Zhuo, J.; Zhang, Y.; Wang, J.; Chen, L.; Yang, Z.; '
        'Chu, C.; Xie, S.; Laird, A.R.; Fox, P.T.; Eickhoff, S.B.; Yu, C.; '
        'Jiang, T. '
        'The Human Brainnetome Atlas: '
        'A New Brain Atlas Based on Connectional Architecture. '
        'Cerebral Cortex, 2016, 26 (8), 3508-3526',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation '
        '[ brainnetome246fs brainnetome246mni ]',
        is_external=True)
    cmdline.add_citation(
        'Glasser, M. F.; Coalson, T. S.; Robinson, E. C.; Hacker, C. D.; '
        'Harwell, J.; Yacoub, E.; Ugurbil, K.; Andersson, J.; Beckmann, C. F.; '
        'Jenkinson, M.; Smith, S. M. & Van Essen, D. C. '
        'A multi-modal parcellation of human cerebral cortex. '
        'Nature, 2016, 536, 171-178',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation hcpmmp1',
        is_external=True)
    cmdline.add_citation(
        'Iglesias, J. E.; Liu, C. Y.; Thompson, P. M. & Tu, Z. '
        'Robust Brain Extraction Across Datasets and '
        'Comparison With Publicly Available Methods. '
        'IEEE Transactions on Medical Imaging, 2011, 30, 1617-1634',
        condition='If performing either preproc-level or participant-level '
        'analysis and ROBEX is used for T1w brain extraction',
        is_external=True)
    cmdline.add_citation(
        'Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. '
        'Multi-tissue constrained spherical deconvolution for improved '
        'analysis of multi-shell diffusion MRI data. '
        'NeuroImage, 2014, 103, 411-426',
        condition='If performing either preproc-level '
        'or participant-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Kellner, E.; Dhital, B.; Kiselev, V. G.; Reisert, M. '
        'Gibbs-ringing artifact removal based on local subvoxel-shifts. '
        'Magnetic Resonance in Medicine, 2006, 76(5), 1574-1581',
        condition='If performing preproc-level analysis',
        is_external=True)
    cmdline.add_citation(
        'Patenaude, B.; Smith, S. M.; Kennedy, D. N. & Jenkinson, M. '
        'A Bayesian model of shape and appearance for '
        'subcortical brain segmentation. '
        'NeuroImage, 2011, 56, 907-922',
        condition='If performing participant-level analysis and '
        f'not using {option_prefix}parcellation '
        '[ brainnetome246fs or brainnetome246mni ]',
        is_external=True)
    cmdline.add_citation(
        'Perry, A.; Wen, W.; Kochan, N. A.; Thalamuthu, A.; Sachdev, P. S.; '
        'Breakspear, M. '
        'The independent influences of age and education on functional brain '
        'networks and cognition in healthy older adults. '
        'Human Brain Mapping, 2017, 38(10), 5094-5114',
        condition='If performing participant-level analysis and '
        f'using {option_prefix}parcellation perry512',
        is_external=True)
    cmdline.add_citation(
        'Smith, S. M. '
        'Fast robust automated brain extraction. '
        'Human Brain Mapping, 2002, 17, 143-155',
        condition='If performing either preproc-level or participant-level '
        'analysis and FSL is used for T1w brain extraction',
        is_external=True)
    cmdline.add_citation(
        'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. '
        'Anatomically-constrained tractography: '
        'Improved diffusion MRI streamlines tractography through '
        'effective use of anatomical information. '
        'NeuroImage, 2012, 62, 1924-1938',
        condition='If performing participant-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. '
        'The effects of SIFT on the reproducibility and biological '
        'accuracy of the structural connectome. '
        'NeuroImage, 2015a, 104, 253-265',
        condition='If performing participant-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. '
        'SIFT2: Enabling dense quantitative assessment of brain white matter '
        'connectivity using streamlines tractography. '
        'NeuroImage, 2015b, 119, 338-351',
        condition='If performing participant-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Smith, R. E.; Raffelt, D.; Tournier, J.-D.; Connelly, A. '
        'Quantitative streamlines tractography: Methods and inter-subject '
        'normalisation. '
        'OHBM Aperture, 2022, doi: 10.52294/ApertureNeuro.2022.2.NEOD9565',
        condition='If performing group-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Tournier, J.-D.; Calamante, F., Gadian, D.G. & Connelly, A. '
        'Direct estimation of the fiber orientation density function from '
        'diffusion-weighted MRI data using spherical deconvolution. '
        'NeuroImage, 2004, 23, 1176-1185',
        condition='If performing either preproc-level or participant-level '
        'analysis',
        is_external=False)
    cmdline.add_citation(
        'Tournier, J.-D.; Calamante, F. & Connelly, A. '
        'Improved probabilistic streamlines tractography by 2nd order '
        'integration over fibre orientation distributions. '
        'Proceedings of the International Society for Magnetic '
        'Resonance in Medicine, 2010, 1670',
        condition='If performing participant-level analysis '
        f'and {option_prefix}streamlines 0 is not set',
        is_external=False)
    cmdline.add_citation(
        'Tustison, N.; Avants, B.; Cook, P.; Zheng, Y.; Egan, A.; '
        'Yushkevich, P. & Gee, J. '
        'N4ITK: Improved N3 Bias Correction. '
        'IEEE Transactions on Medical Imaging, 2010, 29, 1310-1320',
        condition='If performing either preproc-level or participant-level '
        'analysis, and N4 is used for either DWI or T1w bias field correction',
        is_external=True)
    cmdline.add_citation(
        'Tustison, N.; Avants, B. '
        'Explicit B-spline regularization in diffeomorphic image '
        'registration. '
        'Frontiers in Neuroinformatics, 2013, 7, 39',
        condition='If performing participant-level analysis, '
        f'using {option_prefix}parcellation [ aal, aal2, '
        'brainnetome246mni, craddock200, craddock400, perry512, '
        'yeo7mni or yeo17mni ], '
        f'and not using {option_prefix}template_reg fsl',
        is_external=True)
    cmdline.add_citation(
        'Tzourio-Mazoyer, N.; Landeau, B.; Papathanassiou, D.; Crivello, F.; '
        'Etard, O.; Delcroix, N.; Mazoyer, B. & Joliot, M. '
        'Automated Anatomical Labeling of activations in SPM using a '
        'Macroscopic Anatomical Parcellation of the MNI MRI '
        'single-subject brain. '
        'NeuroImage, 15(1), 273-289',
        condition='If performing participant-level analysis '
        f'and using {option_prefix}parcellation [ aal or aal2 ]',
        is_external=True)
    cmdline.add_citation(
        'Veraart, J.; Novikov, D. S.; Christiaens, D.; Ades-aron, B.; '
        'Sijbers, J. & Fieremans, E. '
        'Denoising of diffusion MRI using random matrix theory. '
        'NeuroImage, 2016, 142, 394-406',
        condition='If performing preproc-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Yeh, C.H.; Smith, R.E.; Liang, X.; Calamante, F.; Connelly, A. '
        'Correction for diffusion MRI fibre tracking biases: '
        'The consequences for structural connectomic metrics. '
        'Neuroimage, 2016, 142, 150-162',
        condition='If utilising connectome outputs from '
        'participant-level analysis',
        is_external=False)
    cmdline.add_citation(
        'Yeo, B.T.; Krienen, F.M.; Sepulcre, J.; Sabuncu, M.R.; Lashkari, D.; '
        'Hollinshead, M.; Roffman, J.L.; Smoller, J.W.; Zollei, L.; '
        'Polimeni, J.R.; Fischl, B.; Liu, H. & Buckner, R.L. '
        'The organization of the human cerebral cortex estimated by '
        'intrinsic functional connectivity. '
        'J Neurophysiol, 2011, 106(3), 1125-1165',
        condition='If performing participant-level analysis '
        f'and using {option_prefix}parcellation '
        '[ yeo7fs, yeo7mni, yeo17fs or yeo17mni ]',
        is_external=False)
    cmdline.add_citation(
        'Zalesky, A.; Fornito, A.; Harding, I. H.; Cocchi, L.; Yucel, M.; '
        'Pantelis, C. & Bullmore, E. T. '
        'Whole-brain anatomical networks: Does the choice of nodes matter? '
        'NeuroImage, 2010, 50, 970-983',
        condition='If performing participant-level analysis '
        'and using {option_prefix}parcellation perry512',
        is_external=True)
    cmdline.add_citation(
        'Zhang, Y.; Brady, M. & Smith, S. '
        'Segmentation of brain MR images through a hidden Markov random '
        'field model and the expectation-maximization algorithm. '
        'IEEE Transactions on Medical Imaging, 2001, 20, 45-57',
        condition='If performing participant-level analysis',
        is_external=True)
