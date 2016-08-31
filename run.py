#!/usr/bin/env python

import os, glob, sys
import lib.app, lib.cmdlineParser

from lib.binaryInPath  import binaryInPath
from lib.delFile       import delFile
from lib.errorMessage  import errorMessage
from lib.getFSLSuffix  import getFSLSuffix
from lib.getHeaderInfo import getHeaderInfo
from lib.getUserPath   import getUserPath
from lib.imagesMatch   import imagesMatch
from lib.isWindows     import isWindows
from lib.runCommand    import runCommand


def runSubject (bids_dir, label, output_prefix):
  import lib.app

  output_dir = os.path.join(output_prefix, label);
  if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
  os.makedirs(output_dir)
  os.makedirs(os.path.join(output_dir, 'connectome'))
  os.makedirs(os.path.join(output_dir, 'dwi'))

  fsl_path = os.environ.get('FSLDIR', '')
  if not fsl_path:
    errorMessage('Environment variable FSLDIR is not set; please run appropriate FSL configuration script')

  bet_cmd = 'bet'
  if not binaryInPath(bet_cmd):
    bet_cmd = 'fsl5.0-bet'
    if not binaryInPath(bet_cmd):
      errorMessage('Could not find FSL program bet; please verify FSL install')

  flirt_cmd = 'flirt'
  if not binaryInPath(flirt_cmd):
    flirt_cmd = 'fsl5.0-flirt'
    if not binaryInPath(flirt_cmd):
      errorMessage('Could not find FSL program flirt; please verify FSL install')

  ssroi_cmd = 'standard_space_roi'
  if not binaryInPath(ssroi_cmd):
    ssroi_cmd = 'fsl5.0-standard_space_roi'
    if not binaryInPath(ssroi_cmd):
      errorMessage('Could not find FSL program standard_space_roi; please verify FSL install')

  fsl_suffix = getFSLSuffix()

  if lib.app.args.parc == 'fs_2005' or lib.app.args.parc == 'fs_2009':
    if not os.environ['FREESURFER_HOME']:
      errorMessage('Environment variable FREESURFER_HOME not set; please verify FreeSurfer installation')
    if not binaryInPath('recon_all'):
      errorMessage('Could not find FreeSurfer script recon_all; please verify FreeSurfer installation')

  lib.app.makeTempDir()

  # Need to perform an initial import of JSON data using mrconvert; so let's grab the diffusion gradient table as well
  # If no bvec/bval present, need to go down the directory listing
  # Only try to import JSON file if it's actually present
  grad_prefix = os.path.join(bids_dir, label, 'dwi', label + '_dwi'
  if os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec'):
    grad_prefix = os.path.join(bids_dir, 'dwi')
    if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
      errorMessage('Unable to locate valid diffusion gradient table');
  grad_import_option = ' -fslgrad ' + grad_prefix + '.bvec ' + grad_prefix + '.bval'
  json_path = os.path.join(bids_dir, label, 'dwi', label + '_dwi.json')
  if os.path.isfile(json_path):
    json_import_option = ' -json_import ' + json_path
  else:
    json_import_option = ''
  runCommand('mrconvert ' + os.path.join(bids_dir, label, 'dwi', label + '_dwi.nii.gz')
             + grad_import_option + json_import_option
             + ' ' + os.path.join(lib.app.tempDir, 'input.mif'))

  # Go hunting for reversed phase-encode data
  # TODO Should ideally have compatibility with fieldmap data also
  fmap_dir = os.path.join(bids_dir, label, 'fmap')
  if not os.path.isdir(fmap_dir):
    errorMessage('Subject does not possess fmap data necessary for EPI distortion correction')
  fmap_index = 0;
  while (1):
    prefix = os.path.join(fmap_dir, label + '_dir-' + str(fmap_index))
    if os.path.isfile(prefix + '.nii.gz') and os.path.isfile(prefix + '.json'):
      runCommand('mrconvert ' + prefix + '.nii.gz -json_import ' + prefix + '.json ' + os.path.join(lib.app.tempDir, 'RPE' + str(fmap_index) + '.mif'))
    else:
      break
    fmap_index += 1

  runCommand('mrconvert ' + os.path.join(bids_dir, label, 'anat', label + '_T1w.nii.gz') + ' ' + os.path.join(lib.app.tempDir, 'T1.mif'))

  cwd = os.cwd()
  lib.app.gotoTempDir()

  # Step 1: Denoise
  runCommand('dwidenoise input.mif dwi_denoised.mif')
  delFile('input.mif')

  # Step 2: Distortion correction
  # TODO Need to assess presence of reversed phase-encoding data and act accordingly

  if fmap_index < 2:
    errorMessage('Inadequate number of images in fmap directory for inhomogeneity estimation')
  dwi_pe = getHeaderProperty('dwi_denoised.mif', 'PhaseEncodingDirection')
  if not dwi_pe:
    errorMessage('Phase encoding direction of DWI not defined')
  fmap_pe = [ ]
  for index in range(0, fmap_index):
    pe = getHeaderProperty('RPE' + str(index) + 'mif', 'PhaseEncodingDirection')
    if not pe:
      errorMessage('Field mapping images do not all contain phase encoding information')
    pe = getPEDir(pe)
    if not getPEDir(dwi_pe)[0] == pe[0]:
      errorMessage('Non-collinear phase encoding directions not currently supported')
    fmap_pe.append(pe)
  # Need to detect 'RPE*' volumes with equivalent phase-encoding directions, and concatenate them
  # Either that, or start working on the dwipreproc interface...

  runCommand('dwipreproc dwi_denoised.mif dwi_denoised_preprocessed.mif')
  delFile('dwi_denoised.mif')

  # Step 3: Bias field correction
  runCommand('dwibiascorrect dwi_denoised_preprocessed.mif dwi.mif')
  delFile('dwi_denoised_preprocessed.mif')

  # Step 4: Generate a brain mask for DWI
  runCommand('dwi2mask dwi.mif dwi_mask.mif')

  # Step 5: Perform brain extraction on the T1 image in its original space
  #         (this is necessary for histogram matching prior to registration)
  runCommand('mrconvert T1.mif T1.nii -stride -1,+2,+3')
  mni_mask_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_1mm_brain_mask_dil.nii.gz')
  mni_mask_dilation = 0;
  if os.path.exists (mni_mask_path):
    mni_mask_dilation = 4;
  else:
    mni_mask_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_2mm_brain_mask_dil.nii.gz')
    if os.path.exists (mni_mask_path):
      mni_mask_dilation = 2;
  if mni_mask_dilation:
    runCommand('maskfilter ' + mni_mask_path + ' dilate mni_mask.nii -npass ' + str(mni_mask_dilation))
    runCommand(ssroi_cmd + ' T1.nii T1_preBET' + fsl_suffix + ' -maskMASK mni_mask.nii -roiNONE', False)
  else:
    runCommand(ssroi_cmd + ' T1.nii T1_preBET' + fsl_suffix + ' -b', False)
  if not os.path.exists('T1_preBET' + fsl_suffix):
    warnMessage('FSL command ' + ssroi_cmd + ' appears to have failed; passing T1 directly to BET')
    runCommand('mrconvert T1.nii T1_preBET' + fsl_suffix + ' -stride -1,+2,+3')
  delFile('T1.nii')
  runCommand(bet_cmd + ' T1_preBET' + fsl_suffix + ' T1_BET' + fsl_suffix + ' -f 0.15 -R')

  # Step 6: Generate target image for T1->DWI registration
  runCommand('mrcalc 1 dwi_denoise_preproc.mif -div dwi_mask.mif -mult - | mrhistmatch - T1_BET' + fsl_suffix + ' dwi_pseudoT1.nii -stride -1,+2,+3')

  # Step 7: Perform T1->DWI registration
  #         Since mrregister is currently symmetric, but here we explicitly want an asymmetric
  #         registration (so we don't have to worry about gradient direction reorientation),
  #         for now we'll go with FSL's FLIRT
  #         TODO Switch to least-squares metric, or switch to mrregister
  runCommand('flirt -ref dwi_pseudoT1 -in T1_BET -omat T1_to_DWI_FLIRT.mat -dof 6')
  runCommand('transformconvert T1_to_DWI_FLIRT.mat flirt_import T1_to_DWI_MRtrix.mat')
  delFile('T1_to_DWI_FLIRT.mat')
  runCommand('mrtransform T1.mif T1_registered.mif -linear T1_to_DWI_MRtrix.mat')
  delFile('T1.mif')
  delFile('T1_to_DWI_MRtrix.mat')

  # Step 8: Generate 5TT image for ACT
  runCommand('5ttgen fsl T1_registered.mif 5TT.mif')

  # Step 9: Determine whether we are working with single-shell or multi-shell data
  shells = [ int(round(float(x))) for x in getHeaderInfo('dwi_denoise_preproc.mif', 'shells').split() ]
  multishell = (len(shells) > 2)

  # Step 10: Estimate response function(s) for spherical deconvolution
  if multishell:
    runCommand('dwi2response msmt_5tt dwi_denoise_preproc.mif response_wm.txt response_gm.txt response_csf.txt -mask dwi_mask.mif')
    rf_file_for_scaling = 'response_wm.txt'
  else:
    runCommand('dwi2response tournier dwi_denoise_preproc.mif response.txt -mask dwi_mask.mif')
    rf_file_for_scaling = 'response.txt'

  # Step 12: Perform spherical deconvolution
  #          Use a dilated mask for spherical deconvolution as a 'safety margin' -
  #          ACT should be responsible for stopping streamlines before they reach the edge of the DWI mask
  runCommand('maskfilter dwi_mask.mif dilate dwi_mask_dilated.mif -npass 3')
  if multishell:
    runCommand('dwi2fod msmt_csd dwi_denoise_preproc.mif response_wm.txt FOD_WM.mif response_gm.txt FOD_GM.mif response_csf.txt FOD_CSF.mif -mask dwi_mask_dilated.mif')
  else:
    # Still use the msmt_csd algorithm with single-shell data: Use hard non-negativity constraint
    runCommand('dwi2fod msmt_csd dwi_denoise_preproc.mif response.txt FOD_WM.mif')

  # Step 13: Generate the tractogram
  # TODO Determine the appropriate number of streamlines based on the number of nodes in the parcellation
  # TODO Present this as a command-line option
  runCommand('tckgen FOD_WM.mif tractogram.tck -act 5TT.mif -backtrack -crop_at_gmwmi -cutoff 0.06 -maxlength 250 -number 10M -seed_dynamic FOD_WM.mif')

  # Step 14: Use SIFT2 to determine streamline weights
  fd_scale_gm_option = ''
  if not multishell:
    fd_scale_gm_option = ' -fd_scale_gm'
  runCommand('tcksift2 tractogram.tck FOD_WM.mif weights.csv -act 5TT.mif -out_mu mu.txt' + fod_scale_gm_option + ' -info')

  # Step 15: Generate the grey matter parcellation
  #          The necessary steps here will vary significantly depending on the parcellation scheme selected
  if lib.app.args.parc == 'fs_2005' or lib.app.args.parc == 'fs_2009':

    # Run FreeSurfer pipeline on this subject's T1 image
    # TODO May need to change the FreeSurfer subjects directory in order to have write access in Singularity
    runCommand('T1_registered.mif T1_registered.nii')
    runCommand('recon_all -subjid temp -i T1_registered.nii')
    delFile('T1_registered.nii')
    runCommand('recon_all -subjid temp -all')

    # Grab the relevant parcellation image and target lookup table for conversion
    parc_image = os.path.join(os.environ['FREESURFER_HOME'], 'subjects', 'temp', 'mri')
    lut_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src', 'connectome', 'tables')
    if lib.app.args.parc == 'fs_2005':
      parc_image = os.path.join(parc_image, 'aparc.aseg.mgz')
      lut_file = os.path.join(lut_file, 'fs_default.txt')
    else:
      parc_image = os.path.join(parc_image, 'aparc.a2009.aseg.mgz')
      lut_file = os.path.join(lut_file, 'fs_2009.txt')

    # Perform the index conversion
    runCommand('labelconvert ' + parc_image + ' ' + os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt') + ' ' + lut_file + ' parc_init.mif')
    shutil.rmtree(os.path.join(os.environ['FREESURFER_HOME'], 'subjects', 'temp'))

    # Fix the sub-cortical grey matter parcellations using FSL FIRST
    runCommand('labelsgmfix parc_init.mif T1_registered.mif ' + label_file + ' parc.mif')
    delFile('parc_init.mif')

  # TODO Implement AAL parcellation: Need MNI single-subject T1 image to perform registration (doesn't come with AAL package)

  else:
    errorMessage('Unknown parcellation scheme requested: ' + lib.app.args.parc)

  # Step 16: Generate the connectome
  #          Only provide the standard density-weighted connectome for now
  runCommand('tck2connectome tractogram.tck parc.mif connectome.csv -tck_weights_in weights.csv')
  delFile('weights.csv')

  # Move necessary files to output directory
  shutil.copy('connectome.csv', os.path.join(output_dir, 'connectome', label + '_connectome.csv'))
  runCommand('mrconvert dwi.mif ' + os.path.join(output_dir, 'dwi', label + '_dwi.nii.gz')
             + ' -export_grad_fsl ' + os.path.join(output_dir, 'dwi', label + '_dwi.bvec') + ' ' + os.path.join(output_dir, 'dwi', label + '_dwi.bval')
             + ' -json_export ' + os.path.join(output_dir, 'dwi', label + '_dwi.json'))
  shutil.copy('out_mu.txt', os.path.join(output_dir, 'connectome', label + '_mu.txt'))
  shutil.copy(rf_file_for_scaling, os.path.join(output_dir, 'dwi', label + '_response.txt'))

  # Manually wipe and zero the temp directory (since we might be processing more than one subject)
  os.chdir(cwd)
  shutil.rmtree(lib.app.tempDir)
  lib.app.tempDir = ''

# End of runSubject() function



analysis_choices = [ 'participant', 'group' ]
parcellation_choices = [ 'fs_2005', 'fs_2009' ]

lib.app.author = 'Robert E. Smith (robert.smith@florey.edu.au)'
lib.cmdlineParser.initialise('Generate subject connectomes from raw image data, perform inter-subject connection density normalisation, and perform statistical inference across subjects')
lib.app.parser.add_argument('bids_dir', help='The directory with the input dataset formatted according to the BIDS standard.')
lib.app.parser.add_argument('output_dir', help='The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.')
lib.app.parser.add_argument('analysis_level', help='Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: ' + ', '.join(analysis_choices), choices=analysis_choices)
batch_options = lib.app.parser.add_argument_group('Options specific to the batch processing of subject data')
batch_options.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a comma-separated list.')
connectome_options = lib.app.parser.add_argument_group('Options for setting up the connectome reconstruction')
connectome_options.add_argument('-parc', help='The choice of connectome parcellation scheme. Options are: ' + ', '.join(parcellation_choices), choices=parcellation_choices)
testing_options = lib.app.parser.add_argument_group('Options for testing the run.py script')
lib.app.initialise()

if isWindows():
  errorMessage('Script cannot be run on Windows due to FSL dependency')

subjects_to_analyze = [ ]
# Only run a subset of subjects
if lib.app.args.participant_label:
  index_list = lib.app.args.participant_label.split(' ')
  # TODO May need zero-padding
  subjects_to_analyze = [ 'sub-' + i for i in index_list ]
  for dir in subjects_to_analyze:
    if not os.path.isdir(os.path.join(lib.app.args.bids_dir, dir)):
      print (os.cwd)
      print (lib.app.args.bids_dir)
      print (dir)
      errorMessage('Unable to find directory for subject: ' + dir)
# Run all subjects sequentially
else:
  subject_dirs = glob.glob(os.path.join(lib.app.args.bids_dir, 'sub-*'))
  subjects_to_analyze = subject_dirs
  # subjects_to_analyze = [ dir.split("-")[-1] for dir in subject_dirs ]


# Running participant level
if lib.app.args.analysis_level == "participant":

  for subject_label in subjects_to_analyze:
    runSubject(lib.app.args.bids_dir, subject_label, lib.app.args.output_dir)

# Running group level
elif lib.app.args.analysis_level == "group":

  if lib.app.args.participant_label:
    errorMessage('Do not specify a subgroup of subjects if performing a group analysis')

  pop_template_dir = os.path.join(lib.app.args.output_dir, 'population_template')
  os.makedirs(pop_template_dir)
  os.makedirs(os.path.join(pop_template_dir, 'images'))
  os.makedirs(os.path.join(pop_template_dir, 'masks'))
  os.makedirs(os.path.join(pop_template_dir, 'values'))
  os.makedirs(os.path.join(pop_template_dir, 'warps'))

  # First pass through subject data in group analysis:
  #   - Grab DWI data (written back from single-subject analysis back into BIDS format)
  #   - Generate mask and FA images to be used in populate template generation
  #   - Generate the mean b=0 image that will be used later
  for subject_label in subjects_to_analyze:
    dwi_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_dwi.nii.gz')
    if not os.path.exists(dwi_path):
      errorMessage('Unable to find subject DWI data: ' + dwi_path)
    bvec_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_dwi.bvec')
    bval_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_dwi.bval')
    if not os.path.exists(bvec_path) or not os.path.exists(bval_path):
      errorMessage('Unable to find DWI gradient table for subject: ' + subject_label)
    json_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_dwi.json')
    if not os.path.exists(dwi_path):
      errorMessage('Unable to find DWI JSON file: ' + json_path)
    mask_path = os.path.join(lib.app.args.output_dir, 'population_template', 'masks', subject_label + '.mif')
    grad_import_option = ' -fslgrad ' + bvec_path + ' ' + bval_path
    runCommand('dwi2mask ' + dwi_path + ' ' + mask_path + grad_import_option)
    tensor_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_tensor.mif')
    runCommand('dwi2tensor ' + dwi_path + ' ' + tensor_path + ' -mask ' + mask_path + grad_import_option)
    fa_path = os.path.join(lib.app.args.output_dir, 'population_template', 'images', subject_label + '.mif')
    runCommand('tensor2metric ' + tensor_path + ' ' + fa_path)
    deflFile(tensor_path)
    bzeros_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_bzeros.mif')
    runCommand('dwiextract ' + dwi_path + ' ' + bzeros_path + ' -bzero')
    runCommand('mrmath ' + bzeros_path + ' ' + os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_mean_bzero.mif') + ' -axis 3')
    delFile(bzeros_path)

  # First group-level calculation: Generate the population template
  template_path = os.path.join(pop_template_dir, 'template.mif')
  runCommand('population_template ' + os.path.join(pop_template_dir, 'images') + ' -mask_dir ' + os.path.join(pop_template_dir, 'masks') + ' -warp_dir ' + os.path.join(pop_template_dir, 'warps') + ' ' + template_path + ' -linear_scale 0.25,0.5,1.0,1.0 -nl_scale 0.5,0.75,1.0,1.0,1.0 -nl_niter 5,5,5,5,5')

  # Second pass through subject data in group analysis:
  #   - Warp template FA image back to subject space & threshold to define a WM mask
  #   - Calculate the median subject FA value within this mask
  #   - Store this in a file, and contribute to calculation of the mean of these values across subjects
  #   - Contribute to the group average response function
  mean_median_bzero = 0.0
  mean_RF = [ ]
  for subject_label in subjects_to_analyze:
    warped_template_path = os.path.join(lib.app.args.args.output_dir, subject_label, 'dwi', subject_label + '_template_fa.mif')
    runCommand('mrtransform ' + template_path + ' -warp_full ' + os.path.join(pop_template_dir, 'warps', subject_label + '.mif') + ' ' + warped_template_path + ' -from 2 -template ' + os.path.join(args.output_dir, 'population_template', 'images', subject_label + '.mif'))
    voxel_mask_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_intensity_mask.mif')
    runCommand('mrthreshold ' + warped_template_path + ' ' + voxel_mask_path + ' -abs 0.4')
    delFile(warped_template_path)
    mean_bzero_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_mean_bzero.mif')
    median_bzero = getImageStat(mean_bzero_path, 'median', voxel_mask_path)
    delFile(voxel_mask_path)
    delFile(mean_bzero_path)
    with open(os.path.join(pop_template_dir, 'values', subject_label + '.txt'), 'w') as f:
      f.write (median_bzero)
    mean_median_bzero = mean_median_bzero + float(median_bzero)
    rf_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_response.txt')
    if not os.path.exists(rf_path):
      errorMessage('Unable to find SD response function file: ' + rf_path)
    RF = [ ]
    with open(rf_path, 'r') as f:
      RF.append([ float(v) for v in f.read().split() ])
    RF_lzero = [ line[0] for line in RF ]
    if mean_RF:
      mean_RF = mean_RF + RF_lzero
    else:
      mean_RF = RF_lzero

  # Second group-level calculation:
  #   - Calculate the mean of median b=0 values
  #   - Calculate the mean response function
  mean_median_bzero = mean_median_bzero / len(subjects_to_analyze)
  mean_RF = [ v / len(subjects_to_analyze) for v in mean_RF ]

  # Third pass through subject data in group analysis:
  #   - Scale the connectome strengths:
  #     - Multiply by SIFT proportionality coefficient mu
  #     - Multiply by (mean median b=0) / (subject median b=0)
  #     - Multiply by (subject RF size) / (mean RF size)
  for subject_label in subjects_to_analyze:
    mu_file_path = os.path.join(lib.app.args.output_dir, subject_label, 'connectome', subject_label + '_mu.txt')
    if not os.path.exists(mu_file_path):
      errorMessage('Could not find SIFT proportionality coefficient file: ' + mu_file_path)
    with open(mu_file_path, 'r') as f:
      mu = float(f.read())
    with open(os.path.join(pop_template_dir, 'values', subject_label + '.txt'), 'r') as f:
      median_bzero = float(f.read())
    rf_path = os.path.join(lib.app.args.output_dir, subject_label, 'dwi', subject_label + '_response.txt')
    RF = [ ]
    with open(rf_path, 'r') as f:
      RF.append([ float(v) for v in f.read().split() ])
    RF_lzero = [ line[0] for line in RF ]
    RF_multiplier = 1.0
    for (mean, subj) in zip(mean_RF, RF_lzero):
      RF_multiplier = RF_multiplier * subj / mean
    global_multiplier = mu * (mean_median_bzero / median_bzero) * RF_multiplier

    connectome_path = os.path.join(lib.app.args.output_dir, subject_label, 'connectome', subject_label + '_connectome.csv')
    if not os.path.exists(connectome_path):
      errorMessage('Could not find subject connectome file: ' + connectome_path)
    connectome = [ ]
    with open(connectome_path, 'r') as f:
      connectome.append ( [ float(v) for v in f.read().split() ] )
    with open(os.path.join(lib.app.args.output_dir, subject_label, 'connectome', subject_label + '_connectome_scaled.csv'), 'w') as f:
      for line in connectome:
        f.write( ','.join([ v*global_multiplier for v in line ]) )

  # Third group-level calculation: Perform statistical analysis of connectomes using NBS-TFCE
  # Can't do this without opening the root-directory dataset description file and selecting a contrast to test

  # For now, for the sake of testing the pipeline and doing some form of group data manipulation,
  #   let's just generate the mean connectome
  mean_connectome = [ ]
  for subject_label in subjects_to_analyze:
    path = os.path.join(lib.app.args.output_dir, subject_label, 'connectome', subject_label + '_connectome_scaled.csv')
    connectome = [ ]
    with open(path, 'r') as f:
      connectome.append( [ float(v) for v in f.read().split() ] )
    if mean_connectome:
      for r1,r2 in zip(connectome, mean_connectome):
        r2 = [ c1+c2 for c1,c2 in zip(r1,r2) ]
    else:
      mean_connectome = connectome

  mean_connectome = [ [ v/len(subjects_to_analyze) for v in row ] for row in mean_connectome ]

  with open(os.path.join(pop_template_dir, 'mean_connectome.csv'), 'w') as f:
    for row in mean_connectome:
      f.write(','.join(row))
