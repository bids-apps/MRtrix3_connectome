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
from lib.printMessage  import printMessage
from lib.runCommand    import runCommand
from lib.warnMessage   import warnMessage


__version__ = 'BIDS-App \'MRtrix3_connectome\' version {}'.format(open('/version').read()) if os.path.exists('/version') else 'BIDS-App \'MRtrix3_connectome\' standalone'


def runSubject (bids_dir, label, output_prefix):
  import lib.app, os, shutil
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

  if not binaryInPath('N4BiasFieldCorrection'):
    errorMessage('Could not find ANTs program N4BiasFieldCorrection; please verify ANTs installation')

  if not lib.app.args.parc:
    errorMessage('For participant-level analysis, desired parcellation must be provided using the -parc option')

  parc_image_path = ''
  parc_lut_file = ''
  mrtrix_lut_file = os.path.join(os.path.dirname(os.path.abspath(lib.__file__)), os.pardir, os.pardir, 'src', 'connectome', 'tables')

  if lib.app.args.parc == 'fs_2005' or lib.app.args.parc == 'fs_2009':
    if not os.environ['FREESURFER_HOME']:
      errorMessage('Environment variable FREESURFER_HOME not set; please verify FreeSurfer installation')
    if not binaryInPath('recon-all'):
      errorMessage('Could not find FreeSurfer script recon-all; please verify FreeSurfer installation')
    parc_lut_file = os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')
    if lib.app.args.parc == 'fs_2005':
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'fs_default.txt')
    else:
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'fs_2009.txt')

# TODO Command-line option to provide the path to the parcellation information,
#   for if the script is used outside of the BIDS App container and these data
#   are provided in some other location

  if lib.app.args.parc == 'aal' or lib.app.args.parc == 'aal2':
    mni152_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_1mm.nii.gz')
    if not os.path.isfile(mni152_path):
      errorMessage('Could not find MNI152 template image within FSL installation')
    if lib.app.args.parc == 'aal':
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'aal.txt')
    else:
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'aal2.txt')

  if parc_image_path and not os.path.isfile(parc_image_path):
    errorMessage('Could not find parcellation image (expected location: ' + parc_image_path + ')')
  if not os.path.isfile(parc_lut_file):
    errorMessage('Could not find parcellation lookup table file (expected location: ' + parc_lut_file + ')')
  if not os.path.exists(mrtrix_lut_file):
    errorMessage('Could not find MRtrix3 connectome lookup table file (expected location: ' + mrtrix_lut_file + ')')

  lib.app.makeTempDir()

  # Need to perform an initial import of JSON data using mrconvert; so let's grab the diffusion gradient table as well
  # If no bvec/bval present, need to go down the directory listing
  # Only try to import JSON file if it's actually present
  # TODO May need to concatenate more than one input DWI, since if there's more than one phase-encode
  #   direction in the acquisition they'll need to be split across multiple files
  grad_prefix = os.path.join(bids_dir, label, 'dwi', label + '_dwi')
  if not (os.path.isfile(grad_prefix + '.bval') and os.path.isfile(grad_prefix + '.bvec')):
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
  # TODO If there's an 'IntendedFor' field in the JSON file, and it's NOT the DWI(s) we're using, don't use
  fmap_dir = os.path.join(bids_dir, label, 'fmap')
  if not os.path.isdir(fmap_dir):
    errorMessage('Subject does not possess fmap data necessary for EPI distortion correction')
  fmap_index = 1;
  fmap_image_list = [ ]
  while (1):
    prefix = os.path.join(fmap_dir, label + '_dir-' + str(fmap_index) + '_epi')
    if os.path.isfile(prefix + '.nii.gz') and os.path.isfile(prefix + '.json'):
      runCommand('mrconvert ' + prefix + '.nii.gz -json_import ' + prefix + '.json ' + os.path.join(lib.app.tempDir, 'RPE' + str(fmap_index) + '.mif'))
      fmap_image_list.append('RPE' + str(fmap_index) + '.mif')
    else:
      break
    fmap_index += 1

  runCommand('mrconvert ' + os.path.join(bids_dir, label, 'anat', label + '_T1w.nii.gz') + ' ' + os.path.join(lib.app.tempDir, 'T1.mif'))

  cwd = os.getcwd()
  lib.app.gotoTempDir()

  # Step 1: Denoise
  runCommand('dwidenoise input.mif dwi_denoised.mif')
  delFile('input.mif')

  # Step 2: Distortion correction
  runCommand('mrcat ' + ' '.join(fmap_image_list) + ' fmap_images.mif -axis 3')
  for path in fmap_image_list:
    delFile(path)
  runCommand('dwipreproc dwi_denoised.mif dwi_denoised_preprocessed.mif -rpe_header -topup_images fmap_images.mif')
  delFile('dwi_denoised.mif')
  delFile('fmap_images.mif')

  # Step 3: Bias field correction
  runCommand('dwibiascorrect dwi_denoised_preprocessed.mif dwi.mif -ants')
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
  delFile('T1_preBET' + fsl_suffix)

  # Step 6: Generate target image for T1->DWI registration
  runCommand('dwiextract dwi.mif -bzero - | mrmath - mean - -axis 3 | mrcalc 1 - -div dwi_mask.mif -mult - | mrconvert - - -stride -1,+2,+3 | mrhistmatch - T1_BET' + fsl_suffix + ' dwi_pseudoT1.nii')

  # Step 7: Perform T1->DWI registration
  #         Since mrregister is currently symmetric, but here we explicitly want an asymmetric
  #         registration (so we don't have to worry about gradient direction reorientation),
  #         for now we'll go with FSL's FLIRT
  #         TODO Switch to least-squares metric, or switch to mrregister
  runCommand('flirt -ref dwi_pseudoT1 -in T1_BET -omat T1_to_DWI_FLIRT.mat -dof 6')
  runCommand('transformconvert T1_to_DWI_FLIRT.mat T1_BET' + fsl_suffix + ' dwi_pseudoT1.nii flirt_import T1_to_DWI_MRtrix.mat')
  delFile('T1_to_DWI_FLIRT.mat')
  runCommand('mrtransform T1.mif T1_registered.mif -linear T1_to_DWI_MRtrix.mat')
  delFile('T1.mif')
  delFile('T1_to_DWI_MRtrix.mat')

  # Step 8: Generate 5TT image for ACT
  runCommand('5ttgen fsl T1_registered.mif 5TT.mif')

  # Step 9: Determine whether we are working with single-shell or multi-shell data
  # TODO Detect b=0 shell and remove from list before testing length
  shells = [ int(round(float(x))) for x in getHeaderInfo('dwi.mif', 'shells').split() ]
  multishell = (len(shells) > 2)

  # Step 10: Estimate response function(s) for spherical deconvolution
  if multishell:
    runCommand('dwi2response msmt_5tt dwi.mif response_wm.txt response_gm.txt response_csf.txt -mask dwi_mask.mif')
    rf_file_for_scaling = 'response_wm.txt'
  else:
    runCommand('dwi2response tournier dwi.mif response.txt -mask dwi_mask.mif')
    rf_file_for_scaling = 'response.txt'

  # Step 12: Perform spherical deconvolution
  #          Use a dilated mask for spherical deconvolution as a 'safety margin' -
  #          ACT should be responsible for stopping streamlines before they reach the edge of the DWI mask
  runCommand('maskfilter dwi_mask.mif dilate dwi_mask_dilated.mif -npass 3')
  if multishell:
    runCommand('dwi2fod msmt_csd dwi.mif response_wm.txt FOD_WM.mif response_gm.txt FOD_GM.mif response_csf.txt FOD_CSF.mif -mask dwi_mask_dilated.mif')
  else:
    # Still use the msmt_csd algorithm with single-shell data: Use hard non-negativity constraint
    runCommand('dwiextract dwi.mif - | dwi2fod msmt_csd - response.txt FOD_WM.mif')

  # Step 13: Generate the tractogram
  # TODO Determine the appropriate number of streamlines based on the number of nodes in the parcellation
  num_streamlines = 10000000
  if lib.app.args.streamlines:
    num_streamlines = lib.app.args.streamlines
  runCommand('tckgen FOD_WM.mif tractogram.tck -act 5TT.mif -backtrack -crop_at_gmwmi -cutoff 0.06 -maxlength 250 -number ' + str(num_streamlines) + ' -seed_dynamic FOD_WM.mif')

  # Step 14: Use SIFT2 to determine streamline weights
  fd_scale_gm_option = ''
  if not multishell:
    fd_scale_gm_option = ' -fd_scale_gm'
  runCommand('tcksift2 tractogram.tck FOD_WM.mif weights.csv -act 5TT.mif -out_mu mu.txt' + fd_scale_gm_option + ' -info')

  # Step 15: Generate the grey matter parcellation
  #          The necessary steps here will vary significantly depending on the parcellation scheme selected
  if lib.app.args.parc == 'fs_2005' or lib.app.args.parc == 'fs_2009':

    # Run FreeSurfer pipeline on this subject's T1 image
    runCommand('mrconvert T1_registered.mif T1_registered.nii -stride +1,+2,+3')
    runCommand('recon-all -sd ' + lib.app.tempDir + ' -subjid freesurfer -i T1_registered.nii')
    delFile('T1_registered.nii')
    runCommand('recon-all -sd ' + lib.app.tempDir + ' -subjid freesurfer -all')

    # Grab the relevant parcellation image and target lookup table for conversion
    parc_image_path = os.path.join('freesurfer', 'mri')

    if lib.app.args.parc == 'fs_2005':
      parc_image_path = os.path.join(parc_image_path, 'aparc.aseg.mgz')
    else:
      parc_image_path = os.path.join(parc_image_path, 'aparc.a2009.aseg.mgz')

    # Perform the index conversion
    runCommand('labelconvert ' + parc_image_path + ' ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc_init.mif')
    if not lib.app.args.nocleanup:
      shutil.rmtree('freesurfer')

    # Fix the sub-cortical grey matter parcellations using FSL FIRST
    runCommand('labelsgmfix parc_init.mif T1_registered.mif ' + mrtrix_lut_file + ' parc.mif')
    delFile('parc_init.mif')

  elif lib.app.args.parc == 'aal' or lib.app.args.parc == 'aal2':

    runCommand('mrconvert T1_registered.mif T1_registered.nii')
    # Can use MNI152 image provided with FSL for registration
    runCommand('flirt -ref ' + mni152_path + ' -in T1 -omat T1_to_MNI_FLIRT.mat -dof 12')
    delFile('T1_registered.nii')
    runCommand('transformconvert T1_to_MNI_FLIRT.mat T1.nii ' + mni_path + ' flirt_import T1_to_MNI_MRtrix.mat')
    delFile('T1_to_MNI_FLIRT.mat')
    runCommand('transformcalc T1_to_MNI_MRtrix.mat invert MNI_to_T1_MRtrix.mat')
    delFile('T1_to_MNI_MRtrix.mat')
    runCommand('mrtransform ' + parc_image_path + ' AAL.mif -linear MNI_to_T1_MRtrix.mat -template T1_registered.mif -interp nearest')
    delFile('MNI_to_T1_MRtrix.mat')
    runCommand('labelconvert AAL.mif ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc.mif')
    delFile('AAL.mif')

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



#analysis_choices = [ 'participant', 'group' ]
analysis_choices = [ 'participant' ]
parcellation_choices = [ 'aal', 'aal2', 'fs_2005', 'fs_2009' ]

lib.app.author = 'Robert E. Smith (robert.smith@florey.edu.au)'
lib.cmdlineParser.initialise('Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in MRtrix3')
lib.app.parser.add_argument('bids_dir', help='The directory with the input dataset formatted according to the BIDS standard.')
lib.app.parser.add_argument('output_dir', help='The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.')
lib.app.parser.add_argument('analysis_level', help='Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: ' + ', '.join(analysis_choices), choices=analysis_choices)
lib.app.parser.add_argument('-v', '--version', action='version', version=__version__)
batch_options = lib.app.parser.add_argument_group('Options specific to the batch processing of subject data')
batch_options.add_argument('--participant_label', nargs='+', help='The label(s) of the participant(s) that should be analyzed. The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a space-separated list.')
participant_options = lib.app.parser.add_argument_group('Options that are relevant to participant-level analysis')
participant_options.add_argument('-parc', help='The choice of connectome parcellation scheme. Options are: ' + ', '.join(parcellation_choices), choices=parcellation_choices)
participant_options.add_argument('-streamlines', type=int, help='The number of streamlines to generate for each subject')
# TODO Option to copy particular data files from participant level processing into the output directory
#group_options = lib.app.parser.add_argument_group('Options that are relevant to group-level analysis')
#testing_options = lib.app.parser.add_argument_group('Options for testing the run.py script')
lib.app.parser._option_string_actions['-nthreads'].option_strings = [ '-nthreads', '-n_cpus' ]
lib.app.parser._option_string_actions['-n_cpus'] = lib.app.parser._option_string_actions['-nthreads']
for i in lib.app.parser._actions:
  if i.dest == 'nthreads':
    i.option_strings = [ '-nthreads', '-n_cpus' ]
    break

lib.app.initialise()

if isWindows():
  errorMessage('Script cannot be run on Windows due to FSL dependency')

#runCommand('bids-validator ' + lib.app.args.bids_dir)

subjects_to_analyze = [ ]
# Only run a subset of subjects
if lib.app.args.participant_label:
  subjects_to_analyze = [ 'sub-' + i for i in lib.app.args.participant_label ]
  for subject_dir in subjects_to_analyze:
    if not os.path.isdir(os.path.join(lib.app.args.bids_dir, subject_dir)):
      errorMessage('Unable to find directory for subject: ' + subject_dir)
# Run all subjects sequentially
else:
  subject_dirs = glob.glob(os.path.join(lib.app.args.bids_dir, 'sub-*'))
  # subjects_to_analyze = subject_dirs
  subjects_to_analyze = [ 'sub-' + dir.split("-")[-1] for dir in subject_dirs ]


# Running participant level
if lib.app.args.analysis_level == "participant":

  for subject_label in subjects_to_analyze:
    printMessage('Commencing execution for subject ' + subject_label)
    runSubject(lib.app.args.bids_dir, subject_label, lib.app.args.output_dir)

# Running group level
elif lib.app.args.analysis_level == "group":

  errorMessage('Group-level analysis not yet supported')

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

  # Third group-level calculation: Generate the group mean connectome
  # For any higher-level analysis (e.g. NBSE, computing connectome global measures, etc.),
  #   trying to incorporate such analysis into this particular pipeline script is likely to
  #   overly complicate the interface, and not actually provide much in terms of
  #   convenience / reproducibility guarantees. The primary functionality of this group-level
  #   analysis is therefore to achieve inter-subject connection density normalisation; users
  #   then have the flexibility to subsequently analyse the data however they choose (ideally
  #   based on subject classification data provided with the BIDS-compliant dataset).
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



