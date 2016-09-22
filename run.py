#!/usr/bin/env python


import os, glob, shutil, sys
import lib.app, lib.cmdlineParser

from lib.binaryInPath  import binaryInPath
from lib.delFile       import delFile
from lib.errorMessage  import errorMessage
from lib.getFSLSuffix  import getFSLSuffix
from lib.getHeaderInfo import getHeaderInfo
from lib.getImageStat  import getImageStat
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
  runCommand('mrconvert T1_registered.mif T1_registered.nii -stride +1,+2,+3')
  if lib.app.args.parc == 'fs_2005' or lib.app.args.parc == 'fs_2009':

    # Run FreeSurfer pipeline on this subject's T1 image
    runCommand('recon-all -sd ' + lib.app.tempDir + ' -subjid freesurfer -i T1_registered.nii')
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

    # Can use MNI152 image provided with FSL for registration
    runCommand('flirt -ref ' + mni152_path + ' -in T1_registered.nii -omat T1_to_MNI_FLIRT.mat -dof 12')
    runCommand('transformconvert T1_to_MNI_FLIRT.mat T1_registered.nii ' + mni152_path + ' flirt_import T1_to_MNI_MRtrix.mat')
    
    delFile('T1_to_MNI_FLIRT.mat')
    runCommand('transformcalc T1_to_MNI_MRtrix.mat invert MNI_to_T1_MRtrix.mat')
    delFile('T1_to_MNI_MRtrix.mat')
    runCommand('mrtransform ' + parc_image_path + ' AAL.mif -linear MNI_to_T1_MRtrix.mat -template T1_registered.mif -interp nearest')
    delFile('MNI_to_T1_MRtrix.mat')
    runCommand('labelconvert AAL.mif ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc.mif')
    delFile('AAL.mif')

  else:
    errorMessage('Unknown parcellation scheme requested: ' + lib.app.args.parc)
  delFile('T1_registered.nii')

  # Step 16: Generate the connectome
  #          Only provide the standard density-weighted connectome for now
  runCommand('tck2connectome tractogram.tck parc.mif connectome.csv -tck_weights_in weights.csv')
  delFile('weights.csv')

  # Move necessary files to output directory
  shutil.copy('connectome.csv', os.path.join(output_dir, 'connectome', label + '_connectome.csv'))
  runCommand('mrconvert dwi.mif ' + os.path.join(output_dir, 'dwi', label + '_dwi.nii.gz')
             + ' -export_grad_fsl ' + os.path.join(output_dir, 'dwi', label + '_dwi.bvec') + ' ' + os.path.join(output_dir, 'dwi', label + '_dwi.bval')
             + ' -json_export ' + os.path.join(output_dir, 'dwi', label + '_dwi.json'))
  shutil.copy('mu.txt', os.path.join(output_dir, 'connectome', label + '_mu.txt'))
  shutil.copy(rf_file_for_scaling, os.path.join(output_dir, 'dwi', label + '_response.txt'))

  # Manually wipe and zero the temp directory (since we might be processing more than one subject)
  os.chdir(cwd)
  if lib.app.cleanup:
    printMessage('Deleting temporary directory ' + lib.app.tempDir)
    shutil.rmtree(lib.app.tempDir)
  else:
    printMessage('Contents of temporary directory kept, location: ' + lib.app.tempDir)
  lib.app.tempDir = ''

# End of runSubject() function









# Create runGroup() function, use a temporary directory
# Don't write to the output directory unless the function actually completes
def runGroup(output_dir):
  import lib.app, os, shutil

  # Check presence of all required input files before proceeding
  # Pre-calculate paths of all files since many will be used in more than one location
  class subjectPaths:
    def __init__(self, label):
      self.in_dwi        = os.path.join(output_dir, label, 'dwi', label + '_dwi.nii.gz')
      self.in_bvec       = os.path.join(output_dir, label, 'dwi', label + '_dwi.bvec')
      self.in_bval       = os.path.join(output_dir, label, 'dwi', label + '_dwi.bval')
      self.in_json       = os.path.join(output_dir, label, 'dwi', label + '_dwi.json')
      self.in_rf         = os.path.join(output_dir, label, 'dwi', label + '_response.txt')
      self.in_connectome = os.path.join(output_dir, label, 'connectome', label + '_connectome.csv')
      self.in_mu         = os.path.join(output_dir, label, 'connectome', label + '_mu.txt')

      for path in vars(self).values():
        if not os.path.exists(path):
          errorMessage('Unable to find critical subject data (expected location: ' + path + ')')

      self.temp_mask      = os.path.join('masks',  label + '.mif')
      self.temp_fa        = os.path.join('images', label + '.mif')
      self.temp_bzero     = os.path.join('bzeros', label + '.mif')
      self.temp_warp      = os.path.join('warps',  label + '.mif')
      self.temp_voxels    = os.path.join('voxels', label + '.mif')
      self.out_value      = os.path.join('values', label + '.txt')
      self.out_connectome = os.path.join('connectomes', label + '.csv')

      self.label = label

  subject_list = [ 'sub-' + dir.split("-")[-1] for dir in glob.glob(os.path.join(output_dir, 'sub-*')) ]
  subjects = [ ]
  for label in subject_list:
    subjects.append(subjectPaths(label))

  lib.app.makeTempDir()
  lib.app.gotoTempDir()

  # First pass through subject data in group analysis:
  #   - Grab DWI data (written back from single-subject analysis back into BIDS format)
  #   - Generate mask and FA images to be used in populate template generation
  #   - Generate mean b=0 image for each subject for later use
  os.makedirs('bzeros')
  os.makedirs('images')
  os.makedirs('masks')
  for s in subjects:
    grad_import_option = ' -fslgrad ' + s.in_bvec + ' ' + s.in_bval
    runCommand('dwi2mask ' + s.in_dwi + ' ' + s.temp_mask + grad_import_option)
    runCommand('dwi2tensor ' + s.in_dwi + ' - -mask ' + s.temp_mask + grad_import_option + ' | tensor2metric - -fa ' + s.temp_fa)
    runCommand('dwiextract ' + s.in_dwi + grad_import_option + ' - -bzero | mrmath - mean ' + s.temp_bzero + ' -axis 3')

  # First group-level calculation: Generate the population FA template
  runCommand('population_template images -mask_dir masks -warp_dir warps template.mif -linear_scale 0.25,0.5,1.0,1.0 -nl_scale 0.5,0.75,1.0,1.0,1.0 -nl_niter 5,5,5,5,5')
  if lib.app.cleanup:
    shutil.rmtree('images')
    shutil.rmtree('masks')

  # Second pass through subject data in group analysis:
  #   - Warp template FA image back to subject space & threshold to define a WM mask in subject space
  #   - Calculate the median subject b=0 value within this mask
  #   - Store this in a file, and contribute to calculation of the mean of these values across subjects
  #   - Contribute to the group average response function
  os.makedirs('values')
  os.makedirs('voxels')
  mean_median_bzero = 0.0
  mean_RF = [ ]
  for s in subjects:
    runCommand('mrtransform template.mif -warp_full ' + s.temp_warp + ' - -from 2 -template ' + s.temp_bzero + ' | mrthreshold - ' + s.temp_voxels + ' -abs 0.4')
    median_bzero = getImageStat(s.temp_bzero, 'median', s.temp_voxels)
    delFile(s.temp_bzero)
    delFile(s.temp_voxels)
    delFile(s.temp_warp)
    with open(s.out_value, 'w') as f:
      f.write (median_bzero)
    mean_median_bzero = mean_median_bzero + float(median_bzero)
    RF = [ ]
    with open(s.in_rf, 'r') as f:
      for line in f:
        RF.append([ float(v) for v in line.split() ])
    RF_lzero = [ line[0] for line in RF ]
    if mean_RF:
      mean_RF = mean_RF + RF_lzero
    else:
      mean_RF = RF_lzero
  if lib.app.cleanup:
    shutil.rmtree('bzeros')
    shutil.rmtree('voxels')
    shutil.rmtree('warps')

  # Second group-level calculation:
  #   - Calculate the mean of median b=0 values
  #   - Calculate the mean response function
  # TODO Write these to file?
  # Perhaps a better option would be a text file summary of all inter-subject connection density normalisation parameters per subject, as well as all group mean data
  mean_median_bzero = mean_median_bzero / len(subjects)
  mean_RF = [ v / len(subjects) for v in mean_RF ]

  # Third pass through subject data in group analysis:
  #   - Scale the connectome strengths:
  #     - Multiply by SIFT proportionality coefficient mu
  #     - Multiply by (mean median b=0) / (subject median b=0)
  #     - Multiply by (subject RF size) / (mean RF size)
  #   - Write the result to file
  os.makedirs('connectomes')
  for s in subjects:
    with open(s.in_mu, 'r') as f:
      mu = float(f.read())
    with open(s.out_value, 'r') as f:
      median_bzero = float(f.read())
    RF = [ ]
    with open(s.in_rf, 'r') as f:
      for line in f:
        RF.append([ float(v) for v in line.split() ])
    RF_lzero = [ line[0] for line in RF ]
    RF_multiplier = 1.0
    for (mean, subj) in zip(mean_RF, RF_lzero):
      RF_multiplier = RF_multiplier * subj / mean
    global_multiplier = mu * (mean_median_bzero / median_bzero) * RF_multiplier

    connectome = [ ]
    with open(s.in_connectome, 'r') as f:
      for line in f:
        connectome.append ( [ float(v) for v in line.split() ] )
    with open(s.out_connectome, 'w') as f:
      for line in connectome:
        f.write(' '.join([ str(v*global_multiplier) for v in line ]) + '\n')

  # Third group-level calculation: Generate the group mean connectome
  # For any higher-level analysis (e.g. NBSE, computing connectome global measures, etc.),
  #   trying to incorporate such analysis into this particular pipeline script is likely to
  #   overly complicate the interface, and not actually provide much in terms of
  #   convenience / reproducibility guarantees. The primary functionality of this group-level
  #   analysis is therefore to achieve inter-subject connection density normalisation; users
  #   then have the flexibility to subsequently analyse the data however they choose (ideally
  #   based on subject classification data provided with the BIDS-compliant dataset).
  mean_connectome = [ ]
  for s in subjects:
    connectome = [ ]
    with open(s.out_connectome, 'r') as f:
      for line in f:
        connectome.append( [ float(v) for v in line.split() ] )
    if mean_connectome:
      for r1,r2 in zip(connectome, mean_connectome):
        r2 = [ c1+c2 for c1,c2 in zip(r1,r2) ]
    else:
      mean_connectome = connectome

  mean_connectome = [ [ v/len(subjects) for v in row ] for row in mean_connectome ]

  # Write results of interest back to the output directory
  for s in subjects:
    shutil.copyfile(s.out_connectome, os.path.join(output_dir, s.label, 'connectome', s.label + '_scaled_connectome.csv'))
  with open(os.path.join(output_dir, 'mean_connectome.csv'), 'w') as f:
    for row in mean_connectome:
      f.write(' '.join( [ str(v) for v in row ] ) + '\n')

# End of runGroup() function








analysis_choices = [ 'participant', 'group' ]
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
# TODO Option(s) to copy particular data files from participant level / group level processing into the output directory
#group_options = lib.app.parser.add_argument_group('Options that are relevant to group-level analysis')
#testing_options = lib.app.parser.add_argument_group('Options for testing the run.py script')
# Modify the existing -nthreads option (created in lib.cmdlineParser) to also accept the usage '-n_cpus'
lib.app.parser._option_string_actions['-nthreads'].option_strings = [ '-nthreads', '-n_cpus' ]
lib.app.parser._option_string_actions['-n_cpus'] = lib.app.parser._option_string_actions['-nthreads']
for i in lib.app.parser._actions:
  if i.dest == 'nthreads':
    i.option_strings = [ '-nthreads', '-n_cpus' ]
    break

lib.app.initialise()

if isWindows():
  errorMessage('Script cannot be run on Windows due to FSL dependency')

runCommand('bids-validator ' + lib.app.args.bids_dir)

# Running participant level
if lib.app.args.analysis_level == 'participant':

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
    subjects_to_analyze = [ 'sub-' + dir.split("-")[-1] for dir in subject_dirs ]

  for subject_label in subjects_to_analyze:
    printMessage('Commencing execution for subject ' + subject_label)
    runSubject(lib.app.args.bids_dir, subject_label, os.path.abspath(lib.app.args.output_dir))

# Running group level
elif lib.app.args.analysis_level == 'group':

  if lib.app.args.participant_label:
    errorMessage('Cannot use --participant_label option when performing group analysis')
  runGroup(os.path.abspath(lib.app.args.output_dir))

lib.app.complete()


