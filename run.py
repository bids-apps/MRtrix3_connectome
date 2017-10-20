#!/usr/bin/env python

import glob, json, math, os, shutil
from distutils.spawn import find_executable
from mrtrix3 import app, file, fsl, image, path, run


is_container = os.path.exists('/version')
__version__ = 'BIDS-App \'MRtrix3_connectome\' version {}'.format(open('/version').read()) if is_container else 'BIDS-App \'MRtrix3_connectome\' standalone'
option_prefix = '--' if is_container else '-'


def runSubject(bids_dir, label, output_prefix):

  output_dir = os.path.join(output_prefix, label)
  if os.path.exists(output_dir):
    app.warn('Output directory for subject \'' + label + '\' already exists; contents will be erased when this execution completes')

  fsl_path = os.environ.get('FSLDIR', '')
  if not fsl_path:
    app.error('Environment variable FSLDIR is not set; please run appropriate FSL configuration script')

  flirt_cmd = fsl.exeName('flirt')

  robex_found = find_executable('ROBEX') and find_executable('runROBEX.sh')
  N4_found = find_executable('N4BiasFieldCorrection')

  if robex_found and N4_found:
    app.console('N4BiasFieldCorrection and ROBEX found; will use for bias field correction and brain extraction')
    brain_extraction_cmd = 'runROBEX.sh'
  else:
    if robex_found and not N4_found:
      app.console('N4BiasFieldCorrection not found; will use FSL for bias field correction & brain extraction')
    elif N4_found and not robex_found:
      app.console('ROBEX not found; will use FSL for brain extraction')
    else:
      app.console('N4BiasFieldCorrection and ROBEX not found; will use FSL for brain extraction')
    brain_extraction_cmd = fsl.exeName('fsl_anat')

  dwibiascorrect_algo = '-ants'
  if not N4_found:
    # Can't use findFSLBinary() here, since we want to proceed even if it's not found
    if find_executable('fast') or find_executable('fsl5.0-fast'):
      dwibiascorrect_algo = '-fsl'
      app.console('Could not find ANTs program N4BiasFieldCorrection; '
                  'using FSL FAST for bias field correction')
    else:
      dwibiascorrect_algo = ''
      app.warn('Could not find ANTs program \'N4BiasFieldCorrection\' or FSL program \'fast\'; '
               'will proceed without performing DWI bias field correction')

  if not app.args.parcellation:
    app.error('For participant-level analysis, desired parcellation must be provided using the -parcellation option')

  template_image_path = ''
  parc_image_path = ''
  parc_lut_file = ''
  mrtrix_lut_file = os.path.join(os.path.dirname(os.path.abspath(app.__file__)),
                                 os.pardir,
                                 os.pardir,
                                 'share',
                                 'mrtrix3',
                                 'labelconvert')

  if app.args.parcellation in [ 'fs_2005', 'fs_2009' ]:
    if not 'FREESURFER_HOME' in os.environ:
      app.error('Environment variable FREESURFER_HOME not set; please verify FreeSurfer installation')
    freesurfer_subjects_dir = os.environ['SUBJECTS_DIR'] if 'SUBJECTS_DIR' in os.environ else os.path.join(os.environ['FREESURFER_HOME'], 'subjects')
    if not os.path.isdir(freesurfer_subjects_dir):
      app.error('Could not find FreeSurfer subjects directory (expected location: ' + freesurfer_subjects_dir + ')')
    for subdir in [ 'fsaverage', 'lh.EC_average', 'rh.EC_average' ]:
      if not os.path.isdir(os.path.join(freesurfer_subjects_dir, subdir)):
        app.error('Could not find requisite FreeSurfer subject directory \'' + subdir + '\' (expected location: ' + os.path.join(freesurfer_subjects_dir, subdir) + ')')
    reconall_path = find_executable('recon-all')
    if not reconall_path:
      app.error('Could not find FreeSurfer script recon-all; please verify FreeSurfer installation')
    # Query contents of recon-all script, looking for "-openmp" and "-parallel" occurences
    # Add options to end of recon-all -all call, based on which of these options are available
    #   as well as the value of app.numThreads
    # - In 5.3.0, just the -openmp option is available
    # - In 6.0.0, -openmp needs to be preceded by -parallel
    reconall_multithread_options = []
    if app.numThreads is None or app.numThreads:
      with open(reconall_path, 'r') as f:
        reconall_text = f.read().splitlines()
      for line in reconall_text:
        line = line.strip()
        if line == 'case "-parallel":':
          reconall_multithread_options = [ '-parallel' ] + reconall_multithread_options
        # If number of threads in this script is not being explicitly controlled,
        #   allow recon-all to use its own default number of threads
        elif line == 'case "-openmp":' and app.numThreads is not None:
          reconall_multithread_options.extend([ '-openmp', str(app.numThreads) ])
    if reconall_multithread_options:
      reconall_multithread_options = ' ' + ' '.join(reconall_multithread_options)
    else:
      reconall_multithread_options = ''
    app.var(reconall_multithread_options)
    parc_lut_file = os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')
    if app.args.parcellation == 'fs_2005':
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'fs_default.txt')
    else:
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'fs_a2009s.txt')

  elif app.args.parcellation in [ 'aal', 'aal2' ]:
    template_image_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_1mm.nii.gz')
    if app.args.parcellation == 'aal':
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'aal.txt')
    else:
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'aal2.txt')

  elif app.args.parcellation == 'perry512':
    template_image_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_1mm.nii.gz')
    parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', '512inMNI.nii'))
    parc_lut_file = mrtrix_lut_file = ''

  elif app.args.parcellation == 'sri24':
    template_image_path = os.path.join(os.sep, 'opt', 'sri24','spgr.nii')
    parc_image_path = os.path.join(os.sep, 'opt', 'sri24', 'lpba40.nii')
    parc_lut_file = os.path.join(os.sep, 'opt', 'sri24', 'LPBA40-labels.txt')
    mrtrix_lut_file = os.path.join(mrtrix_lut_file, 'lpba40.txt')

  if template_image_path and not os.path.isfile(template_image_path):
    if getattr(app.args, 'atlas_path', None):
      template_image_path = [ template_image_path, os.path.join(os.path.dirname(app.args.atlas_path), os.path.basename(template_image_path)) ]
      if os.path.isfile(template_image_path[1]):
        template_image_path = template_image_path[1]
      else:
        app.error('Could not find template image (tested locations: ' + str(template_image_path) + ')')
    else:
      app.error('Could not find template image (expected location: ' + template_image_path + ')')

  if parc_image_path and not os.path.isfile(parc_image_path):
    if getattr(app.args, 'atlas_path', None):
      parc_image_path = [ parc_image_path, os.path.join(os.path.dirname(app.args.atlas_path), os.path.basename(parc_image_path)) ]
      if os.path.isfile(parc_image_path[1]):
        parc_image_path = parc_image_path[1]
      else:
        app.error('Could not find parcellation image (tested locations: ' + str(parc_image_path) + ')')
    else:
      app.error('Could not find parcellation image (expected location: ' + parc_image_path + ')')

  if parc_lut_file and not os.path.isfile(parc_lut_file):
    if getattr(app.args, 'atlas_path', None):
      parc_lut_file = [ parc_lut_file, os.path.join(os.path.dirname(app.args.atlas_path), os.path.basename(parc_lut_file)) ]
      if os.path.isfile(parc_lut_file[1]):
        parc_lut_file = parc_lut_file[1]
      else:
        app.error('Could not find parcellation lookup table file (tested locations: ' + str(parc_lut_file) + ')')
    else:
      app.error('Could not find parcellation lookup table file (expected location: ' + parc_lut_file + ')')

  if mrtrix_lut_file and not os.path.exists(mrtrix_lut_file):
    app.error('Could not find MRtrix3 connectome lookup table file (expected location: ' + mrtrix_lut_file + ')')

  app.makeTempDir()

  # Need to perform an initial import of JSON data using mrconvert; so let's grab the diffusion gradient table as well
  # If no bvec/bval present, need to go down the directory listing
  # Only try to import JSON file if it's actually present
  #   direction in the acquisition they'll need to be split across multiple files
  # May need to concatenate more than one input DWI, since if there's more than one phase-encode direction
  #   in the acquired DWIs (i.e. not just those used for estimating the inhomogeneity field), they will
  #   need to be stored as separate NIfTI files in the 'dwi/' directory.
  dwi_image_list = glob.glob(os.path.join(bids_dir, label, 'dwi', label) + '*_dwi.nii*')
  dwi_index = 1
  for entry in dwi_image_list:
    # os.path.split() falls over with .nii.gz extensions; only removes the .gz
    prefix = entry.split(os.extsep)[0]
    if os.path.isfile(prefix + '.bval') and os.path.isfile(prefix + '.bvec'):
      prefix = prefix + '.'
    else:
      prefix = os.path.join(bids_dir, 'dwi')
      if not (os.path.isfile(prefix + 'bval') and os.path.isfile(prefix + 'bvec')):
        app.error('Unable to locate valid diffusion gradient table for image \'' + entry + '\'')
    grad_import_option = ' -fslgrad ' + prefix + 'bvec ' + prefix + 'bval'
    json_path = prefix + 'json'
    if os.path.isfile(json_path):
      json_import_option = ' -json_import ' + json_path
    else:
      json_import_option = ''
    run.command('mrconvert ' + entry + grad_import_option + json_import_option
                + ' ' + path.toTemp('dwi' + str(dwi_index) + '.mif', True))
    dwi_index += 1

  # Go hunting for reversed phase-encode data dedicated to field map estimation
  fmap_image_list = []
  fmap_dir = os.path.join(bids_dir, label, 'fmap')
  fmap_index = 1
  if os.path.isdir(fmap_dir):
    if app.args.preprocessed:
      app.error('fmap/ directory detected for subject \'' + label + '\' despite use of ' + option_prefix + 'preprocessed option')
    fmap_image_list = glob.glob(os.path.join(fmap_dir, label) + '_dir-*_epi.nii*')
    for entry in fmap_image_list:
      prefix = entry.split(os.extsep)[0]
      json_path = prefix + '.json'
      with open(json_path, 'r') as f:
        json_elements = json.load(f)
      if 'IntendedFor' in json_elements and not any(i.endswith(json_elements['IntendedFor']) for i in dwi_image_list):
        app.console('Image \'' + entry + '\' is not intended for use with DWIs; skipping')
        continue
      if os.path.isfile(json_path):
        json_import_option = ' -json_import ' + json_path
        # fmap files will not come with any gradient encoding in the JSON;
        #   therefore we need to add it manually ourselves so that mrcat / mrconvert can
        #   appropriately handle the table once these images are concatenated with the DWIs
        fmap_image_size = image.Header(entry).size()
        fmap_image_num_volumes = 1 if len(fmap_image_size) == 3 else fmap_image_size[3]
        run.command('mrconvert ' + entry + json_import_option +
                    ' -set_property dw_scheme \"' +
                    '\\n'.join(['0,0,1,0'] * fmap_image_num_volumes) +
                    '\" ' +
                    path.toTemp('fmap' + str(fmap_index) + '.mif', True))
        fmap_index += 1
      else:
        app.warn('No corresponding .json file found for image \'' + entry + '\'; skipping')

    fmap_image_list = [ 'fmap' + str(index) + '.mif' for index in range(1, fmap_index) ]
  # If there's no data in fmap/ directory, need to check to see if there's any phase-encoding
  #   contrast within the input DWI(s)
  elif len(dwi_image_list) < 2 and not app.args.preprocessed:
    app.error('Inadequate data for pre-processing of subject \'' + label + '\': No phase-encoding contrast in input DWIs or fmap/ directory')

  dwi_image_list = [ 'dwi' + str(index) + '.mif' for index in range(1, dwi_index) ]

  # Import anatomical image
  run.command('mrconvert ' + os.path.join(bids_dir, label, 'anat', label + '_T1w.nii.gz') + ' ' +
              path.toTemp('T1.mif', True))

  cwd = os.getcwd()
  app.gotoTempDir()

  dwipreproc_se_epi = ''
  dwipreproc_se_epi_option = ''

  # For automated testing, down-sampled images are used. However, this invalidates the requirements of
  #   both MP-PCA denoising and Gibbs ringing removal. In addition, eddy can still take a long time
  #   despite the down-sampling. Therefore, provide images that have been pre-processed to the stage
  #   where it is still only DWI, JSON & bvecs/bvals that need to be provided.
  if app.args.preprocessed:

    if len(dwi_image_list) > 1:
      app.error('If DWIs have been pre-processed, then only a single DWI file should be present')
    app.console('Skipping MP-PCA denoising, Gibbs ringing removal, distortion correction and bias field correction due to use of ' + option_prefix + 'preprocessed option')
    run.function(os.rename, dwi_image_list[0], 'dwi.mif')

  else: # Do initial image pre-processing (denoising, Gibbs ringing removal if available, distortion correction & bias field correction) as normal

    # Concatenate any SE EPI images with the DWIs before denoising (& unringing), then
    #   separate them again after the fact
    # TODO Note however that this may not be possible: the fmap/ images may not lie on the
    #   same grid as the DWIs. When this occurs, cannot apply denoising algorithm to fmap images.
    if fmap_image_list:
      run.command('mrcat ' + ' '.join(fmap_image_list) + ' fmap_cat.mif -axis 3')
      for i in fmap_image_list:
        file.delTemporary(i)
      fmap_header = image.Header('fmap_cat.mif')
      fmap_num_volumes = fmap_header.size()[3]
      dwidenoise_input = 'dwi_fmap_cat.mif'
      run.command('mrcat fmap_cat.mif ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
      file.delTemporary('fmap_cat.mif')
      for i in dwi_image_list:
        file.delTemporary(i)
    else:
      fmap_num_volumes = 0
      # Even if no explicit fmap images, may still need to concatenate multiple DWI inputs
      if len(dwi_image_list) > 1:
        dwidenoise_input = 'dwi_cat.mif'
        run.command('mrcat ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
        for i in dwi_image_list:
          file.delTemporary(i)
      else:
        dwidenoise_input = dwi_image_list[0]


    # Step 1: Denoise
    mrdegibbs_input = os.path.splitext(dwidenoise_input)[0] + '_denoised.mif'
    run.command('dwidenoise ' + dwidenoise_input + ' ' + mrdegibbs_input)
    file.delTemporary(dwidenoise_input)

    # Step 2: Gibbs ringing removal
    mrdegibbs_output = os.path.splitext(mrdegibbs_input)[0] + '_degibbs.mif'
    run.command('mrdegibbs ' + mrdegibbs_input + ' ' + mrdegibbs_output + ' -nshifts 50')
    file.delTemporary(mrdegibbs_input)

    # If fmap images and DWIs have been concatenated, now is the time to split them back apart
    if fmap_num_volumes:
      dwipreproc_input = 'dwipreproc_in.mif'
      dwipreproc_se_epi = 'se_epi.mif'
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_se_epi + ' -coord 3 0:' + str(fmap_num_volumes-1))
      cat_num_volumes = image.Header(mrdegibbs_output).size()[3]
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_input + ' -coord 3 ' + str(fmap_num_volumes) + ':' + str(cat_num_volumes-1))
      file.delTemporary(mrdegibbs_output)
      dwipreproc_se_epi_option = ' -se_epi ' + dwipreproc_se_epi
    else:
      dwipreproc_input = mrdegibbs_output
      dwipreproc_se_epi = None
      dwipreproc_se_epi_option = ''

    # Step 3: Distortion correction
    run.command('dwipreproc ' + dwipreproc_input + ' dwi_preprocessed.mif -rpe_header' + dwipreproc_se_epi_option)
    file.delTemporary(dwipreproc_input)
    if dwipreproc_se_epi:
      file.delTemporary(dwipreproc_se_epi)

    # Step 4: Bias field correction
    if dwibiascorrect_algo:
      run.command('dwibiascorrect dwi_preprocessed.mif dwi.mif ' + dwibiascorrect_algo)
      file.delTemporary('dwi_preprocessed.mif')
    else:
      run.function(shutil.move, 'dwi_preprocessed.mif', 'dwi.mif')

  # No longer branching based on whether or not -preprocessed was specified

  # Step 5: Generate a brain mask for DWI
  #   Also produce a dilated version of the mask for later use
  run.command('dwi2mask dwi.mif dwi_mask.mif')
  run.command('maskfilter dwi_mask.mif dilate dwi_mask_dilated.mif -npass 3')

  # TODO Crop DWIs based on brain mask

  # Step 6: Perform brain extraction and bias field correction on the T1 image
  #         in its original space (this is necessary for histogram matching
  #         prior to registration)
  T1_header = image.Header('T1.mif')
  T1_revert_stride_option = ' -stride ' + ','.join([str(i) for i in T1_header.stride()])
  if brain_extraction_cmd == 'runROBEX.sh':
    # Do a semi-iterative approach here: Get an initial brain mask, use that
    #   mask to estimate a bias field, then re-compute the brain mask
    run.command('mrconvert T1.mif T1.nii -stride +1,+2,+3')
    run.command(brain_extraction_cmd + ' T1.nii T1_initial_brain.nii T1_initial_mask.nii')
    file.delTemporary('T1_initial_brain.nii')
    run.command('N4BiasFieldCorrection -i T1.nii -w T1_initial_mask.nii -o T1_biascorr.nii')
    file.delTemporary('T1.nii')
    file.delTemporary('T1_initial_mask.nii')
    run.command(brain_extraction_cmd + ' T1_biascorr.nii T1_biascorr_brain.nii T1_biascorr_brain_mask.nii')
    file.delTemporary('T1_biascorr.nii')
    run.command('mrconvert T1_biascorr_brain.nii T1_biascorr_brain.mif' + T1_revert_stride_option)
    file.delTemporary('T1_biascorr_brain.nii')
    run.command('mrconvert T1_biascorr_brain_mask.nii T1_mask.mif -datatype bit' + T1_revert_stride_option)
    file.delTemporary('T1_biascorr_brain_mask.nii')
  else:
    run.command('mrconvert T1.mif T1.nii -stride -1,+2,+3')
    run.command(brain_extraction_cmd + ' -i T1.nii --noseg --nosubcortseg')
    file.delTemporary('T1.nii')
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain') + ' T1_biascorr_brain.mif' + T1_revert_stride_option)
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain_mask') + ' T1_mask.mif -datatype bit' + T1_revert_stride_option)
    file.delTemporary('T1.anat')

  # Step 7: Generate target images for T1->DWI registration
  run.command('dwiextract dwi.mif -bzero - | '
              'mrcalc - 0.0 -max - | '
              'mrmath - mean -axis 3 dwi_meanbzero.mif')
  run.command('mrcalc 1 dwi_meanbzero.mif -div dwi_mask.mif -mult - | '
              'mrhistmatch - T1_biascorr_brain.mif dwi_pseudoT1.mif -mask_input dwi_mask.mif -mask_target T1_mask.mif')
  run.command('mrcalc 1 T1_biascorr_brain.mif -div T1_mask.mif -mult - | '
              'mrhistmatch - dwi_meanbzero.mif T1_pseudobzero.mif -mask_input T1_mask.mif -mask_target dwi_mask.mif')

  # Step 8: Perform T1->DWI registration
  #         Note that two registrations are performed: Even though we have a symmetric registration,
  #         generation of the two histogram-matched images means that you will get slightly different
  #         answers depending on which synthesized image & original image you use.
  run.command('mrregister T1_biascorr_brain.mif dwi_pseudoT1.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_T1_to_pseudoT1.txt')
  file.delTemporary('T1_biascorr_brain.mif')
  run.command('mrregister T1_pseudobzero.mif dwi_meanbzero.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_pseudobzero_to_bzero.txt')
  file.delTemporary('dwi_meanbzero.mif')
  run.command('transformcalc rigid_T1_to_pseudoT1.txt rigid_pseudobzero_to_bzero.txt average rigid_T1_to_dwi.txt')
  file.delTemporary('rigid_T1_to_pseudoT1.txt')
  file.delTemporary('rigid_pseudobzero_to_bzero.txt')
  run.command('mrtransform T1.mif T1_registered.mif -linear rigid_T1_to_dwi.txt')
  file.delTemporary('T1.mif')
  # Note: Since we're using a mask from fsl_anat (which crops the FoV), but using it as input to 5ttge fsl
  #   (which is receiving the raw T1), we need to resample in order to have the same dimensions between these two
  run.command('mrtransform T1_mask.mif T1_mask_registered.mif -linear rigid_T1_to_dwi.txt -template T1_registered.mif -interp nearest')
  file.delTemporary('T1_mask.mif')

  # Step 9: Generate 5TT image for ACT
  run.command('5ttgen fsl T1_registered.mif 5TT.mif -mask T1_mask_registered.mif')
  if app.args.output_verbosity > 1:
    run.command('5tt2vis 5TT.mif vis.mif')
  if app.args.output_verbosity <= 1:
    file.delTemporary('T1_mask_registered.mif')

  # Step 10: Estimate response functions for spherical deconvolution
  run.command('dwi2response dhollander dwi.mif response_wm.txt response_gm.txt response_csf.txt -mask dwi_mask.mif')

  # Step 11: Determine whether we are working with single-shell or multi-shell data
  shells = [int(round(float(value))) for value in image.mrinfo('dwi.mif', 'shellvalues').strip().split()]
  multishell = (len(shells) > 2)

  # Step 12: Perform spherical deconvolution
  #          Use a dilated mask for spherical deconvolution as a 'safety margin' -
  #          ACT should be responsible for stopping streamlines before they reach the edge of the DWI mask
  if multishell:
    run.command('dwi2fod msmt_csd dwi.mif response_wm.txt FOD_WM.mif response_gm.txt FOD_GM.mif response_csf.txt FOD_CSF.mif '
                '-mask dwi_mask_dilated.mif -lmax 10,0,0')
    run.command('mrconvert FOD_WM.mif - -coord 3 0 | mrcat FOD_CSF.mif FOD_GM.mif - tissues.mif -axis 3')
  else:
    # Still use the msmt_csd algorithm with single-shell data: Use hard non-negativity constraint
    # Also incorporate the CSF response to provide some fluid attenuation
    run.command('dwi2fod msmt_csd dwi.mif response_wm.txt FOD_WM.mif response_csf.txt FOD_CSF.mif '
                '-mask dwi_mask_dilated.mif -lmax 10,0')
    file.delTemporary('FOD_CSF.mif')

  # Step 13: Generate the grey matter parcellation
  #          The necessary steps here will vary significantly depending on the parcellation scheme selected
  run.command('mrconvert T1_registered.mif T1_registered.nii -stride +1,+2,+3')
  if app.args.parcellation in [ 'fs_2005', 'fs_2009' ]:

    # Since we're instructing recon-all to use a different subject directory, we need to
    #   construct softlinks to a number of directories provided by FreeSurfer that
    #   recon-all will expect to find in the same directory as the overridden subject path
    for subdir in [ 'fsaverage', 'lh.EC_average', 'rh.EC_average' ]:
      run.function(os.symlink, os.path.join(freesurfer_subjects_dir, subdir), subdir)

    # Run FreeSurfer pipeline on this subject's T1 image
    run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -i T1_registered.nii')
    run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -all' + reconall_multithread_options)

    # Grab the relevant parcellation image and target lookup table for conversion
    parc_image_path = os.path.join('freesurfer', 'mri')
    if app.args.parcellation == 'fs_2005':
      parc_image_path = os.path.join(parc_image_path, 'aparc+aseg.mgz')
    else:
      parc_image_path = os.path.join(parc_image_path, 'aparc.a2009s+aseg.mgz')

    # Perform the index conversion
    run.command('labelconvert ' + parc_image_path + ' ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc_init.mif')
    if app.cleanup:
      run.function(shutil.rmtree, 'freesurfer')

    # Fix the sub-cortical grey matter parcellations using FSL FIRST
    run.command('labelsgmfix parc_init.mif T1_registered.mif ' + mrtrix_lut_file + ' parc.mif')
    file.delTemporary('parc_init.mif')

  elif app.args.parcellation in [ 'aal', 'aal2', 'perry512', 'sri24' ]:

    # TODO Currently registration to SRI24 goes well and truly awry
    run.command(flirt_cmd + ' -ref ' + template_image_path + ' -in T1_registered.nii -omat T1_to_template_FLIRT.mat -dof 12')
    run.command('transformconvert T1_to_template_FLIRT.mat T1_registered.nii ' + template_image_path + ' flirt_import T1_to_template_MRtrix.mat')
    file.delTemporary('T1_to_template_FLIRT.mat')
    run.command('transformcalc T1_to_template_MRtrix.mat invert template_to_T1_MRtrix.mat')
    file.delTemporary('T1_to_template_MRtrix.mat')
    run.command('mrtransform ' + parc_image_path + ' atlas_transformed.mif -linear template_to_T1_MRtrix.mat '
                '-template T1_registered.mif -interp nearest')
    file.delTemporary('template_to_T1_MRtrix.mat')
    if parc_lut_file or mrtrix_lut_file:
      assert parc_lut_file and mrtrix_lut_file
      run.command('labelconvert atlas_transformed.mif ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc.mif')
      file.delTemporary('atlas_transformed.mif')
    else: # Not all parcellations need to go through the labelconvert step; they may already be numbered incrementally from 1
      run.function(shutil.move, 'atlas_transformed.mif', 'parc.mif')

  else:
    app.error('Unknown parcellation scheme requested: ' + app.args.parcellation)
  file.delTemporary('T1_registered.nii')
  if app.args.output_verbosity > 3 and mrtrix_lut_file:
    run.command('label2colour parc.mif parcRGB.mif -lut ' + mrtrix_lut_file)

  # Step 14: Generate the tractogram
  # If not manually specified, determine the appropriate number of streamlines based on the number of nodes in the parcellation:
  #   mean edge weight of 1,000 streamlines
  # A smaller FOD amplitude threshold of 0.06 (default 0.1) is used for tracking due to the use of the msmt_csd
  #   algorithm, which imposes a hard rather than soft non-negativity constraint
  num_nodes = int(image.statistic('parc.mif', 'max'))
  num_streamlines = 500 * num_nodes * (num_nodes-1)
  if app.args.streamlines:
    num_streamlines = app.args.streamlines
  run.command('tckgen FOD_WM.mif tractogram.tck -act 5TT.mif -backtrack -crop_at_gmwmi -cutoff 0.06 -maxlength 250 -power 0.33 '
              '-select ' + str(num_streamlines) + ' -seed_dynamic FOD_WM.mif')

  # Step 15: Use SIFT2 to determine streamline weights
  fd_scale_gm_option = ''
  if not multishell:
    fd_scale_gm_option = ' -fd_scale_gm'
  run.command('tcksift2 tractogram.tck FOD_WM.mif weights.csv -act 5TT.mif -out_mu mu.txt' + fd_scale_gm_option)

  # Step 16: Generate a TDI at DWI native resolution, with SIFT mu scaling, and precise mapping
  #   (for comparison to WM ODF l=0 term, to verify that SIFT2 has worked correctly)
  if app.args.output_verbosity > 2:
    with open('mu.txt', 'r') as f:
      mu = float(f.read())
    run.command('tckmap tractogram.tck -tck_weights_in weights.csv -template FOD_WM.mif -precise - | '
                'mrcalc - ' + str(mu) + ' -mult tdi_native.mif')

  # Step 17: Generate a conventional TDI at super-resolution
  #   (mostly just because we can)
  if app.args.output_verbosity > 2:
    run.command('tckmap tractogram.tck -tck_weights_in weights.csv -template vis.mif -vox ' + ','.join([str(value/3.0) for value in image.Header('vis.mif').spacing() ]) + ' -datatype uint16 tdi_highres.mif')

  # Step 18: Generate the connectome
  #          Also get the mean length for each edge; this is the most likely alternative contrast to be useful
  run.command('tck2connectome tractogram.tck parc.mif connectome.csv -tck_weights_in weights.csv -out_assignments assignments.csv')
  run.command('tck2connectome tractogram.tck parc.mif meanlength.csv -tck_weights_in weights.csv -scale_length -stat_edge mean')

  # Step 19: Produce additional data that can be used for visualisation within mrview's connectome toolbar
  if app.args.output_verbosity > 2:
    run.command('connectome2tck tractogram.tck assignments.csv exemplars.tck -tck_weights_in weights.csv -exemplars parc.mif -files single')
    run.command('label2mesh parc.mif nodes.obj')
    run.command('meshfilter nodes.obj smooth nodes_smooth.obj')
    file.delTemporary('nodes.obj')

  # TODO Eventually will want to move certain elements into .json files rather than text files

  # Prepare output path for writing
  app.console('Processing for subject \'' + label + '\' completed; writing results to output directory')
  if os.path.exists(output_dir):
    run.function(shutil.rmtree, output_dir)
  run.function(os.makedirs, output_dir)
  run.function(os.makedirs, os.path.join(output_dir, 'connectome'))
  run.function(os.makedirs, os.path.join(output_dir, 'dwi'))
  run.function(os.makedirs, os.path.join(output_dir, 'tractogram'))
  if app.args.output_verbosity > 1:
    run.function(os.makedirs, os.path.join(output_dir, 'anat'))

  # Copy / convert necessary files to output directory
  run.function(shutil.copy, 'connectome.csv', os.path.join(output_dir, 'connectome', label + '_connectome.csv'))
  run.function(shutil.copy, 'meanlength.csv', os.path.join(output_dir, 'connectome', label + '_meanlength.csv'))
  run.function(shutil.copy, 'mu.txt', os.path.join(output_dir, 'tractogram', label + '_mu.txt'))
  run.function(shutil.copy, 'response_wm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-WM_response.txt'))
  with open(os.path.join(output_dir, 'dwi', label + '_shells.txt'), 'w') as f:
    f.write(' '.join([str(value) for value in shells]))
  if app.args.output_verbosity > 1:
    run.command('mrconvert dwi.mif ' + os.path.join(output_dir, 'dwi', label + '_dwi.nii.gz') + \
                ' -export_grad_fsl ' + os.path.join(output_dir, 'dwi', label + '_dwi.bvec') + ' ' + os.path.join(output_dir, 'dwi', label + '_dwi.bval') + \
                ' -stride +1,+2,+3,+4')
    run.command('mrconvert dwi_mask.mif ' + os.path.join(output_dir, 'dwi', label + '_brainmask.nii.gz') + ' -datatype uint8 -stride +1,+2,+3')
    run.command('mrconvert T1_registered.mif ' + os.path.join(output_dir, 'anat', label + '_T1w.nii.gz') + ' -stride +1,+2,+3')
    run.command('mrconvert T1_mask_registered.mif ' + os.path.join(output_dir, 'anat', label + '_T1w_brainmask.nii.gz') + ' -datatype uint8 -stride +1,+2,+3')
    run.command('mrconvert vis.mif ' + os.path.join(output_dir, 'anat', label + '_tissues3D.nii.gz') + ' -stride +1,+2,+3,+4')
    run.command('mrconvert FOD_WM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-WM_ODF.nii.gz') + ' -stride +1,+2,+3,+4')
    if multishell:
      run.function(shutil.copy, 'response_gm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-GM_response.txt'))
      run.function(shutil.copy, 'response_csf.txt', os.path.join(output_dir, 'dwi', label + '_tissue-CSF_response.txt'))
      run.command('mrconvert FOD_GM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-GM_ODF.nii.gz') + ' -stride +1,+2,+3,+4')
      run.command('mrconvert FOD_CSF.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-CSF_ODF.nii.gz') + ' -stride +1,+2,+3,+4')
      run.command('mrconvert tissues.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-all.nii.gz') + ' -stride +1,+2,+3,+4')
    run.command('mrconvert parc.mif ' + os.path.join(output_dir, 'anat', label + '_parc-' + app.args.parcellation + '_indices.nii.gz') + ' -stride +1,+2,+3')
  if app.args.output_verbosity > 2:
    # Move rather than copying the tractogram just because of its size
    run.function(shutil.move, 'tractogram.tck', os.path.join(output_dir, 'tractogram', label + '_tractogram.tck'))
    run.function(shutil.copy, 'weights.csv', os.path.join(output_dir, 'tractogram', label + '_weights.csv'))
    run.function(shutil.copy, 'assignments.csv', os.path.join(output_dir, 'connectome', label + '_assignments.csv'))
    run.function(shutil.copy, 'exemplars.tck', os.path.join(output_dir, 'connectome', label + '_exemplars.tck'))
    run.function(shutil.copy, 'nodes_smooth.obj', os.path.join(output_dir, 'anat', label + '_parc-' + app.args.parcellation + '.obj'))
    run.command('mrconvert parcRGB.mif ' + os.path.join(output_dir, 'anat', label + '_parc-' + app.args.parcellation +'_colour.nii.gz') + ' -stride +1,+2,+3')
    run.command('mrconvert tdi_native.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-native_tdi.nii.gz') + ' -stride +1,+2,+3')
    run.command('mrconvert tdi_highres.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-highres_tdi.nii.gz') + ' -stride +1,+2,+3')

  # Manually wipe and zero the temp directory (since we might be processing more than one subject)
  os.chdir(cwd)
  if app.cleanup:
    app.console('Deleting temporary directory ' + app.tempDir)
    # Can't use run.function() here; it'll try to write to the log file that resides in the temp directory just deleted
    shutil.rmtree(app.tempDir)
  else:
    app.console('Contents of temporary directory kept; location: ' + app.tempDir)
  app.tempDir = ''

# End of runSubject() function









def runGroup(output_dir):

  # Check presence of all required input files before proceeding
  # Pre-calculate paths of all files since many will be used in more than one location
  class subjectPaths(object):
    def __init__(self, label):
      self.in_dwi = os.path.join(output_dir, label, 'dwi', label + '_dwi.nii.gz')
      self.in_bvec = os.path.join(output_dir, label, 'dwi', label + '_dwi.bvec')
      self.in_bval = os.path.join(output_dir, label, 'dwi', label + '_dwi.bval')
      self.in_rf = os.path.join(output_dir, label, 'dwi', label + '_tissue-WM_response.txt')
      self.in_connectome = os.path.join(output_dir, label, 'connectome', label + '_connectome.csv')
      self.in_mu = os.path.join(output_dir, label, 'tractogram', label + '_mu.txt')

      for entry in vars(self).values():
        if not os.path.exists(entry):
          app.error('Unable to find critical subject data (expected location: ' + entry + ')')

      # Permissible for this to not exist
      self.in_mask = os.path.join(output_dir, label, 'dwi', label + '_brainmask.nii.gz')

      with open(self.in_mu, 'r') as f:
        self.mu = float(f.read())

      self.RF = []
      with open(self.in_rf, 'r') as f:
        for line in f:
          self.RF.append([ float(v) for v in line.split() ])

      self.temp_mask = os.path.join('masks',  label + '.mif')
      self.temp_fa = os.path.join('images', label + '.mif')
      self.temp_bzero = os.path.join('bzeros', label + '.mif')
      self.temp_warp = os.path.join('warps',  label + '.mif')
      self.temp_voxels = os.path.join('voxels', label + '.mif')
      self.median_bzero = 0.0
      self.dwiintensitynorm_factor = 1.0
      self.RF_multiplier = 1.0
      self.global_multiplier = 1.0
      self.temp_connectome = os.path.join('connectomes', label + '.csv')
      self.out_scale_bzero = os.path.join(output_dir, label, 'connectome', label + '_scalefactor_bzero.csv')
      self.out_scale_RF = os.path.join(output_dir, label, 'connectome', label + '_scalefactor_response.csv')
      self.out_connectome = os.path.join(output_dir, label, 'connectome', label + '_connectome_scaled.csv')

      self.label = label

  subject_list = ['sub-' + sub_dir.split("-")[-1] for sub_dir in glob.glob(os.path.join(output_dir, 'sub-*'))]
  if not subject_list:
    app.error('No processed subject data found in output directory for group analysis')
  subjects = []
  for label in subject_list:
    subjects.append(subjectPaths(label))

  app.makeTempDir()
  app.gotoTempDir()

  # First pass through subject data in group analysis:
  #   - Grab DWI data (written back from single-subject analysis back into BIDS format)
  #   - Generate mask and FA images to be used in populate template generation
  #       (though if output_verbosity > 2 a mask is already provided)
  #   - Generate mean b=0 image for each subject for later use
  progress = app.progressBar('Importing and preparing subject data', len(subjects))
  run.function(os.makedirs, 'bzeros')
  run.function(os.makedirs, 'images')
  run.function(os.makedirs, 'masks')
  for s in subjects:
    grad_import_option = ' -fslgrad ' + s.in_bvec + ' ' + s.in_bval
    if os.path.exists(s.in_mask):
      run.command('mrconvert ' + s.in_mask + ' ' + s.temp_mask)
    else:
      run.command('dwi2mask ' + s.in_dwi + ' ' + s.temp_mask + grad_import_option)
    run.command('dwi2tensor ' + s.in_dwi + ' - -mask ' + s.temp_mask + grad_import_option + ' | tensor2metric - -fa ' + s.temp_fa)
    run.command('dwiextract ' + s.in_dwi + grad_import_option + ' - -bzero | mrmath - mean ' + s.temp_bzero + ' -axis 3')
    progress.increment()
  progress.done()

  # First group-level calculation: Generate the population FA template
  app.console('Generating population template for inter-subject intensity normalisation WM mask derivation')
  run.command('population_template images -mask_dir masks -warp_dir warps template.mif '
              '-type rigid_affine_nonlinear -rigid_scale 0.25,0.5,0.8,1.0 -affine_scale 0.7,0.8,1.0,1.0 '
              '-nl_scale 0.5,0.75,1.0,1.0,1.0 -nl_niter 5,5,5,5,5 -linear_no_pause')
  file.delTemporary('images')
  file.delTemporary('masks')

  # Second pass through subject data in group analysis:
  #   - Warp template FA image back to subject space & threshold to define a WM mask in subject space
  #   - Calculate the median subject b=0 value within this mask
  #   - Store this in a file, and contribute to calculation of the mean of these values across subjects
  #   - Contribute to the group average response function
  progress = app.progressBar('Generating group-average response function and intensity normalisation factors', len(subjects)+1)
  run.function(os.makedirs, 'voxels')
  sum_median_bzero = 0.0
  sum_RF = []
  for s in subjects:
    run.command('mrtransform template.mif -warp_full ' + s.temp_warp + ' - -from 2 -template ' + s.temp_bzero + ' | '
                'mrthreshold - ' + s.temp_voxels + ' -abs 0.4')
    s.median_bzero = float(image.statistic(s.temp_bzero, 'median', '-mask ' + s.temp_voxels))
    file.delTemporary(s.temp_bzero)
    file.delTemporary(s.temp_voxels)
    file.delTemporary(s.temp_warp)
    sum_median_bzero += s.median_bzero
    if sum_RF:
      sum_RF = [[a+b for a, b in zip(one, two)] for one, two in zip(sum_RF, s.RF)]
    else:
      sum_RF = s.RF
    progress.increment()
  file.delTemporary('bzeros')
  file.delTemporary('voxels')
  file.delTemporary('warps')
  progress.done()

  # Second group-level calculation:
  #   - Calculate the mean of median b=0 values
  #   - Calculate the mean response function, and extract the l=0 values from it
  mean_median_bzero = sum_median_bzero / len(subjects)
  mean_RF = [[v/len(subjects) for v in line] for line in sum_RF]
  mean_RF_lzero = [line[0] for line in mean_RF]

  # Third pass through subject data in group analysis:
  #   - Scale the connectome strengths:
  #     - Multiply by SIFT proportionality coefficient mu
  #     - Multiply by (mean median b=0) / (subject median b=0)
  #     - Multiply by (subject RF size) / (mean RF size)
  #         (needs to account for multi-shell data)
  #   - Write the result to file
  progress = app.progressBar('Applying normalisation scaling to subject connectomes', len(subjects))
  run.function(os.makedirs, 'connectomes')
  for s in subjects:
    RF_lzero = [line[0] for line in s.RF]
    s.RF_multiplier = 1.0
    for (mean, subj) in zip(mean_RF_lzero, RF_lzero):
      s.RF_multiplier = s.RF_multiplier * subj / mean
    # Don't want to be scaling connectome independently for differences in RF l=0 terms across all shells;
    #   use the geometric mean of the per-shell scale factors
    s.RF_multiplier = math.pow(s.RF_multiplier, 1.0 / len(mean_RF_lzero))

    s.bzero_multiplier = mean_median_bzero / s.median_bzero

    s.global_multiplier = s.mu * s.bzero_multiplier * s.RF_multiplier

    connectome = [ ]
    with open(s.in_connectome, 'r') as f:
      for line in f:
        connectome.append([float(v) for v in line.split()])
    with open(s.temp_connectome, 'w') as f:
      for line in connectome:
        f.write(' '.join([str(v*s.global_multiplier) for v in line]) + '\n')
    progress.increment()
  progress.done()

  # Third group-level calculation: Generate the group mean connectome
  # For any higher-level analysis (e.g. NBSE, computing connectome global measures, etc.),
  #   trying to incorporate such analysis into this particular pipeline script is likely to
  #   overly complicate the interface, and not actually provide much in terms of
  #   convenience / reproducibility guarantees. The primary functionality of this group-level
  #   analysis is therefore to achieve inter-subject connection density normalisation; users
  #   then have the flexibility to subsequently analyse the data however they choose (ideally
  #   based on subject classification data provided with the BIDS-compliant dataset).
  progress = app.progressBar('Calculating group mean connectome', len(subjects)+1)
  mean_connectome = []
  for s in subjects:
    connectome = []
    with open(s.temp_connectome, 'r') as f:
      for line in f:
        connectome.append([float(v) for v in line.split()])
    if mean_connectome:
      mean_connectome = [[c1+c2 for c1, c2 in zip(r1, r2)] for r1, r2 in zip(mean_connectome, connectome)]
    else:
      mean_connectome = connectome
    progress.increment()

  mean_connectome = [[v/len(subjects) for v in row] for row in mean_connectome]
  progress.done()

  # Write results of interest back to the output directory;
  #   both per-subject and group information
  progress = app.progressBar('Writing results to output directory', len(subjects)+2)
  for s in subjects:
    run.function(shutil.copyfile, s.temp_connectome, s.out_connectome)
    with open(s.out_scale_bzero, 'w') as f:
      f.write(str(s.bzero_multiplier))
    with open(s.out_scale_RF, 'w') as f:
      f.write(str(s.RF_multiplier))
    progress.increment()

  with open(os.path.join(output_dir, 'mean_response.txt'), 'w') as f:
    for row in mean_RF:
      f.write(' '.join([str(v) for v in row]) + '\n')
  progress.increment()
  with open(os.path.join(output_dir, 'mean_connectome.csv'), 'w') as f:
    for row in mean_connectome:
      f.write(' '.join([str(v) for v in row]) + '\n')
  progress.done()

# End of runGroup() function








analysis_choices = [ 'participant', 'group' ]
parcellation_choices = [ 'aal', 'aal2', 'fs_2005', 'fs_2009', 'perry512', 'sri24' ]

app.init('Robert E. Smith (robert.smith@florey.edu.au)',
         'Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in MRtrix3')

# If running within a container, erase existing standard options, and fill with only desired options
if is_container:
  for option in reversed(app.cmdline._actions):
    app.cmdline._handle_conflict_resolve(None, [(option.option_strings[0],option)])
  # app.cmdline._action_groups[2] is "Standard options" that was created earlier
  app.cmdline._action_groups[2].add_argument('-d', '--debug', dest='debug', action='store_true', help='In the event of encountering an issue with the script, re-run with this flag set to provide more useful information to the developer')
  app.cmdline._action_groups[2].add_argument('-h', '--help', dest='help', action='store_true', help='Display help information for the script')
  app.cmdline._action_groups[2].add_argument('-n', '--n_cpus', type=int, metavar='number', dest='nthreads', help='Use this number of threads in MRtrix3 multi-threaded applications (0 disables multi-threading)')
  app.cmdline._action_groups[2].add_argument('-v', '--version', action='version', version=__version__)

app.cmdline.add_argument('bids_dir', help='The directory with the input dataset formatted according to the BIDS standard.')
app.cmdline.add_argument('output_dir', help='The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.')
app.cmdline.add_argument('analysis_level', help='Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: ' + ', '.join(analysis_choices), choices=analysis_choices)
batch_options = app.cmdline.add_argument_group('Options specific to the batch processing of subject data')
batch_options.add_argument(option_prefix + 'participant_label', nargs='+', help='The label(s) of the participant(s) that should be analyzed. The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a space-separated list.')
participant_options = app.cmdline.add_argument_group('Options that are relevant to participant-level analysis')
if not is_container:
  participant_options.add_argument(option_prefix + 'atlas_path', help='The path to search for an atlas parcellation (may be necessary when the script is executed outside of the BIDS App container')
participant_options.add_argument(option_prefix + 'output_verbosity', type=int, default=2, help='The verbosity of script output (number from 1 to 3); higher values result in more generated data being included in the output directory')
participant_options.add_argument(option_prefix + 'parcellation', help='The choice of connectome parcellation scheme (compulsory for participant-level analysis). Options are: ' + ', '.join(parcellation_choices), choices=parcellation_choices)
participant_options.add_argument(option_prefix + 'preprocessed', action='store_true', help='Indicate that the subject DWI data have been preprocessed, and hence initial image processing steps will be skipped (also useful for testing)')
participant_options.add_argument(option_prefix + 'streamlines', type=int, help='The number of streamlines to generate for each subject')

app.parse()


if app.isWindows():
  app.error('Script cannot be run on Windows due to FSL dependency')

if find_executable('bids-validator'):
  run.command('bids-validator ' + app.args.bids_dir)
else:
  app.warn('BIDS validator script not installed; proceeding without validation of input data')

# Running participant level
if app.args.analysis_level == 'participant':

  if app.args.output_verbosity < 1 or app.args.output_verbosity > 3:
    app.error('Valid values for ' + option_prefix + 'output_verbosity option are from 1 to 3')

  subjects_to_analyze = [ ]
  # Only run a subset of subjects
  if app.args.participant_label:
    subjects_to_analyze = [ 'sub-' + sub_index for sub_index in app.args.participant_label ]
    for subject_dir in subjects_to_analyze:
      if not os.path.isdir(os.path.join(app.args.bids_dir, subject_dir)):
        app.error('Unable to find directory for subject: ' + subject_dir)
  # Run all subjects sequentially
  else:
    subject_dirs = glob.glob(os.path.join(app.args.bids_dir, 'sub-*'))
    subjects_to_analyze = [ 'sub-' + directory.split("-")[-1] for directory in subject_dirs ]
    if not subjects_to_analyze:
      app.error('Could not find any subjects in BIDS directory')

  for subject_label in subjects_to_analyze:
    app.console('Commencing execution for subject ' + subject_label)
    runSubject(app.args.bids_dir, subject_label, os.path.abspath(app.args.output_dir))

# Running group level
elif app.args.analysis_level == 'group':

  if app.args.participant_label:
    app.error('Cannot use --participant_label option when performing group analysis')
  runGroup(os.path.abspath(app.args.output_dir))

app.complete()
