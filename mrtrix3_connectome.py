#!/usr/bin/env python

import glob, json, math, os, re, shutil, subprocess
from distutils.spawn import find_executable
from mrtrix3 import app, file, fsl, image, path, run


IS_CONTAINER = os.path.exists('/version') and os.path.exists('/mrtrix3_version')
__version__ = 'BIDS-App \'MRtrix3_connectome\' version {}'.format(open('/version').read()) if IS_CONTAINER else 'BIDS-App \'MRtrix3_connectome\' standalone'
OPTION_PREFIX = '--' if IS_CONTAINER else '-'


# If running for multiple subjects in a single command invocation,
#   these initialisation steps need only be performed once
class ParticipantShared(object): #pylint: disable=useless-object-inheritance
  def __init__(self, atlas_path, output_verbosity, parcellation, preprocessed, streamlines, template_reg):

    if not parcellation:
      app.error('For participant-level analysis, desired parcellation must be provided using the ' + OPTION_PREFIX + 'parcellation option')
    self.parcellation = parcellation

    self.output_verbosity = output_verbosity
    self.preprocessed = preprocessed
    self.streamlines = streamlines

    fsl_path = os.environ.get('FSLDIR', '')
    if not fsl_path:
      app.error('Environment variable FSLDIR is not set; please run appropriate FSL configuration script')

    robex_found = find_executable('ROBEX') and find_executable('runROBEX.sh')
    N4_found = find_executable('N4BiasFieldCorrection')

    if robex_found and N4_found:
      app.console('N4BiasFieldCorrection and ROBEX found; will use for bias field correction and brain extraction')
      self.brain_extraction_cmd = 'runROBEX.sh'
    else:
      if robex_found and not N4_found:
        app.console('N4BiasFieldCorrection not found; will use FSL for bias field correction & brain extraction')
      elif N4_found and not robex_found:
        app.console('ROBEX not found; will use FSL for brain extraction')
      else:
        app.console('N4BiasFieldCorrection and ROBEX not found; will use FSL for brain extraction')
      self.brain_extraction_cmd = fsl.exeName('fsl_anat')

    self.dwibiascorrect_algo = '-ants'
    if not N4_found:
      # Can't use fsl.exeName() here, since we want to proceed even if it's not found
      if find_executable('fast') or find_executable('fsl5.0-fast'):
        self.dwibiascorrect_algo = '-fsl'
        app.console('Could not find ANTs program N4BiasFieldCorrection; ' + \
                    'using FSL FAST for bias field correction')
      else:
        self.dwibiascorrect_algo = ''
        app.warn('Could not find ANTs program \'N4BiasFieldCorrection\' or FSL program \'fast\'; ' + \
                 'will proceed without performing initial DWI bias field correction')

    self.eddy_binary = fsl.eddyBinary(True)
    if self.eddy_binary:
      self.eddy_cuda = True
    else:
      self.eddy_binary = fsl.eddyBinary(False)
      self.eddy_cuda = False
    # Using run.command() to query the eddy help page will lead to a warning being
    #   issued due to a non-zero return code
    # Therefore use subprocess directly
    eddy_help_process = subprocess.Popen([ self.eddy_binary, '--help' ], stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    eddy_stderr = eddy_help_process.communicate()[1]
    self.eddy_repol = False
    self.eddy_mporder = False
    self.eddy_mbs = False
    for line in eddy_stderr.decode('utf-8', errors='replace').splitlines():
      line = line.lstrip()
      if line.startswith('--repol'):
        self.eddy_repol = True
      elif line.startswith('--mporder') and self.eddy_cuda:
        self.eddy_mporder = True
      elif line.startswith('--estimate_move_by_susceptibility'):
        self.eddy_mbs = True

    self.do_freesurfer = parcellation in [ 'brainnetome246fs', 'desikan', 'destrieux', 'hcpmmp1', 'yeo7fs', 'yeo17fs' ]
    self.do_mni = parcellation in [ 'aal', 'aal2', 'brainnetome246mni', 'craddock200', 'craddock400', 'perry512', 'yeo7mni', 'yeo17mni' ]
    if parcellation != 'none':
      assert self.do_freesurfer or self.do_mni

    if template_reg:
      if self.do_mni:
        self.template_registration_software = template_reg
      else:
        app.warn('Volumetric template registration not being performed; ' + OPTION_PREFIX + 'template_reg option ignored')
        self.template_registration_software = ''
    else:
      self.template_registration_software = 'ants' if self.do_mni else ''
    if self.template_registration_software == 'ants':
      if not find_executable('ANTS') or not find_executable('WarpImageMultiTransform'):
        app.error('Commands \'ANTS\' and \'WarpImageMultiTransform\' must be present in PATH to use ANTs software for template registration')
    elif self.template_registration_software == 'fsl':
      self.flirt_cmd = fsl.exeName('flirt')
      self.fnirt_cmd = fsl.exeName('fnirt')
      self.invwarp_cmd = fsl.exeName('invwarp')
      self.applywarp_cmd = fsl.exeName('applywarp')
      self.fnirt_config_basename = 'T1_2_MNI152_2mm.cnf'
      self.fnirt_config_path = os.path.join(fsl_path, 'etc', 'flirtsch', self.fnirt_config_basename)
      if not os.path.isfile(self.fnirt_config_path):
        app.error('Unable to find configuration file for FNI FNIRT (expected location: ' + self.fnirt_config_path + ')')

    self.template_image_path = ''
    self.template_mask_path = ''
    self.parc_image_path = ''
    self.parc_lut_file = ''
    self.mrtrix_lut_file = ''

    mrtrix_lut_dir = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(app.__file__)), \
                                                   os.pardir, \
                                                   os.pardir, \
                                                   'share', \
                                                   'mrtrix3', \
                                                   'labelconvert'))

    if self.do_freesurfer:
      self.freesurfer_home = os.environ.get('FREESURFER_HOME', None)
      if not self.freesurfer_home:
        app.error('Environment variable FREESURFER_HOME not set; please verify FreeSurfer installation')
      if not find_executable('recon-all'):
        app.error('Could not find FreeSurfer script recon-all; please verify FreeSurfer installation')
      self.freesurfer_subjects_dir = os.environ['SUBJECTS_DIR'] if 'SUBJECTS_DIR' in os.environ else os.path.join(self.freesurfer_home, 'subjects')
      if not os.path.isdir(self.freesurfer_subjects_dir):
        app.error('Could not find FreeSurfer subjects directory (expected location: ' + self.freesurfer_subjects_dir + ')')
      for subdir in [ 'fsaverage', 'fsaverage5', 'lh.EC_average', 'rh.EC_average' ]:
        if not os.path.isdir(os.path.join(self.freesurfer_subjects_dir, subdir)):
          app.error('Could not find requisite FreeSurfer subject directory \'' + subdir + '\' (expected location: ' + os.path.join(self.freesurfer_subjects_dir, subdir) + ')')
      self.reconall_path = find_executable('recon-all')
      if not self.reconall_path:
        app.error('Could not find FreeSurfer script recon-all; please verify FreeSurfer installation')
      if parcellation in [ 'hcpmmp1', 'yeo7fs', 'yeo17fs' ]:
        if parcellation == 'hcpmmp1':
          self.hcpmmp1_annot_paths = [ os.path.join(self.freesurfer_subjects_dir, 'fsaverage', 'label', hemi + 'h.HCPMMP1.annot') for hemi in [ 'l', 'r' ] ]
          if not all([os.path.isfile(path) for path in self.hcpmmp1_annot_paths]):
            app.error('Could not find necessary annotation labels for applying HCPMMP1 parcellation ' + \
                      '(expected location: ' + os.path.join(self.freesurfer_subjects_dir, 'fsaverage', 'label', '?h.HCPMMP1.annot') + ')')
        else: # yeo7fs, yeo17fs
          self.yeo_annot_paths = [ os.path.join(self.freesurfer_subjects_dir, 'fsaverage5', 'label', hemi + 'h.Yeo2011_' + ('7' if parcellation == 'yeo7fs' else '17') + 'Networks_N1000.split_components.annot') for hemi in [ 'l', 'r' ] ]
          if not all([os.path.isfile(path) for path in self.yeo_annot_paths]):
            app.error('Could not find necessary annotation labels for applying Yeo2011 parcellation ' + \
                      '(expected location: ' + os.path.join(self.freesurfer_subjects_dir, 'fsaverage5', 'label', '?h.Yeo2011_' + ('7' if parcellation == 'yeo7fs' else '17') + 'Networks_N1000.split_components.annot') + ')')
        for cmd in [ 'mri_surf2surf', 'mri_aparc2aseg' ]:
          if not find_executable(cmd):
            app.error('Could not find FreeSurfer command ' + cmd + ' (necessary for applying HCPMMP1 parcellation); please verify FreeSurfer installation')
      elif parcellation == 'brainnetome246fs':
        self.brainnetome_cortex_gcs_paths = [ os.path.join(self.freesurfer_home, 'average', hemi + 'h.BN_Atlas.gcs') for hemi in [ 'l', 'r' ] ]
        if not all([os.path.isfile(path) for path in self.brainnetome_cortex_gcs_paths]):
          app.error('Could not find necessary GCS files for applying Brainnetome cortical parcellation via FreeSurfer ' + \
                    '(expected location: ' + os.path.join(self.freesurfer_home, 'average', 'label', '?h.BN_Atlas.gcs') + ')')
        self.brainnetome_sgm_gca_path = os.path.join(self.freesurfer_home, 'average', 'BN_Atlas_subcortex.gca')
        if not os.path.isfile(self.brainnetome_sgm_gca_path):
          app.error('Could not find necessary GCA file for applying Brainnetome sub-cortical parcellation via FreeSurfer ' + \
                    '(expected location: ' + self.brainnetome_sgm_gca_path + ')')
        for cmd in [ 'mri_label2vol', 'mri_ca_label', 'mris_ca_label' ]:
          if not find_executable(cmd):
            app.error('Could not find FreeSurfer command ' + cmd + ' (necessary for applying Brainnetome parcellation); please verify FreeSurfer installation')

      # Query contents of recon-all script, looking for "-openmp" and "-parallel" occurences
      # Add options to end of recon-all -all call, based on which of these options are available
      #   as well as the value of app.numThreads
      # - In 5.3.0, just the -openmp option is available
      # - In 6.0.0, -openmp needs to be preceded by -parallel
      self.reconall_multithread_options = []
      if app.numThreads is None or app.numThreads:
        with open(self.reconall_path, 'r') as f:
          reconall_text = f.read().splitlines()
        for line in reconall_text:
          line = line.strip()
          if line == 'case "-parallel":':
            self.reconall_multithread_options = [ '-parallel' ] + self.reconall_multithread_options
          # If number of threads in this script is not being explicitly controlled,
          #   allow recon-all to use its own default number of threads
          elif line == 'case "-openmp":' and app.numThreads is not None:
            self.reconall_multithread_options.extend([ '-openmp', str(app.numThreads) ])
      if self.reconall_multithread_options:
        self.reconall_multithread_options = ' ' + ' '.join(self.reconall_multithread_options)
      else:
        self.reconall_multithread_options = ''
      app.debug(self.reconall_multithread_options)

      if parcellation == 'brainnetome246fs':
        self.parc_lut_file = os.path.join(self.freesurfer_home, 'BN_Atlas_246_LUT.txt')
        self.mrtrix_lut_file = ''
      elif parcellation == 'desikan':
        self.parc_lut_file = os.path.join(self.freesurfer_home, 'FreeSurferColorLUT.txt')
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'fs_default.txt')
      elif parcellation == 'destrieux':
        self.parc_lut_file = os.path.join(self.freesurfer_home, 'FreeSurferColorLUT.txt')
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'fs_a2009s.txt')
      elif parcellation == 'hcpmmp1':
        self.parc_lut_file = os.path.join(mrtrix_lut_dir, 'hcpmmp1_original.txt')
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'hcpmmp1_ordered.txt')
      elif parcellation in [ 'yeo7fs', 'yeo17fs' ]:
        self.parc_lut_file = os.path.join(self.freesurfer_home, 'Yeo2011_' + ('7' if parcellation == 'yeo7fs' else '17') + 'networks_Split_Components_LUT.txt')
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'Yeo2011_' + ('7' if parcellation == 'yeo7fs' else '17') + 'N_split.txt')
      else:
        assert False

      # If running in a container environment, and --debug is used (resulting in the
      #   scratch directory being a mounted drive), it's possible that attempting to
      #   construct a softlink may lead to an OSError
      # As such, run a test to determine whether or not it is possible to construct
      #   a softlink within the scratch directory; if it is not possible, revert to
      #   performing deep copies of the relevant FreeSurfer template directories
      self.freesurfer_template_link_function = os.symlink
      try:
        self.freesurfer_template_link_function(self.freesurfer_subjects_dir, 'test_softlink')
        file.delTemporary('test_softlink')
        app.debug('Using softlinks to FreeSurfer template directories')
      except OSError:
        app.debug('Unable to create softlinks; will perform deep copies of FreeSurfer template directories')
        self.freesurfer_template_link_function = shutil.copytree

    elif self.do_mni:
      self.template_image_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_2mm.nii.gz')
      self.template_mask_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_2mm_brain_mask.nii.gz')
      if parcellation == 'aal':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.nii'))
        self.parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.txt'))
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'aal.txt')
      elif parcellation == 'aal2':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.nii'))
        self.parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.txt'))
        self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'aal2.txt')
      elif parcellation == 'brainnetome246mni':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'brainnetome', 'BN_Atlas_246_1mm.nii.gz'))
        self.parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'brainnetome', 'BN_Atlas_246_LUT.txt'))
        self.mrtrix_lut_file = ''
      elif parcellation == 'craddock200':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'ADHD200_parcellate_200.nii.gz'))
        self.parc_lut_file = self.mrtrix_lut_file = ''
      elif parcellation == 'craddock400':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'ADHD200_parcellate_400.nii.gz'))
        self.parc_lut_file = self.mrtrix_lut_file = ''
      elif parcellation == 'perry512':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', '512inMNI.nii'))
        self.parc_lut_file = self.mrtrix_lut_file = ''
      elif parcellation == 'yeo7mni':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'Yeo2011', 'Yeo2011_7Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz'))
        self.parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'Yeo2011', '7Networks_ColorLUT_freeview.txt'))
        self.mrtrix_lut_file = ''
      elif parcellation == 'yeo17mni':
        self.parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'Yeo2011', 'Yeo2011_17Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz'))
        self.parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'Yeo2011', '17Networks_ColorLUT_freeview.txt'))
        self.mrtrix_lut_file = ''
      else:
        assert False

    def findAtlasFile(filepath, description):
      if not filepath:
        return ''
      if os.path.isfile(filepath):
        return filepath
      if not atlas_path:
        app.error('Could not find ' + description + ' (expected location: ' + filepath + ')')
      newpath = os.path.join(os.path.dirname(atlas_path), os.path.basename(filepath))
      if os.path.isfile(newpath):
        return newpath
      app.error('Could not find ' + description + ' (tested locations: \'' + filepath + '\', \'' + newpath + '\')')
      return None

    self.template_image_path = findAtlasFile(self.template_image_path, 'template image')
    self.template_mask_path = findAtlasFile(self.template_mask_path, 'template brain mask image')
    self.parc_image_path = findAtlasFile(self.parc_image_path, 'parcellation image')
    self.parc_lut_file = findAtlasFile(self.parc_lut_file, 'parcellation lookup table file')

    if self.mrtrix_lut_file and not os.path.exists(self.mrtrix_lut_file):
      app.error('Could not find MRtrix3 connectome lookup table file (expected location: ' + self.mrtrix_lut_file + ')')







def runSubject(bids_dir, label, shared, output_prefix):

  output_dir = os.path.join(output_prefix, label)
  if os.path.exists(output_dir):
    app.warn('Output directory for subject \'' + label + '\' already exists; contents will be erased when this execution completes')

  app.makeTempDir()

  # Need to perform an initial import of JSON data using mrconvert; so let's grab the diffusion gradient table as well
  # If no bvec/bval present, need to go down the directory listing
  # Only try to import JSON file if it's actually present
  #   direction in the acquisition they'll need to be split across multiple files
  # May need to concatenate more than one input DWI, since if there's more than one phase-encode direction
  #   in the acquired DWIs (i.e. not just those used for estimating the inhomogeneity field), they will
  #   need to be stored as separate NIfTI files in the 'dwi/' directory.
  app.console('Importing DWI data into temporary directory')
  dwi_image_list = glob.glob(os.path.join(bids_dir, label, 'dwi', label + '*_dwi.nii*'))
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
    json_import_option = ''
    if os.path.isfile(json_path):
      json_import_option = ' -json_import ' + json_path
    elif not shared.preprocessed:
      app.error('No sidecar JSON file found for image \'' + entry + '\'; cannot proceed with DWI preprocessing without this information')
    run.command('mrconvert ' + entry + grad_import_option + json_import_option
                + ' ' + path.toTemp('dwi' + str(dwi_index) + '.mif', True))
    dwi_index += 1

  # Go hunting for reversed phase-encode data dedicated to field map estimation
  fmap_image_list = []
  fmap_dir = os.path.join(bids_dir, label, 'fmap')
  fmap_index = 1
  if os.path.isdir(fmap_dir):
    if shared.preprocessed:
      app.error('fmap/ directory detected for subject \'' + label + '\' despite use of ' + OPTION_PREFIX + 'preprocessed option; this directory should not be present if DWIs have already been pre-processed')
    app.console('Importing fmap data into temporary directory')
    fmap_image_list = glob.glob(os.path.join(fmap_dir, label + '*_dir-*_epi.nii*'))
    for entry in fmap_image_list:
      prefix = entry.split(os.extsep)[0]
      json_path = prefix + '.json'
      with open(json_path, 'r') as f:
        json_elements = json.load(f)
      if 'IntendedFor' in json_elements and not any(i.endswith(json_elements['IntendedFor']) for i in dwi_image_list):
        app.console('Image \'' + entry + '\' is not intended for use with DWIs; skipping')
        continue
      if not os.path.isfile(json_path):
        app.error('No sidecar JSON file found for image \'' + entry + '\'')
      # fmap files will not come with any gradient encoding in the JSON;
      #   therefore we need to add it manually ourselves so that mrcat / mrconvert can
      #   appropriately handle the table once these images are concatenated with the DWIs
      fmap_image_size = image.Header(entry).size()
      fmap_image_num_volumes = 1 if len(fmap_image_size) == 3 else fmap_image_size[3]
      fmap_dwscheme_file = 'fmap' + str(fmap_index) + '.b'
      with open(path.toTemp(fmap_dwscheme_file, False), 'w') as f:
        for _ in range(0, fmap_image_num_volumes):
          f.write('0,0,1,0\n')
      run.command('mrconvert ' + entry + ' -json_import ' + json_path + ' ' + \
                  '-grad ' + path.toTemp(fmap_dwscheme_file, True) + ' ' + \
                  path.toTemp('fmap' + str(fmap_index) + '.mif', True))
      file.delTemporary(fmap_dwscheme_file)
      fmap_index += 1

    fmap_image_list = [ 'fmap' + str(index) + '.mif' for index in range(1, fmap_index) ]
  # If there's no data in fmap/ directory, need to check to see if there's any phase-encoding
  #   contrast within the input DWI(s)
  elif len(dwi_image_list) < 2 and not shared.preprocessed:
    app.error('Inadequate data for pre-processing of subject \'' + label + '\': No phase-encoding contrast in input DWIs, and no fmap/ directory, so EPI distortion correction cannot be performed')

  dwi_image_list = [ 'dwi' + str(index) + '.mif' for index in range(1, dwi_index) ]

  # Import anatomical image
  app.console('Importing T1 image into temporary directory')
  t1w_image_list = glob.glob(os.path.join(bids_dir, label, 'anat', label + '*_T1w.nii*'))
  if len(t1w_image_list) > 1:
    app.error('More than one T1-weighted image found for subject ' + label + '; script not yet compatible with this')
  elif not t1w_image_list:
    app.error('No T1-weighted image found for subject ' + label)
  run.command('mrconvert ' + t1w_image_list[0] + ' ' + path.toTemp('T1.mif', True))

  cwd = os.getcwd()
  app.gotoTempDir()

  dwipreproc_se_epi = ''
  dwipreproc_se_epi_option = ''

  # For automated testing, down-sampled images are used. However, this invalidates the requirements of
  #   both MP-PCA denoising and Gibbs ringing removal. In addition, eddy can still take a long time
  #   despite the down-sampling. Therefore, provide images that have been pre-processed to the stage
  #   where it is still only DWI, JSON & bvecs/bvals that need to be provided.
  if shared.preprocessed:

    if len(dwi_image_list) > 1:
      app.error('If DWIs have been pre-processed, then only a single DWI file should be present')
    app.console('Skipping MP-PCA denoising, Gibbs ringing removal, distortion correction and bias field correction due to use of ' + OPTION_PREFIX + 'preprocessed option')
    run.function(os.rename, dwi_image_list[0], 'dwi.mif')

  else: # Do initial image pre-processing (denoising, Gibbs ringing removal if available, distortion correction & bias field correction) as normal

    # Concatenate any SE EPI images with the DWIs before denoising (& unringing), then
    #   separate them again after the fact
    # TODO Note however that this may not be possible: the fmap/ images may not lie on the
    #   same grid as the DWIs. When this occurs, cannot apply denoising algorithm to fmap images.
    if fmap_image_list:
      app.console('Concatenating DWI and fmap data for combined pre-processing')
      dwidenoise_input = 'dwi_fmap_cat.mif'
      if len(fmap_image_list) > 1:
        run.command('mrcat ' + ' '.join(fmap_image_list) + ' fmap_cat.mif -axis 3')
        fmap_header = image.Header('fmap_cat.mif')
        fmap_num_volumes = fmap_header.size()[3]
        run.command('mrcat fmap_cat.mif ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
        file.delTemporary('fmap_cat.mif')
      else:
        fmap_header = image.Header(fmap_image_list[0])
        fmap_num_volumes = 1 if len(fmap_header.size()) < 4 else fmap_header.size()[3]
        run.command('mrcat ' + fmap_image_list[0] + ' ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
      for i in fmap_image_list:
        file.delTemporary(i)
      for i in dwi_image_list:
        file.delTemporary(i)
    else:
      fmap_num_volumes = 0
      # Even if no explicit fmap images, may still need to concatenate multiple DWI inputs
      if len(dwi_image_list) > 1:
        app.console('Concatenating input DWI series')
        dwidenoise_input = 'dwi_cat.mif'
        run.command('mrcat ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
        for i in dwi_image_list:
          file.delTemporary(i)
      else:
        dwidenoise_input = dwi_image_list[0]


    # Step 1: Denoise
    app.console('Performing MP-PCA denoising of DWI' + (' and fmap' if fmap_num_volumes else '') + ' data')
    mrdegibbs_input = os.path.splitext(dwidenoise_input)[0] + '_denoised.mif'
    run.command('dwidenoise ' + dwidenoise_input + ' ' + mrdegibbs_input)
    file.delTemporary(dwidenoise_input)

    # Step 2: Gibbs ringing removal
    app.console('Performing Gibbs ringing removal for DWI' + (' and fmap' if fmap_num_volumes else '') + ' data')
    mrdegibbs_output = os.path.splitext(mrdegibbs_input)[0] + '_degibbs.mif'
    run.command('mrdegibbs ' + mrdegibbs_input + ' ' + mrdegibbs_output + ' -nshifts 50')
    file.delTemporary(mrdegibbs_input)

    # If fmap images and DWIs have been concatenated, now is the time to split them back apart
    if fmap_num_volumes:
      app.console('Separating DWIs and fmap images from concatenated series')
      dwipreproc_input = 'dwipreproc_in.mif'
      dwipreproc_se_epi = 'se_epi.mif'
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_se_epi + ' -coord 3 0:' + str(fmap_num_volumes-1))
      cat_num_volumes = image.Header(mrdegibbs_output).size()[3]
      app.var(cat_num_volumes)
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_input + ' -coord 3 ' + str(fmap_num_volumes) + ':' + str(cat_num_volumes-1))
      file.delTemporary(mrdegibbs_output)
      dwipreproc_se_epi_option = ' -se_epi ' + dwipreproc_se_epi + ' -align_seepi'
    else:
      dwipreproc_input = mrdegibbs_output
      dwipreproc_se_epi = None
      dwipreproc_se_epi_option = ''

    # Step 3: Distortion correction
    app.console('Performing various geometric corrections of DWIs')
    dwipreproc_input_header = image.Header(dwipreproc_input)
    have_slice_timing = 'SliceTiming' in dwipreproc_input_header.keyval()
    app.var(have_slice_timing)
    mb_factor = int(dwipreproc_input_header.keyval().get('MultibandAccelerationFactor', '1'))
    app.var(mb_factor)
    if 'SliceDirection' in dwipreproc_input_header.keyval():
      slice_direction_code = dwipreproc_input_header.keyval()['SliceDirection']
      if 'i' in slice_direction_code:
        num_slices = dwipreproc_input_header.size()[0]
      elif 'j' in slice_direction_code:
        num_slices = dwipreproc_input_header.size()[1]
      elif 'k' in slice_direction_code:
        num_slices = dwipreproc_input_header.size()[2]
      else:
        num_slices = dwipreproc_input_header.size()[2]
        app.warn('Error reading BIDS field \'SliceDirection\' (value: \'' + slice_direction_code + '\'); assuming third axis')
    else:
      num_slices = dwipreproc_input_header.size()[2]
    app.var(num_slices)
    mporder = 1 + num_slices/(mb_factor*4)
    app.var(mporder)

    eddy_binary = fsl.eddyBinary(True)
    if eddy_binary:
      eddy_cuda = True
    else:
      eddy_binary = fsl.eddyBinary(False)
      eddy_cuda = False
    app.var(eddy_binary, eddy_cuda)

    eddy_options = []
    if shared.eddy_repol:
      eddy_options.append('--repol')
    if shared.eddy_mporder and have_slice_timing:
      eddy_options.append('--mporder=' + str(mporder))
    if shared.eddy_mbs:
      eddy_options.append('--estimate_move_by_susceptibility')
    dwipreproc_eddy_option = ' -eddy_options \" ' + ' '.join(eddy_options) + '\"' if eddy_options else ''
    run.function(os.makedirs, 'eddyqc')
    run.command('dwipreproc ' + dwipreproc_input + ' dwi_preprocessed.mif -rpe_header' + dwipreproc_se_epi_option + dwipreproc_eddy_option + ' -eddyqc_text eddyqc/')
    file.delTemporary(dwipreproc_input)
    if dwipreproc_se_epi:
      file.delTemporary(dwipreproc_se_epi)

    # Step 4: Bias field correction
    if shared.dwibiascorrect_algo:
      app.console('Performing initial B1 bias field correction of DWIs')
      run.command('dwibiascorrect dwi_preprocessed.mif dwi.mif ' + shared.dwibiascorrect_algo)
      file.delTemporary('dwi_preprocessed.mif')
    else:
      run.function(shutil.move, 'dwi_preprocessed.mif', 'dwi.mif')

  # No longer branching based on whether or not -preprocessed was specified

  # Step 5: Generate a brain mask for DWI
  #   Also produce a dilated version of the mask for later use
  app.console('Estimating a brain mask for DWIs')
  run.command('dwi2mask dwi.mif dwi_mask.mif')
  # TODO Crop DWIs based on brain mask
  run.command('maskfilter dwi_mask.mif dilate dwi_mask_dilated.mif -npass 3')

  # Step 6: Generate an FA image
  # This is required for group-level analysis
  app.console('Generating FA image for group-level analysis')
  run.command('dwi2tensor dwi.mif -mask dwi_mask.mif - | tensor2metric - -fa fa.mif')

  # Step 7: Estimate response functions for spherical deconvolution
  app.console('Estimating tissue response functions for spherical deconvolution')
  run.command('dwi2response dhollander dwi.mif response_wm.txt response_gm.txt response_csf.txt -mask dwi_mask.mif')

  # Determine whether we are working with single-shell or multi-shell data
  bvalues = [int(round(float(value))) for value in image.mrinfo('dwi.mif', 'shell_bvalues').strip().split()]
  multishell = (len(bvalues) > 2)

  # Step 8: Perform spherical deconvolution
  #          Use a dilated mask for spherical deconvolution as a 'safety margin' -
  #          ACT should be responsible for stopping streamlines before they reach the edge of the DWI mask
  app.console('Estimating' + ('' if multishell else ' white matter') + ' Fibre Orientation Distribution' + ('s' if multishell else ''))
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

  # Step 9: Perform brain extraction and bias field correction on the T1 image
  #         in its original space (this is necessary for histogram matching
  #         prior to registration)
  app.console('Performing brain extraction and B1 bias field correction of T1 image')
  T1_header = image.Header('T1.mif')
  T1_revert_strides_option = ' -strides ' + ','.join([str(i) for i in T1_header.strides()])
  if shared.brain_extraction_cmd == 'runROBEX.sh':
    # Do a semi-iterative approach here: Get an initial brain mask, use that
    #   mask to estimate a bias field, then re-compute the brain mask
    run.command('mrconvert T1.mif T1.nii -strides +1,+2,+3')
    run.command(shared.brain_extraction_cmd + ' T1.nii T1_initial_brain.nii T1_initial_mask.nii')
    file.delTemporary('T1_initial_brain.nii')
    run.command('N4BiasFieldCorrection -i T1.nii -w T1_initial_mask.nii -o T1_biascorr.nii')
    file.delTemporary('T1.nii')
    file.delTemporary('T1_initial_mask.nii')
    run.command(shared.brain_extraction_cmd + ' T1_biascorr.nii T1_biascorr_brain.nii T1_biascorr_brain_mask.nii')
    file.delTemporary('T1_biascorr.nii')
    run.command('mrconvert T1_biascorr_brain.nii T1_biascorr_brain.mif' + T1_revert_strides_option)
    file.delTemporary('T1_biascorr_brain.nii')
    run.command('mrconvert T1_biascorr_brain_mask.nii T1_mask.mif -datatype bit' + T1_revert_strides_option)
    file.delTemporary('T1_biascorr_brain_mask.nii')
  else:
    run.command('mrconvert T1.mif T1.nii -strides -1,+2,+3')
    run.command(shared.brain_extraction_cmd + ' -i T1.nii --noseg --nosubcortseg')
    file.delTemporary('T1.nii')
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain') + ' T1_biascorr_brain.mif' + T1_revert_strides_option)
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain_mask') + ' T1_mask.mif -datatype bit' + T1_revert_strides_option)
    file.delTemporary('T1.anat')

  # Step 10: Generate target images for T1->DWI registration
  app.console('Generating contrast-matched images for inter-modal registration between DWIs and T1')
  run.command('dwiextract dwi.mif -bzero - | '
              'mrcalc - 0.0 -max - | '
              'mrmath - mean -axis 3 dwi_meanbzero.mif')
  run.command('mrcalc 1 dwi_meanbzero.mif -div dwi_mask.mif -mult - | '
              'mrhistmatch nonlinear - T1_biascorr_brain.mif dwi_pseudoT1.mif -mask_input dwi_mask.mif -mask_target T1_mask.mif')
  run.command('mrcalc 1 T1_biascorr_brain.mif -div T1_mask.mif -mult - | '
              'mrhistmatch nonlinear - dwi_meanbzero.mif T1_pseudobzero.mif -mask_input T1_mask.mif -mask_target dwi_mask.mif')

  # Step 11: Perform T1->DWI registration
  #          Note that two registrations are performed: Even though we have a symmetric registration,
  #          generation of the two histogram-matched images means that you will get slightly different
  #          answers depending on which synthesized image & original image you use
  app.console('Performing registration between DWIs and T1')
  run.command('mrregister T1_biascorr_brain.mif dwi_pseudoT1.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_T1_to_pseudoT1.txt')
  file.delTemporary('T1_biascorr_brain.mif')
  run.command('mrregister T1_pseudobzero.mif dwi_meanbzero.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_pseudobzero_to_bzero.txt')
  run.command('transformcalc rigid_T1_to_pseudoT1.txt rigid_pseudobzero_to_bzero.txt average rigid_T1_to_dwi.txt')
  file.delTemporary('rigid_T1_to_pseudoT1.txt')
  file.delTemporary('rigid_pseudobzero_to_bzero.txt')
  run.command('mrtransform T1.mif T1_registered.mif -linear rigid_T1_to_dwi.txt')
  file.delTemporary('T1.mif')
  # Note: Since we're using a mask from fsl_anat (which crops the FoV), but using it as input to 5ttge fsl
  #   (which is receiving the raw T1), we need to resample in order to have the same dimensions between these two
  run.command('mrtransform T1_mask.mif T1_mask_registered.mif -linear rigid_T1_to_dwi.txt -template T1_registered.mif -interp nearest -datatype bit')
  file.delTemporary('T1_mask.mif')

  # Step 12: Generate 5TT image for ACT
  app.console('Generating five-tissue-type (5TT) image for Anatomically-Constrained Tractography (ACT)')
  run.command('5ttgen fsl T1_registered.mif 5TT.mif -mask T1_mask_registered.mif')
  if shared.output_verbosity > 1:
    run.command('5tt2vis 5TT.mif vis.mif')

  # Step 13: Generate the grey matter parcellation
  #          The necessary steps here will vary significantly depending on the parcellation scheme selected

  if shared.do_freesurfer:

    app.console('Getting grey matter parcellation in subject space using FreeSurfer')

    run.command('mrconvert T1_registered.mif T1_registered.nii -strides +1,+2,+3')

    # Since we're instructing recon-all to use a different subject directory, we need to
    #   construct softlinks to a number of directories provided by FreeSurfer that
    #   recon-all will expect to find in the same directory as the overridden subject path
    subdirs = [ 'fsaverage', 'lh.EC_average', 'rh.EC_average' ]
    if shared.parcellation in [ 'yeo7fs', 'yeo17fs' ]:
      subdirs.append('fsaverage5')
    for subdir in subdirs:
      run.function(shared.freesurfer_template_link_function, os.path.join(shared.freesurfer_subjects_dir, subdir), subdir)

    # Run FreeSurfer pipeline on this subject's T1 image
    run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -i T1_registered.nii')
    run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -all' + shared.reconall_multithread_options)

    # Grab the relevant parcellation image and target lookup table for conversion
    parc_image_path = os.path.join('freesurfer', 'mri')
    if shared.parcellation == 'desikan':
      parc_image_path = os.path.join(parc_image_path, 'aparc+aseg.mgz')
    elif shared.parcellation == 'destrieux':
      parc_image_path = os.path.join(parc_image_path, 'aparc.a2009s+aseg.mgz')
    else:
      # Non-standard parcellations are not applied as part of the recon-all command;
      #   need to explicitly map them to the subject
      # This requires SUBJECTS_DIR to be set; commands don't have a corresponding -sd option like recon-all
      run._env['SUBJECTS_DIR'] = app.tempDir
      if shared.parcellation == 'brainnetome246fs':
        for index, hemi in enumerate([ 'l', 'r' ]):
          run.command('mris_ca_label -l ' + os.path.join('freesurfer', 'label', hemi + 'h.cortex.label') + ' ' + \
                      'freesurfer ' + hemi + 'h ' + \
                      os.path.join('freesurfer', 'surf', hemi + 'h.sphere.reg') + ' ' + \
                      shared.brainnetome_cortex_gcs_paths[index] + ' ' + \
                      os.path.join('freesurfer', 'label', hemi + 'h.BN_Atlas.annot'))
          run.command('mri_label2vol --annot ' + os.path.join('freesurfer', 'label', hemi + 'h.BN_Atlas.annot') + ' ' + \
                      '--temp ' + os.path.join('freesurfer', 'mri', 'brain.mgz') + ' ' + \
                      '--o ' + os.path.join('freesurfer', 'mri', hemi + 'h.BN_Atlas.mgz') + ' ' + \
                      '--subject freesurfer --hemi ' + hemi + 'h --identity --proj frac 0 1 .1')
        run.command('mri_ca_label ' + os.path.join('freesurfer', 'mri', 'brain.mgz') + ' ' + \
                    os.path.join('freesurfer', 'mri', 'transforms', 'talairach.m3z') + ' ' + \
                    shared.brainnetome_sgm_gca_path + ' ' + \
                    os.path.join('freesurfer', 'mri', 'BN_Atlas_subcortex.mgz'))
        parc_image_path = os.path.join(parc_image_path, 'aparc.BN_Atlas+aseg.mgz')
        # Need to deal with prospect of overlapping mask labels
        # - Any overlap between the two hemisphere ribbons = set to zero
        # - Any overlap between cortex and sub-cortical = retain cortex
        run.command('mrcalc ' + ' '.join([ os.path.join('freesurfer', 'mri', hemi + 'h.BN_Atlas.mgz') for hemi in [ 'l', 'r' ] ]) + ' ' + \
                    '-mult cortex_overlap.mif -datatype bit')
        run.command('mrcalc ' + ' '.join([ os.path.join('freesurfer', 'mri', hemi + 'h.BN_Atlas.mgz') for hemi in [ 'l', 'r' ] ]) + ' ' + \
                    '-add ' + os.path.join('freesurfer', 'mri', 'BN_Atlas_subcortex.mgz') + ' -mult sgm_overlap.mif -datatype bit')
        run.command('mrcalc ' + ' '.join([ os.path.join('freesurfer', 'mri', hemi + 'h.BN_Atlas.mgz') for hemi in [ 'l', 'r' ] ]) + ' ' + \
                    '-add 1.0 cortex_overlap.mif -sub -mult ' + \
                    os.path.join('freesurfer', 'mri', 'BN_Atlas_subcortex.mgz') + ' 1.0 sgm_overlap.mif -sub -mult ' + \
                    '-add ' + parc_image_path)
        file.delTemporary('cortex_overlap.mif')
        file.delTemporary('sgm_overlap.mif')

        # TODO:
        # - Add cerebellar hemispheres in explicitly
        # - Remove these data from the resulting connectome matrix explicitly in Python
        #
        # - Use 274 combined image to obtain cerebellar subdivisions?
        #   Will require executing both FreeSurfer and MNI registration
        #   Will also require generation of a new LUT

      elif shared.parcellation == 'hcpmmp1':
        parc_image_path = os.path.join(parc_image_path, 'aparc.HCPMMP1+aseg.mgz')
        for index, hemi in enumerate([ 'l', 'r' ]):
          run.command('mri_surf2surf --srcsubject fsaverage --trgsubject freesurfer --hemi ' + hemi + 'h ' + \
                      '--sval-annot ' + shared.hcpmmp1_annot_paths[index] + ' ' + \
                      '--tval ' + os.path.join('freesurfer', 'label', hemi + 'h.HCPMMP1.annot'))
        run.command('mri_aparc2aseg --s freesurfer --old-ribbon --annot HCPMMP1 --o ' + parc_image_path)
      elif shared.parcellation in [ 'yeo7fs', 'yeo17fs' ]:
        num = '7' if shared.parcellation == 'yeo7fs' else '17'
        parc_image_path = os.path.join(parc_image_path, 'aparc.Yeo' + num + '+aseg.mgz')
        for index, hemi in enumerate([ 'l', 'r' ]):
          run.command('mri_surf2surf --srcsubject fsaverage5 --trgsubject freesurfer --hemi ' + hemi + 'h ' + \
                      '--sval-annot ' + shared.yeo_annot_paths[index] + ' ' + \
                      '--tval ' + os.path.join('freesurfer', 'label', hemi + 'h.Yeo' + num + '.annot'))
        run.command('mri_aparc2aseg --s freesurfer --old-ribbon --annot Yeo' + num + ' --o ' + parc_image_path)
      else:
        assert False

    if shared.mrtrix_lut_file:
      # If necessary:
      # Perform the index conversion
      run.command('labelconvert ' + parc_image_path + ' ' + shared.parc_lut_file + ' ' + shared.mrtrix_lut_file + ' parc_init.mif')
      # Fix the sub-cortical grey matter parcellations using FSL FIRST
      run.command('labelsgmfix parc_init.mif T1_registered.mif ' + shared.mrtrix_lut_file + ' parc.mif')
      file.delTemporary('parc_init.mif')
    else:
      # Non-standard sub-cortical parcellation; labelsgmfix not applicable
      run.command('mrconvert ' + parc_image_path + ' parc.mif -datatype uint32')
    if app.cleanup:
      run.function(shutil.rmtree, 'freesurfer')


  elif shared.do_mni:

    app.console('Registering to MNI template and transforming grey matter parcellation back to subject space')

    # Use non-dilated brain masks for performing histogram matching & linear registration
    run.command('mrhistmatch linear T1_registered.mif -mask_input T1_mask_registered.mif ' + shared.template_image_path + ' -mask_target ' + shared.template_mask_path + ' - | ' \
                'mrconvert - T1_registered_histmatch.nii -strides -1,+2,+3')

    assert shared.template_registration_software
    if shared.template_registration_software == 'ants':

      # Use ANTs SyN for registration to template
      # From Klein et al., NeuroImage 2009:
      # Warp: ANTS 3 -m PR[htargeti.nii, hsourcei.nii, 1, 2] -o houtput transformi.nii -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 -use-Histogram-Matching
      # Reslice: WarpImageMultiTransform 3 hlabeled sourcei.nii houtput labelsi.nii -R htargeti.nii htransformiWarp.nii htransformiAffine.txt --use-NN
      run.command('ANTS 3 -m PR[' + shared.template_image_path + ', T1_registered_histmatch.nii, 1, 2] -o ANTS -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching')
      transformed_atlas_path = 'atlas_transformed.nii'
      run.command('WarpImageMultiTransform 3 ' + shared.parc_image_path + ' ' + transformed_atlas_path + ' -R T1_registered_histmatch.nii -i ANTSAffine.txt ANTSInverseWarp.nii --use-NN')
      file.delTemporary(glob.glob('ANTSWarp.nii*')[0])
      file.delTemporary(glob.glob('ANTSInverseWarp.nii*')[0])
      file.delTemporary('ANTSAffine.txt')

    elif shared.template_registration_software == 'fsl':

      # Subject T1, brain masked; for flirt -in
      run.command('mrcalc T1_registered_histmatch.nii T1_mask_registered.mif -mult T1_registered_histmatch_masked.nii')
      # Template T1, brain masked; for flirt -ref
      run.command('mrcalc ' + shared.template_image_path + ' ' + shared.template_mask_path + ' -mult template_masked.nii')
      # Now have data required to run flirt
      run.command(shared.flirt_cmd + ' -ref template_masked.nii -in T1_registered_histmatch_masked.nii -omat T1_to_template.mat -dof 12 -cost leastsq')
      file.delTemporary('T1_registered_histmatch_masked.nii')
      file.delTemporary('template_masked.nii')

      # Use dilated brain masks for non-linear registration to mitigate mask edge effects
      # Subject T1, unmasked; for fnirt --in
      run.command('mrconvert T1_registered.mif T1_registered.nii -strides -1,+2,+3')
      # Subject brain mask, dilated; for fnirt --inmask
      run.command('maskfilter T1_mask_registered.mif dilate - -npass 3 | ' \
                  'mrconvert - T1_mask_registered_dilated.nii -strides -1,+2,+3')
      # Template brain mask, dilated; for fnirt --refmask
      run.command('maskfilter ' + shared.template_mask_path + ' dilate template_mask_dilated.nii -npass 3')
      # Now have data required to run fnirt
      run.command(shared.fnirt_cmd + ' --config=' + shared.fnirt_config_basename + ' --ref=' + shared.template_image_path + ' --in=T1_registered_histmatch.nii ' \
                  '--aff=T1_to_template.mat --refmask=template_mask_dilated.nii --inmask=T1_mask_registered_dilated.nii ' \
                  '--cout=T1_to_template_warpcoef.nii')
      file.delTemporary('T1_mask_registered_dilated.nii')
      file.delTemporary('template_mask_dilated.nii')
      file.delTemporary('T1_to_template.mat')
      fnirt_warp_subject2template_path = fsl.findImage('T1_to_template_warpcoef')

      # Use result of registration to transform atlas parcellation to subject space
      run.command(shared.invwarp_cmd + ' --ref=T1_registered.nii --warp=' + fnirt_warp_subject2template_path + ' --out=template_to_T1_warpcoef.nii')
      file.delTemporary(fnirt_warp_subject2template_path)
      fnirt_warp_template2subject_path = fsl.findImage('template_to_T1_warpcoef')
      run.command(shared.applywarp_cmd + ' --ref=T1_registered.nii --in=' + shared.parc_image_path + ' --warp=' + fnirt_warp_template2subject_path + ' --out=atlas_transformed.nii --interp=nn')
      file.delTemporary(fnirt_warp_template2subject_path)
      transformed_atlas_path = fsl.findImage('atlas_transformed')

    file.delTemporary('T1_registered_histmatch.nii')

    if shared.parc_lut_file and shared.mrtrix_lut_file:
      run.command('labelconvert ' + transformed_atlas_path + ' ' + shared.parc_lut_file + ' ' + shared.mrtrix_lut_file + ' parc.mif')
    else: # Not all parcellations need to go through the labelconvert step; they may already be numbered incrementally from 1
      run.command('mrconvert ' + transformed_atlas_path + ' parc.mif -strides T1_registered.mif')
    file.delTemporary(transformed_atlas_path)

  file.delTemporary('T1_registered.nii')

  if shared.output_verbosity > 2:
    if shared.mrtrix_lut_file:
      label2colour_lut_option = ' -lut ' + shared.mrtrix_lut_file
    elif shared.parc_lut_file:
      label2colour_lut_option = ' -lut ' + shared.parc_lut_file
    else: # Use random colouring if no LUT available, but still generate the image
      label2colour_lut_option = ''
    run.command('label2colour parc.mif parcRGB.mif' + label2colour_lut_option)

  # If no parcellation is requested, it is still possible to generate a whole-brain tractogram by
  #   explicitly providing the -streamlines option
  num_streamlines = None
  if shared.streamlines:
    num_streamlines = shared.streamlines
  elif shared.parcellation != 'none':
    # If not manually specified, determine the appropriate number of streamlines based on the number of nodes in the parcellation:
    #   mean edge weight of 1,000 streamlines
    num_nodes = int(image.statistic('parc.mif', 'max'))
    num_streamlines = 500 * num_nodes * (num_nodes-1)
  if num_streamlines:

    # Step 14: Generate the tractogram
    app.console('Performing whole-brain fibre-tracking')
    tractogram_filepath = 'tractogram.tck'
    run.command('tckgen FOD_WM.mif ' + tractogram_filepath + ' -act 5TT.mif -backtrack -crop_at_gmwmi -maxlength 250 -power 0.33 ' \
                '-select ' + str(num_streamlines) + ' -seed_dynamic FOD_WM.mif')

    # Step 15: Use SIFT2 to determine streamline weights
    app.console('Running the SIFT2 algorithm to assign weights to individual streamlines')
    fd_scale_gm_option = ''
    if not multishell:
      fd_scale_gm_option = ' -fd_scale_gm'
    # If SIFT2 fails, reduce number of streamlines and try again
    while num_streamlines:
      run.command('tcksift2 ' + tractogram_filepath + ' FOD_WM.mif weights.csv -act 5TT.mif -out_mu mu.txt' + fd_scale_gm_option, False)
      if os.path.isfile('weights.csv'):
        break
      app.warn('SIFT2 failed, likely due to running out of RAM; reducing number of streamlines and trying again')
      num_streamlines = int(num_streamlines // 2)
      new_tractogram_filepath = 'tractogram_' + str(num_streamlines) + '.tck'
      run.command('tckedit ' + tractogram_filepath + ' ' + new_tractogram_filepath + ' -number ' + str(num_streamlines))
      file.delTemporary(tractogram_filepath)
      tractogram_filepath = new_tractogram_filepath
    if not num_streamlines:
      app.error('Unable to run SIFT2 algorithm for any number of streamlines')


    if shared.output_verbosity > 2:
      # Generate TDIs:
      # - A TDI at DWI native resolution, with SIFT mu scaling, and precise mapping
      #   (for comparison to WM ODF l=0 term, to verify that SIFT2 has worked correctly)
      app.console('Producing Track Density Images (TDIs)')
      with open('mu.txt', 'r') as f:
        mu = float(f.read())
      run.command('tckmap ' + tractogram_filepath + ' -tck_weights_in weights.csv -template FOD_WM.mif -precise - | ' \
                'mrcalc - ' + str(mu) + ' -mult tdi_native.mif')
      # - Conventional TDI at super-resolution (mostly just because we can)
      run.command('tckmap ' + tractogram_filepath + ' -tck_weights_in weights.csv -template vis.mif -vox ' + ','.join([str(value/3.0) for value in image.Header('vis.mif').spacing() ]) + ' -datatype uint16 tdi_highres.mif')


  if shared.parcellation != 'none':
    # Step 16: Generate the connectome
    #          Also get the mean length for each edge; this is the most likely alternative contrast to be useful
    app.console('Combining whole-brain tractogram with grey matter parcellation to produce the connectome')
    # TODO Does FreeSurfer have comparable gaps between parcels to the MNI volumetric parcellations?
    assignment_option = ' -assignment_radial_search 5' if shared.parcellation in [ 'yeo7mni', 'yeo17mni' ] else ''
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif connectome.csv -tck_weights_in weights.csv -out_assignments assignments.csv' + assignment_option)
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif meanlength.csv -tck_weights_in weights.csv -scale_length -stat_edge mean' + assignment_option)

    if shared.output_verbosity > 2:
      # Produce additional data that can be used for visualisation within mrview's connectome toolbar
      app.console('Generating geometric data for enhanced connectome visualisation')
      run.command('connectome2tck ' + tractogram_filepath + ' assignments.csv exemplars.tck -tck_weights_in weights.csv -exemplars parc.mif -files single')
      run.command('label2mesh parc.mif nodes.obj')
      run.command('meshfilter nodes.obj smooth nodes_smooth.obj')
      file.delTemporary('nodes.obj')



  # Prepare output path for writing
  app.console('Processing for subject \'' + label + '\' completed; writing results to output directory')
  if os.path.exists(output_dir):
    run.function(shutil.rmtree, output_dir)
  run.function(os.makedirs, output_dir)
  run.function(os.makedirs, os.path.join(output_dir, 'connectome'))
  run.function(os.makedirs, os.path.join(output_dir, 'dwi'))
  run.function(os.makedirs, os.path.join(output_dir, 'tractogram'))
  if shared.output_verbosity > 1:
    run.function(os.makedirs, os.path.join(output_dir, 'anat'))

  parc_string = '_parc-' + shared.parcellation

  # Generate a copy of the lookup table file:
  #   - Use the post-labelconvert file if it's used; otherwise, if the atlas
  #     itself comes with a lookup table that didn't require conversion, write that;
  #   - In the group directory rather than the subject directory;
  #   - If it doesn't already exist.
  lut_export_file = shared.mrtrix_lut_file if shared.mrtrix_lut_file else shared.parc_lut_file
  if lut_export_file:
    lut_export_path = os.path.join(output_prefix, parc_string[1:] + '_lookup' + os.path.splitext(lut_export_file)[1])
    try:
      shutil.copy(lut_export_file, lut_export_path)
    except OSError:
      pass

  # Copy / convert necessary files to output directory
  if shared.parcellation != 'none':
    run.function(shutil.copy, 'connectome.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_level-participant_connectome.csv'))
  if num_streamlines:
    run.function(shutil.copy, 'mu.txt', os.path.join(output_dir, 'tractogram', label + '_mu.txt'))
  run.command('mrconvert dwi_meanbzero.mif ' + os.path.join(output_dir, 'dwi', label + '_meanbzero.nii.gz') + ' -strides +1,+2,+3')
  run.command('mrconvert fa.mif ' + os.path.join(output_dir, 'dwi', label + '_model-tensor_fa.nii.gz') + ' -strides +1,+2,+3')
  run.function(shutil.copy, 'response_wm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-WM_response.txt'))
  with open(os.path.join(output_dir, 'dwi', label + '_bvalues.txt'), 'w') as f:
    f.write(' '.join([str(value) for value in bvalues]))
  if shared.output_verbosity > 1:
    if shared.parcellation != 'none':
      run.command('mrconvert parc.mif ' + os.path.join(output_dir, 'anat', label + parc_string + '_indices.nii.gz') + ' -strides +1,+2,+3')
      run.function(shutil.copy, 'meanlength.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_meanlength.csv'))
    run.command('mrconvert dwi.mif ' + os.path.join(output_dir, 'dwi', label + '_dwi.nii.gz') + \
                ' -export_grad_fsl ' + os.path.join(output_dir, 'dwi', label + '_dwi.bvec') + ' ' + os.path.join(output_dir, 'dwi', label + '_dwi.bval') + \
                ' -strides +1,+2,+3,+4')
    run.command('mrconvert dwi_mask.mif ' + os.path.join(output_dir, 'dwi', label + '_brainmask.nii.gz') + ' -datatype uint8 -strides +1,+2,+3')
    run.command('mrconvert T1_registered.mif ' + os.path.join(output_dir, 'anat', label + '_T1w.nii.gz') + ' -strides +1,+2,+3')
    run.command('mrconvert T1_mask_registered.mif ' + os.path.join(output_dir, 'anat', label + '_brainmask.nii.gz') + ' -datatype uint8 -strides +1,+2,+3')
    run.command('mrconvert 5TT.mif ' + os.path.join(output_dir, 'anat', label + '_5TT.nii.gz') + ' -strides +1,+2,+3,+4')
    run.command('mrconvert vis.mif ' + os.path.join(output_dir, 'anat', label + '_tissues3D.nii.gz') + ' -strides +1,+2,+3')
    run.command('mrconvert FOD_WM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-WM_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
    if multishell:
      run.function(shutil.copy, 'response_gm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-GM_response.txt'))
      run.function(shutil.copy, 'response_csf.txt', os.path.join(output_dir, 'dwi', label + '_tissue-CSF_response.txt'))
      run.command('mrconvert FOD_GM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-GM_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
      run.command('mrconvert FOD_CSF.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-CSF_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
      run.command('mrconvert tissues.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-all.nii.gz') + ' -strides +1,+2,+3,+4')
    if not shared.preprocessed:
      run.function(shutil.copytree, 'eddyqc', os.path.join(output_dir, 'dwi', 'eddyqc'))
  if shared.output_verbosity > 2:
    if num_streamlines:
      # Move rather than copying the tractogram just because of its size
      run.function(shutil.move, tractogram_filepath, os.path.join(output_dir, 'tractogram', label + '_tractogram.tck'))
      run.function(shutil.copy, 'weights.csv', os.path.join(output_dir, 'tractogram', label + '_weights.csv'))
      run.command('mrconvert tdi_native.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-native_tdi.nii.gz') + ' -strides +1,+2,+3')
      run.command('mrconvert tdi_highres.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-highres_tdi.nii.gz') + ' -strides +1,+2,+3')
    if shared.parcellation != 'none':
      run.function(shutil.copy, 'assignments.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_assignments.csv'))
      run.function(shutil.copy, 'exemplars.tck', os.path.join(output_dir, 'connectome', label + parc_string + '_exemplars.tck'))
      run.function(shutil.copy, 'nodes_smooth.obj', os.path.join(output_dir, 'anat', label + parc_string + '.obj'))
      run.command('mrconvert parcRGB.mif ' + os.path.join(output_dir, 'anat', label + parc_string +'_colour.nii.gz') + ' -strides +1,+2,+3')

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
      self.in_bzero = os.path.join(output_dir, label, 'dwi', label + '_meanbzero.nii.gz')
      self.in_fa = os.path.join(output_dir, label, 'dwi', label + '_model-tensor_fa.nii.gz')
      self.in_rf = os.path.join(output_dir, label, 'dwi', label + '_tissue-WM_response.txt')
      connectome_files = glob.glob(os.path.join(output_dir, label, 'connectome', label + '_parc-*_level-participant_connectome.csv'))
      if not connectome_files:
        app.error('No participant-level connectome file found for subject \'' + label + '\'')
      elif len(connectome_files) > 1:
        app.error('Connectomes from multiple parcellations detected for subject \'' + label + '\'; this is not yet supported')
      self.in_connectome = connectome_files[0]
      self.in_mu = os.path.join(output_dir, label, 'tractogram', label + '_mu.txt')

      for entry in vars(self).values():
        if not os.path.exists(entry):
          app.error('Unable to find critical subject data (expected location: ' + entry + ')')

      self.parcellation = re.findall('(?<=_parc-)[a-zA-Z0-9]*', os.path.basename(self.in_connectome))[0]

      # Permissible for this to not exist
      self.in_mask = os.path.join(output_dir, label, 'dwi', label + '_brainmask.nii.gz')

      with open(self.in_mu, 'r') as f:
        self.mu = float(f.read())

      self.RF = []
      with open(self.in_rf, 'r') as f:
        for line in f:
          self.RF.append([ float(v) for v in line.split() ])

      self.temp_mask = os.path.join('masks', label + '.nii.gz')
      self.link_fa = os.path.join('images', label + '.nii.gz')
      self.link_bzero = os.path.join('bzeros', label + '.nii.gz')
      self.temp_warp = os.path.join('warps', label + '.mif')
      self.temp_voxels = os.path.join('voxels', label + '.mif')
      self.median_bzero = 0.0
      self.dwiintensitynorm_factor = 1.0
      self.RF_multiplier = 1.0
      self.global_multiplier = 1.0
      self.temp_connectome = os.path.join('connectomes', label + '.csv')
      self.out_scale_intensity = os.path.join(output_dir, label, 'connectome', label + '_factor-intensity_multiplier.txt')
      self.out_scale_RF = os.path.join(output_dir, label, 'connectome', label + '_factor-response_multiplier.txt')
      self.out_connectome = os.path.join(output_dir, label, 'connectome', os.path.basename(self.in_connectome).replace('_level-participant', '_level-group'))

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
  #   Generate mask and FA image directories to be used in populate template generation.
  #   If output_verbosity >= 2 then a mask is already provided; if not, then one can be
  #     quickly calculated from the mean b=0 image, which must be provided
  progress = app.progressBar('Importing and preparing subject data', len(subjects))
  run.function(os.makedirs, 'bzeros')
  run.function(os.makedirs, 'images')
  run.function(os.makedirs, 'masks')
  for s in subjects:
    if os.path.exists(s.in_mask):
      run.function(os.symlink, s.in_mask, s.temp_mask)
    else:
      run.command('mrcalc ' + s.in_bzero + ' 0.0 -gt ' + s.temp_mask)
    run.function(os.symlink, s.in_fa, s.link_fa)
    run.function(os.symlink, s.in_bzero, s.link_bzero)
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
    run.command('mrtransform template.mif -warp_full ' + s.temp_warp + ' - -from 2 -template ' + s.link_bzero + ' | '
                'mrthreshold - ' + s.temp_voxels + ' -abs 0.4')
    s.median_bzero = float(image.statistic(s.link_bzero, 'median', '-mask ' + s.temp_voxels))
    file.delTemporary(s.link_bzero)
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
  # Can only do this if the parcellation is identical across subjects;
  #   this needs to be explicitly checked
  parcellation = subjects[0].parcellation
  consistent_parcellation = all(s.parcellation == parcellation for s in subjects)
  if consistent_parcellation:
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
  else:
    app.warn('Different parcellations across subjects; cannot calculate a group mean connectome')

  # Write results of interest back to the output directory;
  #   both per-subject and group information
  progress = app.progressBar('Writing results to output directory', len(subjects)+2)
  for s in subjects:
    run.function(shutil.copyfile, s.temp_connectome, s.out_connectome)
    with open(s.out_scale_intensity, 'w') as f:
      f.write(str(s.bzero_multiplier))
    with open(s.out_scale_RF, 'w') as f:
      f.write(str(s.RF_multiplier))
    progress.increment()

  with open(os.path.join(output_dir, 'tissue-WM_response.txt'), 'w') as f:
    for row in mean_RF:
      f.write(' '.join([str(v) for v in row]) + '\n')
  progress.increment()
  if consistent_parcellation:
    with open(os.path.join(output_dir, 'parc-' + parcellation + '_connectome.csv'), 'w') as f:
      for row in mean_connectome:
        f.write(' '.join([str(v) for v in row]) + '\n')
  progress.done()

# End of runGroup() function






# TODO Add Lausanne FreeSurfer atlas parcellations

analysis_choices = [ 'participant', 'group' ]
parcellation_choices = [ 'aal', 'aal2', 'brainnetome246fs', 'brainnetome246mni', 'craddock200', 'craddock400', 'desikan', 'destrieux', 'hcpmmp1', 'none', 'perry512', 'yeo7fs', 'yeo7mni', 'yeo17fs', 'yeo17mni' ]
registration_choices = [ 'ants', 'fsl' ]

app.init('Robert E. Smith (robert.smith@florey.edu.au)',
         'Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in MRtrix3')

# If running within a container, erase existing standard options, and fill with only desired options
if IS_CONTAINER:
  for option in reversed(app.cmdline._actions):
    app.cmdline._handle_conflict_resolve(None, [(option.option_strings[0],option)])
  # app.cmdline._action_groups[2] is "Standard options" that was created earlier
  app.cmdline._action_groups[2].add_argument('-d', '--debug', dest='debug', action='store_true', help='In the event of encountering an issue with the script, re-run with this flag set to provide more useful information to the developer')
  app.cmdline._action_groups[2].add_argument('-h', '--help', dest='help', action='store_true', help='Display help information for the script')
  app.cmdline._action_groups[2].add_argument('-n', '--n_cpus', type=int, metavar='number', dest='nthreads', help='Use this number of threads in MRtrix3 multi-threaded applications (0 disables multi-threading)')
  app.cmdline._action_groups[2].add_argument('-s', '--skip-bids-validator', dest='skipbidsvalidator', action='store_true', help='Skip BIDS validation')
  app.cmdline._action_groups[2].add_argument('-v', '--version', action='version', version=__version__)
else:
  app.cmdline._action_groups[2].add_argument(OPTION_PREFIX + 'skip-bids-validator', dest='skipbidsvalidator', action='store_true', help='Skip BIDS validation')
app.cmdline.add_argument('bids_dir', help='The directory with the input dataset formatted according to the BIDS standard.')
app.cmdline.add_argument('output_dir', help='The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.')
app.cmdline.add_argument('analysis_level', help='Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: ' + ', '.join(analysis_choices), choices=analysis_choices)
batch_options = app.cmdline.add_argument_group('Options specific to the batch processing of subject data')
batch_options.add_argument(OPTION_PREFIX + 'participant_label', nargs='+', help='The label(s) of the participant(s) that should be analyzed. The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a space-separated list.')
participant_options = app.cmdline.add_argument_group('Options that are relevant to participant-level analysis')
if not IS_CONTAINER:
  participant_options.add_argument(OPTION_PREFIX + 'atlas_path', metavar='path', help='The filesystem path in which to search for an atlas parcellation (may be necessary if not installed in the same location as they are placed within the BIDS App container)')
participant_options.add_argument(OPTION_PREFIX + 'output_verbosity', type=int, default=2, help='The verbosity of script output (number from 1 to 3); higher values result in more generated data being included in the output directory')
participant_options.add_argument(OPTION_PREFIX + 'parcellation', help='The choice of connectome parcellation scheme (compulsory for participant-level analysis). Options are: ' + ', '.join(parcellation_choices), choices=parcellation_choices)
participant_options.add_argument(OPTION_PREFIX + 'preprocessed', action='store_true', help='Indicate that the subject DWI data have been preprocessed, and hence initial image processing steps will be skipped')
participant_options.add_argument(OPTION_PREFIX + 'streamlines', type=int, default=0, help='The number of streamlines to generate for each subject (will be determined heuristically if not explicitly set)')
participant_options.add_argument(OPTION_PREFIX + 'template_reg', metavar='software', help='The choice of registration software for mapping subject to template space. Options are: ' + ', '.join(registration_choices), choices=registration_choices)

app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'template_reg fsl', 'Andersson, J. L. R.; Jenkinson, M. & Smith, S. Non-linear registration, aka spatial normalisation. FMRIB technical report, 2010, TR07JA2', True)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078', True)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Andersson, J. L. R.; Graham, M. S.; Zsoldos, E. & Sotiropoulos, S. N. Incorporating outlier detection and replacement into a non-parametric framework for movement and distortion correction of diffusion MR images. NeuroImage, 2016, 141, 556-572', True)
if not IS_CONTAINER:
  app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H.; Filippini, N. & Bastiani, M. Towards a comprehensive framework for movement and distortion correction of diffusion MR images: Within volume movement. NeuroImage, 2017, 152, 450-466', True)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Andersson, J. L. R.; Graham,, M. S.; Drobnjak, I.; Zhang, H. & Campbell, J. Susceptibility-induced distortion that varies due to motion: Correction in diffusion MR without acquiring additional data. NeuroImage, 2018, 171, 277-295', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ aal, aal2, brainnetome246mni, craddock200, craddock400, perry512, yeo7mni or yeo17mni ], and not using ' + OPTION_PREFIX + 'template_reg fsl', 'Avants, B. B.; Epstein, C. L.; Grossman, M.; Gee, J. C. Symmetric diffeomorphic image registration with cross-correlation: Evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 2008, 12, 26-41', True)
app.cmdline.addCitation('', 'Bhushan, C.; Haldar, J. P.; Choi, S.; Joshi, A. A.; Shattuck, D. W. & Leahy, R. M. Co-registration and distortion correction of diffusion and anatomical images based on inverse contrast normalization. NeuroImage, 2015, 115, 269-280', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ craddock200 or craddock400 ]', 'Craddock, R. C.; James, G. A.; Holtzheimer, P. E.; Hu, X. P.; Mayberg, H. S. A whole brain fMRI atlas generated via spatially constrained spectral clustering. Human Brain Mapping, 2012, 33(8), 1914-1928', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ brainnetome246fs, desikan, destrieux, hcpmmp1, yeo7fs or yeo17fs ]', 'Dale, A. M.; Fischl, B. & Sereno, M. I. Cortical Surface-Based Analysis: I. Segmentation and Surface Reconstruction. NeuroImage, 1999, 9, 179-194', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation desikan', 'Desikan, R. S.; Segonne, F.; Fischl, B.; Quinn, B. T.; Dickerson, B. C.; Blacker, D.; Buckner, R. L.; Dale, A. M.; Maguire, R. P.; Hyman, B. T.; Albert, M. S. & Killiany, R. J. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest NeuroImage, 2006, 31, 968-980', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation destrieux', 'Destrieux, C.; Fischl, B.; Dale, A. & Halgren, E. Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature NeuroImage, 2010, 53, 1-15', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ brainnetome246fs brainnetome246mni ]', 'Fan, L.; Li, H.; Zhuo, J.; Zhang, Y.; Wang, J.; Chen, L.; Yang, Z.; Chu, C.; Xie, S.; Laird, A.R.; Fox, P.T.; Eickhoff, S.B.; Yu, C.; Jiang, T. The Human Brainnetome Atlas: A New Brain Atlas Based on Connectional Architecture. Cerebral Cortex, 2016, 26 (8), 3508-3526', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation hcpmmp1', 'Glasser, M. F.; Coalson, T. S.; Robinson, E. C.; Hacker, C. D.; Harwell, J.; Yacoub, E.; Ugurbil, K.; Andersson, J.; Beckmann, C. F.; Jenkinson, M.; Smith, S. M. & Van Essen, D. C. A multi-modal parcellation of human cerebral cortex. Nature, 2016, 536, 171-178', True)
app.cmdline.addCitation('' if IS_CONTAINER else 'If using ROBEX for brain extraction', 'Iglesias, J. E.; Liu, C. Y.; Thompson, P. M. & Tu, Z. Robust Brain Extraction Across Datasets and Comparison With Publicly Available Methods. IEEE Transactions on Medical Imaging, 2011, 30, 1617-1634', True)
app.cmdline.addCitation('', 'Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data NeuroImage, 2014, 103, 411-426', False)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Kellner, E.; Dhital, B.; Kiselev, V. G.; Reisert, M. Gibbs-ringing artifact removal based on local subvoxel-shifts. Magnetic Resonance in Medicine, 2006, 76(5), 1574-1581', True)
app.cmdline.addCitation('If not using ' + OPTION_PREFIX + 'parcellation [ brainnetome246fs or brainnetome246mni ]', 'Patenaude, B.; Smith, S. M.; Kennedy, D. N. & Jenkinson, M. A Bayesian model of shape and appearance for subcortical brain segmentation. NeuroImage, 2011, 56, 907-922', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation perry512', 'Perry, A.; Wen, W.; Kochan, N. A.; Thalamuthu, A.; Sachdev, P. S.; Breakspear, M. The independent influences of age and education on functional brain networks and cognition in healthy older adults. Human Brain Mapping, 2017, 38(10), 5094-5114', True)
if not IS_CONTAINER:
  app.cmdline.addCitation('If not using ROBEX for brain extraction', 'Smith, S. M. Fast robust automated brain extraction. Human Brain Mapping, 2002, 17, 143-155', True)
app.cmdline.addCitation('', 'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. Anatomically-constrained tractography: Improved diffusion MRI streamlines tractography through effective use of anatomical information. NeuroImage, 2012, 62, 1924-1938', False)
app.cmdline.addCitation('', 'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015a, 104, 253-265', False)
app.cmdline.addCitation('', 'Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography. NeuroImage, 2015b, 119, 338-351', False)
app.cmdline.addCitation('', 'Tournier, J.-D.; Calamante, F., Gadian, D.G. & Connelly, A. Direct estimation of the fiber orientation density function from diffusion-weighted MRI data using spherical deconvolution. NeuroImage,     2004, 23, 1176-1185', False)
app.cmdline.addCitation('', 'Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670', False)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Tustison, N.; Avants, B.; Cook, P.; Zheng, Y.; Egan, A.; Yushkevich, P. & Gee, J. N4ITK: Improved N3 Bias Correction. IEEE Transactions on Medical Imaging, 2010, 29, 1310-1320', True)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ aal or aal2 ]', 'Tzourio-Mazoyer, N.; Landeau, B.; Papathanassiou, D.; Crivello, F.; Etard, O.; Delcroix, N.; Mazoyer, B. & Joliot, M. Automated Anatomical Labeling of activations in SPM using a Macroscopic Anatomical Parcellation of the MNI MRI single-subject brain. NeuroImage, 15(1), 273-289', True)
app.cmdline.addCitation('If ' + OPTION_PREFIX + 'preprocessed is not used', 'Veraart, J.; Novikov, D. S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406', False)
app.cmdline.addCitation('', 'Yeh, C.H.; Smith, R.E.; Liang, X.; Calamante, F.; Connelly, A. Correction for diffusion MRI fibre tracking biases: The consequences for structural connectomic metrics. Neuroimage, 2016, 142, 150-162', False)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation [ yeo7fs, yeo7mni, yeo17fs or yeo17mni ]', 'Yeo, B.T.; Krienen, F.M.; Sepulcre, J.; Sabuncu, M.R.; Lashkari, D.; Hollinshead, M.; Roffman, J.L.; Smoller, J.W.; Zollei, L.; Polimeni, J.R.; Fischl, B.; Liu, H. & Buckner, R.L. The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol, 2011, 106(3), 1125-1165', False)
app.cmdline.addCitation('If using ' + OPTION_PREFIX + 'parcellation perry512', 'Zalesky, A.; Fornito, A.; Harding, I. H.; Cocchi, L.; Yucel, M.; Pantelis, C. & Bullmore, E. T. Whole-brain anatomical networks: Does the choice of nodes matter? NeuroImage, 2010, 50, 970-983', True)
app.cmdline.addCitation('', 'Zhang, Y.; Brady, M. & Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Transactions on Medical Imaging, 2001, 20, 45-57', True)

app.parse()


# If running within a container, and the --debug option has been provided,
#   modify the interlly-stored MRtrix3 configuration contents, so that any
#   temporary directories will be constructed within the mounted output
#   directory, and therefore temporary directory contents will not be lost
#   upon container instance destruction if the script fails at any point.
if IS_CONTAINER and app.args.debug and 'ScriptTmpDir' not in app.config:
  app.config['ScriptTmpDir'] = os.path.abspath(app.args.output_dir)


if app.isWindows():
  app.error('Script cannot be run on Windows due to FSL dependency')

if app.args.skipbidsvalidator:
  app.console('Skipping BIDS validation based on user request')
elif find_executable('bids-validator'):
  run.command('bids-validator ' + app.args.bids_dir)
else:
  app.warn('BIDS validator script not installed; proceeding without validation of input data')

# Running participant level
if app.args.analysis_level == 'participant':

  if app.args.output_verbosity < 1 or app.args.output_verbosity > 3:
    app.error('Valid values for ' + OPTION_PREFIX + 'output_verbosity option are from 1 to 3')

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

  participant_shared = ParticipantShared(getattr(app.args, 'atlas_path', None),
                                         app.args.output_verbosity,
                                         app.args.parcellation,
                                         app.args.preprocessed,
                                         app.args.streamlines,
                                         app.args.template_reg)
  for subject_label in subjects_to_analyze:
    app.console('Commencing execution for subject ' + subject_label)
    runSubject(app.args.bids_dir, subject_label, participant_shared, os.path.abspath(app.args.output_dir))

# Running group level
elif app.args.analysis_level == 'group':

  if app.args.participant_label:
    app.error('Cannot use --participant_label option when performing group analysis')
  runGroup(os.path.abspath(app.args.output_dir))

app.complete()
