#!/usr/bin/env python

import glob
import json
import math
import os
import re
import shutil
from collections import namedtuple
from distutils.spawn import find_executable
import mrtrix3 #pylint: disable=import-error
from mrtrix3 import CONFIG, MRtrixError #pylint: disable=import-error
from mrtrix3 import app, fsl, image, matrix, path, run, utils #pylint: disable=import-error


IS_CONTAINER = os.path.exists('/version') \
               and os.path.exists('/mrtrix3_version')
__version__ = 'BIDS-App \'MRtrix3_connectome\' version {}' \
                  .format(open('/version').read()) \
              if IS_CONTAINER \
              else 'BIDS-App \'MRtrix3_connectome\' standalone'
OPTION_PREFIX = '--' if IS_CONTAINER else '-'

OUT_DWI_JSON_DATA = {'SkullStripped': False}
OUT_5TT_JSON_DATA = {'LabelMap': ['CGM', 'SGM', 'WM', 'CSF', 'Path']}




class T1wShared(object): #pylint: disable=useless-object-inheritance
    def __init__(self):
        try:
            self.fsl_anat_cmd = fsl.exe_name(find_executable('fsl_anat'))
        except MRtrixError:
            self.fsl_anat_cmd = None
        if find_executable('ROBEX'):
            self.robex_cmd = find_executable('runROBEX.sh')
        else:
            self.robex_cmd = None
        self.N4_cmd = find_executable('N4BiasFieldCorrection')

        if not self.fsl_anat_cmd and not self.robex_cmd:
            app.warn('No commands for T1 image processing found; '
                     'will be unable to proceed for any session without '
                     'pre-processed T1-weighted data')




class PreprocShared(object): #pylint: disable=useless-object-inheritance
    def __init__(self):
        fsl_path = os.environ.get('FSLDIR', '')
        if not fsl_path:
            raise MRtrixError(
                'Environment variable FSLDIR is not set; '
                'please run appropriate FSL configuration script')

        self.t1w_shared = T1wShared()

        self.dwibiascorrect_algo = 'ants'
        if not self.t1w_shared.N4_cmd:
            self.dwibiascorrect_algo = None
            app.warn('Could not find ANTs program "N4BiasFieldCorrection"; '
                     'will proceed without performing initial b=0 - based '
                     'DWI bias field correction')

        def get_eddy_help(binary_name):
            try:
                return run.command([binary_name, '--help'], show=False).stderr
            except run.MRtrixCmdError as eddy_except:
                return eddy_except.stderr

        self.eddy_binary = fsl.eddy_binary(True)
        if self.eddy_binary:
            self.eddy_cuda = True
            eddy_help = get_eddy_help(self.eddy_binary)
            if 'error while loading shared libraries' in eddy_help:
                app.warn('CUDA version of FSL "eddy" present on system, '
                         'but does not execute successfully; OpenMP version '
                         'will instead be used')
                self.eddy_binary = None
                self.eddy_cuda = False
                eddy_help = ''
        if not self.eddy_binary:
            self.eddy_binary = fsl.eddy_binary(False)
            if not self.eddy_binary:
                raise MRtrixError('Could not find FSL program "eddy"')
            self.eddy_cuda = False
            eddy_help = get_eddy_help(self.eddy_binary)

        app.debug('Eddy binary: ' + str(self.eddy_binary))
        app.debug('Eddy is CUDA version: ' + str(self.eddy_cuda))

        self.eddy_repol = False
        self.eddy_mporder = False
        self.eddy_mbs = False
        for line in eddy_help.splitlines():
            line = line.lstrip()
            if line.startswith('--repol'):
                self.eddy_repol = True
            elif line.startswith('--mporder') and self.eddy_cuda:
                self.eddy_mporder = True
            elif line.startswith('--estimate_move_by_susceptibility'):
                self.eddy_mbs = True
# End of PreprocShared() class







class ParticipantShared(object): #pylint: disable=useless-object-inheritance
    def __init__(self, atlas_path, parcellation,
                 streamlines, template_reg):

        if not parcellation:
            raise MRtrixError(
                'For participant-level analysis, '
                'desired parcellation must be provided using the '
                + OPTION_PREFIX + 'parcellation option')
        self.parcellation = parcellation

        self.streamlines = streamlines

        fsl_path = os.environ.get('FSLDIR', '')
        if not fsl_path:
            raise MRtrixError(
                'Environment variable FSLDIR is not set; '
                'please run appropriate FSL configuration script')

        self.t1w_shared = T1wShared()

        self.do_freesurfer = parcellation in ['brainnetome246fs',
                                              'desikan',
                                              'destrieux',
                                              'hcpmmp1',
                                              'yeo7fs',
                                              'yeo17fs']
        self.do_mni = parcellation in ['aal',
                                       'aal2',
                                       'brainnetome246mni',
                                       'craddock200',
                                       'craddock400',
                                       'perry512',
                                       'yeo7mni',
                                       'yeo17mni']
        if parcellation != 'none':
            assert self.do_freesurfer or self.do_mni

        if template_reg:
            if self.do_mni:
                self.template_registration_software = template_reg
            else:
                app.warn('Volumetric template registration '
                         'not being performed; '
                         + OPTION_PREFIX + 'template_reg option ignored')
                self.template_registration_software = ''
        else:
            self.template_registration_software = 'ants' if self.do_mni else ''
        if self.template_registration_software == 'ants':
            if not find_executable('ANTS') \
                    or not find_executable('WarpImageMultiTransform'):
                raise MRtrixError(
                    'Commands \'ANTS\' and \'WarpImageMultiTransform\' '
                    'must be present in PATH to use '
                    'ANTs software for template registration')
        elif self.template_registration_software == 'fsl':
            self.flirt_cmd = fsl.exe_name('flirt')
            self.fnirt_cmd = fsl.exe_name('fnirt')
            self.invwarp_cmd = fsl.exe_name('invwarp')
            self.applywarp_cmd = fsl.exe_name('applywarp')
            self.fnirt_config_basename = 'T1_2_MNI152_2mm.cnf'
            self.fnirt_config_path = os.path.join(fsl_path,
                                                  'etc',
                                                  'flirtsch',
                                                  self.fnirt_config_basename)
            if not os.path.isfile(self.fnirt_config_path):
                raise MRtrixError(
                    'Unable to find configuration file for FNI FNIRT '
                    + '(expected location: '
                    + self.fnirt_config_path + ')')

        self.template_image_path = ''
        self.template_mask_path = ''
        self.parc_image_path = ''
        self.parc_lut_file = ''
        self.mrtrix_lut_file = ''

        mrtrix_lut_dir = os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(app.__file__)),
                os.pardir,
                os.pardir,
                'share',
                'mrtrix3',
                'labelconvert'))

        if self.do_freesurfer:
            self.freesurfer_home = os.environ.get('FREESURFER_HOME', None)
            if not self.freesurfer_home:
                raise MRtrixError(
                    'Environment variable FREESURFER_HOME not set; '
                    'please verify FreeSurfer installation')
            if not find_executable('recon-all'):
                raise MRtrixError(
                    'Could not find FreeSurfer script "recon-all"; '
                    'please verify FreeSurfer installation')
            self.freesurfer_subjects_dir = os.environ['SUBJECTS_DIR'] \
                                           if 'SUBJECTS_DIR' in os.environ \
                                           else os.path.join(
                                               self.freesurfer_home,
                                               'subjects')
            if not os.path.isdir(self.freesurfer_subjects_dir):
                raise MRtrixError(
                    'Could not find FreeSurfer subjects directory '
                    '(expected location: '
                    + self.freesurfer_subjects_dir + ')')
            for subdir in ['fsaverage',
                           'fsaverage5',
                           'lh.EC_average',
                           'rh.EC_average']:
                if not os.path.isdir(os.path.join(self.freesurfer_subjects_dir,
                                                  subdir)):
                    raise MRtrixError(
                        'Could not find requisite FreeSurfer subject '
                        'directory \'' + subdir + '\' '
                        '(expected location: '
                        + os.path.join(self.freesurfer_subjects_dir,
                                       subdir) + ')')
            self.reconall_path = find_executable('recon-all')
            if not self.reconall_path:
                raise MRtrixError(
                    'Could not find FreeSurfer script "recon-all"; '
                    'please verify FreeSurfer installation')
            if parcellation in ['hcpmmp1', 'yeo7fs', 'yeo17fs']:
                if parcellation == 'hcpmmp1':

                    def hcpmmp_annot_path(hemi):
                        return os.path.join(self.freesurfer_subjects_dir,
                                            'fsaverage',
                                            'label',
                                            hemi + 'h.HCPMMP1.annot')

                    self.hcpmmp1_annot_paths = [hcpmmp_annot_path(hemi)
                                                for hemi in ['l', 'r']]
                    if not all([os.path.isfile(path) \
                                for path in self.hcpmmp1_annot_paths]):
                        raise MRtrixError(
                            'Could not find necessary annotation labels '
                            'for applying HCPMMP1 parcellation '
                            '(expected location: '
                            + hcpmmp_annot_path('?') + ')')
                else: # yeo7fs, yeo17fs

                    def yeo_annot_path(hemi):
                        return os.path.join(
                            self.freesurfer_subjects_dir,
                            'fsaverage5',
                            'label',
                            hemi + 'h.Yeo2011_'
                            + ('7' if parcellation == 'yeo7fs' else '17')
                            + 'Networks_N1000.split_components.annot')

                    self.yeo_annot_paths = [yeo_annot_path(hemi) \
                                            for hemi in ['l', 'r']]
                    if not all([os.path.isfile(path) \
                               for path in self.yeo_annot_paths]):
                        raise MRtrixError(
                            'Could not find necessary annotation labels '
                            'for applying Yeo2011 parcellation '
                            '(expected location: '
                            + yeo_annot_path('?') + ')')
                for cmd in ['mri_surf2surf', 'mri_aparc2aseg']:
                    if not find_executable(cmd):
                        raise MRtrixError(
                            'Could not find FreeSurfer command '
                            + cmd + ' '
                            '(necessary for applying '
                            'HCPMMP1 parcellation); '
                            'please verify FreeSurfer installation')
            elif parcellation == 'brainnetome246fs':

                def brainnetome_gcs_path(hemi):
                    return os.path.join(self.freesurfer_home,
                                        'average',
                                        hemi + 'h.BN_Atlas.gcs')

                self.brainnetome_cortex_gcs_paths = [
                    brainnetome_gcs_path(hemi)
                    for hemi in ['l', 'r']]
                if not all([os.path.isfile(path)
                            for path in self.brainnetome_cortex_gcs_paths]):
                    raise MRtrixError(
                        'Could not find necessary GCS files for '
                        'applying Brainnetome cortical parcellation via '
                        'FreeSurfer (expected location: '
                        + brainnetome_gcs_path('?') + ')')
                self.brainnetome_sgm_gca_path = \
                    os.path.join(self.freesurfer_home,
                                 'average',
                                 'BN_Atlas_subcortex.gca')
                if not os.path.isfile(self.brainnetome_sgm_gca_path):
                    raise MRtrixError(
                        'Could not find necessary GCA file for applying '
                        'Brainnetome sub-cortical parcellation '
                        'via FreeSurfer (expected location: '
                        + self.brainnetome_sgm_gca_path + ')')
                for cmd in ['mri_label2vol',
                            'mri_ca_label',
                            'mris_ca_label']:
                    if not find_executable(cmd):
                        raise MRtrixError(
                            'Could not find FreeSurfer command '
                            + cmd + ' '
                            '(necessary for applying '
                            'Brainnetome parcellation); '
                            'please verify FreeSurfer installation')

            # Query contents of recon-all script,
            #   looking for "-openmp" and "-parallel" occurences
            # Add options to end of recon-all -all call,
            #   based on which of these options are available
            #   as well as the value of app.numThreads
            # - In 5.3.0, just the -openmp option is available
            # - In 6.0.0, -openmp needs to be preceded by -parallel
            self.reconall_multithread_options = []
            if app.NUM_THREADS is None or app.NUM_THREADS > 1:
                with open(self.reconall_path, 'r') as f:
                    reconall_text = f.read().splitlines()
                for line in reconall_text:
                    line = line.strip()
                    if line == 'case "-parallel":':
                        self.reconall_multithread_options = \
                            ['-parallel'] + self.reconall_multithread_options
                    # If number of threads in this script is not being
                    #   explicitly controlled, allow recon-all to use
                    #   its own default number of threads
                    elif line == 'case "-openmp":' \
                            and app.NUM_THREADS is not None:
                        self.reconall_multithread_options.extend(
                            ['-openmp', str(app.NUM_THREADS)])
            if self.reconall_multithread_options:
                self.reconall_multithread_options = \
                    ' ' + ' '.join(self.reconall_multithread_options)
            else:
                self.reconall_multithread_options = ''
            app.debug(self.reconall_multithread_options)

            if parcellation == 'brainnetome246fs':
                self.parc_lut_file = os.path.join(self.freesurfer_home,
                                                  'BN_Atlas_246_LUT.txt')
                self.mrtrix_lut_file = ''
            elif parcellation == 'desikan':
                self.parc_lut_file = os.path.join(self.freesurfer_home,
                                                  'FreeSurferColorLUT.txt')
                self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir,
                                                    'fs_default.txt')
            elif parcellation == 'destrieux':
                self.parc_lut_file = os.path.join(self.freesurfer_home,
                                                  'FreeSurferColorLUT.txt')
                self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir,
                                                    'fs_a2009s.txt')
            elif parcellation == 'hcpmmp1':
                self.parc_lut_file = os.path.join(mrtrix_lut_dir,
                                                  'hcpmmp1_original.txt')
                self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir,
                                                    'hcpmmp1_ordered.txt')
            elif parcellation in ['yeo7fs', 'yeo17fs']:
                self.parc_lut_file = \
                    os.path.join(self.freesurfer_home,
                                 'Yeo2011_'
                                 + ('7' if parcellation == 'yeo7fs' else '17')
                                 + 'networks_Split_Components_LUT.txt')
                self.mrtrix_lut_file = \
                    os.path.join(mrtrix_lut_dir,
                                 'Yeo2011_'
                                 + ('7' if parcellation == 'yeo7fs' else '17')
                                 + 'N_split.txt')
            else:
                assert False

            # If running in a container environment, and --debug is used
            #   (resulting in the scratch directory being a mounted drive),
            #   it's possible that attempting to construct a softlink may
            #   lead to an OSError
            # As such, run a test to determine whether or not it is
            #   possible to construct a softlink within the scratch
            #   directory; if it is not possible, revert to performing
            #   deep copies of the relevant FreeSurfer template directories
            self.freesurfer_template_link_function = os.symlink
            try:
                self.freesurfer_template_link_function(
                    self.freesurfer_subjects_dir,
                    'test_softlink')
                app.cleanup('test_softlink')
                app.debug('Using softlinks to FreeSurfer template directories')
            except OSError:
                app.debug('Unable to create softlinks; '
                          'will perform deep copies of FreeSurfer '
                          'template directories')
                self.freesurfer_template_link_function = shutil.copytree

        elif self.do_mni:
            self.template_image_path = \
                os.path.join(fsl_path,
                             'data',
                             'standard',
                             'MNI152_T1_2mm.nii.gz')
            self.template_mask_path = \
                os.path.join(fsl_path,
                             'data',
                             'standard',
                             'MNI152_T1_2mm_brain_mask.nii.gz')
            if parcellation == 'aal':
                self.parc_image_path = \
                    os.path.abspath(os.path.join(os.sep,
                                                 'opt',
                                                 'aal',
                                                 'ROI_MNI_V4.nii'))
                self.parc_lut_file = \
                    os.path.abspath(os.path.join(os.sep,
                                                 'opt',
                                                 'aal',
                                                 'ROI_MNI_V4.txt'))
                self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir,
                                                    'aal.txt')
            elif parcellation == 'aal2':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'aal',
                                     'ROI_MNI_V5.nii'))
                self.parc_lut_file = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'aal',
                                     'ROI_MNI_V5.txt'))
                self.mrtrix_lut_file = os.path.join(mrtrix_lut_dir,
                                                    'aal2.txt')
            elif parcellation == 'brainnetome246mni':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'brainnetome',
                                     'BN_Atlas_246_1mm.nii.gz'))
                self.parc_lut_file = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'brainnetome',
                                     'BN_Atlas_246_LUT.txt'))
                self.mrtrix_lut_file = ''
            elif parcellation == 'craddock200':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'ADHD200_parcellate_200.nii.gz'))
                self.parc_lut_file = ''
                self.mrtrix_lut_file = ''
            elif parcellation == 'craddock400':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'ADHD200_parcellate_400.nii.gz'))
                self.parc_lut_file = ''
                self.mrtrix_lut_file = ''
            elif parcellation == 'perry512':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     '512inMNI.nii'))
                self.parc_lut_file = ''
                self.mrtrix_lut_file = ''
            elif parcellation == 'yeo7mni':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'Yeo2011',
                                     'Yeo2011_7Networks_N1000.split_components'
                                     + '.FSL_MNI152_1mm.nii.gz'))
                self.parc_lut_file = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'Yeo2011',
                                     '7Networks_ColorLUT_freeview.txt'))
                self.mrtrix_lut_file = ''
            elif parcellation == 'yeo17mni':
                self.parc_image_path = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'Yeo2011',
                                     'Yeo2011_17Networks_N1000'
                                     + '.split_components'
                                     + '.FSL_MNI152_1mm.nii.gz'))
                self.parc_lut_file = \
                    os.path.abspath(
                        os.path.join(os.sep,
                                     'opt',
                                     'Yeo2011',
                                     '17Networks_ColorLUT_freeview.txt'))
                self.mrtrix_lut_file = ''
            else:
                assert False

        def find_atlas_file(filepath, description):
            if not filepath:
                return ''
            if os.path.isfile(filepath):
                return filepath
            if not atlas_path:
                raise MRtrixError('Could not find ' + description + ' '
                                  '(expected location: ' + filepath + ')')
            newpath = os.path.join(os.path.dirname(atlas_path),
                                   os.path.basename(filepath))
            if os.path.isfile(newpath):
                return newpath
            raise MRtrixError('Could not find ' + description + ' '
                              '(tested locations: '
                              '\'' + filepath + '\', '
                              '\'' + newpath + '\')')

        self.template_image_path = \
            find_atlas_file(self.template_image_path,
                            'template image')
        self.template_mask_path = \
            find_atlas_file(self.template_mask_path,
                            'template brain mask image')
        self.parc_image_path = \
            find_atlas_file(self.parc_image_path,
                            'parcellation image')
        self.parc_lut_file = \
            find_atlas_file(self.parc_lut_file,
                            'parcellation lookup table file')

        if self.mrtrix_lut_file and not os.path.exists(self.mrtrix_lut_file):
            raise MRtrixError(
                'Could not find MRtrix3 connectome lookup table file '
                '(expected location: ' + self.mrtrix_lut_file + ')')
# End of ParticipantShared() class










# Regardless of the source of T1-weighted image information,
#   scratch directory will contain at completion of this function:
# - Either:
#   - T1.mif
#     or
#   - T1_premasked.mif
#     , depending on software used (full unmasked T1-weighted image data
#     may not be available)
# - T1_mask.mif
#
# TODO Think after all that this does need to export a greater
#   amount of information with respect to what was done & how it was derived
def get_t1w_preproc_images(import_path,
                           session,
                           t1w_shared,
                           t1w_preproc):

    session_label = '_'.join(session)
    preproc_image_path = None
    preproc_image_is_masked = None
    preproc_mask_path = None
    raw_image_path = None

    if t1w_preproc:

        # Multiple possibilities for how such data may have been provided:
        # - Raw path to the image itself
        # - Path to anat/ directory
        # - Path to subject directory within BIDS Derivatives dataset
        # - Path to BIDS Derivatives dataset
        if os.path.isfile(t1w_preproc):
            preproc_image_path = t1w_preproc
        else:
            expected_image_basename = session_label + '*_T1w.nii*'
            for candidate in [
                    os.path.join(t1w_preproc,
                                 expected_image_basename),
                    os.path.join(t1w_preproc,
                                 'anat',
                                 expected_image_basename),
                    os.path.join(os.path.join(t1w_preproc, *session),
                                 'anat',
                                 expected_image_basename)]:

                glob_result = glob.glob(candidate)
                if glob_result:
                    if len(glob_result) == 1:
                        preproc_image_path = glob_result[0]
                        break
                    glob_refined_result = \
                        [item for item in glob_result \
                            if not '_space-' in item]
                    if len(glob_refined_result) == 1:
                        preproc_image_path = glob_refined_result[0]
                        break
                    raise MRtrixError('Unable to unambiguously select pre-'
                                      'processed T1-weighted image due to '
                                      'multiple candidates in location "'
                                      + candidate
                                      + '": '
                                      + ';'.join(glob_result))

            if preproc_image_path is None:
                raise MRtrixError('No pre-processed T1w image found from '
                                  'specified path "' + t1w_preproc + '"')

    else:

        # Look inside of import_path to see if there is a pre-processed
        #   T1w image already there
        glob_result = glob.glob(os.path.join(os.path.join(import_path,
                                                          *session),
                                             'anat',
                                             session_label
                                             + '*_desc-preproc*_T1w.nii*'))
        if glob_result:
            if len(glob_result) == 1:
                preproc_image_path = glob_result[0]
            else:
                raise MRtrixError('Multiple pre-processed T1w images found in '
                                  + 'import directory "'
                                  + import_path
                                  + '": '
                                  + ';'.join(glob_result))

    # Same checks regardless of whether the existing pre-processed image
    #   comes from the output directory or a user-specified location
    if preproc_image_path:

        if '_desc-preproc' not in preproc_image_path:
            raise MRtrixError('Selected T1-weighted image "'
                              + preproc_image_path
                              + '" not flagged as pre-processed')

        # Check to see if there's a JSON file along with the T1-weighted
        #   image; if they is, parse it to find out whether or not the
        #   pre-processed image has been brain-extracted
        expected_json_path = preproc_image_path \
                             .rstrip('.gz') \
                             .rstrip('.nii') \
                             + '.json'
        try:
            with open(expected_json_path, 'r') as t1_json_file:
                t1_json_data = json.load(t1_json_file)
            preproc_image_is_masked = t1_json_data.get('SkullStripped', None)
        except IOError:
            pass
        if preproc_image_is_masked is None:
            # Try to assess whether or not skull-stripping has occurred
            #   based on the prevalence of NaNs or zero values
            # - Obtain mask that contains:
            #   - All voxels with non-finite value
            #   and:
            #   - All voxels with a value of zero
            # - Feed to mrstats, extracting the mean
            # - If this is > 25% of the image, it's skull-stripped
            frac_voxels_outside_mask = \
                float(run.command('mrcalc '
                                  + preproc_image_path
                                  + ' 0 -eq 1 '
                                  + preproc_image_path
                                  + ' -finite -sub -add - '
                                  + '| '
                                  + ' mrstats - -output mean').stdout)
            preproc_image_is_masked = \
                frac_voxels_outside_mask > 0.25
            app.warn('No sidecar information for pre-processed '
                     + 'T1-weighted image "'
                     + preproc_image_path
                     + '" regarding skull-stripping; '
                     + 'image has been inferred to '
                     + ('be' if preproc_image_is_masked else 'not be')
                     + ' pre-masked based on image data ('
                     + str(int(round(100.0 * frac_voxels_outside_mask)))
                     + '% of voxels contain no data)')

        # Copy pre-procesed T1-weighted image into scratch directory
        run.command('mrconvert '
                    + preproc_image_path
                    + ' '
                    + path.to_scratch('T1_masked.mif' \
                                      if preproc_image_is_masked \
                                      else 'T1.mif'))

        # If we have been provided with a pre-processed T1-weighted image
        #   (regardless of where it has come from), check to see if there
        #   is a corresponding mask image
        preproc_mask_path = preproc_image_path \
                            .replace('_desc-preproc', '_desc-brain') \
                            .replace('_T1w.nii', '_mask.nii')
        if os.path.isfile(preproc_mask_path):
            run.command('mrconvert '
                        + preproc_mask_path
                        + ' '
                        + path.to_scratch('T1_mask.mif')
                        + ' -datatype bit')
        elif preproc_image_is_masked:
            run.command('mrcalc '
                        + preproc_image_path
                        + ' 0 -gt '
                        + path.to_scratch('T1_mask.mif')
                        + ' -datatype bit')
            # No pre-existing mask image, but we also don't want to
            #   run our own brain extraction
            preproc_mask_path = ''
        else:
            app.console('No brain mask image found alongside '
                        'pre-processed T1-weighted image "'
                        + preproc_image_path
                        + '"; will generate one manually')
            preproc_mask_path = None

    else:

        # Check input path for raw un-processed T1w image
        glob_result = glob.glob(os.path.join(os.path.join(import_path,
                                                          *session),
                                             'anat',
                                             session_label + '*_T1w.nii*'))
        if not glob_result:
            raise MRtrixError('No raw or pre-processed T1-weighted images '
                              + 'could be found in input directory "'
                              + import_path
                              + '" for session '
                              + session_label)
        if len(glob_result) > 1:
            raise MRtrixError('Multiple raw T1w images found in '
                              + 'input directory "'
                              + import_path
                              + '" for session '
                              + session_label
                              + ': '
                              + ';'.join(glob_result))
        raw_image_path = glob_result[0]

    # Do we need to do any pre-processing of our own at all?
    if preproc_mask_path is None:

        app.console('Performing requisite processing of '
                    'T1-weighted data')
        cwd = os.getcwd()
        run.function(os.makedirs,
                     path.to_scratch('t1w_preproc'))
        run.function(os.chdir, path.to_scratch('t1w_preproc'))

        if preproc_image_path:

            if t1w_shared.robex_cmd:
                app.console('Using ROBEX for brain extraction for session '
                            + session_label
                            + ', operating on existing pre-processed '
                            + 'T1-weighted image')
            elif t1w_shared.fsl_anat_path:
                app.console('Using fsl_anat for brain extraction for session '
                            + session_label
                            + ' (due to ROBEX not being installed)'
                            + ', operating on existing pre-processed '
                            + 'T1-weighted image')
            else:
                raise MRtrixError('Unable to continue processing for session '
                                  + session_label
                                  + ': no pre-processed T1-weighted image mask '
                                  + 'available / provided, and no appropriate '
                                  + 'brain masking software installed')

            run.command('mrconvert '
                        + preproc_image_path
                        + ' T1.nii -strides '
                        + ('+1,+2,+3' if t1w_shared.robex_cmd else '-1,+2,+3'))

            if t1w_shared.robex_cmd:
                run.command(t1w_shared.robex_cmd
                            + ' T1.nii T1_brain.nii T1_mask.nii')
                run.command('mrconvert T1_mask.nii '
                            + path.to_scratch('T1_mask.mif')
                            + ' -datatype bit')
            elif t1w_shared.fsl_anat_cmd:
                run.command(t1w_shared.fsl_anat_cmd
                            + ' -i T1.nii --noseg --nosubcortseg --nobias')
                run.command('mrconvert '
                            + fsl.find_image(
                                os.path.join('T1.anat',
                                             'T1_brain_mask'))
                            + ' '
                            + path.to_scratch('T1_mask.mif')
                            + ' -datatype bit')
            else:
                assert False

        else:

            # No pre-processed T1-weighted image available:
            #   do everything based on the raw T1-weighted image
            if t1w_shared.robex_cmd and t1w_shared.N4_cmd:
                app.console('No pre-processed T1-weighted image '
                            + 'found for session '
                            + session_label
                            + '; will use ROBEX and N4 for '
                            + 'iterative brain extraction and bias field '
                            + 'correction from raw T1-weighted image input')
            elif t1w_shared.fsl_anat_cmd:
                app.console('No pre-processed T1-weighted image '
                            + 'found for session '
                            + session_label
                            + '; will use fsl_anat for brain extraction and '
                            + 'bias field correction from raw T1-weighted '
                            + 'image input'
                            + (''
                               if t1w_shared.robex_cmd
                               else ' (ROBEX not installed)')
                            + (''
                               if t1w_shared.N4_cmd
                               else ' (N4 not installed)'))
            else:
                raise MRtrixError('Cannot complete processing for session '
                                  + session_label
                                  + ': no pre-processed T1-weighted image '
                                  + 'available, and software tools for '
                                  + 'processing raw T1-weighted image '
                                  + 'not installed')

            run.command('mrconvert '
                        + raw_image_path
                        + ' T1.nii -strides '
                        + ('+1,+2,+3'
                           if t1w_shared.robex_cmd and t1w_shared.N4_cmd
                           else '-1,+2,+3'))

            if t1w_shared.robex_cmd and t1w_shared.N4_cmd:

                # Do a semi-iterative approach here:
                #   Get an initial brain mask, use that mask to estimate a
                #   bias field, then re-compute the brain mask
                # TODO Consider making this fully iterative, just like the
                #   approach in preproc with dwi2mask and mtnormalise
                run.command(t1w_shared.robex_cmd
                            + ' T1.nii T1_initial_brain.nii'
                            + ' T1_initial_mask.nii')
                app.cleanup('T1_initial_brain.nii')
                run.command(t1w_shared.N4_cmd
                            + ' -i T1.nii'
                            + ' -w T1_initial_mask.nii'
                            + ' -o T1_biascorr.nii')
                app.cleanup('T1.nii')
                app.cleanup('T1_initial_mask.nii')
                run.command(t1w_shared.robex_cmd
                            + ' T1_biascorr.nii T1_biascorr_brain.nii'
                            + ' T1_biascorr_brain_mask.nii')
                app.cleanup('T1_biascorr_brain.nii')
                run.command('mrconvert T1_biascorr.nii '
                            + path.to_scratch('T1.mif'))
                app.cleanup('T1_biascorr.nii')
                run.command('mrconvert T1_biascorr_brain_mask.nii '
                            + path.to_scratch('T1_mask.mif')
                            + ' -datatype bit')
                app.cleanup('T1_biascorr_brain_mask.nii')

            elif t1w_shared.fsl_anat_cmd:

                run.command(t1w_shared.fsl_anat_cmd
                            + ' -i T1.nii --noseg --nosubcortseg')
                app.cleanup('T1.nii')
                run.command('mrconvert '
                            + fsl.find_image(
                                os.path.join('T1.anat',
                                             'T1_biascorr'))
                            + ' '
                            + path.to_scratch('T1_premasked.mif'))
                run.command('mrconvert '
                            + fsl.find_image(
                                os.path.join('T1.anat',
                                             'T1_biascorr_brain_mask'))
                            + ' '
                            + path.to_scratch('T1_mask.mif')
                            + ' -datatype bit')
                app.cleanup('T1.anat')

            else:
                assert False

        run.function(os.chdir, cwd)
        app.cleanup(path.to_scratch('t1w_preproc'))

# Completed function get_t1w_preproc()















def run_preproc(bids_dir, session, shared,
                t1w_preproc_path, output_verbosity, output_app_dir):

    session_label = '_'.join(session)
    output_subdir = os.path.join(output_app_dir, 'preproc', *session)
    if os.path.exists(output_subdir):
        app.warn('Output directory for session "' + session_label + '" '
                 'already exists; all contents will be erased when this '
                 'execution completes')

    app.make_scratch_dir()

    # Need to perform an initial import of JSON data using mrconvert;
    #   so let's grab the diffusion gradient table as well
    # If no bvec/bval present, need to go down the directory listing
    # Only try to import JSON file if it's actually present
    #   direction in the acquisition they'll need to be split
    #   across multiple files
    # May need to concatenate more than one input DWI, since if there's
    #   more than one phase-encode direction in the acquired DWIs
    #   (i.e. not just those used for estimating the inhomogeneity field),
    #   they will need to be stored as separate NIfTI files in the
    #   'dwi/' directory.
    app.console('Importing DWI data into scratch directory')
    in_dwi_path = os.path.join(os.path.join(bids_dir, *session),
                               'dwi',
                               '*_dwi.nii*')
    in_dwi_image_list = sorted(glob.glob(in_dwi_path))
    if not in_dwi_image_list:
        raise MRtrixError('No DWI data found for session \''
                          + session_label
                          + '\' (search location: ' + in_dwi_path)
    dwi_index = 0

    re_is_complex = re.compile(r'_part-(mag|phase)_')

    for entry in in_dwi_image_list:
        # Is this one image in a magnitude-phase pair?
        is_complex = re_is_complex.search(os.path.basename(entry))
        if is_complex:
            matching_text = is_complex.group(1)
            complex_part = matching_text.strip('_').split('-')[-1]
            assert complex_part in ['mag', 'phase']
            if complex_part == 'mag':
                # Find corresponding phase image
                in_phase_image = entry.replace('_part-mag_', '_part-phase_')
                if in_phase_image not in in_dwi_image_list:
                    raise MRtrixError(
                        'Image '
                        + entry
                        + ' does not have corresponding phase image')
                # Check if phase image is stored in radians
                phase_stats = image.statistics(in_phase_image, allvolumes=True)
                if abs(2.0*math.pi - (phase_stats.max - phase_stats.min)) \
                    > 0.01:
                    app.warn('Phase image '
                             + in_phase_image
                             + ' is not stored in radian units '
                             + '(values from '
                             + str(phase_stats.min)
                             + ' to '
                             + str(phase_stats.max)
                             + '); data will be rescaled automatically')
                    # Are the values stored as integers? If so, assume that
                    #   _one greater_ than the max phase corresponds to 2pi
                    add_to_max_phase = 1 \
                        if phase_stats.min.is_integer() \
                        and phase_stats.max.is_integer() \
                        else 0
                    phase_rescale_factor = 2.0 * math.pi / \
                                           (phase_stats.max
                                            + add_to_max_phase
                                            - phase_stats.min)
                else:
                    phase_rescale_factor = None

            else:
                # Make sure we also have the corresponding magnitude image
                if entry.replace('_part-phase_',
                                 '_part-mag_') not in in_dwi_image_list:
                    raise MRtrixError(
                        'Image '
                        + entry
                        + ' does not have corresponding mag image')
                # Do nothing for the second image in the pair
                continue
        else:

            in_phase_image = None

        # Find sidecar files
        # dcm2bids will have separate bvecs / bvals / json for the
        #   magnitude and phase images;
        # for proper BIDS compliance, these sidecar files will be stored
        #   without the "_part-[mag|phase]" part

        # os.path.split() falls over with .nii.gz extensions;
        #   only removes the .gz
        dwi_prefix = entry.split(os.extsep)[0] + '.'

        sidecar_prefixes = [dwi_prefix,
                            dwi_prefix.replace('_part-mag_', '_'),
                            os.path.join(bids_dir, dwi_prefix),
                            os.path.join(
                                bids_dir,
                                dwi_prefix.replace('_part-mag_', '_'))]
        entry_bval = None
        entry_bvec = None
        entry_json = None
        for prefix in sidecar_prefixes:
            if not entry_bvec and \
                os.path.isfile(prefix + 'bval') and \
                os.path.isfile(prefix + 'bvec'):
                entry_bval = prefix + 'bval'
                entry_bvec = prefix + 'bvec'
            if not entry_json and \
                os.path.isfile(prefix + 'json'):
                entry_json = prefix + 'json'
        if not entry_bvec:
            raise MRtrixError(
                'Unable to locate valid diffusion gradient table '
                'for image \'' + entry + '\'')
        if not entry_json:
            raise MRtrixError(
                'Unable to locate valid JSON sidecar file '
                'for image \'' + entry + '\'')

        grad_import_option = ' -fslgrad ' + entry_bvec + ' ' + entry_bval
        json_import_option = ' -json_import ' + entry_json

        # Import the data
        dwi_index += 1
        if in_phase_image:
            run.command('mrconvert '
                        + entry
                        + grad_import_option
                        + json_import_option
                        + ' - '
                        + '| '
                        + 'mrcalc '
                        + '- '
                        + in_phase_image + ' '
                        + (str(phase_stats.min) + ' -sub '
                           + str(phase_rescale_factor) + ' -mult ' \
                           if phase_rescale_factor else '')
                        + '-polar '
                        + path.to_scratch('dwi' + str(dwi_index) + '.mif'))
        else:
            run.command('mrconvert '
                        + entry + ' '
                        + path.to_scratch('dwi' + str(dwi_index) + '.mif')
                        + grad_import_option
                        + json_import_option)

    dwi_image_list = ['dwi' + str(index) + '.mif'
                      for index in range(1, dwi_index+1)]

    if len(dwi_image_list) > 1:
        dwi_first_header = image.Header(
            path.to_scratch(dwi_image_list[0], False))
        for i in dwi_image_list[1:]:
            if not image.match(dwi_first_header,
                               path.to_scratch(i, False),
                               up_to_dim=3):
                raise MRtrixError(
                    'DWI series not defined on same image grid; '
                    'script not yet capable of handling such data')

    # Go hunting for reversed phase-encode data
    #   dedicated to field map estimation
    in_fmap_image_list = []
    fmap_dir = os.path.join(os.path.join(bids_dir, *session), 'fmap')
    fmap_index = 0
    fmap_image_list = []
    if os.path.isdir(fmap_dir):
        app.console('Importing fmap data into scratch directory')
        in_fmap_image_list = sorted(
            glob.glob(os.path.join(fmap_dir, '*_dir-*_epi.nii*')))
        for entry in in_fmap_image_list:
            prefix = os.path.splitext(entry.rstrip('.gz'))[0]
            json_path = prefix + '.json'
            try:
                with open(json_path, 'r') as f:
                    json_elements = json.load(f)
            except OSError:
                app.warn('No JSON file found for image "'
                         + entry
                         + '"; not importing')
                continue
            if 'IntendedFor' in json_elements:
                if isinstance(json_elements['IntendedFor'], list) and \
                    not any(any(i.endswith(target) for i in in_dwi_image_list)
                            for target in json_elements['IntendedFor']):
                    app.console('Image \'' + entry + '\' is not intended '
                                'for use with DWIs; skipping')
                    continue
                if not any(i.endswith(json_elements['IntendedFor'])
                           for i in in_dwi_image_list):
                    app.console('Image \'' + entry + '\' is not intended '
                                'for use with DWIs; skipping')
                    continue
            if not os.path.isfile(json_path):
                raise MRtrixError('No sidecar JSON file found '
                                  'for image \'' + entry + '\'')
            # fmap files will not come with any gradient encoding in the JSON;
            #   therefore we need to add it manually ourselves so that
            #   mrcat / mrconvert can appropriately handle the table once
            #   these images are concatenated with the DWIs
            fmap_index += 1
            fmap_image_size = image.Header(entry).size()
            fmap_image_num_volumes = \
                1 if len(fmap_image_size) == 3 else fmap_image_size[3]
            fmap_dwscheme_file = 'fmap' + str(fmap_index) + '.b'
            with open(path.to_scratch(fmap_dwscheme_file, False), 'w') as f:
                for _ in range(0, fmap_image_num_volumes):
                    f.write('0,0,1,0\n')
            run.command('mrconvert '
                        + entry + ' '
                        + path.to_scratch('fmap' + str(fmap_index) + '.mif',
                                          True)
                        + ' -json_import ' + json_path
                        + ' -grad ' + path.to_scratch(fmap_dwscheme_file,
                                                      True))
            app.cleanup(fmap_dwscheme_file)

        fmap_image_list = ['fmap' + str(index) + '.mif'
                           for index in range(1, fmap_index+1)]

    # No need to explicitly check whether fmap/ images are defined
    #   on a common image grid; these we are happy to resample

    # If there's no usable data in fmap/ directory,
    #   need to check to see if there's any phase-encoding
    #   contrast within the input DWI(s)
    if not fmap_image_list and len(dwi_image_list) < 2:
        raise MRtrixError('Inadequate data for pre-processing of session '
                          '\"' + session_label + '": '
                          'No phase-encoding contrast in input DWIs, '
                          'and no fmap/ directory, '
                          'so EPI distortion correction cannot be performed')

    # Get T1-weighted image data
    #   (could be generated from raw data, or grabbed from a
    #   user-specified path source;
    #   don't look in output directory for preproc)
    get_t1w_preproc_images(bids_dir,
                           session,
                           shared.t1w_shared,
                           t1w_preproc_path)
    T1_is_premasked = os.path.isfile(path.to_scratch('T1_premasked.mif',
                                                     False))
    T1_image = 'T1_premasked.mif' if T1_is_premasked else 'T1.mif'

    cwd = os.getcwd()
    app.goto_scratch_dir()

    dwifslpreproc_se_epi = ''
    dwifslpreproc_se_epi_option = ''

    # For DWI data, denoise each individually, then concatenate;
    #  for fmap/ data, just concatenate and run mrdegibbs

    if len(dwi_image_list) == 1:
        run.function(os.rename, dwi_image_list[0], 'dwi.mif', show=False)
        dwi_image_list[0] = 'dwi.mif'

    # Step 1: Denoise
    app.console('Denoising DWI data')
    for entry in dwi_image_list:
        run.command('dwidenoise ' + entry + ' '
                    + os.path.splitext(entry)[0] + '_denoised.mif')
        app.cleanup(entry)
    dwi_image_list = [os.path.splitext(entry)[0] + '_denoised.mif'
                      for entry in dwi_image_list]

    # If data are complex, take the magnitude
    new_dwi_image_list = []
    for entry in dwi_image_list:
        if image.Header(entry).datatype().startswith('CFloat'):
            mag_entry = os.path.splitext(entry)[0] + '_mag.mif'
            run.command('mrcalc ' + entry + ' -abs ' + mag_entry)
            app.cleanup(entry)
            new_dwi_image_list.append(mag_entry)
        else:
            new_dwi_image_list.append(entry)
    dwi_image_list = new_dwi_image_list

    # Step 2: Gibbs ringing removal
    app.console('Performing Gibbs ringing removal for DWI'
                + (' and fmap ' if fmap_image_list else ' ')
                + 'data')
    for i in dwi_image_list:
        run.command('mrdegibbs ' + i + ' '
                    + os.path.splitext(i)[0] + '_degibbs.mif'
                    + ' -nshifts 50')
        app.cleanup(i)
    dwi_image_list = [os.path.splitext(i)[0] + '_degibbs.mif'
                      for i in dwi_image_list]
    for i in fmap_image_list:
        run.command('mrdegibbs ' + i + ' '
                    + os.path.splitext(i)[0] + '_degibbs.mif'
                    + ' -nshifts 50')
        app.cleanup(i)
    fmap_image_list = [os.path.splitext(i)[0] + '_degibbs.mif'
                       for i in fmap_image_list]

    # We need to concatenate the DWI and fmap/ data (separately)
    #   before they can be fed into dwifslpreproc
    if len(dwi_image_list) == 1:
        dwifslpreproc_input = dwi_image_list[0]
    else:
        app.console('Concatenating DWI data')
        dwifslpreproc_input = 'dwifslpreproc_in.mif'
        run.command('mrcat ' + ' '.join(dwi_image_list) + ' '
                    + dwifslpreproc_input + ' -axis 3')

    # Some decisions regarding pre-processing depend on whether or not
    #   a twice-refocused sequence has been used: a single-refocused
    #   sequence may have residual eddy current distortions in b=0
    #   volumes
    dwifslpreproc_input_header = image.Header(dwifslpreproc_input)
    monopolar = "DiffusionScheme" in dwifslpreproc_input_header.keyval() \
                and dwifslpreproc_input_header \
                    .keyval()["DiffusionScheme"] == "Monopolar"

    if not fmap_image_list:
        dwifslpreproc_se_epi = ''
        dwifslpreproc_se_epi_option = ''

        # If no images in fmap/ directory, but DWIs are monopolar, then
        #   don't want to let dwifslpreproc automatically grab all of the
        #   b=0 volumes and just use those; instead, grab, for each
        #   input DWI series, just the b=0 volumes at the start of the
        #   series
        if len(dwi_image_list) > 1 and monopolar:
            try:
                bzero_image_list = []
                for dwi_image in dwi_image_list:
                    first_nonzero_volume = \
                        min([int(indices[0]) for indices in
                             image.mrinfo(dwi_image, 'shell_indices')
                             .split(' ')[1:]])
                    if not first_nonzero_volume:
                        raise MRtrixError('First DWI volume is not b=0; '
                                          'cannot utilise b=0 volumes '
                                          'prior to DWI volumes only')
                    bzero_image = os.path.splitext(dwi_image)[0] \
                                  + '_bzero.mif'
                    run.command(['mrconvert',
                                 dwi_image,
                                 bzero_image,
                                 '-coord',
                                 '3',
                                 ','.join(str(i) for i in
                                          range(0, first_nonzero_volume))])
                    bzero_image_list.append(bzero_image)
                dwifslpreproc_se_epi = 'dwifslpreproc_seepi.mif'
                run.command(['mrcat',
                             bzero_image_list,
                             dwifslpreproc_se_epi,
                             '-axis',
                             '3'])
                dwifslpreproc_se_epi_option = ' -se_epi ' \
                                                + dwifslpreproc_se_epi
            except MRtrixError:
                dwifslpreproc_se_epi = ''
                dwifslpreproc_se_epi_option = ''
                app.warn('DWIs detected as using monopolar diffusion '
                         'sensitisation, but error encountered in extracting '
                         'pre-DWI b=0 volumes; topup field estimate may be '
                         'affected by eddy current distortions in b=0 '
                         'volumes')

    elif len(fmap_image_list) == 1:
        dwifslpreproc_se_epi = fmap_image_list[0]
        dwifslpreproc_se_epi_option = ' -se_epi ' + dwifslpreproc_se_epi \
                                    + ' -align_seepi'
    else:
        # Can we do a straight concatenation?
        # If not, we need to resample all images onto a common voxel grid;
        #   the DWI grid makes most sense, since dwifslpreproc will
        #   perform that resampling itself otherwise
        dwifslpreproc_se_epi = 'dwifslpreproc_seepi.mif'
        try:
            run.command('mrcat ' + ' '.join(fmap_image_list) + ' '
                        + dwifslpreproc_se_epi + ' -axis 3')
            app.cleanup(fmap_image_list)
        except run.MRtrixCmdError:
            app.console('Unable to concatenate fmap/ data directly; '
                        'resampling images onto DWI voxel grid')
            for i in fmap_image_list:
                run.command('mrtransform ' + i
                            + ' -template ' + dwifslpreproc_input + ' '
                            + os.path.splitext(i)[0] + '_regrid.mif'
                            + ' -interp sinc')
                app.cleanup(i)
            fmap_image_list = [os.path.splitext(i)[0] + '_regrid.mif'
                               for i in fmap_image_list]
            run.command('mrcat ' + ' '.join(fmap_image_list) + ' '
                        + dwifslpreproc_se_epi + ' -axis 3')
            app.cleanup(fmap_image_list)
        dwifslpreproc_se_epi_option = ' -se_epi ' + dwifslpreproc_se_epi \
                                    + ' -align_seepi'
    if len(dwi_image_list) > 1:
        app.cleanup(dwi_image_list)

    # Step 3: Distortion correction
    app.console('Performing various geometric corrections of DWIs')
    dwifslpreproc_input_header = image.Header(dwifslpreproc_input)
    have_slice_timing = 'SliceTiming' in dwifslpreproc_input_header.keyval()
    app.debug('Have slice timing: ' + str(have_slice_timing))
    mb_factor = int(dwifslpreproc_input_header.keyval()
                    .get('MultibandAccelerationFactor', '1'))
    app.debug('Multiband factor: ' + str(mb_factor))
    if 'SliceDirection' in dwifslpreproc_input_header.keyval():
        slice_direction_code = \
            dwifslpreproc_input_header.keyval()['SliceDirection']
        if 'i' in slice_direction_code:
            num_slices = dwifslpreproc_input_header.size()[0]
        elif 'j' in slice_direction_code:
            num_slices = dwifslpreproc_input_header.size()[1]
        elif 'k' in slice_direction_code:
            num_slices = dwifslpreproc_input_header.size()[2]
        else:
            num_slices = dwifslpreproc_input_header.size()[2]
            app.warn('Error reading BIDS field \'SliceDirection\' '
                     '(value: \'' + slice_direction_code + '\'); '
                     'assuming third axis')
    else:
        num_slices = dwifslpreproc_input_header.size()[2]
    app.debug('Number of slices: ' + str(num_slices))
    mporder = 1 + int(math.ceil(num_slices/(mb_factor*4)))
    app.debug('MPorder: ' + str(mporder))

    eddy_options = []
    if shared.eddy_repol:
        eddy_options.append('--repol')
    if shared.eddy_mporder and have_slice_timing:
        eddy_options.append('--mporder=' + str(mporder))
    # Disabled: is only necessary for cohorts with very large motion,
    #    and significantly increases execution time
    #if shared.eddy_mbs:
    #    eddy_options.append('--estimate_move_by_susceptibility')
    #
    # High b-value monopolar data still has eddy current distortions
    #   in b=0 images
    # This appears to result in processing failure too regularly
    # Error messages include the following:
    # - terminate called after throwing an instance of
    #   'NEWMAT::SingularException'
    # - matrix multiplication: problem with matrix inverse;
    #   suggest to use solve() instead
    #   EDDY:::  ECScanClasses.cpp:::  void EDDY::ECScanManager::
    #   SeparateFieldOffsetFromMovement(EDDY::ScanType, EDDY::OffsetModel):
    #   Exception thrown
    # - eddy: msg=ECScanManager::set_slice_to_vol_reference:
    #   ref index out of bounds
    #if monopolar:
    #    eddy_options.append('--b0_flm=linear')
    #

    shell_asymmetries = \
        [float(value) for value in
         run.command('dirstat ' + dwifslpreproc_input + ' -output asym',
                     show=False)[0]
         .splitlines()]
    app.debug('Shell asymmetries: ' + str(shell_asymmetries))
    if any(value > 0.1 for value in shell_asymmetries):
        app.console('Utilising eddy linear second-level model due to poor '
                    'distribution of diffusion gradient direction polarities')
        eddy_options.append('--slm=linear')

    run.function(os.makedirs, 'eddyqc', show=False)
    dwifslpreproc_output = 'dwifslpreproc_out.mif' \
                           if dwifslpreproc_input == 'dwifslpreproc_in.mif' \
                           else (os.path.splitext(dwifslpreproc_input)[0]
                                 + '_preproc.mif')

    eddy_olnstd_value = 4.0 # The internal eddy default
    eddy_olnstd_option = []

    while not os.path.isfile(dwifslpreproc_output):

        try:

            # If dwifslpreproc fails due to:
            # EDDY:::  DoVolumeToVolumeRegistration: Unable to find volume
            #          with no outliers in shell 0 with b-value=549.375
            # , want to progressively increase the outlier rejection threshold
            #   until this error no longer occurs.
            eddy_all_options = eddy_options + eddy_olnstd_option
            dwifslpreproc_eddy_option = \
                ' -eddy_options " ' \
                + ' '.join(eddy_all_options) + '"' if eddy_all_options else ''

            run.command('dwifslpreproc '
                        + dwifslpreproc_input
                        + ' '
                        + dwifslpreproc_output
                        + dwifslpreproc_se_epi_option
                        + dwifslpreproc_eddy_option
                        + ' -rpe_header -eddyqc_text eddyqc/'
                        + ('' if app.DO_CLEANUP else
                           ' -scratch ' + app.SCRATCH_DIR + ' -nocleanup'))

        except run.MRtrixCmdError as e_dwifslpreproc:
            if any(item in str(e_dwifslpreproc) for item in [
                    'msg=ECScanManager::set_slice_to_vol_reference: ' \
                    'ref index out of bounds',
                    'Unable to find volume with no outliers']):
                eddy_olnstd_value += 0.5
                eddy_olnstd_option = ['--ol_nstd=' + str(eddy_olnstd_value)]
                app.warn('FSL eddy failed due to outlier rejection; '
                         're-running with increased threshold')
            else:
                raise

    app.cleanup(dwifslpreproc_input)
    app.cleanup(dwifslpreproc_se_epi)

    # Step 4: b=0-based bias field correction
    dwi_image = 'dwi_init.mif'
    if shared.dwibiascorrect_algo:
        app.console('Performing initial B1 bias field correction of DWIs')
        run.command('dwibiascorrect '
                    + shared.dwibiascorrect_algo
                    + ' '
                    + dwifslpreproc_output
                    + ' '
                    + dwi_image)
        app.cleanup(dwifslpreproc_output)
    else:
        run.function(shutil.move, dwifslpreproc_output, dwi_image)

    # Determine whether we are working with single-shell or multi-shell data
    bvalues = [
        int(round(float(value)))
        for value in image.mrinfo(dwi_image, 'shell_bvalues') \
                                 .strip().split()]
    multishell = (len(bvalues) > 2)

    # Step 5: Initial DWI brain mask
    dwi_mask_image = 'dwi_mask_init.mif'
    app.console('Performing intial DWI brain masking')
    run.command('dwi2mask ' + dwi_image + ' ' + dwi_mask_image)

    # Step 6: Combined RF estimation / CSD / mtnormalise / mask revision
    # DWI brain masking may be inaccurate due to residual bias field.
    #   We want to perform:
    #     - Response function estimation;
    #     - Multi-tissue CSD (with a lower lmax for speed);
    #     - Mtnormalise to remove any bias field;
    #     - Re-calculation of brain mask;
    #   in an iterative fashion, as all steps may influence the others.
    class Tissue(object): #pylint: disable=useless-object-inheritance
        def __init__(self, name, index):
            self.name = name
            iter_string = '_iter' + str(index)
            self.rf = 'response_' + name + iter_string + '.txt'
            self.fod_init = 'FODinit_' + name + iter_string + '.mif'
            self.fod_norm = 'FODnorm_' + name + iter_string + '.mif'

    app.console('Commencing iterative DWI bias field correction and '
                'brain masking')
    iteration = 0
    step = 'initialisation'
    dice_coefficient = 0.0
    def msg():
        return 'Iteration {0}; {1} step; previous Dice coefficient {2}' \
               .format(iteration, step, dice_coefficient)
    progress = app.ProgressBar(msg)
    for iteration in range(0, 10):
        iter_string = '_iter' + str(iteration+1)

        tissues = [Tissue('WM', iteration),
                   Tissue('GM', iteration),
                   Tissue('CSF', iteration)]

        step = 'dwi2response'
        progress.increment()
        run.command('dwi2response dhollander '
                    + dwi_image
                    + ' -mask '
                    + dwi_mask_image
                    + ' '
                    + ' '.join(tissue.rf for tissue in tissues))


        # Remove GM if we can't deal with it
        lmaxes = '4,0,0'
        if not multishell:
            app.cleanup(tissues[1].rf)
            tissues = tissues[::2]
            lmaxes = '4,0'

        step = 'dwi2fod'
        progress.increment()
        run.command('dwi2fod msmt_csd '
                    + dwi_image
                    + ' -mask ' + dwi_mask_image
                    + ' -lmax ' + lmaxes
                    + ' '
                    + ' '.join(tissue.rf + ' ' + tissue.fod_init
                               for tissue in tissues))

        step = 'mtnormalise'
        progress.increment()
        field_path = 'field' + iter_string + '.mif'
        factors_path = 'factors' + iter_string + '.txt'
        run.command('maskfilter ' + dwi_mask_image + ' erode - |'
                    + ' mtnormalise -mask -'
                    + ' -check_norm ' + field_path
                    + ' -check_factors ' + factors_path
                    + ' '
                    + ' '.join(tissue.fod_init + ' ' + tissue.fod_norm
                               for tissue in tissues))
        app.cleanup([tissue.fod_init for tissue in tissues])
        app.cleanup([tissue.fod_norm for tissue in tissues])

        # Apply both estimated bias field, and appropiate
        #   scaling factor, to DWIs
        step = 'mrcalc_dwi'
        progress.increment()
        csf_rf = matrix.load_matrix(tissues[-1].rf)
        csf_rf_bzero_lzero = csf_rf[0][0]
        app.cleanup([tissue.rf for tissue in tissues])
        balance_factors = matrix.load_vector(factors_path)
        csf_balance_factor = balance_factors[-1]
        app.cleanup(factors_path)
        scale_multiplier = (1000.0 * math.sqrt(4.0*math.pi)) / \
                           (csf_rf_bzero_lzero / csf_balance_factor)
        new_dwi_image = 'dwi' + iter_string + '.mif'
        run.command('mrcalc ' + dwi_image + ' '
                    + field_path + ' -div '
                    + str(scale_multiplier) + ' -mult '
                    + new_dwi_image)
        app.cleanup(field_path)
        app.cleanup(dwi_image)
        dwi_image = new_dwi_image

        step = 'dwi2mask'
        progress.increment()
        new_dwi_mask_image = 'dwi_mask' + iter_string + '.mif'
        run.command('dwi2mask ' + dwi_image + ' ' + new_dwi_mask_image)

        # Compare input and output masks
        step = 'mrcalc_mask'
        dwi_old_mask_count = image.statistics(dwi_mask_image,
                                              mask=dwi_mask_image).count
        dwi_new_mask_count = image.statistics(new_dwi_mask_image,
                                              mask=new_dwi_mask_image).count
        app.debug('Old mask: ' + str(dwi_old_mask_count))
        app.debug('New mask: ' + str(dwi_new_mask_count))
        dwi_mask_overlap_image = 'dwi_mask_overlap' + iter_string + '.mif'
        run.command('mrcalc '
                    + dwi_mask_image
                    + ' '
                    + new_dwi_mask_image
                    + ' -mult '
                    + dwi_mask_overlap_image)
        app.cleanup(dwi_mask_image)
        dwi_mask_image = new_dwi_mask_image
        mask_overlap_count = image.statistics(dwi_mask_overlap_image,
                                              mask=dwi_mask_overlap_image).count
        app.debug('Mask overlap: ' + str(mask_overlap_count))
        dice_coefficient = 2.0 * mask_overlap_count / \
                           (dwi_old_mask_count + dwi_new_mask_count)
        app.debug('Dice coefficient: ' + str(dice_coefficient))
        if dice_coefficient > (1.0 - 1e-3):
            progress.done()
            app.console('Exiting iterative loop due to mask convergence')
            break


    # Step 7: Crop images to reduce storage space
    #   (but leave some padding on the sides)
    dwi_cropped_image = 'dwi_crop.mif'
    dwi_cropped_mask_image = 'mask_crop.mif'
    run.command('mrgrid ' + dwi_image + ' crop ' + dwi_cropped_image
                + ' -mask ' + dwi_mask_image + ' -uniform -3')
    app.cleanup(dwi_image)
    dwi_image = dwi_cropped_image
    run.command('mrgrid ' + dwi_mask_image + ' crop ' + dwi_cropped_mask_image
                + ' -mask ' + dwi_mask_image + ' -uniform -3')
    app.cleanup(dwi_mask_image)
    dwi_mask_image = dwi_cropped_mask_image

    # Step 8: Generate target images for T1->DWI registration
    app.console('Generating contrast-matched images for '
                'inter-modal registration between DWIs and T1')
    run.command('dwiextract ' + dwi_image + ' -bzero - | '
                'mrcalc - 0.0 -max - | '
                'mrmath - mean -axis 3 dwi_meanbzero.mif')
    run.command('mrcalc 1 dwi_meanbzero.mif -div '
                + dwi_mask_image
                + ' -mult -'
                + ' | '
                + 'mrhistmatch nonlinear - '
                + T1_image
                + ' dwi_pseudoT1.mif'
                + ' -mask_input '
                + dwi_mask_image
                + ' -mask_target T1_mask.mif')
    run.command('mrcalc 1 '
                + T1_image
                + ' -div'
                + ('' if T1_is_premasked else (' T1_mask.mif -mult'))
                + ' - | '
                + 'mrhistmatch nonlinear'
                + ' - dwi_meanbzero.mif T1_pseudobzero.mif'
                + ' -mask_input T1_mask.mif'
                + ' -mask_target ' + dwi_mask_image)

    # Step 9: Perform DWI->T1 registration
    #   Note that two registrations are performed:
    #   Even though we have a symmetric registration, generation of the
    #   two histogram-matched images means that you will get slightly
    #   different answers depending on which synthesized image &
    #   original image you use
    app.console('Performing registration between DWIs and T1')
    transform_pT1_T1 = 'rigid_pseudoT1_to_T1.txt'
    transform_b0_pb0 = 'rigid_bzero_to_pseudobzero.txt'
    run.command('mrregister dwi_pseudoT1.mif '
                + T1_image
                + ' -type rigid'
                + ' -mask1 ' + dwi_mask_image
                + ' -mask2 T1_mask.mif'
                + ' -rigid ' + transform_pT1_T1)
    run.command('mrregister dwi_meanbzero.mif T1_pseudobzero.mif'
                + ' -type rigid '
                + ' -mask1 ' + dwi_mask_image
                + ' -mask2 T1_mask.mif'
                + ' -rigid ' + transform_b0_pb0)
    app.cleanup('dwi_meanbzero.mif')

    # Step 10: Perform DWI->T1 transformation
    # In this scenario, we're going to transform the DWI data to the T1
    #   rather than the other way around, since the T1 is more likely to
    #   be used as a common reference across multiple analysis pipelines,
    #   and we're transforming DWIs rather than FODs
    transform_average = 'rigid_dwi_to_T1.txt'
    run.command('transformcalc '
                + transform_pT1_T1
                + ' '
                + transform_b0_pb0
                + ' average '
                + transform_average)
    app.cleanup(transform_pT1_T1)
    app.cleanup(transform_b0_pb0)
    transformed_dwi_image = os.path.splitext(dwi_image)[0] \
                            + '_transform.mif'
    transformed_dwi_mask_image = os.path.splitext(dwi_mask_image)[0] \
                                 + '_transform.mif'
    run.command('mrtransform '
                + dwi_image
                + ' '
                + transformed_dwi_image
                + ' -linear '
                + transform_average)
    app.cleanup(dwi_image)
    dwi_image = transformed_dwi_image
    run.command('mrtransform '
                + dwi_mask_image
                + ' '
                + transformed_dwi_mask_image
                + ' -linear '
                + transform_average)
    app.cleanup(dwi_mask_image)
    app.cleanup(transform_average)
    dwi_mask_image = transformed_dwi_mask_image



    # Processing completed; export
    app.console('Processing completed for session "'
                + session_label
                + '"; writing results to output directory')
    if os.path.exists(output_subdir):
        run.function(shutil.rmtree, output_subdir)
    run.function(os.makedirs, output_subdir)
    run.function(os.makedirs, os.path.join(output_subdir, 'anat'))
    run.function(os.makedirs, os.path.join(output_subdir, 'dwi'))
    run.command('mrconvert '
                + dwi_image
                + ' '
                + os.path.join(output_subdir,
                               'dwi',
                               session_label + '_desc-preproc_dwi.nii.gz')
                + ' -export_grad_fsl '
                + os.path.join(output_subdir,
                               'dwi',
                               session_label + '_desc-preproc_dwi.bvec')
                + ' '
                + os.path.join(output_subdir,
                               'dwi',
                               session_label + '_desc-preproc_dwi.bval')
                + ' -strides +1,+2,+3,+4')
    with open(os.path.join(output_subdir,
                           'dwi',
                           session_label + '_desc-preproc_dwi.json'),
              'w') as out_dwi_json_file:
        json.dump(OUT_DWI_JSON_DATA, out_dwi_json_file)
    run.command('mrconvert '
                + dwi_mask_image
                + ' '
                + os.path.join(output_subdir,
                               'dwi',
                               session_label + '_desc-brain_mask.nii.gz')
                + ' -datatype uint8'
                + ' -strides +1,+2,+3')
    # Even if EddyQC software is not installed, a directory is still
    #   generated containing some eddy outputs
    run.function(shutil.copytree,
                 'eddyqc',
                 os.path.join(output_subdir,
                              'dwi',
                              'eddyqc'))

    run.command('mrconvert '
                + T1_image
                + ' '
                + os.path.join(output_subdir,
                               'anat',
                               session_label + '_desc-preproc_T1w.nii.gz')
                + ' -strides +1,+2,+3')
    T1_json_data = {"SkullStripped": T1_is_premasked}
    with open(os.path.join(output_subdir,
                           'anat',
                           session_label + '_desc-preproc_T1w.json'),
              'w') as T1_json_file:
        json.dump(T1_json_data, T1_json_file)
    run.command('mrconvert T1_mask.mif '
                + os.path.join(output_subdir,
                               'anat',
                               session_label + '_desc-brain_mask.nii.gz')
                + ' -datatype uint8'
                + ' -strides +1,+2,+3')

    # Manually wipe and zero the scratch directory
    #   (since we might be processing more than one subject)
    os.chdir(cwd)
    if app.DO_CLEANUP:
        app.console('Deleting scratch directory ' + app.SCRATCH_DIR)
        # Can't use run.function() here;
        #   it'll try to write to the log file that resides
        #   in the scratch directory just deleted
        app.cleanup(app.SCRATCH_DIR)
    elif output_verbosity == 4:
        app.console('Copying scratch directory to output location')
        run.function(shutil.copytree,
                     app.SCRATCH_DIR,
                     os.path.join(output_subdir, 'scratch'))
    else:
        app.console('Contents of scratch directory kept; '
                    'location: ' + app.SCRATCH_DIR)
    app.SCRATCH_DIR = ''

# End of run_preproc() function




















def run_participant(bids_dir, session, shared,
                    t1w_preproc_path, output_verbosity, output_app_dir):

    session_label = '_'.join(session)
    output_analysis_level_path = os.path.join(output_app_dir, 'participant')
    output_subdir = os.path.join(output_analysis_level_path, *session)

    if os.path.exists(output_subdir):
        app.warn('Output directory for session "' + session_label + '" '
                 'already exists; all contents will be erased when this '
                 'execution completes')

    app.make_scratch_dir()

    # Check paths of individual output files before script completion
    #   by building a database of what files are to be written to output
    parc_string = '_desc-' + shared.parcellation
    OutputItem = \
        namedtuple(
            'OutputItem',
            'is_image min_verbosity needs_multishell options path')
    output_items = {
        'response_wm.txt': \
            OutputItem(False, 1, False, None,
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-WM_response.txt')),
        'T1_mask.mif': \
            OutputItem(True, 1, False,
                       '-strides +1,+2,+3 -datatype uint8',
                       os.path.join('anat',
                                    session_label
                                    + '_desc-brain_mask.nii.gz')),
        '5TT.mif': \
            OutputItem(True, 2, False,
                       '-strides +1,+2,+3,+4',
                       os.path.join('anat',
                                    session_label
                                    + '_desc-5tt_probseg.nii.gz')),
        '5TT.json': \
            OutputItem(False, 2, False, None,
                       os.path.join('anat',
                                    session_label
                                    + '_desc-5tt_probseg.json')),
        'vis.mif': \
            OutputItem(True, 2, False, '-strides +1,+2,+3',
                       os.path.join('anat',
                                    session_label
                                    + '_desc-vis_probseg.nii.gz')),
        'FOD_WM.mif': \
            OutputItem(True, 2, False, '-strides +1,+2,+3,+4',
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-WM_ODF.nii.gz')),
        'response_gm.txt': \
            OutputItem(False, 2, True, None,
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-GM_response.txt')),
        'response_csf.txt': \
            OutputItem(False, 2, True, None,
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-CSF_response.txt')),
        'FOD_GM.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-GM_ODF.nii.gz')),
        'FOD_CSF.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-CSF_ODF.nii.gz')),
        'tissues.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       os.path.join('dwi',
                                    session_label
                                    + '_tissue-all_probseg.nii.gz'))
    }

    if shared.parcellation != 'none':
        output_items['connectome.csv'] = \
            OutputItem(False, 1, False, None,
                       os.path.join('connectome',
                                    session_label
                                    + parc_string
                                    + '_connectome.csv'))
        output_items['mu.txt'] = \
            OutputItem(False, 1, False, None,
                       os.path.join('tractogram',
                                    session_label + '_mu.txt'))
        output_items['parc.mif'] = \
            OutputItem(True, 2, False, '-strides +1,+2,+3',
                       os.path.join('anat',
                                    session_label
                                    + parc_string
                                    + '_dseg.nii.gz'))
        output_items['meanlength.csv'] = \
            OutputItem(False, 2, False, None,
                       os.path.join('connectome',
                                    session_label
                                    + parc_string
                                    + '_meanlength.csv'))
        output_items['assignments.csv'] = \
            OutputItem(False, 3, False, None,
                       os.path.join('connectome',
                                    session_label
                                    + parc_string
                                    + '_assignments.csv'))
        output_items['nodes_smooth.obj'] = \
            OutputItem(False, 3, False, None,
                       os.path.join('anat',
                                    session_label
                                    + parc_string
                                    + '_dseg.obj'))
        output_items['exemplars.tck'] = \
            OutputItem(False, 3, False, None,
                       os.path.join('connectome',
                                    session_label
                                    + parc_string
                                    + '_exemplars.tck'))
        output_items['parcRGB.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3,+4',
                       os.path.join('anat',
                                    session_label
                                    + parc_string
                                    + '_desc-rgb_dseg.nii.gz'))

    if shared.streamlines or shared.parcellation != 'none':
        output_items['tractogram.tck'] = \
            OutputItem(False, 3, False, None,
                       os.path.join('tractogram',
                                    session_label + '_tractogram.tck'))
        output_items['weights.csv'] = \
            OutputItem(False, 3, False, None,
                       os.path.join('tractogram',
                                    session_label + '_weights.csv'))
        output_items['tdi_dwi.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       os.path.join('tractogram',
                                    session_label
                                    + '_space-dwi_tdi.nii.gz'))
        output_items['tdi_t1.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       os.path.join('tractogram',
                                    session_label
                                    + '_space-T1w_tdi.nii.gz'))
        output_items['tdi_hires.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       os.path.join('tractogram',
                                    session_label
                                    + '_space-superres_tdi.nii.gz'))

    subdirs_to_make = ['tractogram']
    if shared.parcellation != 'none':
        subdirs_to_make.insert(0, 'connectome')

    app.make_scratch_dir()


    def do_import(import_path):
        in_dwi_path = os.path.join(import_path,
                                   'dwi',
                                   '*_dwi.nii*')
        in_dwi_image_list = glob.glob(in_dwi_path)
        if len(in_dwi_image_list) > 1:
            raise MRtrixError('To run participant-level analysis, '
                              + 'input directory should contain only one '
                              + 'DWI image file; session "'
                              + session_label
                              + '" loaded from "'
                              + import_path
                              + '" contains '
                              + str(len(in_dwi_image_list)))
        in_dwi_path = in_dwi_image_list[0]
        if not '_desc-preproc_' in in_dwi_path:
            raise MRtrixError('Input DWI image "'
                              + in_dwi_path
                              + '" loaded from "'
                              + import_path
                              + '" not flagged as pre-processed data')
        in_dwi_path_prefix = in_dwi_path.split(os.extsep)[0]
        # Don't look for bvec / bval in a lower directory in this case
        in_bvec_path = in_dwi_path_prefix + '.bvec'
        in_bval_path = in_dwi_path_prefix + '.bval'
        if not os.path.isfile(in_bvec_path) \
            or not os.path.isfile(in_bval_path):
            raise MRtrixError('Did not find bvec / bval pair '
                              + 'corresponding to image '
                              + in_dwi_path
                              + ' (expected locations: "'
                              + in_bvec_path
                              + '" "'
                              + in_bval_path
                              + '")')
        # JSON isn't compulsory in this case
        in_dwi_json_path = in_dwi_path_prefix + 'json'
        in_dwi_json_import_option = ' -json_import ' + in_dwi_json_path \
                                    if os.path.isfile(in_dwi_json_path) \
                                    else ''
        # Is there a mask present?
        in_dwi_mask_image_list = \
            glob.glob(os.path.join(output_subdir,
                                   'dwi',
                                   '*_desc-brain*_mask.nii*'))
        if len(in_dwi_mask_image_list) > 1:
            raise MRtrixError('More than one DWI mask found for session "'
                              + session_label
                              + '"')
        in_dwi_mask_path = in_dwi_mask_image_list[0] \
                        if in_dwi_mask_image_list \
                        else None
        if not in_dwi_mask_path:
            output_items['dwi_mask.mif'] = \
                OutputItem(True, 1, False,
                           '-strides +1,+2,+3 -datatype uint8',
                           os.path.join('dwi',
                                        session_label
                                        + '_desc-brain_mask.nii.gz'))

        app.console('Importing pre-processed data into scratch directory')

        run.command('mrconvert '
                    + in_dwi_path
                    + ' '
                    + path.to_scratch('dwi.mif')
                    + ' -fslgrad ' + in_bvec_path + ' ' + in_bval_path
                    + in_dwi_json_import_option
                    + ' -strides 0,0,0,1')

        if in_dwi_mask_path:
            run.command('mrconvert '
                        + in_dwi_mask_path
                        + ' '
                        + path.to_scratch('dwi_mask.mif')
                        + ' -datatype bit')

        get_t1w_preproc_images(import_path,
                               session,
                               shared.t1w_shared,
                               t1w_preproc_path)
        T1_is_premasked = os.path.isfile(path.to_scratch('T1_premasked.mif',
                                                         False))
        if shared.do_freesurfer and T1_is_premasked:
            raise MRtrixError('Cannot execute FreeSurfer for obtaining '
                              'parcellation: input T1-weighted image is '
                              'already skull-stripped')
    # End of do_import() function


    # We first make an attempt at loading all requisite data from
    #   "bids_dir" (since the user may have used that path to request
    #   that the pre-processed data be utilised from some path other than
    #   "mrtrix3_connectome/preproc/"); if that doesn't work, we wipe the
    #   scratch directory and try again based on the latter
    try:
        do_import(bids_dir)
    except MRtrixError as e_frombids:
        for item in os.listdir(app.SCRATCH_DIR):
            os.remove(os.path.join(app.SCRATCH_DIR, item))
        try:
            do_import(os.path.join(output_app_dir, 'preproc'))
        except MRtrixError as e_fromoutput:
            err = 'Unable to import requisite pre-processed data from ' \
                  'either specified input directory or MRtrix3_connectome ' \
                  'output directory\n\n'
            err += 'Error when loading from "' + bids_dir + '":\n'
            err += str(e_frombids) + '\n'
            err += 'Error when loading from "' + output_app_dir + '":\n'
            err += str(e_fromoutput) + '\n'
            raise MRtrixError(err)


    cwd = os.getcwd()
    app.goto_scratch_dir()


    # T1-weighted data are always written to output directory regardless;
    #   output paths can only be constructed now
    T1_is_premasked = os.path.isfile(path.to_scratch('T1_premasked.mif',
                                                     False))
    T1_image = 'T1_premasked.mif' if T1_is_premasked else 'T1.mif'
    output_items[T1_image] = \
        OutputItem(True, 1, False, ' -strides +1,+2,+3',
                   os.path.join(output_subdir,
                                'anat',
                                session_label + '_desc-preproc_T1w.nii.gz'))
    T1_json_path = os.path.splitext(T1_image)[0] + '.json'
    output_items[T1_json_path] = \
        OutputItem(False, 1, False, None,
                   os.path.join(output_subdir,
                                'anat',
                                session_label + '_desc-preproc_T1w.json'))
    T1_json_data = {"SkullStripped": T1_is_premasked}

    # Before we can begin: Are there any data we require
    #   that were not imported from the output directory?
    if not os.path.isfile('dwi_mask.mif'):
        app.console('Generating DWI brain mask '
                    '(was not already present in pre-processing directory)')
        run.command('dwi2mask dwi.mif dwi_mask.mif')

    # Step 1: Estimate response functions for spherical deconvolution
    app.console('Estimating tissue response functions for '
                'spherical deconvolution')
    run.command('dwi2response dhollander dwi.mif '
                'response_wm.txt response_gm.txt response_csf.txt '
                '-mask dwi_mask.mif')

    # Determine whether we are working with single-shell or multi-shell data
    bvalues = [
        int(round(float(value)))
        for value in image.mrinfo('dwi.mif', 'shell_bvalues') \
                                 .strip().split()]
    multishell = (len(bvalues) > 2)

    # Step 2: Perform spherical deconvolution
    #   Don't even use a processing mask:
    #     ACT should be responsible for stopping streamlines before they
    #     reach the edge of the DWI mask
    #   Also means that any subsequent manual use of the FOD
    #     images can't possibly be detrimentally affected by
    #     bad masking

    app.console('Estimating '
                + (' multi-tissue ODF images'
                   if multishell
                   else 'Fibre Orientation Distribution image'))
    # TODO Update to use similar code to preproc?
    # Would have consequences for group-level analysis...
    if multishell:
        run.command('dwi2fod msmt_csd dwi.mif '
                    'response_wm.txt FOD_WM.mif '
                    'response_gm.txt FOD_GM.mif '
                    'response_csf.txt FOD_CSF.mif '
                    '-lmax 10,0,0')
        run.command('mrconvert FOD_WM.mif - -coord 3 0 | '
                    'mrcat FOD_CSF.mif FOD_GM.mif - tissues.mif -axis 3')
    else:
        # Still use the msmt_csd algorithm with single-shell data:
        #   Use hard non-negativity constraint
        # Also incorporate the CSF response to provide some fluid attenuation
        run.command('dwi2fod msmt_csd dwi.mif '
                    'response_wm.txt FOD_WM.mif '
                    'response_csf.txt FOD_CSF.mif '
                    '-lmax 10,0')
        app.cleanup('FOD_CSF.mif')

    # Step 3: Generate 5TT image for ACT
    # Use T1 brain mask generated from elsewhere:
    #   don't particularly trust the raw "bet" call inside
    #   5ttgen fsl
    app.console('Generating five-tissue-type (5TT) image for '
                'Anatomically-Constrained Tractography (ACT)')
    run.command('5ttgen fsl '
                + T1_image
                + ' 5TT.mif'
                + (' -premasked' \
                   if T1_is_premasked \
                   else (' -mask T1_mask.mif')))
    if output_verbosity > 1:
        with open('5TT.json', 'w') as out_5tt_json_file:
            json.dump(OUT_5TT_JSON_DATA, out_5tt_json_file)
        run.command('5tt2vis 5TT.mif vis.mif')

    # Step 4: Generate the grey matter parcellation
    #   The necessary steps here will vary significantly depending on
    #   the parcellation scheme selected
    if shared.do_freesurfer:
        app.console('Getting grey matter parcellation in '
                    'subject space using FreeSurfer')

        # Since we're instructing recon-all to use a different subject
        #   directory, we need to construct softlinks to a number of
        #   directories provided by FreeSurfer that recon-all will
        #   expect to find in the same directory as the overridden
        #   subject path
        subdirs = ['fsaverage', 'lh.EC_average', 'rh.EC_average']
        if shared.parcellation in ['yeo7fs', 'yeo17fs']:
            subdirs.append('fsaverage5')
        for subdir in subdirs:
            run.function(shared.freesurfer_template_link_function,
                         os.path.join(shared.freesurfer_subjects_dir, subdir),
                         subdir)

        # Run FreeSurfer pipeline on this subject's T1 image
        run.command('recon-all -sd ' + app.SCRATCH_DIR + ' -subjid freesurfer '
                    '-i T1_raw.nii')
        run.command('recon-all -sd ' + app.SCRATCH_DIR + ' -subjid freesurfer '
                    '-all' + shared.reconall_multithread_options)

        # Grab the relevant parcellation image and
        #   target lookup table for conversion
        parc_image_path = os.path.join('freesurfer', 'mri')
        if shared.parcellation == 'desikan':
            parc_image_path = os.path.join(parc_image_path,
                                           'aparc+aseg.mgz')
        elif shared.parcellation == 'destrieux':
            parc_image_path = os.path.join(parc_image_path,
                                           'aparc.a2009s+aseg.mgz')
        else:
            # Non-standard parcellations are not applied as part of
            #   the recon-all command; need to explicitly map them to
            #   the subject
            # This requires SUBJECTS_DIR to be set;
            #   commands don't have a corresponding -sd option like recon-all
            env = run.shared.env
            env['SUBJECTS_DIR'] = app.SCRATCH_DIR
            if shared.parcellation == 'brainnetome246fs':
                for index, hemi in enumerate(['l', 'r']):
                    run.command(
                        'mris_ca_label'
                        + ' -l ' + os.path.join('freesurfer',
                                                'label',
                                                hemi + 'h.cortex.label')
                        + ' freesurfer ' + hemi + 'h '
                        + os.path.join('freesurfer',
                                       'surf',
                                       hemi + 'h.sphere.reg')
                        + ' '
                        + shared.brainnetome_cortex_gcs_paths[index]
                        + ' '
                        + os.path.join('freesurfer',
                                       'label',
                                       hemi + 'h.BN_Atlas.annot'),
                        env=env)
                    run.command(
                        'mri_label2vol'
                        + ' --annot ' + os.path.join('freesurfer',
                                                     'label',
                                                     hemi + 'h.BN_Atlas.annot')
                        + ' --temp ' + os.path.join('freesurfer',
                                                    'mri',
                                                    'brain.mgz')
                        + ' --o ' + os.path.join('freesurfer',
                                                 'mri',
                                                 hemi + 'h.BN_Atlas.mgz')
                        + ' --subject freesurfer'
                        + ' --hemi ' + hemi + 'h'
                        + ' --identity'
                        + ' --proj frac 0 1 .1',
                        env=env)
                run.command(
                    'mri_ca_label '
                    + os.path.join('freesurfer',
                                   'mri',
                                   'brain.mgz')
                    + ' '
                    + os.path.join('freesurfer',
                                   'mri',
                                   'transforms',
                                   'talairach.m3z')
                    + ' '
                    + shared.brainnetome_sgm_gca_path
                    + ' '
                    + os.path.join('freesurfer',
                                   'mri',
                                   'BN_Atlas_subcortex.mgz'),
                    env=env)
                parc_image_path = os.path.join(parc_image_path,
                                               'aparc.BN_Atlas+aseg.mgz')
                # Need to deal with prospect of overlapping mask labels
                # - Any overlap between the two hemisphere ribbons
                #   = set to zero
                # - Any overlap between cortex and sub-cortical
                #   = retain cortex
                run.command('mrcalc '
                            + ' '.join([os.path.join('freesurfer',
                                                     'mri',
                                                     hemi + 'h.BN_Atlas.mgz')
                                        for hemi in ['l', 'r']])
                            + ' -mult cortex_overlap.mif'
                            + ' -datatype bit')
                run.command('mrcalc '
                            + ' '.join([os.path.join('freesurfer',
                                                     'mri',
                                                     hemi + 'h.BN_Atlas.mgz')
                                        for hemi in ['l', 'r']])
                            + ' -add '
                            + os.path.join('freesurfer',
                                           'mri',
                                           'BN_Atlas_subcortex.mgz')
                            + ' -mult sgm_overlap.mif'
                            + ' -datatype bit')
                run.command('mrcalc '
                            + ' '.join([os.path.join('freesurfer',
                                                     'mri',
                                                     hemi + 'h.BN_Atlas.mgz')
                                        for hemi in ['l', 'r']])
                            + ' -add 1.0 cortex_overlap.mif -sub -mult '
                            + os.path.join('freesurfer',
                                           'mri',
                                           'BN_Atlas_subcortex.mgz')
                            + ' 1.0 sgm_overlap.mif -sub -mult -add '
                            + parc_image_path)
                app.cleanup('cortex_overlap.mif')
                app.cleanup('sgm_overlap.mif')

            elif shared.parcellation == 'hcpmmp1':
                parc_image_path = os.path.join(parc_image_path,
                                               'aparc.HCPMMP1+aseg.mgz')
                for index, hemi in enumerate(['l', 'r']):
                    run.command('mri_surf2surf '
                                '--srcsubject fsaverage '
                                '--trgsubject freesurfer '
                                '--hemi ' + hemi + 'h '
                                '--sval-annot '
                                + shared.hcpmmp1_annot_paths[index]
                                + ' --tval '
                                + os.path.join('freesurfer',
                                               'label',
                                               hemi + 'h.HCPMMP1.annot'),
                                env=env)
                run.command('mri_aparc2aseg '
                            '--s freesurfer '
                            '--old-ribbon '
                            '--annot HCPMMP1 '
                            '--o ' + parc_image_path,
                            env=env)
            elif shared.parcellation in ['yeo7fs', 'yeo17fs']:
                num = '7' if shared.parcellation == 'yeo7fs' else '17'
                parc_image_path = os.path.join(parc_image_path,
                                               'aparc.Yeo' + num + '+aseg.mgz')
                for index, hemi in enumerate(['l', 'r']):
                    run.command('mri_surf2surf '
                                '--srcsubject fsaverage5 '
                                '--trgsubject freesurfer '
                                '--hemi ' + hemi + 'h '
                                '--sval-annot ' + shared.yeo_annot_paths[index]
                                + ' --tval '
                                + os.path.join('freesurfer',
                                               'label',
                                               hemi + 'h.Yeo'
                                               + num + '.annot'),
                                env=env)
                run.command('mri_aparc2aseg '
                            '--s freesurfer '
                            '--old-ribbon '
                            '--annot Yeo' + num + ' '
                            '--o ' + parc_image_path,
                            env=env)
            else:
                assert False

        if shared.mrtrix_lut_file:
            # If necessary:
            # Perform the index conversion
            run.command('labelconvert ' + parc_image_path + ' '
                        + shared.parc_lut_file + ' '
                        + shared.mrtrix_lut_file
                        + ' parc_init.mif')
            # Fix the sub-cortical grey matter parcellations using FSL FIRST
            run.command('labelsgmfix parc_init.mif T1_raw.nii '
                        + shared.mrtrix_lut_file
                        + ' parc.mif')
            app.cleanup('T1_raw.nii')
            app.cleanup('parc_init.mif')
        else:
            # Non-standard sub-cortical parcellation;
            #   labelsgmfix not applicable
            run.command('mrconvert ' + parc_image_path + ' parc.mif '
                        '-datatype uint32')
        app.cleanup('freesurfer')


    elif shared.do_mni:
        app.console('Registering to MNI template and transforming grey '
                    'matter parcellation back to subject space')

        # Use non-dilated brain masks for performing
        #   histogram matching & linear registration
        T1_histmatched_path = 'T1_histmatch.nii'
        run.command('mrhistmatch linear '
                    + T1_image
                    + ' '
                    + shared.template_image_path
                    + ' -mask_input T1_mask.mif'
                    + ' -mask_target ' + shared.template_mask_path
                    + ' - |'
                    + ' mrconvert - '
                    + T1_histmatched_path
                    + ' -strides '
                    + ('-1,+2,+3' \
                       if shared.template_registration_software == 'fsl' \
                       else '+1,+2,+3'))

        assert shared.template_registration_software
        if shared.template_registration_software == 'ants':

            # Use ANTs SyN for registration to template
            # From Klein et al., NeuroImage 2009:
            run.command('ANTS 3 '
                        + '-m PR['
                        + shared.template_image_path
                        + ', '
                        + T1_histmatched_path
                        + ', 1, 2]'
                        + ' -o ANTS'
                        + ' -r Gauss[2,0]'
                        + ' -t SyN[0.5]'
                        + ' -i 30x99x11'
                        + ' --use-Histogram-Matching')
            transformed_atlas_path = 'atlas_transformed.nii'
            run.command('WarpImageMultiTransform 3 '
                        + shared.parc_image_path
                        + ' '
                        + transformed_atlas_path
                        + ' -R '
                        + T1_histmatched_path
                        + ' -i ANTSAffine.txt ANTSInverseWarp.nii'
                        + ' --use-NN')
            app.cleanup(glob.glob('ANTSWarp.nii*'))
            app.cleanup(glob.glob('ANTSInverseWarp.nii*'))
            app.cleanup('ANTSAffine.txt')

        elif shared.template_registration_software == 'fsl':

            # Subject T1, brain masked; for flirt -in
            if T1_is_premasked:
                flirt_in_path = T1_histmatched_path
            else:
                flirt_in_path = \
                    os.path.splitext(T1_histmatched_path)[0] \
                    + '_masked.nii'
                run.command('mrcalc '
                            + T1_histmatched_path
                            + ' T1_mask.mif -mult '
                            + flirt_in_path)
            # Template T1, brain masked; for flirt -ref
            flirt_ref_path = 'template_masked.nii'
            run.command('mrcalc '
                        + shared.template_image_path
                        + ' '
                        + shared.template_mask_path
                        + ' -mult '
                        + flirt_ref_path
                        + ' -strides -1,+2,+3')
            # Now have data required to run flirt
            run.command(shared.flirt_cmd
                        + ' -ref ' + flirt_ref_path
                        + ' -in ' + flirt_in_path
                        + ' -omat T1_to_template.mat'
                        + ' -dof 12'
                        + ' -cost leastsq')
            if not T1_is_premasked:
                app.cleanup(flirt_in_path)
            app.cleanup(flirt_ref_path)

            # If possible, use dilated brain masks for non-linear
            #   registration to mitigate mask edge effects;
            #   if T1-weighted image is premasked, can't do this
            fnirt_in_path = T1_histmatched_path
            fnirt_ref_path = shared.template_image_path
            if T1_is_premasked:
                fnirt_in_mask_path = 'T1_mask.nii'
                run.command('mrconvert T1_mask.mif '
                            + fnirt_in_mask_path
                            + ' -strides -1,+2,+3')
                fnirt_ref_mask_path = shared.template_mask_path
            else:
                fnirt_in_mask_path = 'T1_mask_dilated.nii'
                run.command('maskfilter T1_mask.mif dilate -'
                            + ' -npass 3'
                            + ' |'
                            + ' mrconvert - '
                            + fnirt_in_mask_path
                            + ' -strides -1,+2,+3')
                fnirt_ref_mask_path = 'template_mask_dilated.nii'
                run.command('maskfilter '
                            + shared.template_mask_path
                            + ' dilate '
                            + fnirt_ref_mask_path
                            + '-npass 3')

            run.command(shared.fnirt_cmd
                        + ' --config=' + shared.fnirt_config_basename
                        + ' --ref=' + fnirt_ref_path
                        + ' --in=' + fnirt_in_path
                        + ' --aff=T1_to_template.mat'
                        + ' --refmask=' + fnirt_ref_mask_path
                        + ' --inmask=' + fnirt_in_mask_path
                        + ' --cout=T1_to_template_warpcoef.nii')
            app.cleanup(fnirt_in_mask_path)
            if not T1_is_premasked:
                app.cleanup(fnirt_ref_mask_path)
            app.cleanup('T1_to_template.mat')
            fnirt_warp_subject2template_path = \
                fsl.find_image('T1_to_template_warpcoef')

            # Use result of registration to transform atlas
            #   parcellation to subject space
            run.command(shared.invwarp_cmd
                        + ' --ref=' + T1_histmatched_path
                        + ' --warp=' + fnirt_warp_subject2template_path
                        + ' --out=template_to_T1_warpcoef.nii')
            app.cleanup(fnirt_warp_subject2template_path)
            fnirt_warp_template2subject_path = \
                fsl.find_image('template_to_T1_warpcoef')
            run.command(shared.applywarp_cmd
                        + ' --ref=' + T1_histmatched_path
                        + ' --in=' + shared.parc_image_path
                        + ' --warp=' + fnirt_warp_template2subject_path
                        + ' --out=atlas_transformed.nii'
                        + ' --interp=nn')
            app.cleanup(fnirt_warp_template2subject_path)
            transformed_atlas_path = fsl.find_image('atlas_transformed')

        app.cleanup(T1_histmatched_path)

        if shared.parc_lut_file and shared.mrtrix_lut_file:
            run.command(['labelconvert',
                         transformed_atlas_path,
                         shared.parc_lut_file,
                         shared.mrtrix_lut_file,
                         'parc.mif'])
        else:
            # Not all parcellations need to go through the labelconvert step;
            #   they may already be numbered incrementally from 1
            run.command(['mrconvert',
                         transformed_atlas_path,
                         'parc.mif'])
        app.cleanup(transformed_atlas_path)


    if output_verbosity > 2 and shared.parcellation != 'none':
        if shared.mrtrix_lut_file:
            label2colour_lut_option = ' -lut ' + shared.mrtrix_lut_file
        elif shared.parc_lut_file:
            label2colour_lut_option = ' -lut ' + shared.parc_lut_file
        else:
            # Use random colouring if no LUT available, but
            #   still generate the image
            label2colour_lut_option = ''
        run.command('label2colour parc.mif parcRGB.mif'
                    + label2colour_lut_option)

    # If no parcellation is requested, it is still possible to
    #   generate a whole-brain tractogram by explicitly providing
    #   the -streamlines option
    num_streamlines = None
    if shared.streamlines:
        num_streamlines = shared.streamlines
    elif shared.parcellation != 'none':
        # If not manually specified, determine the appropriate
        #   number of streamlines based on the number of nodes
        #   in the parcellation:
        #   mean edge weight of 1,000 streamlines
        num_nodes = int(image.statistics('parc.mif').max)
        num_streamlines = 500 * num_nodes * (num_nodes-1)
    if num_streamlines:

        # Step 5: Generate the tractogram
        app.console('Performing whole-brain fibre-tracking')
        tractogram_filepath = 'tractogram_' + str(num_streamlines) + '.tck'
        run.command('tckgen FOD_WM.mif ' + tractogram_filepath + ' '
                    '-act 5TT.mif -backtrack -crop_at_gmwmi '
                    '-maxlength 250 '
                    '-power 0.33 '
                    '-select ' + str(num_streamlines) + ' '
                    '-seed_dynamic FOD_WM.mif')

        # Step 6: Use SIFT2 to determine streamline weights
        app.console('Running the SIFT2 algorithm to assign '
                    'weights to individual streamlines')
        fd_scale_gm_option = ''
        if not multishell:
            fd_scale_gm_option = ' -fd_scale_gm'
        # If SIFT2 fails, reduce number of streamlines and try again
        while num_streamlines:
            try:
                run.command('tcksift2 '
                            + tractogram_filepath
                            + ' FOD_WM.mif weights.csv '
                            + '-act 5TT.mif '
                            + '-out_mu mu.txt'
                            + fd_scale_gm_option)
                break
            except run.MRtrixCmdError:
                app.warn('SIFT2 failed, likely due to running out of RAM; '
                         'reducing number of streamlines and trying again')
                num_streamlines = int(num_streamlines // 2)
                new_tractogram_filepath = \
                    'tractogram_' + str(num_streamlines) + '.tck'
                run.command('tckedit '
                            + tractogram_filepath + ' '
                            + new_tractogram_filepath
                            + ' -number ' + str(num_streamlines))
                app.cleanup(tractogram_filepath)
                tractogram_filepath = new_tractogram_filepath
        if not num_streamlines:
            raise MRtrixError('Unable to run SIFT2 algorithm for '
                              'any number of streamlines')
        run.function(shutil.move, tractogram_filepath, 'tractogram.tck')
        tractogram_filepath = 'tractogram.tck'


        if output_verbosity > 2:
            # Generate TDIs:
            # - A TDI at DWI native resolution, with SIFT mu scaling,
            #   and precise mapping
            #     (for comparison to WM ODF l=0 term, to
            #     verify that SIFT2 has worked correctly)
            app.console('Producing Track Density Images (TDIs)')
            with open('mu.txt', 'r') as f:
                mu = float(f.read())
            # In the space of the DWI image
            run.command('tckmap tractogram.tck -'
                        ' -tck_weights_in weights.csv'
                        ' -template FOD_WM.mif'
                        ' -precise'
                        ' | '
                        'mrcalc - ' + str(mu) + ' -mult tdi_dwi.mif')
            # In the space of the T1-weighted image
            run.command('tckmap tractogram.tck -'
                        + ' -tck_weights_in weights.csv'
                        + ' -template ' + T1_image
                        + ' -precise'
                        + ' | '
                        + 'mrcalc - ' + str(mu) + ' -mult tdi_T1.mif')
            # - Conventional TDI at super-resolution
            #   (mostly just because we can)
            run.command('tckmap tractogram.tck tdi_hires.mif'
                        ' -tck_weights_in weights.csv'
                        ' -vox 0.25'
                        ' -datatype uint16')


    if shared.parcellation != 'none':
        # Step 7: Generate the connectome
        #   Also get the mean length for each edge;
        #   this is the most likely alternative contrast to be useful
        app.console('Combining whole-brain tractogram with grey matter '
                    'parcellation to produce the connectome')
        assignment_option = \
            ' -assignment_radial_search 5' \
            if shared.parcellation in ['yeo7mni', 'yeo17mni'] \
            else ''
        run.command('tck2connectome tractogram.tck'
                    ' parc.mif connectome.csv'
                    ' -tck_weights_in weights.csv'
                    ' -out_assignments assignments.csv'
                    + assignment_option)
        run.command('tck2connectome tractogram.tck'
                    ' parc.mif meanlength.csv'
                    ' -tck_weights_in weights.csv'
                    ' -scale_length'
                    ' -stat_edge mean'
                    + assignment_option)

        if output_verbosity > 2:
            # Produce additional data that can be used for
            #   visualisation within mrview's connectome toolbar
            app.console('Generating geometric data for '
                        'enhanced connectome visualisation')
            run.command('connectome2tck '
                        + tractogram_filepath
                        + ' assignments.csv exemplars.tck'
                        + ' -tck_weights_in weights.csv'
                        + ' -exemplars parc.mif'
                        + ' -files single')
            run.command('label2mesh parc.mif nodes.obj')
            run.command('meshfilter nodes.obj smooth nodes_smooth.obj')
            app.cleanup('nodes.obj')



    # Prepare output path for writing
    app.console('Processing for session "' + session_label
                + '" completed; writing results to output directory')
    for subdir in subdirs_to_make:
        full_subdir_path = os.path.join(output_subdir, subdir)
        if os.path.exists(full_subdir_path):
            run.function(shutil.rmtree, full_subdir_path)
        run.function(os.makedirs, full_subdir_path)
    if not os.path.isdir(os.path.join(output_subdir, 'anat')):
        run.function(os.makedirs, os.path.join(output_subdir, 'anat'))

    # Generate a copy of the lookup table file:
    #   - Use the post-labelconvert file if it's used;
    #     otherwise, if the atlas itself comes with a lookup table
    #     that didn't require conversion, write that;
    #   - In the group directory rather than the subject directory;
    #   - If it doesn't already exist.
    lut_export_file = shared.mrtrix_lut_file \
                      if shared.mrtrix_lut_file \
                      else shared.parc_lut_file
    if lut_export_file:
        lut_export_path = \
            os.path.join(output_analysis_level_path,
                         parc_string[1:] + '_lookup'
                         + os.path.splitext(lut_export_file)[1])
        try:
            shutil.copy(lut_export_file, lut_export_path)
        except OSError:
            pass

    # Copy / convert necessary files to output directory
    for scratch_file, output_item in output_items.items():
        if output_verbosity >= output_item.min_verbosity \
                and (multishell or not output_item.needs_multishell):
            full_output_path = os.path.join(output_subdir, output_item.path)
            if output_item.is_image:
                run.command('mrconvert '
                            + scratch_file
                            + ' '
                            + full_output_path
                            + (' ' + output_item.options
                               if output_item.options
                               else '')
                            + ' -clear_property comments',
                            force=os.path.exists(full_output_path))
            else:
                run.function(shutil.copyfile,
                             scratch_file,
                             full_output_path)

    with open(T1_json_path, 'w') as T1_json_file:
        json.dump(T1_json_data, T1_json_file)

    # Manually wipe and zero the scratch directory
    #   (since we might be processing more than one subject)
    os.chdir(cwd)
    if app.DO_CLEANUP:
        app.console('Deleting scratch directory ' + app.SCRATCH_DIR)
        # Can't use run.function() here;
        #   it'll try to write to the log file that resides
        #   in the scratch directory just deleted
        app.cleanup(app.SCRATCH_DIR)
    elif output_verbosity == 4:
        app.console('Copying scratch directory to output location')
        run.function(shutil.copytree,
                     app.SCRATCH_DIR,
                     os.path.join(output_subdir, 'scratch'))
    else:
        app.console('Contents of scratch directory kept; '
                    'location: ' + app.SCRATCH_DIR)
    app.SCRATCH_DIR = ''

# End of run_participant() function















GROUP_BRAINMASKS_DIR = 'brainmasks'
GROUP_BZEROS_DIR = 'bzeros'
GROUP_CONNECTOMES_DIR = 'connectomes'
GROUP_FA_DIR = 'fa'
GROUP_RESPONSES_DIR = 'responses'
GROUP_WMVOXELS_DIR = 'wmvoxels'
GROUP_WARPS_DIR = 'warps'

def run_group(bids_dir, output_verbosity, output_app_dir):

    participant_dir = os.path.join(output_app_dir, 'participant')
    group_dir = os.path.join(output_app_dir, 'group')

    # Participant-level analysis no longer generates FA and mean b=0 images
    # These really should not be that expensive to compute in series,
    #   and will keep the output directory cleaner

    # Check presence of all required input files before proceeding
    # Pre-calculate paths of all files since many will be used in
    #   more than one location
    class SessionPaths(object):
        def __init__(self, session):
            session_label = '_'.join(session)
            participant_root = os.path.join(group_dir, *session)
            group_root = os.path.join(group_dir, *session)
            # Get input DWI path here rather than in function
            in_dwi_image_list = glob.glob(os.path.join(participant_root,
                                                       'dwi',
                                                       '*_dwi.nii*'))
            if not in_dwi_image_list:
                raise MRtrixError('No DWI data found for session "'
                                  + session_label
                                  + '" output')
            if len(in_dwi_image_list) > 1:
                raise MRtrixError('More than one DWI mage found in session "'
                                  + session_label
                                  + '" output')
            self.in_dwi = in_dwi_image_list[0]
            if not '_desc-preproc_' in self.in_dwi:
                raise MRtrixError('DWI image in output directory for session "'
                                  + session_label
                                  + '" not flagged as pre-processed')
            in_dwi_prefix = self.in_dwi.split(os.extsep)[0]
            self.in_bvec = in_dwi_prefix + '.bvec'
            self.in_bval = in_dwi_prefix + '.bval'
            self.bvalues = [float(value) for value in \
                    run.command('mrinfo '
                                + self.in_dwi
                                + ' -fslgrad '
                                + self.in_bvec
                                + ' '
                                + self.in_bval
                                + ' -shell_bvalues').stdout.split()]

            self.in_rf = os.path.join(participant_root,
                                      'dwi',
                                      session_label
                                      + '_tissue-WM_response.txt')
            connectome_files = \
                glob.glob(
                    os.path.join(participant_root,
                                 'connectome',
                                 session_label
                                 + '_desc-*'
                                 + '_connectome.csv'))
            if not connectome_files:
                raise MRtrixError('No participant-level connectome file '
                                  'found for session "'
                                  + session_label
                                  + '"')
            if len(connectome_files) > 1:
                raise MRtrixError('Connectomes from multiple parcellations '
                                  'detected for session "'
                                  + session_label
                                  + '"; this is not yet supported')
            self.in_connectome = connectome_files[0]
            self.in_mu = os.path.join(participant_root,
                                      'tractogram',
                                      session_label + '_mu.txt')

            for entry in vars(self).values():
                if not os.path.exists(entry):
                    raise MRtrixError('Unable to find critical data '
                                      'for session "' + session_label + '"'
                                      '(expected location: ' + entry + ')')

            self.grad_import_option = ' -fslgrad ' \
                                      + self.in_bvec \
                                      + ' ' \
                                      + self.in_bval

            self.parcellation = \
                re.findall('(?<=_desc-)[a-zA-Z0-9]*',
                           os.path.basename(self.in_connectome))[0]

            # Permissible for this to not exist
            self.in_mask = os.path.join(participant_root,
                                        'dwi',
                                        session_label
                                        + '_desc-brain_mask.nii.gz')

            self.mu = matrix.load_vector(self.in_mu)[0]
            self.RF = matrix.load_matrix(self.in_rf)

            self.temp_mask = os.path.join(GROUP_BRAINMASKS_DIR,
                                          session_label + '.mif')
            self.temp_fa = os.path.join(GROUP_FA_DIR,
                                        session_label + '.mif')
            self.temp_bzero = os.path.join(GROUP_BZEROS_DIR,
                                           session_label + '.mif')
            self.temp_warp = os.path.join(GROUP_WARPS_DIR,
                                          session_label + '.mif')
            self.temp_voxels = os.path.join(GROUP_WMVOXELS_DIR,
                                            session_label + '.mif')
            self.temp_rf = os.path.join(GROUP_RESPONSES_DIR,
                                        session_label + '.txt')
            self.median_bzero = 0.0
            self.dwiintensitynorm_factor = 1.0
            self.RF_multiplier = 1.0
            self.volume_multiplier = 1.0
            for spacing in image.Header(self.in_dwi).spacing()[0:3]:
                self.volume_multiplier *= spacing
            self.global_multiplier = 1.0
            self.temp_connectome = os.path.join(GROUP_CONNECTOMES_DIR,
                                                session_label + '.csv')
            self.out_dir = group_root
            self.out_scale_intensity = \
                os.path.join(group_root,
                             'connectome',
                             session_label
                             + '_factor-intensity'
                             + '_multiplier.txt')
            self.out_scale_RF = \
                os.path.join(group_root,
                             'connectome',
                             session_label
                             + '_factor-response'
                             + '_multiplier.txt')
            self.out_scale_volume = \
                os.path.join(group_root,
                             'connectome',
                             session_label
                             + '_factor-volume'
                             + '_multiplier.txt')
            self.out_connectome = \
                os.path.join(group_root,
                             'connectome',
                             os.path.basename(self.in_connectome))

            self.session_label = session_label

    session_list = get_sessions(participant_dir)
    if not session_list:
        raise MRtrixError(
            'No processed session data found in output directory '
            'directory \'' + participant_dir + '\' for group analysis')
    if len(session_list) == 1:
        app.warn('Only one session present in participant directory; '
                 'some group-level analysis steps will be skipped')
    if os.path.exists(group_dir):
        app.warn('Output directory for group-level analysis '
                 'already exists; all contents will be erased when this '
                 'execution completes')


    bids_session_list = get_sessions(bids_dir)
    not_processed = [session for session in bids_session_list \
                     if session not in session_list]
    if not_processed:
        app.warn(str(len(not_processed)) + ' session'
                 + ('s' if len(not_processed) > 1 else '')
                 + ' present in BIDS directory '
                 + ('have' if len(not_processed) > 1 else 'has')
                 + ' not yet undergone participant-level processing: '
                 + ', '.join('_'.join(session) for session in not_processed))

    sessions = []
    for session in session_list:
        sessions.append(SessionPaths(session))



    # Connectome-based calculations can only be performed if the
    #   parcellation is consistent across all sessions
    parcellation = sessions[0].parcellation
    consistent_parcellation = \
        all(s.parcellation == parcellation for s in sessions)
    out_connectome_path = os.path.join(group_dir,
                                       'desc-'
                                       + parcellation
                                       + '_connectome.csv') \
                          if consistent_parcellation \
                          else None

    app.make_scratch_dir()
    app.goto_scratch_dir()

    # Before proceeding, compile session b-values and make sure that:
    #   - the number of shells is equivalent across sessions
    #   - the b-values don't vary too much within those shells across sessions
    if not all(len(session.bvalues) == len(sessions[0].bvalues)
               for session in sessions):
        raise MRtrixError('Not all sessions DWI data contain the same '
                          'number of b-value shells')
    all_bvalues = [[session.bvalues[index] for session in sessions]
                   for index in range(0, len(sessions[0].bvalues))]
    for shell in all_bvalues:
        shell_mean = sum(shell) / len(shell)
        if max([max(shell)-shell_mean, shell_mean-min(shell)]) > 50.0:
            raise MRtrixError('Excessive deviation of b-values: '
                              + 'mean across subjects b='
                              + str(shell_mean)
                              + '; '
                              + 'range '
                              + str(min(shell))
                              + '-'
                              + str(max(shell)))


    # First pass through subject data in group analysis:
    #   Generate mask and FA image directories to be used in
    #   population template generation.
    #   If output_verbosity >= 2 then a mask is already provided;
    #   if not, then one can be quickly calculated from the
    #   mean b=0 image, which must be provided
    progress = app.ProgressBar('Importing and preparing session data',
                               len(sessions))
    run.function(os.makedirs, GROUP_BRAINMASKS_DIR)
    run.function(os.makedirs, GROUP_BZEROS_DIR)
    run.function(os.makedirs, GROUP_FA_DIR)
    run.function(os.makedirs, GROUP_RESPONSES_DIR)
    for s in sessions:
        # We need three images for each session:
        # - Brain mask: Convert if present, otherwise generate from DWI
        # - Mean b=0 image (for scaling): Generate from DWI
        # - FA image (for registration): Generate from DWI
        if os.path.exists(s.in_mask):
            run.command('mrconvert '
                        + s.in_mask
                        + ' '
                        + s.temp_mask
                        + ' -datatype bit')
        else:
            run.command('dwi2mask '
                        + s.in_dwi
                        + s.grad_import_option
                        + ' '
                        + s.temp_mask)
        run.command('dwiextract '
                    + s.in_dwi
                    + s.grad_import_option
                    + ' -bzero - |'
                    + ' mrmath - mean -axis 3 '
                    + s.temp_bzero)
        run.command('dwi2tensor '
                    + s.in_dwi
                    + s.grad_import_option
                    + ' -mask '
                    + s.temp_mask
                    + ' - |'
                    + ' tensor2metric - -fa '
                    + s.temp_fa
                    + ' -mask '
                    + s.temp_mask)
        run.function(shutil.copy,
                     s.in_rf,
                     s.temp_rf)
        progress.increment()
    progress.done()

    # First group-level calculation:
    # Generate the population FA template
    if len(session_list) == 1:
        app.console('Duplicating single-subject FA image as '
                    'population template image')
        run.function(shutil.copyfile,
                     session_list[0].temp_fa,
                     'template.mif')
    else:
        app.console('Generating population template for '
                    'intensity normalisation WM mask derivation')
        run.command('population_template '
                    + GROUP_FA_DIR
                    + ' template.mif'
                    + ' -mask_dir '
                    + GROUP_BRAINMASKS_DIR
                    + ' -warp_dir '
                    + GROUP_WARPS_DIR
                    + ' -type rigid_affine_nonlinear'
                    + ' -rigid_scale 0.25,0.5,0.8,1.0'
                    + ' -affine_scale 0.7,0.8,1.0,1.0'
                    + ' -nl_scale 0.5,0.75,1.0,1.0,1.0'
                    + ' -nl_niter 5,5,5,5,5'
                    + ' -linear_no_pause')
    app.cleanup(GROUP_FA_DIR)
    app.cleanup(GROUP_BRAINMASKS_DIR)

    # Generate the group average response function
    if len(session_list) == 1:
        app.console('Duplicating single-subject WM response function as '
                    'group-average response function')
        run.function(shutil.copyfile,
                     session_list[0].temp_temp_rf,
                     'response.txt')
    else:
        app.console('Calculating group-average WM response function')
        run.command(['responsemean',
                     [s.temp_rf for s in sessions],
                     'response.txt'])
    app.cleanup(GROUP_RESPONSES_DIR)
    mean_RF = matrix.load_matrix('response.txt')
    mean_RF_lzero = [line[0] for line in mean_RF]

    # Second pass through subject data in group analysis:
    #     - Warp template FA image back to subject space &
    #       threshold to define a WM mask in subject space
    #     - Calculate the median subject b=0 value within this mask
    #     - Store this in a file, and contribute to calculation of the
    #       mean of these values across subjects
    #     - Contribute to the group average response function
    progress = app.ProgressBar('Generating group-average response function '
                               'and intensity normalisation factors',
                               len(sessions)+1)
    run.function(os.makedirs, GROUP_WMVOXELS_DIR)
    sum_median_bzero = 0.0
    for s in sessions:
        run.command('mrtransform template.mif '
                    '-warp_full ' + s.temp_warp + ' '
                    '-from 2 '
                    '-template ' + s.temp_bzero + ' '
                    '- | '
                    'mrthreshold - ' + s.temp_voxels + ' -abs 0.4')
        s.median_bzero = image.statistics(s.temp_bzero,
                                          mask=s.temp_voxels).median
        app.cleanup(s.temp_bzero)
        app.cleanup(s.temp_voxels)
        app.cleanup(s.temp_warp)
        sum_median_bzero += s.median_bzero
        progress.increment()
    app.cleanup(GROUP_BZEROS_DIR)
    app.cleanup(GROUP_WMVOXELS_DIR)
    app.cleanup(GROUP_WARPS_DIR)
    app.cleanup('template.mif')
    progress.done()

    # Second group-level calculation:
    # - Calculate the mean of median b=0 values
    mean_median_bzero = sum_median_bzero / len(sessions)

    # Third pass through session data in group analysis:
    # - Scale the connectome strengths:
    #   - Multiply by SIFT proportionality coefficient mu
    #   - Multiply by (mean median b=0) / (subject median b=0)
    #   - Multiply by (subject RF size) / (mean RF size)
    #     (needs to account for multi-shell data)
    #   - Multiply by voxel volume
    # - Write the result to file
    progress = app.ProgressBar('Applying normalisation scaling to '
                               'subject connectomes',
                               len(sessions))
    run.function(os.makedirs, GROUP_CONNECTOMES_DIR)
    for s in sessions:
        RF_lzero = [line[0] for line in s.RF]
        s.RF_multiplier = 1.0
        for (mean, subj) in zip(mean_RF_lzero, RF_lzero):
            s.RF_multiplier = s.RF_multiplier * subj / mean
        # Don't want to be scaling connectome independently for
        #   differences in RF l=0 terms across all shells;
        #   use the geometric mean of the per-shell scale factors
        s.RF_multiplier = math.pow(s.RF_multiplier, 1.0 / len(mean_RF_lzero))

        s.bzero_multiplier = mean_median_bzero / s.median_bzero

        s.global_multiplier = s.mu \
                              * s.bzero_multiplier \
                              * s.RF_multiplier \
                              * s.volume_multiplier

        connectome = matrix.load_matrix(s.in_connectome)
        temp_connectome = [[v*s.global_multiplier for v in line]
                           for line in connectome]
        matrix.save_matrix(s.temp_connectome, temp_connectome)
        progress.increment()
    progress.done()

    # Third group-level calculation: Generate the group mean connectome
    # Can only do this if the parcellation is identical across subjects;
    #     this needs to be explicitly checked
    if consistent_parcellation:
        progress = app.ProgressBar('Calculating group mean connectome',
                                   len(sessions)+1)
        # TODO Calculate geometric rather than arithmetic mean
        # Requires setting a minimum connectivity value per edge;
        #   this should be equivalent to 1 streamline prior to
        #   application of the multiplier
        mean_connectome = []
        for s in sessions:
            connectome = matrix.load_matrix(s.temp_connectome)
            if mean_connectome:
                mean_connectome = [[c1+c2 for c1, c2 in zip(r1, r2)]
                                   for r1, r2 in zip(mean_connectome,
                                                     connectome)]
            else:
                mean_connectome = connectome
            progress.increment()

        mean_connectome = [[v/len(sessions) for v in row]
                           for row in mean_connectome]
        progress.done()
    else:
        app.warn('Different parcellations across sessions; '
                 'cannot calculate a group mean connectome')

    # Write results of interest back to the output directory;
    #     both per-subject and group information
    progress = app.ProgressBar('Writing results to output directory',
                               len(sessions)+2)
    if os.path.exists(group_dir):
        run.function(shutil.rmtree, group_dir)
    run.function(os.makedirs, group_dir)
    for s in sessions:
        run.function(os.makedirs, s.out_dir)
        run.function(shutil.copyfile,
                     s.temp_connectome,
                     s.out_connectome)
        matrix.save_vector(s.out_scale_intensity,
                           [s.bzero_multiplier],
                           force=IS_CONTAINER)
        matrix.save_vector(s.out_scale_RF,
                           [s.RF_multiplier],
                           force=IS_CONTAINER)
        matrix.save_vector(s.out_scale_volume,
                           [s.volume_multiplier],
                           force=IS_CONTAINER)
        progress.increment()
    app.cleanup(GROUP_CONNECTOMES_DIR)

    matrix.save_matrix(os.path.join(group_dir, 'tissue-WM_response.txt'),
                       mean_RF,
                       force=IS_CONTAINER)
    progress.increment()
    if consistent_parcellation:
        matrix.save_matrix(out_connectome_path,
                           mean_connectome,
                           force=IS_CONTAINER)
    progress.done()

    # For group-level analysis, function is only executed once, so
    #   no need to bypass the default scratch cleanup
    # Only exception is if we want to capture the whole scratch directory
    #   in the output path
    if output_verbosity == 4:
        app.console('Copying scratch directory to output location')
        run.function(shutil.copytree,
                     app.SCRATCH_DIR,
                     os.path.join(group_dir, 'scratch'))

# End of run_group() function










# Examine the contents of a directory (whether the raw BIDS dataset or a
#   derivatives directory), and return a list of sessions.
# This list may optionally be filtered based on the use of batch processing
#   command-line options; e.g. resticting the participant or session IDs.
def get_sessions(root_dir, **kwargs):

    participant_labels = kwargs.pop('participant_label', None)
    session_labels = kwargs.pop('session_label', None)
    if kwargs:
        raise TypeError(
            'Unsupported keyword arguments passed to get_session(): '
            + str(kwargs))

    # Perform a recursive search through the BIDS dataset directory,
    #   looking for anything that resembles a BIDS session
    # For any sub-directory that itself contains directories "anat/" and "dwi/",
    #   store the list of sub-directories required to navigate to that point
    # This becomes the list of feasible processing targets for any level
    #   of analysis
    # From there:
    #   - "--participant_label" can be used to remove entries from the list
    #   - Other options can be added to restrict processing targets;
    #     e.g. "--session_label" to remove based on "ses-*/"
    all_sessions = []
    for dir_name, subdir_list, _ in os.walk(root_dir):
        if 'anat' in subdir_list and 'dwi' in subdir_list:
            all_sessions.append(os.path.relpath(dir_name, start=root_dir))
            del subdir_list
    all_sessions = sorted(all_sessions)
    app.debug(str(all_sessions))

    result = []

    # Need to alert user if they have nominated a particular participant /
    #   session label, and no such data were found in the input dataset
    sub_found = {label: False for label in participant_labels} \
                if participant_labels \
                else {}
    ses_found = {label: False for label in session_labels} \
                if session_labels \
                else {}

    # Define worker function for applying the --participant_label and
    #   --session_label restrictions
    def find_and_flag(ses, prefix, labels, found):
        for dirname in ses:
            if dirname.startswith(prefix):
                present = False
                for label in labels:
                    if label == dirname[len(prefix):]:
                        found[label] = True
                        present = True
                        break
                if not present:
                    return False
        return True

    for session in all_sessions:
        session = os.path.normpath(session).split(os.sep)
        process = True
        if participant_labels:
            if not find_and_flag(session,
                                 'sub-',
                                 participant_labels,
                                 sub_found):
                process = False
        if session_labels:
            if not find_and_flag(session,
                                 'ses-',
                                 session_labels,
                                 ses_found):
                process = False
        if process:
            result.append(session)

    if not result:
        raise MRtrixError('No sessions were selected for processing')

    app.console(str(len(all_sessions)) + ' total sessions found '
                + 'in directory \'' + root_dir + '\'; '
                + ('all'
                   if len(result) == len(all_sessions)
                   else str(len(result)))
                + ' will be processed')
    sub_not_found = [key for key, value in sub_found.items() if not value]
    if sub_not_found:
        app.warn(str(len(sub_not_found)) + ' nominated participant label'
                 + ('s were' if len(sub_not_found) > 1 else ' was')
                 + ' not found in input dataset: '
                 + ', '.join(sub_not_found))
    ses_not_found = [key for key, value in ses_found.items() if not value]
    if ses_not_found:
        app.warn(str(len(ses_not_found)) + ' nominated session label'
                 + ('s were' if len(ses_not_found) > 1 else ' was')
                 + ' not found in input dataset: '
                 + ', '.join(ses_not_found))

    app.debug(str(result))
    return result















ANALYSIS_CHOICES = ['preproc', 'participant', 'group']

PARCELLATION_CHOICES = ['aal',
                        'aal2',
                        'brainnetome246fs',
                        'brainnetome246mni',
                        'craddock200',
                        'craddock400',
                        'desikan',
                        'destrieux',
                        'hcpmmp1',
                        'none',
                        'perry512',
                        'yeo7fs',
                        'yeo7mni',
                        'yeo17fs',
                        'yeo17mni']

REGISTRATION_CHOICES = ['ants', 'fsl']



def usage(cmdline): #pylint: disable=unused-variable
    cmdline.set_author('Robert E. Smith (robert.smith@florey.edu.au)')
    cmdline.set_synopsis(
        'Generate structural connectomes based on diffusion-weighted '
        'and T1-weighted image data using state-of-the-art reconstruction '
        'tools, particularly those provided in MRtrix3')

    cmdline.set_copyright(
        '''Copyright (c) 2016-2020 The Florey Institute of Neuroscience
and Mental Health.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Covered Software is provided under this License on an "as is"
basis, without warranty of any kind, either expressed, implied, or
statutory, including, without limitation, warranties that the
Covered Software is free of defects, merchantable, fit for a
particular purpose or non-infringing.
See the Mozilla Public License v. 2.0 for more details.''')

    # If running within a container, erase existing standard options, and
    #   fill with only desired options
    if IS_CONTAINER:
        # pylint: disable=protected-access
        for option in reversed(cmdline._actions):
            cmdline._handle_conflict_resolve(
                None, [(option.option_strings[0], option)])
        # cmdline._action_groups[2] is "Standard options"
        #   that was created earlier by the API
        cmdline._action_groups[2].add_argument(
            '-d', '--debug',
            dest='debug',
            action='store_true',
            help='In the event of encountering an issue with the script, '
                 're-run with this flag set to provide more useful '
                 'information to the developer')
        cmdline._action_groups[2].add_argument(
            '-h', '--help',
            dest='help',
            action='store_true',
            help='Display help information for the script')
        cmdline._action_groups[2].add_argument(
            '-n', '--n_cpus',
            type=int,
            metavar='number',
            dest='nthreads',
            help='Use this number of threads in MRtrix3 '
                 'multi-threaded applications '
                 '(0 disables multi-threading)')
        cmdline._action_groups[2].add_argument(
            '-scratch', '--scratch',
            dest='scratch',
            help='Set location for script scratch directory')
        cmdline._action_groups[2].add_argument(
            '-skip', '--skip-bids-validator',
            dest='skipbidsvalidator',
            action='store_true',
            help='Skip BIDS validation')
        cmdline._action_groups[2].add_argument(
            '-v', '--version',
            action='version',
            version=__version__)
    else:
        cmdline._action_groups[2].add_argument( # pylint: disable=protected-access
            OPTION_PREFIX + 'skip-bids-validator',
            dest='skipbidsvalidator',
            action='store_true',
            help='Skip BIDS validation')

    cmdline.add_description(
        'While preproc-level analysis only requires data within the '
        'BIDS directory, participant-level analysis requires that the '
        'output directory be pre-populated with the results from '
        'preproc-level processing; similarly, group-level analysis '
        'requires that the output directory be pre-populated with the '
        'results from participant-level analysis.')
    cmdline.add_description(
        'The operations performed by each of the three levels of analysis '
        'are as follows:')
    cmdline.add_description(
        '"preproc": '
        'DWI: Denoising; Gibbs ringing removal; motion, '
        'eddy current and EPI distortion correction and outlier detection & '
        'replacement; brain masking, bias field correction and intensity '
        'normalisation; rigid-body registration & transformation to '
        'T1-weighted image. '
        'T1-weighted image: bias field correction; brain masking.')
    cmdline.add_description(
        '"participant": '
        'DWI: Response function estimation; FOD estimation. '
        'T1-weighted image (if ' + OPTION_PREFIX + 'parcellation '
        'is not none): '
        'Tissue segmentation; grey matter parcellation. '
        'Combined (if ' + OPTION_PREFIX + 'parcellation is not none, '
        'or ' + OPTION_PREFIX + 'streamlines is provided): '
        'Whole-brain streamlines tractography; SIFT2; '
        'connectome construction.')
    cmdline.add_description(
        '"group": '
        'Generation of FA-based population template; '
        'warping of template-based white matter mask to subject spaces; '
        'calculation of group mean white matter response function; '
        'scaling of connectomes based on white matter b=0 intensity, '
        'response function used during participant-level analysis, and '
        'SIFT model proportioinality coefficient; '
        'generation of group mean connectome.')
    cmdline.add_description(
        'The label(s) provided to the '
        + OPTION_PREFIX + 'participant_label and '
        + OPTION_PREFIX + 'session_label options '
        + 'correspond(s) to sub-<participant_label> and '
        + 'ses-<session_label> from the BIDS spec (so they do _not_ '
        + 'include "sub-" or "ses-"). Multiple participants / sessions '
        + 'can be specified with a space-separated list.')
    cmdline.add_description(
        'For both preproc-level and participant-level analyses, if no '
        'specific participants or sessions are nominated by the user '
        '(or the user explicitly specifies multiple participants / '
        'sessions), the script will process each of these in series. '
        'It is additionally possible for the user to invoke multiple '
        'instances of this script in order to process multiple subjects '
        'at once in parallel, ensuring that no single participant / '
        'session is being processed in parallel, and that preproc-level '
        'output data are written fully before commencing participant-level '
        'analysis.')
    cmdline.add_description(
        'The ' + OPTION_PREFIX + 'output_verbosity option principally '
        'affects the participant-level analysis, modulating how many '
        'derivative files are written to the output directory. Permitted '
        'values are from 1 to 4: 1 writes only those files requisite for '
        'group-level analysis; 2 additionally writes files typically '
        'useful for post-hoc analysis (the default); 3 additionally '
        'generates files for enhanced connectome visualisation and copies '
        'the entire whole-brain tractogram; 4 additionally generates a '
        'full copy of the script scratch directory (with all intermediate '
        'files retained) to the output directory (and this applies to '
        'all analysis levels)')
    if not IS_CONTAINER:
        cmdline.add_description(
            'If running participant-level analysis using the script as a '
            'standalone tool rather than inside the provided container, '
            'data pertaining to atlas parcellations can no longer be '
            'guaranteed to be stored at a specific location on the '
            'filesystem. In this case, the user will most likely need to '
            'manually specify the location where the corresponding '
            'parcellation is stored using the -atlas_path option.')

    cmdline.add_argument(
        'bids_dir',
        help='The directory with the input dataset formatted '
             'according to the BIDS standard.')
    cmdline.add_argument(
        'output_dir',
        help='The directory where the output files should be stored.')
    cmdline.add_argument(
        'analysis_level',
        help='Level of analysis that will be performed; '
             'options are: ' + ', '.join(ANALYSIS_CHOICES) + '.',
        choices=ANALYSIS_CHOICES)

    cmdline.add_argument(
        OPTION_PREFIX + 'output_verbosity',
        type=int,
        default=2,
        help='The verbosity of script output (number from 1 to 4).')

    batch_options = cmdline.add_argument_group(
        'Options specific to the batch processing of participant data')
    batch_options.add_argument(
        OPTION_PREFIX + 'participant_label',
        nargs='+',
        help='The label(s) of the participant(s) that should be analyzed.')
    batch_options.add_argument(
        OPTION_PREFIX + 'session_label',
        nargs='+',
        help='The session(s) within each participant that should be analyzed.')

    preproc_participant_options = \
        cmdline.add_argument_group(
            'Options that are relevant to both preproc-level and '
            'participant-level analyses')
    preproc_participant_options.add_argument(
        OPTION_PREFIX + 't1w_preproc',
        metavar='path',
        help='Provide a path by which pre-processed T1-weighted image data '
             'may be found for the processed participant(s) / session(s)')

    participant_options = \
        cmdline.add_argument_group(
            'Options that are relevant to participant-level analysis')
    if not IS_CONTAINER:
        participant_options.add_argument(
            OPTION_PREFIX + 'atlas_path',
            metavar='path',
            help='The filesystem path in which to search for atlas '
                 'parcellation files.')
    participant_options.add_argument(
        OPTION_PREFIX + 'parcellation',
        help='The choice of connectome parcellation scheme '
             '(compulsory for participant-level analysis); '
             'options are: ' + ', '.join(PARCELLATION_CHOICES) + '.',
        choices=PARCELLATION_CHOICES)
    participant_options.add_argument(
        OPTION_PREFIX + 'streamlines',
        type=int,
        default=0,
        help='The number of streamlines to generate for each subject '
             '(will be determined heuristically if not explicitly set).')
    participant_options.add_argument(
        OPTION_PREFIX + 'template_reg',
        metavar='software',
        help='The choice of registration software for mapping subject to '
             'template space; '
             'options are: ' + ', '.join(REGISTRATION_CHOICES) + '.',
        choices=REGISTRATION_CHOICES)

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
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'template_reg fsl',
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
    if not IS_CONTAINER:
        cmdline.add_citation(
            'Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H.; '
            'Filippini, N. & Bastiani, M. '
            'Towards a comprehensive framework for movement and distortion '
            'correction of diffusion MR images: Within volume movement. '
            'NeuroImage, 2017, 152, 450-466',
            condition='If performing preproc-level analysis',
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
        + 'using ' + OPTION_PREFIX + 'parcellation [ aal, aal2, '
        + 'brainnetome246mni, craddock200, craddock400, perry512, '
        + 'yeo7mni or yeo17mni ], '
        + 'and not using ' + OPTION_PREFIX + 'template_reg fsl',
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
        + 'using ' + OPTION_PREFIX + 'parcellation '
        + '[ craddock200 or craddock400 ]',
        is_external=True)
    cmdline.add_citation(
        'Dale, A. M.; Fischl, B. & Sereno, M. I. '
        'Cortical Surface-Based Analysis: '
        'I. Segmentation and Surface Reconstruction. '
        'NeuroImage, 1999, 9, 179-194',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation '
        + '[ brainnetome246fs, desikan, destrieux, hcpmmp1, '
        + 'yeo7fs or yeo17fs ]',
        is_external=True)
    cmdline.add_citation(
        'Desikan, R. S.; Segonne, F.; Fischl, B.; Quinn, B. T.; '
        'Dickerson, B. C.; Blacker, D.; Buckner, R. L.; Dale, A. M.; '
        'Maguire, R. P.; Hyman, B. T.; Albert, M. S. & Killiany, R. J. '
        'An automated labeling system for subdividing the human cerebral '
        'cortex on MRI scans into gyral based regions of interest. '
        'NeuroImage, 2006, 31, 968-980',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation desikan',
        is_external=True)
    cmdline.add_citation(
        'Destrieux, C.; Fischl, B.; Dale, A. & Halgren, E. '
        'Automatic parcellation of human cortical gyri and sulci using '
        'standard anatomical nomenclature. '
        'NeuroImage, 2010, 53, 1-15',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation destrieux',
        is_external=True)
    cmdline.add_citation(
        'Fan, L.; Li, H.; Zhuo, J.; Zhang, Y.; Wang, J.; Chen, L.; Yang, Z.; '
        'Chu, C.; Xie, S.; Laird, A.R.; Fox, P.T.; Eickhoff, S.B.; Yu, C.; '
        'Jiang, T. '
        'The Human Brainnetome Atlas: '
        'A New Brain Atlas Based on Connectional Architecture. '
        'Cerebral Cortex, 2016, 26 (8), 3508-3526',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation '
        + '[ brainnetome246fs brainnetome246mni ]',
        is_external=True)
    cmdline.add_citation(
        'Glasser, M. F.; Coalson, T. S.; Robinson, E. C.; Hacker, C. D.; '
        'Harwell, J.; Yacoub, E.; Ugurbil, K.; Andersson, J.; Beckmann, C. F.; '
        'Jenkinson, M.; Smith, S. M. & Van Essen, D. C. '
        'A multi-modal parcellation of human cerebral cortex. '
        'Nature, 2016, 536, 171-178',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation hcpmmp1',
        is_external=True)
    cmdline.add_citation(
        'Iglesias, J. E.; Liu, C. Y.; Thompson, P. M. & Tu, Z. '
        'Robust Brain Extraction Across Datasets and '
        'Comparison With Publicly Available Methods. '
        'IEEE Transactions on Medical Imaging, 2011, 30, 1617-1634',
        condition='If performing either preproc-level or participant-level '
        'analysis and ROBEX is used for T1 brain extraction',
        is_external=True)
    cmdline.add_citation(
        'Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. '
        'Multi-tissue constrained spherical deconvolution for improved '
        'analysis of multi-shell diffusion MRI data. '
        'NeuroImage, 2014, 103, 411-426',
        condition='If performing either preproc-level or participant-level '
        'analysis',
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
        + 'not using ' + OPTION_PREFIX + 'parcellation '
        + '[ brainnetome246fs or brainnetome246mni ]',
        is_external=True)
    cmdline.add_citation(
        'Perry, A.; Wen, W.; Kochan, N. A.; Thalamuthu, A.; Sachdev, P. S.; '
        'Breakspear, M. '
        'The independent influences of age and education on functional brain '
        'networks and cognition in healthy older adults. '
        'Human Brain Mapping, 2017, 38(10), 5094-5114',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation perry512',
        is_external=True)
    if not IS_CONTAINER:
        cmdline.add_citation(
            'Smith, S. M. '
            'Fast robust automated brain extraction. '
            'Human Brain Mapping, 2002, 17, 143-155',
            condition='If performing either preproc-level or participant-level '
            'analysis and FSL is used for T1 brain extraction',
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
        condition='If performing participant-level analysis and '
        + OPTION_PREFIX + 'streamlines 0 is not set',
        is_external=False)
    cmdline.add_citation(
        'Tournier, J.-D.; Smith, R. E.; Raffelt, D. A.; Tabbara, '
        'R.; Dhollander, T.; Pietsch, M; Christiaens, D.; Jeurissen, '
        'B.; Y, C.-H.; Connelly, A.;'
        'MRtrix3: A fast, flexible and open software framework for '
        'medical image processing and visualisation. '
        'NeuroImage, 2019, 202, 116137',
        is_external=False)
    cmdline.add_citation(
        'Tustison, N.; Avants, B.; Cook, P.; Zheng, Y.; Egan, A.; '
        'Yushkevich, P. & Gee, J. '
        'N4ITK: Improved N3 Bias Correction. '
        'IEEE Transactions on Medical Imaging, 2010, 29, 1310-1320',
        condition='If performing either preproc-level or participant-level '
        'analysis, and N4 is used for either DWI or T1 bias field correction',
        is_external=True)
    cmdline.add_citation(
        'Tzourio-Mazoyer, N.; Landeau, B.; Papathanassiou, D.; Crivello, F.; '
        'Etard, O.; Delcroix, N.; Mazoyer, B. & Joliot, M. '
        'Automated Anatomical Labeling of activations in SPM using a '
        'Macroscopic Anatomical Parcellation of the MNI MRI '
        'single-subject brain. '
        'NeuroImage, 15(1), 273-289',
        condition='If performing participant-level analysis and '
        'using ' + OPTION_PREFIX + 'parcellation [ aal or aal2 ]',
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
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation '
        + '[ yeo7fs, yeo7mni, yeo17fs or yeo17mni ]',
        is_external=False)
    cmdline.add_citation(
        'Zalesky, A.; Fornito, A.; Harding, I. H.; Cocchi, L.; Yucel, M.; '
        'Pantelis, C. & Bullmore, E. T. '
        'Whole-brain anatomical networks: Does the choice of nodes matter? '
        'NeuroImage, 2010, 50, 970-983',
        condition='If performing participant-level analysis and '
        + 'using ' + OPTION_PREFIX + 'parcellation perry512',
        is_external=True)
    cmdline.add_citation(
        'Zhang, Y.; Brady, M. & Smith, S. '
        'Segmentation of brain MR images through a hidden Markov random '
        'field model and the expectation-maximization algorithm. '
        'IEEE Transactions on Medical Imaging, 2001, 20, 45-57',
        condition='If performing participant-level analysis',
        is_external=True)






def execute(): #pylint: disable=unused-variable

    # If running within a container, and the --debug option has been
    #     provided, modify the interlly-stored MRtrix3 configuration
    #     contents, so that any temporary directories will be constructed
    #     within the mounted output directory, and therefore temporary
    #     directory contents will not be lost upon container instance
    #     destruction if the script fails at any point.
    if IS_CONTAINER and app.ARGS.debug:
        app.DO_CLEANUP = False
        if 'ScriptScratchDir' not in CONFIG and not app.ARGS.scratch:
            CONFIG['ScriptScratchDir'] = os.path.abspath(app.ARGS.output_dir)

    if utils.is_windows():
        raise MRtrixError(
            'Script cannot be run on Windows due to FSL dependency')

    if app.ARGS.skipbidsvalidator:
        app.console('Skipping BIDS validation based on user request')
    elif find_executable('bids-validator'):
        run.command('bids-validator ' + app.ARGS.bids_dir)
    else:
        app.warn('BIDS validator script not installed; '
                 'proceeding without validation of input data')

    if app.ARGS.output_verbosity < 1 or app.ARGS.output_verbosity > 4:
        raise MRtrixError('Valid values for '
                          + OPTION_PREFIX + 'output_verbosity '
                          + 'option are from 1 to 4')

    # At output verbosity level 4 we retain all data and move the
    #   scratch directory to the output
    if app.ARGS.output_verbosity == 4:
        app.DO_CLEANUP = False

    # Locate root directory of BIDS App output based on user input
    #   (i.e. location of directory "mrtrix3_connectome", which itself
    #   contains sub-directories for each analysis level)
    # Note: abspath() removes any trailing path separator
    output_abspath = os.path.abspath(app.ARGS.output_dir)
    output_basename = os.path.basename(output_abspath)
    output_dirname = os.path.dirname(output_abspath)
    # Basename of output path is the analysis level being requested
    if output_basename.lower() == app.ARGS.analysis_level:
        output_app_path = os.path.dirname(output_dirname)
        if os.path.basename(os.path.dirname(output_app_path)) \
                           .lower() != 'mrtrix3_connectome':
            raise MRtrixError('Output directory structure malformed')
    # Basename of output path is "mrtrix3_connectome"
    elif output_basename.lower() == 'mrtrix3_connectome':
        output_app_path = output_dirname
    # Basename of output path is unknown
    else:
        output_app_path = os.path.join(output_abspath, 'MRtrix3_connectome')
    if not os.path.isdir(output_app_path):
        run.function(os.makedirs, output_app_path)


    if app.ARGS.analysis_level in ['preproc', 'participant']:
        sessions_to_analyze = get_sessions(
            app.ARGS.bids_dir,
            participant_label=app.ARGS.participant_label,
            session_label=app.ARGS.session_label)
        t1w_preproc_path = os.path.abspath(app.ARGS.t1w_preproc) \
                           if app.ARGS.t1w_preproc \
                           else None

    if app.ARGS.analysis_level == 'preproc':

        preproc_shared = PreprocShared()

        for session_to_process in sessions_to_analyze:
            app.console('Commencing execution for session: \''
                        + '_'.join(session_to_process) + '\'')
            run_preproc(os.path.abspath(app.ARGS.bids_dir),
                        session_to_process,
                        preproc_shared,
                        t1w_preproc_path,
                        app.ARGS.output_verbosity,
                        output_app_path)

    if app.ARGS.analysis_level == 'participant':

        participant_shared = \
            ParticipantShared(getattr(app.ARGS, 'atlas_path', None),
                              app.ARGS.parcellation,
                              app.ARGS.streamlines,
                              app.ARGS.template_reg)

        for session_to_process in sessions_to_analyze:
            app.console('Commencing execution for session: \''
                        + '_'.join(session_to_process) + '\'')
            run_participant(os.path.abspath(app.ARGS.bids_dir),
                            session_to_process,
                            participant_shared,
                            t1w_preproc_path,
                            app.ARGS.output_verbosity,
                            output_app_path)

    elif app.ARGS.analysis_level == 'group':

        if app.ARGS.participant_label:
            raise MRtrixError('Cannot use '
                              + OPTION_PREFIX + 'participant_label option '
                              + 'when performing group-level analysis')
        if app.ARGS.session_label:
            raise MRtrixError('Cannot use '
                              + OPTION_PREFIX + 'session_label option '
                              + 'when performing group-level analysis')

        run_group(os.path.abspath(app.ARGS.bids_dir),
                  app.ARGS.output_verbosity,
                  output_app_path)



mrtrix3.execute()
