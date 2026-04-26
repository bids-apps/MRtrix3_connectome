import os
import pathlib
import shutil
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import fsl
from .. import OPTION_PREFIX
from ..anat.shared import Shared as T1wShared

class Shared(object): #pylint: disable=useless-object-inheritance
    def __init__(self, atlas_path, parcellation,
                 streamlines, template_reg):

        if not parcellation:
            raise MRtrixError(
                'For participant-level analysis, '
                'desired parcellation must be provided using the '
                f'{OPTION_PREFIX}parcellation option')
        self.parcellation = parcellation

        self.streamlines = streamlines

        fsl_path = os.environ.get('FSLDIR', None)
        if fsl_path:
            fsl_path = pathlib.Path(fsl_path)
            if not fsl_path.is_dir():
                raise MRtrixError(
                    'Environment variable FSLDIR does not point'
                    ' to an existing directory')
        else:
            raise MRtrixError(
                'Environment variable FSLDIR is not set; '
                'please run appropriate FSL configuration script')

        self.dwi2mask_algo = 'synthstrip'
        if not shutil.which('mri_synthstrip'):
            self.dwi2mask_algo = 'legacy'
            app.warn('FreeSurfer command mri_synthstrip not present;'
                     ' legacy dwi2mask algorithm will be used if necessary'
                     ' (ie. if pre-processed data do not provide a brain mask)')

        # No GDC data provided;
        #   won't be doing any T1w pre-processing
        self.t1w_shared = T1wShared({})

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
                         f'{OPTION_PREFIX}template_reg option ignored')
                self.template_registration_software = ''
        else:
            self.template_registration_software = 'ants' if self.do_mni else ''
        if self.template_registration_software == 'ants':
            if not shutil.which('antsRegistration') \
                    or not shutil.which('antsApplyTransforms'):
                raise MRtrixError(
                    'Commands \'antsRegistration\' and \'antsApplyTransforms\' '
                    'must be present in PATH to use '
                    'ANTs software for template registration')
        elif self.template_registration_software == 'fsl':
            self.flirt_cmd = fsl.exe_name('flirt')
            self.fnirt_cmd = fsl.exe_name('fnirt')
            self.invwarp_cmd = fsl.exe_name('invwarp')
            self.applywarp_cmd = fsl.exe_name('applywarp')
            self.fnirt_config_basename = 'T1_2_MNI152_2mm.cnf'
            self.fnirt_config_path = fsl_path / \
                                     'etc' / \
                                     'flirtsch' / \
                                     self.fnirt_config_basename
            if not self.fnirt_config_path.is_file():
                raise MRtrixError(
                    'Unable to find configuration file for FNI FNIRT '
                    f'(expected location: {self.fnirt_config_path})')

        self.template_image_path = None
        self.template_mask_path = None
        self.parc_image_path = None
        self.parc_lut_file = None
        self.mrtrix_lut_file = None

        # TODO This might need to change if MRtrix3 is installed
        #   rather than accessing the build directory
        mrtrix_lut_dir = pathlib.Path(app.__file__).parents[2] / \
                         'share' / \
                         'mrtrix3' / \
                         'labelconvert'

        if self.do_freesurfer:
            self.freesurfer_home = os.environ.get('FREESURFER_HOME', None)
            if not self.freesurfer_home:
                raise MRtrixError(
                    'Environment variable FREESURFER_HOME not set; '
                    'please verify FreeSurfer installation')
            self.freesurfer_home = pathlib.Path(self.freesurfer_home)
            if not self.freesurfer_home.is_dir():
                raise MRtrixError(
                    'Environment variable FREESURFER_HOME'
                    ' does not point to an existing directory'
                    f' ({self.freesurfer_home})'
                )
            if not shutil.which('recon-all'):
                raise MRtrixError(
                    'Could not find FreeSurfer script "recon-all"; '
                    'please verify FreeSurfer installation')
            self.freesurfer_subjects_dir = pathlib.Path(os.environ['SUBJECTS_DIR']) \
                                           if 'SUBJECTS_DIR' in os.environ \
                                           else (self.freesurfer_home / 'subjects')
            if not self.freesurfer_subjects_dir.is_dir():
                raise MRtrixError(
                    'Could not find FreeSurfer subjects directory '
                    f'(expected location: {self.freesurfer_subjects_dir})')
            for subdir in ['fsaverage',
                           'fsaverage5',
                           'lh.EC_average',
                           'rh.EC_average']:
                if not (self.freesurfer_subjects_dir / subdir).is_dir():
                    raise MRtrixError(
                        'Could not find requisite FreeSurfer subject '
                        f'directory "{subdir}" '
                        f'(expected location: {self.freesurfer_subjects_dir / subdir})')
            self.reconall_path = shutil.which('recon-all')
            if not self.reconall_path:
                raise MRtrixError(
                    'Could not find FreeSurfer script "recon-all"; '
                    'please verify FreeSurfer installation')
            if parcellation in ['hcpmmp1', 'yeo7fs', 'yeo17fs']:
                if parcellation == 'hcpmmp1':

                    def hcpmmp_annot_path(hemi):
                        return self.freesurfer_subjects_dir / \
                               'fsaverage' / \
                               'label' / \
                               f'{hemi}h.HCPMMP1.annot'

                    self.hcpmmp1_annot_paths = [hcpmmp_annot_path(hemi)
                                                for hemi in ['l', 'r']]
                    if not all(path.is_file() for path in self.hcpmmp1_annot_paths):
                        raise MRtrixError(
                            'Could not find necessary annotation labels '
                            'for applying HCPMMP1 parcellation '
                            f'(expected location: {hcpmmp_annot_path("?")})')
                else: # yeo7fs, yeo17fs

                    def yeo_annot_path(hemi):
                        return self.freesurfer_subjects_dir / \
                            'fsaverage5' / \
                            'label' / \
                            f'{hemi}h.Yeo2011_' \
                            f'{"7" if parcellation == "yeo7fs" else "17"}' \
                            'Networks_N1000.split_components.annot'

                    self.yeo_annot_paths = [yeo_annot_path(hemi) \
                                            for hemi in ['l', 'r']]
                    if not all(path.is_file() for path in self.yeo_annot_paths):
                        raise MRtrixError(
                            'Could not find necessary annotation labels '
                            'for applying Yeo2011 parcellation '
                            f'(expected location: {yeo_annot_path("?")})')
                for cmd in ['mri_surf2surf', 'mri_aparc2aseg']:
                    if not shutil.which(cmd):
                        raise MRtrixError(
                            f'Could not find FreeSurfer command {cmd} '
                            '(necessary for applying HCPMMP1 parcellation); '
                            'please verify FreeSurfer installation')
            elif parcellation == 'brainnetome246fs':

                def brainnetome_gcs_path(hemi):
                    return self.freesurfer_home / 'average' / f'{hemi}h.BN_Atlas.gcs'

                self.brainnetome_cortex_gcs_paths = [
                    brainnetome_gcs_path(hemi)
                    for hemi in ['l', 'r']]
                if not all(path.is_file() for path in self.brainnetome_cortex_gcs_paths):
                    raise MRtrixError(
                        'Could not find necessary GCS files for applying '
                        'Brainnetome cortical parcellation via FreeSurfer '
                        f'(expected location: {brainnetome_gcs_path("?")})')
                self.brainnetome_sgm_gca_path = \
                    self.freesurfer_home / 'average' / 'BN_Atlas_subcortex.gca'
                if not self.brainnetome_sgm_gca_path.is_file():
                    raise MRtrixError(
                        'Could not find necessary GCA file for applying '
                        'Brainnetome sub-cortical parcellation '
                        'via FreeSurfer '
                        f'(expected location: {self.brainnetome_sgm_gca_path})')
                for cmd in ['mri_label2vol',
                            'mri_ca_label',
                            'mris_ca_label']:
                    if not shutil.which(cmd):
                        raise MRtrixError(
                            f'Could not find FreeSurfer command {cmd} '
                            '(necessary for applying Brainnetome parcellation); '
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
                with open(self.reconall_path, 'r', encoding='utf-8') as f:
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
            app.debug(self.reconall_multithread_options)

            if parcellation == 'brainnetome246fs':
                self.parc_lut_file = self.freesurfer_home / 'BN_Atlas_246_LUT.txt'
                self.mrtrix_lut_file = None
            elif parcellation == 'desikan':
                self.parc_lut_file = self.freesurfer_home / 'FreeSurferColorLUT.txt'
                self.mrtrix_lut_file = mrtrix_lut_dir / 'fs_default.txt'
            elif parcellation == 'destrieux':
                self.parc_lut_file = self.freesurfer_home / 'FreeSurferColorLUT.txt'
                self.mrtrix_lut_file = mrtrix_lut_dir / 'fs_a2009s.txt'
            elif parcellation == 'hcpmmp1':
                self.parc_lut_file = mrtrix_lut_dir / 'hcpmmp1_original.txt'
                self.mrtrix_lut_file = mrtrix_lut_dir / 'hcpmmp1_ordered.txt'
            elif parcellation in ['yeo7fs', 'yeo17fs']:
                self.parc_lut_file = \
                    self.freesurfer_home / \
                    'Yeo2011_' \
                    f'{"7" if parcellation == "yeo7fs" else "17"}' \
                    'networks_Split_Components_LUT.txt'
                self.mrtrix_lut_file = \
                    mrtrix_lut_dir / \
                    'Yeo2011_' \
                    f'{"7" if parcellation == "yeo7fs" else "17"}' \
                    'N_split.txt'
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
                os.remove('test_softlink')
                app.debug('Using softlinks to FreeSurfer template directories')
            except OSError:
                app.debug('Unable to create softlinks; '
                          'will perform deep copies of FreeSurfer '
                          'template directories')
                self.freesurfer_template_link_function = shutil.copytree

        elif self.do_mni:
            self.template_image_path = \
                fsl_path / 'data' / 'standard' / 'MNI152_T1_2mm.nii.gz'
            self.template_mask_path = \
                fsl_path / 'data' / 'standard' / 'MNI152_T1_2mm_brain_mask.nii.gz'
            if parcellation == 'aal':
                self.parc_image_path = pathlib.Path(pathlib.Path().root,
                                                    'opt',
                                                    'aal',
                                                    'ROI_MNI_V4.nii')
                self.parc_lut_file = self.parc_image_path.with_suffix('.txt')
                self.mrtrix_lut_file = mrtrix_lut_dir / 'aal.txt'
            elif parcellation == 'aal2':
                self.parc_image_path = pathlib.Path(pathlib.Path().root,
                                                    'opt',
                                                    'aal',
                                                    'ROI_MNI_V5.nii')
                self.parc_lut_file = self.parc_image_path.with_suffix('.txt')
                self.mrtrix_lut_file = mrtrix_lut_dir / 'aal2.txt'
            elif parcellation == 'brainnetome246mni':
                self.parc_image_path = \
                    pathlib.Path(pathlib.Path().root,
                                 'opt',
                                 'brainnetome',
                                 'BNA_MPM_thr25_1.25mm.nii.gz')
                self.parc_lut_file = \
                    pathlib.Path(pathlib.Path().root,
                                 'opt',
                                 'brainnetome',
                                 'BN_Atlas_246_LUT.txt')
                self.mrtrix_lut_file = None
            elif parcellation == 'craddock200':
                self.parc_image_path = \
                    pathlib.Path(pathlib.Path().root,
                                 'opt',
                                 'ADHD200_parcellate_200.nii.gz')
                self.parc_lut_file = None
                self.mrtrix_lut_file = None
            elif parcellation == 'craddock400':
                self.parc_image_path = \
                    pathlib.Path(pathlib.Path().root,
                                 'opt',
                                 'ADHD200_parcellate_400.nii.gz')
                self.parc_lut_file = None
                self.mrtrix_lut_file = None
            elif parcellation == 'perry512':
                self.parc_image_path = \
                    pathlib.Path(pathlib.Path().root,
                                 'opt',
                                 '512inMNI.nii')
                self.parc_lut_file = None
                self.mrtrix_lut_file = None
            elif parcellation == 'yeo7mni':
                self.parc_image_path = \
                    pathlib.Path(
                        pathlib.Path().root,
                        'opt',
                        'Yeo2011',
                        'Yeo2011_7Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz')
                self.parc_lut_file = \
                    pathlib.Path(
                        pathlib.Path().root,
                        'opt',
                        'Yeo2011',
                        '7Networks_ColorLUT_freeview.txt')
                self.mrtrix_lut_file = None
            elif parcellation == 'yeo17mni':
                self.parc_image_path = \
                    pathlib.Path(
                        pathlib.Path().root,
                        'opt',
                        'Yeo2011',
                        'Yeo2011_17Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz')
                self.parc_lut_file = \
                    pathlib.Path(
                        pathlib.Path().root,
                        'opt',
                        'Yeo2011',
                        '17Networks_ColorLUT_freeview.txt')
                self.mrtrix_lut_file = None
            else:
                assert False

        def find_atlas_file(filepath, description):
            if not filepath:
                return None
            if filepath.is_file():
                return filepath
            if not atlas_path:
                raise MRtrixError(f'Could not find {description} '
                                  f'(expected location: {filepath})')
            newpath = atlas_path.parent / filepath.name
            if newpath.is_file():
                return newpath
            raise MRtrixError(f'Could not find {description} '
                              f'(tested locations: "{filepath}", '
                              f'"{newpath}")')

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

        if self.mrtrix_lut_file and not self.mrtrix_lut_file.is_file():
            raise MRtrixError(
                'Could not find MRtrix3 connectome lookup table file '
                f'(expected location: {self.mrtrix_lut_file})')
