import os
import shutil
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import fsl
from mrtrix3 import run
from ..anat.shared import Shared as T1wShared

class Shared(object): #pylint: disable=useless-object-inheritance
    def __init__(self, gdc_dir, concat_denoise, eddy_cubicflm, eddy_mbs):
        self.gdc_dir = gdc_dir
        self.concat_denoise = concat_denoise
        self.eddy_cubicflm = eddy_cubicflm
        self.eddy_mbs = eddy_mbs

        self.gdc_images = {}
        if self.gdc_dir is not None:
            gdc_dir_images = self.gdc_dir.glob('*.*')
            for item in gdc_dir_images:
                scanner_name = item.name.split('.')[0]
                if scanner_name in self.gdc_images:
                    raise MRtrixError(
                        f'Duplicate images for scanner "{scanner_name}" '
                        f'in GDC directory {self.gdc_dir}')
                self.gdc_images[item.name.split('.')[0]] = item
            app.debug(f'{len(self.gdc_images)} gradient non-linearity warp '
                      f'field images found in directory {self.gdc_dir}:'
                      f'{self.gdc_images}')

        fsl_path = os.environ.get('FSLDIR', '')
        if not fsl_path:
            raise MRtrixError(
                'Environment variable FSLDIR is not set; '
                'please run appropriate FSL configuration script')

        self.t1w_shared = T1wShared(self.gdc_images)

        dwidenoise2_cmd = shutil.which('dwidenoise2')
        if dwidenoise2_cmd:
            self.dwidenoise_cmd = [dwidenoise2_cmd]
        else:
            app.warn('dwidenoise2 command not available; '
                     'original dwidenoise implementation will be used')
            self.dwidenoise_cmd = 'dwidenoise'

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

        self.dwibiascorrect_algo = 'ants'
        if not self.t1w_shared.n4_cmd:
            self.dwibiascorrect_algo = None
            app.warn('Could not find ANTs program "N4BiasFieldCorrection"; '
                     'will proceed without performing initial b=0 - based '
                     'DWI bias field correction')

        self.dwi2mask_algo = 'synthstrip'
        if not shutil.which('mri_synthstrip'):
            self.dwi2mask_algo = 'legacy'
            app.warn('FreeSurfer command mri_synthstrip not present;'
                     ' legacy dwi2mask algorithm will be used')
