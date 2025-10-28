import shutil
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import fsl

class Shared(object): #pylint: disable=useless-object-inheritance
    def __init__(self, gdc_images):
        self.gdc_images = gdc_images
        # TODO If synthstrip is available,
        #   prioritise using that for brain extraction
        try:
            self.fsl_anat_cmd = shutil.which(fsl.exe_name('fsl_anat'))
        except MRtrixError:
            self.fsl_anat_cmd = None
        robex_cmd = shutil.which('ROBEX')
        self.robex_cmd = robex_cmd if robex_cmd else shutil.which('runROBEX.sh')
        self.n4_cmd = shutil.which('N4BiasFieldCorrection')

        if not self.fsl_anat_cmd and not self.robex_cmd:
            app.warn('No commands for T1w image processing found; '
                     'command can only proceed if either '
                     'existing pre-processed T1w image data can be found, '
                     'or if doing preproc-level analysis '
                     'where registration of DWI to T1-weighted image data '
                     'will be excluded from processing')
