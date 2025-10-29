import json
import math
import os
import pathlib
import re
import shutil
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import image
from mrtrix3 import matrix
from mrtrix3 import run
from ..anat.get import get_t1w_preproc_images

# Seem that for problematic data, running more than two iterations may
#   cause divergence from the ideal mask; therefore cap at two iterations
DWIBIASNORMMASK_MAX_ITERS = 2

# Use a threshold on the balanced tissue sum image
#   as a replacement of dwi2mask within the iterative
#   bias field correction / brain masking loop in preprpoc
TISSUESUM_THRESHOLD = 0.5 / math.sqrt(4.0 * math.pi)

OUT_DWI_JSON_DATA = {'SkullStripped': False}

# TODO Further split code across multiple files

def run_preproc(bids_dir, session, shared,
                t1w_preproc_path, concat_denoise,
                output_verbosity, output_app_dir):

    session_label = '_'.join(session)
    output_subdir = pathlib.Path(output_app_dir, 'MRtrix3_connectome-preproc', *session)
    if output_subdir.exists():
        app.warn(f'Preproc-level output directory for session "{session_label}" already exists; '
                 'all contents will be erased when this execution completes')

    cwd = pathlib.Path.cwd()
    app.activate_scratch_dir()

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
    in_dwi_path = pathlib.Path(bids_dir, *session, 'dwi')
    in_dwi_image_list = sorted(in_dwi_path.glob('*_dwi.nii*'))
    if not in_dwi_image_list:
        raise MRtrixError(f'No DWI data found for session {session_label} '
                          f'(search location: {in_dwi_path}')
    dwi_index = 0

    re_is_complex = re.compile(r'_part-(mag|phase)_')

    # TODO Check the contents of each JSON:
    # - Check to see if gradient non-linearity distortion correction
    #     has already been applied prior;
    #     if it has, then denoising and Gibbs ringing removal
    #     should be omitted from the pipeline
    #   (For now, don't attempt to handle this on a per-image basis;
    #     just make sure whatever the result is is consistent
    #     across all input images for the session)
    # - If no gradient non-linearity distortion correction
    #     has yet been applied,
    #     and the user has specified a directory
    #     containing gradient non-linearity warp fields,
    #     then load the name of the scanner,
    #     and check to see if a warp field is available
    scanner_name = None
    gdc_already_applied = None
    gdc_to_be_applied = None

    for entry in in_dwi_image_list:
        phase_rescale_factor = None
        # Is this one image in a magnitude-phase pair?
        is_complex = re_is_complex.search(entry.name)
        if is_complex:
            matching_text = is_complex.group(1)
            complex_part = matching_text.strip('_').split('-')[-1]
            assert complex_part in ['mag', 'phase']
            if complex_part == 'mag':
                # Find corresponding phase image
                in_phase_image = entry.replace('_part-mag_', '_part-phase_')
                if in_phase_image not in in_dwi_image_list:
                    raise MRtrixError(
                        f'Image {entry} does not have corresponding phase image')
                # Check if phase image is stored in radians
                phase_stats = image.statistics(in_phase_image, allvolumes=True)
                if abs(2.0*math.pi - (phase_stats.max - phase_stats.min)) \
                    > 0.01:
                    app.console(f'Phase image {in_phase_image} is not stored in radian units '
                                f'(values from {phase_stats.min}  to {phase_stats.max}); '
                                'data will be rescaled automatically')
                    # Are the values stored as integers? If so, assume that
                    #   taking the maximum phase value observed in the image
                    #   intensities, incrementing it by one, and having it
                    #   offset + scaled based on the header properties, would
                    #   result in a phase that is 2pi greater than the minimum
                    # It would be better to do this rescaling based on testing
                    #   whether the image datatype is an integer; but there
                    #   does not yet exist a module for interpreting MRtrix3
                    #   data types
                    phase_header = image.Header(in_phase_image)
                    add_to_max_phase = phase_header.intensity_scale() \
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
                        f'Image {entry} does not have corresponding mag image')
                # Do nothing for the second image in the pair
                continue
        else:

            in_phase_image = None

        # Find sidecar files
        # dcm2bids will have separate bvecs / bvals / json for the
        #   magnitude and phase images;
        # for proper BIDS compliance, these sidecar files will be stored
        #   without the "_part-[mag|phase]" part

        dwi_prefix = pathlib.Path(entry)
        while dwi_prefix.suffix:
            dwi_prefix = dwi_prefix.with_suffix('')
        sidecar_prefixes = [dwi_prefix,
                            dwi_prefix.parent / dwi_prefix.name.replace('_part-mag_', '_'),
                            bids_dir / dwi_prefix,
                            bids_dir / \
                                dwi_prefix.parent / \
                                dwi_prefix.name.replace('_part-mag_', '_')]
        entry_bval = None
        entry_bvec = None
        entry_json = None
        for prefix in sidecar_prefixes:
            candidate_bval = prefix.with_suffix('.bval')
            candidate_bvec = prefix.with_suffix('.bvec')
            candidate_json = prefix.with_suffix('.json')
            if not entry_bvec and \
                candidate_bval.is_file() and \
                candidate_bvec.is_file():
                entry_bval = candidate_bval
                entry_bvec = candidate_bvec
            if not entry_json and \
                candidate_json.is_file():
                entry_json = candidate_json
        if not entry_bvec:
            raise MRtrixError(
                'Unable to locate valid diffusion gradient table '
                f'for image "{entry}"')
        if not entry_json:
            raise MRtrixError(
                'Unable to locate valid JSON sidecar file '
                f'for image "{entry}"')

        grad_import_option = ['-fslgrad', entry_bvec, entry_bval]
        json_import_option = ['-json_import', entry_json]

        # Import the data
        dwi_index += 1
        if in_phase_image:
            run.command(['mrconvert', entry, '-']
                        + grad_import_option
                        + json_import_option
                        + ['|', 'mrcalc', '-', in_phase_image]
                        + ([str(phase_rescale_factor), '-mult'] \
                           if phase_rescale_factor \
                           else [])
                        + ['-polar', f'dwi{dwi_index}.mif'])
        else:
            run.command(['mrconvert', entry, f'dwi{dwi_index}.mif']
                        + grad_import_option
                        + json_import_option)

        # Load the JSON explicitly as we need to check coupla things
        with open(entry_json, 'r', encoding='utf-8') as json_file:
            json_data = json.load(json_file)

        if any(item in json_data for item in ('ImageType', 'ImageTypeText')):
            # Load whichever field is present, and convert into list-of-strings if necessary
            # TODO Does this ever come in a form other than a list-of-strings?
            image_type = json_data['ImageTypeText'] \
                if 'ImageTypeText' in json_data \
                else json_data['ImageType']
            this_gdc_already_applied = any(item in image_type for item in ('DIS2D', 'DIS3D'))
            if gdc_already_applied is None:
                gdc_already_applied = this_gdc_already_applied
            elif this_gdc_already_applied != gdc_already_applied:
                raise MRtrixError('Inconsistency in prior application '
                                  'of gradient non-linearity distortion correction '
                                  f'for DWI data for session {session_label}')

        if 'ManufacturersModelName' in json_data:
            if scanner_name is None:
                scanner_name = json_data['ManufacturersModelName']
                if shared.gdc_dir is not None and gdc_to_be_applied is None:
                    if scanner_name in shared.gdc_images:
                        gdc_to_be_applied = True
                    else:
                        app.warn('No gradient non-linearity warp file '
                                 f'for scanner "{scanner_name}" '
                                 f'found in GDC directory {shared.gdc_dir}; '
                                 'no gradient non-linearity distortion correction to be applied '
                                 f'for session {session_label}')
                        gdc_to_be_applied = False
            elif json_data['ManufacturersModelName'] != scanner_name:
                raise MRtrixError('Inconsistent "ManufacturersModelName" value in input DWI data')
        elif scanner_name is not None:
            raise MRtrixError('Inconsistent appearance of "ManufacturersModelName" in input DWI data')
        else:
            gdc_to_be_applied = False


    if gdc_already_applied:
        app.warn('Gradient non-linearity distortion correction already applied to input DWI data; '
                 'some pre-processing steps will need to be omitted accordingly')
        gdc_to_be_applied = False
    if gdc_to_be_applied:
        app.console(f'Scanner model "{scanner_name}" matches file "{shared.gdc_images[scanner_name]}",'
                    ' and no prior application of gradient non-linearity distortion correction indicated in image metadata;'
                    ' correction to be applied after dwifslpreproc')

    dwi_image_list = [pathlib.Path(f'dwi{index}.mif') for index in range(1, dwi_index+1)]

    # Go hunting for reversed phase-encode data
    #   dedicated to field map estimation
    in_fmap_image_list = []
    fmap_dir = pathlib.Path(bids_dir, *session, 'fmap')
    fmap_index = 0
    fmap_image_list = []
    if fmap_dir.is_dir():
        app.console('Importing fmap data into scratch directory')
        in_fmap_image_list = sorted(fmap_dir.glob('*_dir-*_epi.nii*'))
        for entry in in_fmap_image_list:
            prefix = pathlib.Path(entry)
            while prefix.suffix:
                prefix = prefix.with_suffix('')
            json_path = prefix.with_suffix('.json')
            try:
                with open(json_path, 'r', encoding='utf-8') as json_file:
                    json_elements = json.load(json_file)
            except OSError:
                app.warn(f'No JSON file found for image "{entry}"; not importing')
                continue
            if 'IntendedFor' in json_elements:
                if isinstance(json_elements['IntendedFor'], list):
                    if not any(any(str(i).endswith(target) for i in in_dwi_image_list)
                            for target in json_elements['IntendedFor']):
                        app.console(f'Image "{entry}" is not intended '
                                    'for use with DWIs; skipping')
                        continue
                elif not any(str(i).endswith(json_elements['IntendedFor'])
                           for i in in_dwi_image_list):
                    app.console(f'Image "{entry}" is not intended '
                                'for use with DWIs; skipping')
                    continue
            # Verify that the presence / absence
            #   of prior gradient non-linearity distortion correction
            #   is consistent between DWIs and fmap/ data
            if any(item in json_elements for item in ('ImageType', 'ImageTypeText')) \
                and gdc_already_applied is not None:
                image_type = json_elements['ImageType'] \
                             if 'ImageType' in json_elements \
                             else json_elements['ImageTypeText']
                if any(item in image_type for item in ('DIS2D', 'DIS3D')) != gdc_already_applied:
                    raise MRtrixError('Inconsistency in prior application '
                                      'of gradient non-linearity distortion correction '
                                      'between dwi/ and fmap/ data '
                                      f'for session {session_label}')
            # fmap files will not come with any gradient encoding in the JSON;
            #   therefore we need to add it manually ourselves so that
            #   mrcat / mrconvert can appropriately handle the table once
            #   these images are concatenated with the DWIs
            fmap_index += 1
            fmap_image_size = image.Header(entry).size()
            fmap_image_num_volumes = \
                1 if len(fmap_image_size) == 3 else fmap_image_size[3]
            fmap_dwscheme_file = pathlib.Path(f'fmap{fmap_index}.b')
            with open(fmap_dwscheme_file, 'w', encoding='utf-8') as f:
                for _ in range(0, fmap_image_num_volumes):
                    f.write('0,0,1,0\n')
            run.command(['mrconvert',
                         entry,
                         f'fmap{fmap_index}.mif',
                         '-json_import', json_path,
                         '-grad', fmap_dwscheme_file])
            app.cleanup(fmap_dwscheme_file)

        fmap_image_list = [pathlib.Path(f'fmap{index}.mif') for index in range(1, fmap_index+1)]

    # No need to explicitly check whether fmap/ images are defined
    #   on a common image grid; these we are happy to resample

    # If there's no usable data in fmap/ directory,
    #   need to check to see if there's any phase-encoding
    #   contrast within the input DWI(s)
    if not fmap_image_list and len(dwi_image_list) < 2:
        raise MRtrixError('Inadequate data for pre-processing of session '
                          f'"{session_label}": '
                          'No phase-encoding contrast in input DWIs, '
                          'and no fmap/ directory, '
                          'so EPI distortion correction cannot be performed')

    # Get T1-weighted image data
    #   (could be generated from raw data, or grabbed from a
    #   user-specified path source)
    get_t1w_preproc_images(bids_dir,
                           session,
                           shared.t1w_shared,
                           t1w_preproc_path)
    # TODO This will no longer be compulsory:
    #   if these data are not available,
    #   then issue a warning to the user,
    #   and exclude the relevant registration & export
    if pathlib.Path('T1w_premasked.mif').is_file():
        t1w_image = pathlib.Path('T1w_premasked.mif')
        t1w_is_premasked = True
    elif pathlib.Path('T1w.mif').is_file():
        t1w_image = pathlib.Path('T1w.mif')
        t1w_is_premasked = False
    else:
        app.warn(f'Pre-processing of session "{session_label}" '
                 'proceeding in the absence of any T1-weighted image; '
                 'DWI data will be pre-processed, '
                 'but no alignment to T1-weighted image will occur, '
                 'and participant-level analysis will not be applicable')
        t1w_image = None
        t1w_is_premasked = None

    dwifslpreproc_se_epi = ''
    dwifslpreproc_se_epi_option = ''

    if len(dwi_image_list) == 1:
        run.function(os.rename, dwi_image_list[0], 'dwi.mif', show=False)
        dwi_image_list[0] = pathlib.Path('dwi.mif')

    # First two steps of pre-processing are not applicable
    #   if the input data have undergone gradient non-linearity distortion correction:
    #   that requires interpolation of image intensities:
    #   which violates the assumptions of these algorithms
    if gdc_already_applied:
        app.warn('Skipping MPPCA denoising and Gibbs ringing removal '
                 'due to prior application of gradient non-linearity distortion correction '
                 'to the input DWI data')
    else:

        if concat_denoise == 'before' and len(dwi_image_list) > 1:
            # TODO We need to determine first whether these images can be trivially concatenated:
            #   if execution of dwicat would result in resampling,
            #   and that would then preclude the application of denoising / Gibbs ringing removal,
            #   then we should override this and do the concatenation after these steps
            # How do we tell if the concatenation will not involve resampling?
            # Think we should just run the process and then query the contents of stderr
            app.console('Concatenating DWI series prior to denoising')
            new_dwi_image = pathlib.Path('dwi_cat.mif')
            dwicat_stderr = run.command(['dwicat',
                                         list(map(str, dwi_image_list)),
                                         new_dwi_image]).stderr
            if 'data will be resampled onto a new average header space' in dwicat_stderr:
                app.warn('DWI series could not be concatenated without involving resampling, '
                         'which would invalidate denoising and Gibbs ringing removal; '
                         'will instead denoise separately, and concatenate after')
                os.remove(new_dwi_image)
            else:
                app.cleanup(dwi_image_list)
                dwi_image_list = [new_dwi_image]

        # Step 1: Denoise
        app.console('Denoising DWI data')
        for entry in dwi_image_list:
            run.command([shared.dwidenoise_cmd,
                         entry,
                         f'{entry.with_suffix("")}_denoised.mif'])
            app.cleanup(entry)
        dwi_image_list = [pathlib.Path(f'{entry.with_suffix("")}_denoised.mif') \
                          for entry in dwi_image_list]

        # If data are complex, take the magnitude
        new_dwi_image_list = []
        for entry in dwi_image_list:
            if image.Header(entry).datatype().startswith('CFloat'):
                mag_entry = pathlib.Path(f'{entry.with_suffix("")}_mag.mif')
                run.command(f'mrcalc {entry} -abs {mag_entry}')
                app.cleanup(entry)
                new_dwi_image_list.append(mag_entry)
            else:
                new_dwi_image_list.append(entry)
        dwi_image_list = new_dwi_image_list

        # Step 2: Gibbs ringing removal
        # TODO Can newer implementations use complex data
        #   for 2D Gibbs ringing removal?
        app.console('Performing Gibbs ringing removal for DWI'
                    f'{" and fmap" if fmap_image_list else ""} data')
        for i in dwi_image_list:
            run.command(f'mrdegibbs {i} {i.with_suffix("")}_degibbs.mif'
                        ' -nshifts 50')
            app.cleanup(i)
        dwi_image_list = [pathlib.Path(f'{i.with_suffix("")}_degibbs.mif') \
                          for i in dwi_image_list]
        for i in fmap_image_list:
            run.command(f'mrdegibbs {i} {i.with_suffix("")}_degibbs.mif'
                        ' -nshifts 50')
            app.cleanup(i)
        fmap_image_list = [pathlib.Path(f'{i.with_suffix("")}_degibbs.mif') \
                           for i in fmap_image_list]

    # We need to concatenate the DWI and fmap/ data (separately)
    #   before they can be fed into dwifslpreproc
    if len(dwi_image_list) > 1 or fmap_image_list:
        app.console(f'{"Concatenating" if len(dwi_image_list) > 1 else "Preparing"}'
                    f' DWI{" and fmap" if fmap_image_list else ""} data'
                    ' onto common voxel grid')
    if len(dwi_image_list) == 1:
        dwifslpreproc_input = dwi_image_list[0]
    else:
        dwifslpreproc_input = pathlib.Path('dwifslpreproc_in.mif')
        run.command(['dwicat', list(map(str, dwi_image_list)), dwifslpreproc_input])

    # Some decisions regarding pre-processing depend on whether or not
    #   a twice-refocused sequence has been used: a single-refocused
    #   sequence may have residual eddy current distortions in b=0
    #   volumes
    dwifslpreproc_input_header = image.Header(dwifslpreproc_input)
    monopolar = 'DiffusionScheme' in dwifslpreproc_input_header.keyval() \
                and dwifslpreproc_input_header \
                    .keyval()['DiffusionScheme'] == 'Monopolar'

    if fmap_image_list:

        fmap_transformed_image_list = []
        for item in fmap_image_list:
            affine_transform_filepath = pathlib.Path(f'{item.with_suffix("")}2dwi_affine.txt')
            rigid_transform_filepath = pathlib.Path(f'{item.with_suffix("")}2dwi_rigid.txt')
            fmap_transformed_filepath = pathlib.Path(f'{item.with_suffix("")}_transformed.mif')
            run.command(['transformcalc',
                         item,
                         dwifslpreproc_input,
                         'header',
                         affine_transform_filepath])
            run.command(['transformcalc',
                         affine_transform_filepath,
                         'rigid',
                         rigid_transform_filepath])
            run.command(['mrtransform',
                         item,
                         '-linear', rigid_transform_filepath,
                         '-reorient_fod', 'no',
                         fmap_transformed_filepath])
            fmap_transformed_image_list.append(fmap_transformed_filepath)
            app.cleanup(affine_transform_filepath)
            app.cleanup(rigid_transform_filepath)
        app.cleanup(fmap_image_list)
        if len(fmap_transformed_image_list) == 1:
            dwifslpreproc_se_epi = fmap_transformed_image_list[0]
        else:
            dwifslpreproc_se_epi = pathlib.Path('dwifslpreproc_seepi.mif')
            try:
                run.command(['mrcat',
                             list(map(str, fmap_transformed_image_list)),
                             dwifslpreproc_se_epi,
                             '-axis', '3'])
            except run.MRtrixCmdError:
                app.warn('Unable to rigidly align fmap/ images to DWI voxel grid; '
                         'performing explicit interpolation')
                fmap_resampled_image_list = []
                for item in fmap_transformed_image_list:
                    fmap_resampled_image_path = pathlib.Path(f'{item.with_suffix("")}_resampled.mif')
                    run.command(['mrtransform',
                                 item,
                                 '-template', dwifslpreproc_input,
                                 '-interp', 'sinc',
                                 fmap_resampled_image_path])
                    fmap_resampled_image_list.append(fmap_resampled_image_path)
                app.cleanup(fmap_transformed_image_list)
                run.command(['mrcat',
                             list(map(str, fmap_resampled_image_list)),
                             dwifslpreproc_se_epi,
                             '-axis', '3'])
                app.cleanup(fmap_resampled_image_list)
        dwifslpreproc_se_epi_option = ['-se_epi', dwifslpreproc_se_epi,
                                       '-align_seepi']

    else: # No fmap/ images

        dwifslpreproc_se_epi = None
        dwifslpreproc_se_epi_option = []

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
                    bzero_image = pathlib.Path(f'{dwi_image.with_suffix("")}_bzero.mif')
                    run.command(['mrconvert',
                                 dwi_image,
                                 bzero_image,
                                 '-coord',
                                 '3',
                                 ','.join(str(i) for i in
                                          range(0, first_nonzero_volume))])
                    bzero_image_list.append(bzero_image)
                dwifslpreproc_se_epi = pathlib.Path('dwifslpreproc_seepi.mif')
                run.command(['mrcat',
                             bzero_image_list,
                             dwifslpreproc_se_epi,
                             '-axis', '3'])
                app.cleanup(bzero_image_list)
                dwifslpreproc_se_epi_option = ['-se_epi', str(dwifslpreproc_se_epi)]
            except MRtrixError:
                dwifslpreproc_se_epi = None
                dwifslpreproc_se_epi_option = None
                app.warn('DWIs detected as using monopolar diffusion sensitisation,'
                         ' but error encountered in extracting pre-DWI b=0 volumes;'
                         ' topup field estimate may be affected by eddy current distortions'
                         ' in b=0 volumes')

    # If only one image, this is fed directly to dwifslpreproc as-is
    if len(dwi_image_list) > 1:
        app.cleanup(dwi_image_list)

    # Step 3: Distortion correction
    app.console('Performing various geometric corrections of DWIs')
    dwifslpreproc_input_header = image.Header(dwifslpreproc_input)
    have_slice_timing = 'SliceTiming' in dwifslpreproc_input_header.keyval()
    app.debug(f'Have slice timing: {have_slice_timing}')
    mb_factor = int(dwifslpreproc_input_header.keyval()
                    .get('MultibandAccelerationFactor', '1'))
    app.debug(f'Multiband factor: {mb_factor}')
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
            app.warn('Error reading BIDS field "SliceDirection" '
                     f'(value: "{slice_direction_code}"); '
                     'assuming third axis')
    else:
        num_slices = dwifslpreproc_input_header.size()[2]
    app.debug(f'Number of slices: {num_slices}')
    mporder = 1 + int(math.ceil(num_slices/(mb_factor*4)))
    app.debug(f'MPorder: {mporder}')

    eddy_options = ['--flm=cubic'] if shared.eddy_cubicflm else []
    if shared.eddy_repol:
        eddy_options.append('--repol')
    if shared.eddy_mporder and have_slice_timing:
        eddy_options.append('--mporder=' + str(mporder))
    if shared.eddy_mbs:
       eddy_options.append('--estimate_move_by_susceptibility')
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
    # TODO Investigate whether this is still the case
    #if monopolar:
    #    eddy_options.append('--b0_flm=linear')

    # Make sure that MRtrix3 identifies the concatenated data as shelled
    # (If FSL eddy complains, we can force it to comply;
    #   but it's essential for many aspects of future processing that MRtrix3
    #   consider the comprehensive set of data to be shelled)
    try:
        run.command(['mrinfo', dwifslpreproc_input, '-shell_bvalues'], show=False)
    except run.MRtrixCmdError as exc:
        raise MRtrixError('Combined DWI data are not classified as shelled') from exc

    shell_asymmetries = \
        [float(value) for value in
         run.command(f'dirstat {dwifslpreproc_input} -output asym',
                     show=False)[0]
         .splitlines()]
    app.debug(f'Shell asymmetries: {shell_asymmetries}')
    if any(value > 0.1 for value in shell_asymmetries):
        app.console('Utilising eddy linear second-level model due to poor '
                    'distribution of diffusion gradient direction polarities')
        eddy_options.append('--slm=linear')

    run.function(os.makedirs, 'eddyqc', show=False)
    dwifslpreproc_output = pathlib.Path(
                               'dwifslpreproc_out.mif' \
                               if dwifslpreproc_input == 'dwifslpreproc_in.mif' \
                               else (f'{dwifslpreproc_input.with_suffix("")}_preproc.mif'))

    eddy_olnstd_value = 4.0 # The internal eddy default
    eddy_olnstd_option = []
    eddy_force_shelled_option = []

    while not dwifslpreproc_output.is_file():

        try:

            # If dwifslpreproc fails due to:
            # EDDY:::  DoVolumeToVolumeRegistration: Unable to find volume
            #          with no outliers in shell 0 with b-value=549.375
            # , want to progressively increase the outlier rejection threshold
            #   until this error no longer occurs.
            eddy_all_options = eddy_options \
                + eddy_olnstd_option \
                + eddy_force_shelled_option
            dwifslpreproc_eddy_options = \
                ['-eddy_options',
                 ' '.join(eddy_all_options)] \
                if eddy_all_options \
                else []

            run.command(['dwifslpreproc',
                         dwifslpreproc_input,
                         dwifslpreproc_output,
                         '-rpe_header',
                         '-eddyqc_text', 'eddyqc/']
                        + dwifslpreproc_se_epi_option
                        + dwifslpreproc_eddy_options
                        + ([] \
                           if app.DO_CLEANUP \
                           else ['-scratch', app.SCRATCH_DIR, '-nocleanup']))

        except run.MRtrixCmdError as e_dwifslpreproc:
            if any(item in str(e_dwifslpreproc) for item in [
                    'msg=ECScanManager::set_slice_to_vol_reference: ' \
                    'ref index out of bounds',
                    'Unable to find volume with no outliers']):
                eddy_olnstd_value += 0.5
                eddy_olnstd_option = [f'--ol_nstd={eddy_olnstd_value}']
                app.warn('FSL eddy failed due to outlier rejection; '
                         're-running with increased threshold')
            elif 'Data not shelled' in str(e_dwifslpreproc):
                eddy_force_shelled_option = ['--data_is_shelled']
                app.warn('FSL eddy failed due to reporting DWI data as not '
                         'being shelled; despite MRtrix3 classifying as '
                         'shelled; re-running with --data_is_shelled option')
            elif not eddy_force_shelled_option:
                eddy_force_shelled_option = ['--data_is_shelled']
                app.warn('FSL eddy failed with unrecognised error; '
                         're-running with --data_is_shelled option '
                         'in case it resolves issue')
            else:
                raise
            shutil.rmtree('eddyqc/')

    app.cleanup(dwifslpreproc_input)
    app.cleanup(dwifslpreproc_se_epi)

    # TODO Step 4: Gradient non-linearity distortion correction
    # This needs to happen before computing mask or registration to the T1w image
    # TODO Longer-term, could compose the gradient non-linearity warp field
    #   with the DWI->T1w transformation
    #   to resample the DWI data on a voxel grid aligned with the T1w image grid
    #   with only a single interpolation step
    # TODO Would like to change where / how this is applied:
    #   - If only doing volume motion correction,
    #     this should be applied after Gibbs ringing removal before eddy;
    #   - If doing slice-to-volume motion correction,
    #     would like to apply in-plane 2D distortion correction
    #     after Gibbs ringing removal before eddy,
    #     and apply the distortion correction along the slice axis after eddy
    #   To achieve this would require augmentation of probably warpconvert
    #     to have the ability to project a warp field to identity along specific directions
    if gdc_to_be_applied:
        assert scanner_name is not None
        assert scanner_name in shared.gdc_images
        app.console('Applying gradient non-linearity distortion correction')
        dwi_gdc_image = pathlib.Path('dwi_gdc.mif')
        run.command(['mrtransform', dwifslpreproc_output, dwi_gdc_image,
                    '-template', dwifslpreproc_output,
                    '-warp', shared.gdc_images[scanner_name],
                    '-interp', 'cubic',
                    '-modulate', 'jac'])
        dwi_image = dwi_gdc_image
    else:
        app.console('Skipping gradient non-linearity distortion correction')
        dwi_image = dwifslpreproc_output

    # Step 5:
    #   Combined bias field correction,
    #   intensity normalisation,
    #   and brain mask derivation
    app.console('Running simultaneous bias field correction, '
                'intensity normalisation, '
                'and DWI brain mask derivation '
                'via dwibiasnormmask command')
    dwi_biasnorm_image = pathlib.Path('dwi_biasnorm.mif')
    dwi_mask_image = pathlib.Path('dwi_mask.mif')
    # Note that:
    # 1. The first of these results in the synthstrip command
    #    being utilised in the initial dwi2mask call
    #    to derive an initial mask prior to the first iteration
    # 2. The second of these results in the synthstrip command
    #    being run with the ODF sum image as input
    #    during the iterative process
    mask_algo_options = ['-config', 'Dwi2maskAlgo', 'synthstrip',
                         '-mask_algo', 'synthstrip'] \
                        if shared.dwi2mask_algo == 'synthstrip' \
                        else []
    run.command(['dwibiasnormmask',
                 dwi_image,
                 dwi_biasnorm_image,
                 dwi_mask_image]
                + mask_algo_options)
    app.cleanup(dwi_image)
    dwi_image = dwi_biasnorm_image

    # Step 6: Crop images to reduce storage space
    #   (but leave some padding on the sides)
    dwi_cropped_image = pathlib.Path('dwi_crop.mif')
    dwi_cropped_mask_image = pathlib.Path('mask_crop.mif')
    run.command(f'mrgrid {dwi_image} crop {dwi_cropped_image} '
                f'-mask {dwi_mask_image} -uniform -3')
    app.cleanup(dwi_image)
    dwi_image = dwi_cropped_image
    run.command(f'mrgrid {dwi_mask_image} crop {dwi_cropped_mask_image} '
                f'-mask {dwi_mask_image} -uniform -3')
    app.cleanup(dwi_mask_image)
    dwi_mask_image = dwi_cropped_mask_image

    # Step 7: DWI->T1 registration
    if t1w_image:

        # Step 7.1: Generate target images for T1w->DWI registration
        app.console('Generating contrast-matched images for '
                    'inter-modal registration between DWIs and T1w')
        run.command(f'dwiextract {dwi_image} -bzero - | '
                    'mrcalc - 0.0 -max - | '
                    'mrmath - mean -axis 3 dwi_meanbzero.mif')
        run.command(f'mrcalc 1 dwi_meanbzero.mif -div {dwi_mask_image} -mult - | '
                    f'mrhistmatch nonlinear - {t1w_image} dwi_pseudoT1w.mif '
                    f'-mask_input {dwi_mask_image} '
                    '-mask_target T1w_mask.mif')
        run.command(f'mrcalc 1 {t1w_image} -div '
                    f'{"" if t1w_is_premasked else "T1w_mask.mif -mult "}- | '
                    'mrhistmatch nonlinear - dwi_meanbzero.mif T1w_pseudobzero.mif '
                    '-mask_input T1w_mask.mif '
                    f'-mask_target {dwi_mask_image}')

        # Step 7.2: Perform DWI->T1w registration
        #   Note that two registrations are performed:
        #   Even though we have a symmetric registration, generation of the
        #   two histogram-matched images means that you will get slightly
        #   different answers depending on which synthesized image &
        #   original image you use
        app.console('Performing registration between DWIs and T1w')
        transform_pt1w_t1w = pathlib.Path('rigid_pseudoT1w_to_T1w.txt')
        transform_b0_pb0 = pathlib.Path('rigid_bzero_to_pseudobzero.txt')
        run.command(f'mrregister dwi_pseudoT1w.mif {t1w_image}'
                    ' -type rigid'
                    f' -mask1 {dwi_mask_image}'
                    ' -mask2 T1w_mask.mif'
                    f' -rigid {transform_pt1w_t1w}')
        run.command('mrregister dwi_meanbzero.mif T1w_pseudobzero.mif'
                    ' -type rigid'
                    f' -mask1 {dwi_mask_image}'
                    ' -mask2 T1w_mask.mif'
                    f' -rigid {transform_b0_pb0}')
        app.cleanup('dwi_meanbzero.mif')

        # Step 7.3: Perform DWI->T1w transformation
        # In this scenario, we're going to transform the DWI data to the T1w
        #   rather than the other way around, since the T1w is more likely to
        #   be used as a common reference across multiple analysis pipelines,
        #   and we're transforming DWIs rather than FODs
        transform_average = pathlib.Path('rigid_dwi_to_T1w.txt')
        run.command(['transformcalc',
                     transform_pt1w_t1w,
                     transform_b0_pb0,
                     'average',
                     transform_average])
        app.cleanup(transform_pt1w_t1w)
        app.cleanup(transform_b0_pb0)
        transformed_dwi_image = pathlib.Path(f'{dwi_image.with_suffix("")}_transform.mif')
        transformed_dwi_mask_image = pathlib.Path(f'{dwi_mask_image.with_suffix("")}_transform.mif')
        run.command(['mrtransform',
                     dwi_image,
                     transformed_dwi_image,
                     '-linear', transform_average])
        app.cleanup(dwi_image)
        dwi_image = transformed_dwi_image
        run.command(['mrtransform',
                     dwi_mask_image,
                     transformed_dwi_mask_image,
                     '-linear', transform_average])
        app.cleanup(dwi_mask_image)
        app.cleanup(transform_average)
        dwi_mask_image = transformed_dwi_mask_image

    # Processing completed; export
    app.console(f'Processing completed for session "{session_label}"; '
                'writing results to output directory')
    if output_subdir.exists():
        run.function(shutil.rmtree, output_subdir)
    run.function(os.makedirs, output_subdir)
    if t1w_image:
        run.function(os.makedirs, output_subdir / 'anat')
    run.function(os.makedirs, output_subdir / 'dwi')
    run.command(['mrconvert',
                 dwi_image,
                 output_subdir / 'dwi' / f'{session_label}_desc-preproc_dwi.nii.gz',
                '-export_grad_fsl',
                 output_subdir / 'dwi' / f'{session_label}_desc-preproc_dwi.bvec',
                 output_subdir / 'dwi' / f'{session_label}_desc-preproc_dwi.bval',
                '-strides', '+1,+2,+3,+4'])
    with open(output_subdir / 'dwi' / f'{session_label}_desc-preproc_dwi.json',
              'w',
              encoding='utf-8') as out_dwi_json_file:
        json.dump(OUT_DWI_JSON_DATA, out_dwi_json_file)
    run.command(['mrconvert',
                 dwi_mask_image,
                 output_subdir / 'dwi' / f'{session_label}_desc-brain_mask.nii.gz',
                 '-datatype', 'uint8',
                 '-strides', '+1,+2,+3'])
    # Even if EddyQC software is not installed, a directory is still
    #   generated containing some eddy outputs
    run.function(shutil.copytree,
                 'eddyqc',
                 output_subdir / 'dwi' / 'eddyqc')
    if t1w_image:
        run.command(['mrconvert',
                     t1w_image,
                     output_subdir / 'anat' / f'{session_label}_desc-preproc_T1w.nii.gz',
                     '-strides', '+1,+2,+3'])
        t1w_json_data = {"SkullStripped": t1w_is_premasked}
        with open(output_subdir / 'anat' / f'{session_label}_desc-preproc_T1w.json',
                  'w',
                  encoding='utf-8') as t1w_json_file:
            json.dump(t1w_json_data, t1w_json_file)
        run.command(['mrconvert',
                     'T1w_mask.mif',
                     output_subdir / 'anat' / f'{session_label}_desc-brain_mask.nii.gz',
                     '-datatype', 'uint8',
                     '-strides', '+1,+2,+3'])

    # Manually wipe and zero the scratch directory
    #   (since we might be processing more than one subject)
    os.chdir(cwd)
    if app.DO_CLEANUP:
        app.console(f'Deleting scratch directory {app.SCRATCH_DIR}')
        # Can't use run.function() here;
        #   it'll try to write to the log file that resides
        #   in the scratch directory just deleted
        app.cleanup(app.SCRATCH_DIR)
    elif output_verbosity == 4:
        app.console('Copying scratch directory to output location')
        run.function(shutil.copytree,
                     app.SCRATCH_DIR,
                     output_subdir / 'scratch')
    else:
        app.console('Contents of scratch directory kept; '
                    f'location: {app.SCRATCH_DIR}')
    app.SCRATCH_DIR = None
