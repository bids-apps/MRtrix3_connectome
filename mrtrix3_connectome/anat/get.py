import json
import pathlib
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import fsl
from mrtrix3 import run

# Regardless of the source of T1-weighted image information,
#   scratch directory will contain at completion of this function:
# - Either:
#   - T1w.mif
#     or
#   - T1w_premasked.mif
#   , depending on software used
#   (full unmasked T1-weighted image data may not be available)
# - T1w_mask.mif
# TODO Proposing modifications where this is not necessarily guaranteed to be the case;
#   DWI pre-processing could proceed omitting the inter-modal registration step,
#   though a big fat warning would need to be issued
#   that the participant-level analysis is not immediately applicable
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
        if t1w_preproc.is_file():
            preproc_image_path = t1w_preproc
        else:
            expected_image_basename = f'{session_label}*_T1w.nii*'
            for candidate in [
                    pathlib.Path(''),
                    pathlib.Path('anat'),
                    pathlib.Path(*session, 'anat')
            ]:

                glob_result = list((t1w_preproc / candidate).glob(expected_image_basename))
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
                    raise MRtrixError(
                        'Unable to unambiguously select '
                        'pre-processed T1-weighted image '
                        'due to multiple candidates '
                        f'in location "{t1w_preproc / candidate}": '
                        f'{";".join(glob_result)}')

            if preproc_image_path is None:
                raise MRtrixError(
                    'No pre-processed T1w image found '
                    f'from specified path "{t1w_preproc}" '
                    f'for session {session_label}')

    else:

        # Look inside of import_path to see if there is a pre-processed
        #   T1w image already there
        glob_result = list(pathlib.Path(import_path, *session, 'anat')
                           .glob(f'{session_label}*_desc-preproc*_T1w.nii*'))
        if glob_result:
            if len(glob_result) == 1:
                preproc_image_path = glob_result[0]
            else:
                raise MRtrixError(
                    'Multiple pre-processed T1w images found '
                    f'in import directory "{import_path}": '
                    f'{";".join(glob_result)}')

    # Same checks regardless of whether the existing pre-processed image
    #   comes from the output directory or a user-specified location
    if preproc_image_path:

        if '_desc-preproc' not in preproc_image_path.name:
            raise MRtrixError(
                f'Selected T1-weighted image "{preproc_image_path}" '
                'not flagged as pre-processed')

        # Check to see if there's a JSON file along with the T1-weighted
        #   image; if they is, parse it to find out whether or not the
        #   pre-processed image has been brain-extracted
        expected_json_path = pathlib.Path(preproc_image_path)
        while expected_json_path.suffix:
            expected_json_path = expected_json_path.with_suffix('')
        expected_json_path = expected_json_path.with_suffix('.json')
        try:
            with open(expected_json_path, 'r', encoding='utf-8') as t1_json_file:
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
                float(run.command(f'mrcalc {preproc_image_path} 0 -eq 1 '
                                  f'{preproc_image_path} -finite -sub -add - | '
                                  ' mrstats - -output mean').stdout)
            preproc_image_is_masked = frac_voxels_outside_mask > 0.25
            app.warn('No sidecar information for pre-processed '
                     f'T1-weighted image "{preproc_image_path}" regarding skull-stripping; '
                     'image has been inferred'
                     f' to {"be" if preproc_image_is_masked else "not be"} '
                     'pre-masked based on image data '
                     f'({int(round(100.0 * frac_voxels_outside_mask))}% of voxels'
                     ' contain no data)')

        # Copy pre-procesed T1-weighted image into scratch directory
        run.command(['mrconvert',
                     preproc_image_path,
                     'T1w_premasked.mif' if preproc_image_is_masked else 'T1w.mif'])

        # If we have been provided with a pre-processed T1-weighted image
        #   (regardless of where it has come from),
        #   check to see if there is a corresponding mask image
        preproc_mask_path = preproc_image_path.parent / \
                            preproc_image_path.name \
                            .replace('_desc-preproc', '_desc-brain') \
                            .replace('_T1w.nii', '_mask.nii')
        if preproc_mask_path.is_file():
            run.command(['mrconvert',
                         preproc_mask_path,
                         'T1w_mask.mif',
                         '-datatype', 'bit'])
        elif preproc_image_is_masked:
            run.command(f'mrcalc {preproc_image_path} 0 -gt T1w_mask.mif '
                        '-datatype bit')
            # No pre-existing mask image, but we also don't want to
            #   run our own brain extraction
            preproc_mask_path = ''
        else:
            app.console('No brain mask image found alongside '
                        f'pre-processed T1-weighted image "{preproc_image_path}"; '
                        'will generate one manually')
            preproc_mask_path = None

    else:

        # Check input path for raw un-processed T1w image
        glob_result = list(pathlib.Path(import_path, *session, 'anat')
                           .glob(f'{session_label}*_T1w.nii*'))
        if not glob_result:
            raise MRtrixError('No raw or pre-processed T1-weighted images '
                              f'could be found in input directory "{import_path}" '
                              f'for session {session_label}')
        if len(glob_result) > 1:
            raise MRtrixError('Multiple raw T1w images found in '
                              f'input directory "{import_path}" '
                              f'for session {session_label}: '
                              f'{";".join(glob_result)}')
        raw_image_path = glob_result[0]

    # Do we need to do any pre-processing of our own at all?
    if preproc_mask_path is None:

        app.console('Performing requisite processing of T1-weighted data')

        # A pre-processed T1-weighted image is present,
        #   it's just the mask that is absent
        if preproc_image_path:

            if t1w_shared.synthstrip_cmd:
                app.console(f'Using SynthStrip for brain extraction for session {session_label}, '
                            'operating on existing pre-processed T1-weighted image')
            if t1w_shared.robex_cmd:
                app.console(f'Using ROBEX for brain extraction for session {session_label}, '
                            'operating on existing pre-processed T1-weighted image')
            elif t1w_shared.fsl_anat_path:
                app.console(f'Using fsl_anat for brain extraction for session {session_label} '
                            '(due to ROBEX not being installed), '
                            'operating on existing pre-processed T1-weighted image')
            else:
                raise MRtrixError(f'Unable to continue processing for session {session_label}: '
                                  'no pre-processed T1-weighted image mask available / provided, '
                                  'and no appropriate brain masking software installed')

            run.command(['mrconvert',
                        preproc_image_path,
                        'T1w.nii',
                        '-strides',
                        '+1,+2,+3' \
                            if t1w_shared.synthstrip_cmd or t1w_shared.robex_cmd \
                            else '-1,+2,+3'])

            if t1w_shared.synthstrip_cmd:
                run.command(f'{t1w_shared.synthstrip_cmd} -i T1w.nii -m T1w_mask.nii')
                run.command(['mrconvert',
                             'T1w_mask.nii',
                             'T1w_mask.mif',
                             '-datatype', 'bit'])
            elif t1w_shared.robex_cmd:
                run.command(f'{t1w_shared.robex_cmd} T1w.nii T1w_brain.nii T1w_mask.nii')
                run.command(['mrconvert',
                             'T1w_mask.nii',
                             'T1w_mask.mif',
                             '-datatype', 'bit'])
            elif t1w_shared.fsl_anat_cmd:
                run.command(f'{t1w_shared.fsl_anat_cmd} -i T1w.nii'
                            ' --noseg --nosubcortseg --nobias')
                run.command(['mrconvert',
                             fsl.find_image(pathlib.PurePath('T1w.anat', 'T1_brain_mask')),
                             'T1w_mask.mif',
                             '-datatype', 'bit'])
            else:
                assert False

        # No pre-processed T1-weighted image available:
        #   do everything based on the raw T1-weighted image
        else:

            # If we're doing pre-processing of the T1w from scratch,
            #   check to see if GDC has already been applied to the raw data,
            #   and if not, whether we have the right warp field available;
            #   if we do, then we'll perform GDC before anything else
            gdc_already_applied = None
            gdc_to_be_applied = None
            raw_json_path = pathlib.Path(raw_image_path)
            while raw_json_path.suffix:
                raw_json_path = raw_json_path.with_suffix('')
            raw_json_path = raw_json_path.with_suffix('.json')
            if raw_json_path.is_file():
                with open(raw_json_path, 'r', encoding='utf-8') as json_file:
                    json_data = json.load(json_file)
                if any(item in json_data for item in ('ImageType', 'ImageTypeText')):
                    image_type = json_data['ImageType'] \
                        if 'ImageType' in json_data \
                        else 'ImageTypeText'
                    gdc_already_applied = any(item in image_type for item in ('DIS2D', 'DIS3D'))
                    if not gdc_already_applied:
                        # Only now in the scenario where we would like to be applying GDC ourselves
                        #   will we see if we can find the corresponding warp field
                        gdc_to_be_applied = 'ManufacturersModelName' in json_data \
                            and json_data['ManufacturersModelName'] in t1w_shared.gdc_images

                else:
                    app.warn('Unable to establish whether GDC has already been applied '
                             f'to raw T1-weighted image "{raw_image_path}" '
                             'from contents of sidecar JSON file; '
                             'will not apply this correction')
                gdc_to_be_applied = False

            else:
                app.warn('Unable to establish whether GDC has already been applied '
                         f'to raw T1-weighted image "{raw_image_path}" '
                         'due to absence of sidecar JSON file; '
                         'will not apply this correction')
                gdc_to_be_applied = False

            if (t1w_shared.synthstrip_cmd or t1w_shared.robex_cmd) and t1w_shared.n4_cmd:
                app.console('No pre-processed T1-weighted image '
                            f'found for session {session_label}; '
                            f'will use {"SynthStrip" if t1w_shared.synthstrip_cmd else "ROBEX"} '
                            'and N4 for iterative brain extraction and bias field '
                            'correction from raw T1-weighted image input')
            elif t1w_shared.fsl_anat_cmd:
                app.console('No pre-processed T1-weighted image '
                            f'found for session {session_label}; '
                            'will use fsl_anat for brain extraction and '
                            'bias field correction from raw T1-weighted image input'
                            '(other softwares not installed)"')
            else:
                # TODO Make this acceptable in preproc-level analysis with adequate warning
                return
                #raise MRtrixError(f'Cannot complete processing for session {session_label}:'
                #                  'no pre-processed T1-weighted image available,'
                #                  ' and software tools for processing raw T1-weighted image'
                #                  ' not installed')

            if gdc_to_be_applied:
                run.command(['mrtransform', raw_image_path, 'T1w_raw.nii',
                             '-template', raw_image_path,
                             '-warp', t1w_shared.gdc_images[json_data['ManufacturersModelName']],
                             '-interp', 'linear',
                             '-modulate', 'jac',
                             '-strides', '+1,+2,+3'])
            else:
                run.command(['mrconvert', raw_image_path, 'T1w_raw.nii',
                             '-strides', '+1,+2,+3'])

            if (t1w_shared.synthstrip_cmd or t1w_shared.robex_cmd) and t1w_shared.n4_cmd:

                # Do a semi-iterative approach here:
                #   Get an initial brain mask, use that mask to estimate a
                #   bias field, then re-compute the brain mask
                # TODO Consider making this fully iterative, just like the
                #   approach in preproc with dwi2mask and mtnormalise
                if t1w_shared.synthstrip_cmd:
                    run.command([t1w_shared.synthstrip_cmd,
                                 '-i', 'T1w_raw.nii',
                                 '-m', 'T1w_initial_mask.nii'])
                else:
                    run.command([t1w_shared.robex_cmd,
                                'T1w_raw.nii',
                                'T1w_initial_brain.nii',
                                'T1w_initial_mask.nii'])
                    app.cleanup('T1w_initial_brain.nii')
                run.command([t1w_shared.n4_cmd,
                             '-i', 'T1w_raw.nii',
                             '-w', 'T1w_initial_mask.nii',
                             '-o', 'T1w_biascorr.nii'])
                app.cleanup('T1w_initial_mask.nii')
                if t1w_shared.synthstrip_cmd:
                    run.command([t1w_shared.synthstrip_cmd,
                                 '-i', 'T1w_biascorr.nii',
                                 '-m', 'T1w_biascorr_mask.nii'])
                else:
                    run.command([t1w_shared.robex_cmd,
                                'T1w_biascorr.nii',
                                'T1w_biascorr_brain.nii',
                                'T1w_biascorr_mask.nii'])
                    app.cleanup('T1w_biascorr_brain.nii')
                run.command(['mrconvert',
                             'T1w_biascorr.nii',
                             'T1w.mif'])
                app.cleanup('T1w_biascorr.nii')
                run.command(['mrconvert',
                             'T1w_biascorr_mask.nii',
                             'T1w_mask.mif',
                             '-datatype', 'bit'])
                app.cleanup('T1w_biascorr_mask.nii')

            else:
                assert t1w_shared.fsl_anat_cmd

                run.command(f'{t1w_shared.fsl_anat_cmd} -i T1w_raw.nii --noseg --nosubcortseg')
                run.command(['mrconvert',
                             fsl.find_image(pathlib.PurePath('T1w_raw.anat',
                                                             'T1_biascorr')),
                             'T1w_premasked.mif'])
                run.command(['mrconvert',
                             fsl.find_image(pathlib.PurePath('T1w_raw.anat',
                                                             'T1_biascorr_brain_mask')),
                             'T1w_mask.mif',
                             '-datatype', 'bit'])
                app.cleanup('T1w_raw.anat')
