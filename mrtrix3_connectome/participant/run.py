import glob
import json
import os
import pathlib
import shutil
from collections import namedtuple
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import fsl
from mrtrix3 import image
from mrtrix3 import run
from ..anat.get import get_t1w_preproc_images

OUT_5TT_JSON_DATA = {'LabelMap': ['CGM', 'SGM', 'WM', 'CSF', 'Path']}

def run_participant(bids_dir, session, shared,
                    t1w_preproc_path, output_verbosity, output_app_dir):

    session_label = '_'.join(session)
    output_analysis_level_path = output_app_dir / 'MRtrix3_connectome-participant'
    output_subdir = pathlib.Path(output_analysis_level_path, *session)

    if output_subdir.exists():
        app.warn('Participant-level output directory'
                 f' for session "{session_label}" already exists;'
                 ' all contents will be erased when this execution completes')

    # Check paths of individual output files before script completion
    #   by building a database of what files are to be written to output
    parc_string = f'_desc-{shared.parcellation}'
    OutputItem = \
        namedtuple(
            'OutputItem',
            'is_image min_verbosity needs_multishell options path')
    output_items = {
        'response_wm.txt': \
            OutputItem(False, 1, False, None,
                       pathlib.PurePath('dwi', f'{session_label}_tissue-WM_response.txt')),
        'T1w_mask.mif': \
            OutputItem(True, 1, False,
                       '-strides +1,+2,+3 -datatype uint8',
                       pathlib.PurePath('anat', f'{session_label}_desc-brain_mask.nii.gz')),
        '5TT.mif': \
            OutputItem(True, 2, False,
                       '-strides +1,+2,+3,+4',
                       pathlib.PurePath('anat', f'{session_label}_desc-5tt_probseg.nii.gz')),
        '5TT.json': \
            OutputItem(False, 2, False, None,
                       pathlib.PurePath('anat', f'{session_label}_desc-5tt_probseg.json')),
        'vis.mif': \
            OutputItem(True, 2, False, '-strides +1,+2,+3',
                       pathlib.PurePath('anat', f'{session_label}_desc-vis_probseg.nii.gz')),
        'FOD_WM.mif': \
            OutputItem(True, 2, False, '-strides +1,+2,+3,+4',
                       pathlib.PurePath('dwi', f'{session_label}_tissue-WM_ODF.nii.gz')),
        'response_gm.txt': \
            OutputItem(False, 2, True, None,
                       pathlib.PurePath('dwi', f'{session_label}_tissue-GM_response.txt')),
        'response_csf.txt': \
            OutputItem(False, 2, True, None,
                       pathlib.PurePath('dwi', f'{session_label}_tissue-CSF_response.txt')),
        'FOD_GM.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       pathlib.PurePath('dwi', f'{session_label}_tissue-GM_ODF.nii.gz')),
        'FOD_CSF.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       pathlib.PurePath('dwi', f'{session_label}_tissue-CSF_ODF.nii.gz')),
        'tissues.mif': \
            OutputItem(True, 2, True, '-strides +1,+2,+3,+4',
                       pathlib.PurePath('dwi', f'{session_label}_tissue-all_probseg.nii.gz'))
    }

    if shared.parcellation != 'none':
        output_items['connectome.csv'] = \
            OutputItem(False, 1, False, None,
                       pathlib.PurePath('connectome',
                                        f'{session_label}{parc_string}_connectome.csv'))
        output_items['mu.txt'] = \
            OutputItem(False, 1, False, None,
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_mu.txt'))
        output_items['parc.mif'] = \
            OutputItem(True, 2, False, '-strides +1,+2,+3',
                       pathlib.PurePath('anat',
                                        f'{session_label}{parc_string}_dseg.nii.gz'))
        output_items['meanlength.csv'] = \
            OutputItem(False, 2, False, None,
                       pathlib.PurePath('connectome',
                                        f'{session_label}{parc_string}_meanlength.csv'))
        output_items['assignments.csv'] = \
            OutputItem(False, 3, False, None,
                       pathlib.PurePath('connectome',
                                        f'{session_label}{parc_string}_assignments.csv'))
        output_items['nodes_smooth.obj'] = \
            OutputItem(False, 3, False, None,
                       pathlib.PurePath('anat',
                                        f'{session_label}{parc_string}_dseg.obj'))
        output_items['exemplars.tck'] = \
            OutputItem(False, 3, False, None,
                       pathlib.PurePath('connectome',
                                        f'{session_label}{parc_string}_exemplars.tck'))
        output_items['parcRGB.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3,+4',
                       pathlib.PurePath('anat',
                                        f'{session_label}{parc_string}_desc-rgb_dseg.nii.gz'))

    if shared.streamlines or shared.parcellation != 'none':
        output_items['tractogram.tck'] = \
            OutputItem(False, 3, False, None,
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_tractogram.tck'))
        output_items['weights.csv'] = \
            OutputItem(False, 3, False, None,
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_weights.csv'))
        output_items['tdi_dwi.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_space-dwi_tdi.nii.gz'))
        output_items['tdi_T1w.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_space-T1w_tdi.nii.gz'))
        output_items['tdi_hires.mif'] = \
            OutputItem(True, 3, False, '-strides +1,+2,+3',
                       pathlib.PurePath('tractogram',
                                        f'{session_label}_space-superres_tdi.nii.gz'))


    def do_import(import_path):
        in_dwi_image_list = list(pathlib.Path(import_path, *session, 'dwi').glob('*_dwi.nii*'))
        if not in_dwi_image_list:
            raise MRtrixError(f'No DWIs found for session "{session_label}"')
        if len(in_dwi_image_list) > 1:
            raise MRtrixError('To run participant-level analysis, '
                              'input directory should contain only one DWI image file; '
                              f'session "{session_label}" contains {len(in_dwi_image_list)}')
        in_dwi_path = in_dwi_image_list[0]
        if not '_desc-preproc_' in in_dwi_path:
            raise MRtrixError(f'Input DWI image "{in_dwi_path}" not flagged as pre-processed data')
        in_dwi_path_prefix = pathlib.Path(in_dwi_path)
        while in_dwi_path_prefix.suffix:
            in_dwi_path_prefix = in_dwi_path_prefix.with_suffix('')
        # Don't look for bvec / bval in a lower directory in this case
        in_bvec_path = in_dwi_path_prefix.with_suffix('.bvec')
        in_bval_path = in_dwi_path_prefix.with_suffix('.bval')
        if not in_bvec_path.is_file() \
            or not in_bval_path.is_file():
            raise MRtrixError('Did not find bvec / bval pair '
                              f'corresponding to image {in_dwi_path} '
                              f'(expected locations: "{in_bvec_path}" "{in_bval_path}")')
        # JSON isn't compulsory in this case
        in_dwi_json_path = in_dwi_path_prefix.with_suffix('.json')
        in_dwi_json_import_option = ['-json_import', in_dwi_json_path] \
                                    if in_dwi_json_path.is_file() \
                                    else []
        # Is there a mask present?
        in_dwi_mask_image_list = list(pathlib.Path(import_path, *session, 'dwi')
                                      .glob('*_desc-brain*_mask.nii*'))
        if len(in_dwi_mask_image_list) > 1:
            raise MRtrixError(f'More than one DWI mask found for session "{session_label}"')
        in_dwi_mask_path = in_dwi_mask_image_list[0] \
                           if in_dwi_mask_image_list \
                           else None
        if not in_dwi_mask_path:
            output_items['dwi_mask.mif'] = \
                OutputItem(True, 1, False,
                           '-strides +1,+2,+3 -datatype uint8',
                           pathlib.PurePath('dwi', f'{session_label}_desc-brain_mask.nii.gz'))

        app.console('Importing pre-processed data into scratch directory')

        run.command(['mrconvert',
                     in_dwi_path,
                     'dwi.mif',
                     '-fslgrad', in_bvec_path, in_bval_path,
                     '-strides', '0,0,0,1']
                    + in_dwi_json_import_option)

        if in_dwi_mask_path:
            run.command(['mrconvert',
                         in_dwi_mask_path,
                         'dwi_mask.mif',
                         '-datatype', 'bit'])

        get_t1w_preproc_images(import_path,
                               session,
                               shared.t1w_shared,
                               t1w_preproc_path)
        if pathlib.Path('T1w_premasked.mif').is_file():
            t1w_is_premasked = True
        elif pathlib.Path('T1w.mif').is_file():
            t1w_is_premasked = False
        else:
            raise MRtrixError('No T1-weighted image found '
                              f'for session {session_label}; '
                              'cannot perform participant-level analysis')
        if shared.do_freesurfer and t1w_is_premasked:
            raise MRtrixError('Cannot execute FreeSurfer for obtaining parcellation:'
                              ' input T1-weighted image is already skull-stripped')
    # End of do_import() function


    # We first make an attempt at loading all requisite data from
    #   "bids_dir" (since the user may have used that path to request
    #   that the pre-processed data be utilised from some path other than
    #   "mrtrix3_connectome-preproc/"); if that doesn't work, we wipe the
    #   scratch directory and try again based on the latter
    cwd = pathlib.Path.cwd()
    app.activate_scratch_dir()
    try:
        do_import(bids_dir)
    except MRtrixError as e_frombids:
        for item in app.SCRATCH_DIR.iterdir():
            os.remove(item)
        try:
            preproc_dir = output_app_dir / 'MRtrix3_connectome-preproc'
            do_import(preproc_dir)
        except MRtrixError as e_fromoutput:
            err = 'Unable to import requisite pre-processed data' \
                  ' from either specified input directory' \
                  '  or MRtrix3_connectome output directory\n'
            err += f'Error when attempting load from "{bids_dir}":\n'
            err += str(e_frombids) + '\n'
            err += f'Error when attempting load from "{preproc_dir}":\n'
            err += str(e_fromoutput)
            raise MRtrixError(err) # pylint: disable=raise-missing-from


    # T1-weighted data are always written to output directory regardless;
    #   output paths can only be constructed now
    t1w_image = pathlib.Path('T1w_premasked.mif')
    t1w_is_premasked = True
    if not t1w_image.is_file():
        t1w_image = pathlib.Path('T1w.mif')
        assert t1w_image.is_file()
        t1w_is_premasked = False
    output_items[str(t1w_image)] = \
        OutputItem(True, 1, False, ' -strides +1,+2,+3',
                   pathlib.PurePath('anat', f'{session_label}_desc-preproc_T1w.nii.gz'))
    t1w_json_data = {"SkullStripped": t1w_is_premasked}
    t1w_json_path = pathlib.Path(t1w_image)
    while t1w_json_path.suffix:
        t1w_json_path = t1w_json_path.with_suffix('')
    t1w_json_path = t1w_json_path.with_suffix('.json')
    with open(t1w_json_path, 'w', encoding='utf-8') as t1w_json_file:
        json.dump(t1w_json_data, t1w_json_file)
    output_items[str(t1w_json_path)] = \
        OutputItem(False, 1, False, None,
                   pathlib.PurePath('anat', f'{session_label}_desc-preproc_T1w.json'))

    # Before we can begin: Are there any data we require
    #   that were not imported from the output directory?
    if not pathlib.Path('dwi_mask.mif').is_file():
        app.console('Generating DWI brain mask '
                    '(was not already present in pre-processing directory)')
        run.command(f'dwi2mask {shared.dwi2mask_algo} dwi.mif dwi_mask.mif')

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
    multishell = len(bvalues) > 2

    # Step 2: Perform spherical deconvolution
    #   Don't even use a processing mask:
    #     ACT should be responsible for stopping streamlines before they
    #     reach the edge of the DWI mask
    #   Also means that any subsequent manual use of the FOD
    #     images can't possibly be detrimentally affected by
    #     bad masking
    deconvolution_msg = 'multi-tissue ODF images' \
                        if multishell \
                        else 'Fibre Orientation Distribution image'
    app.console(f'Estimating {deconvolution_msg}')
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
    # Use T1w brain mask generated from elsewhere:
    #   don't particularly trust the raw "bet" call inside
    #   5ttgen fsl
    app.console('Generating five-tissue-type (5TT) image for '
                'Anatomically-Constrained Tractography (ACT)')
    run.command(['5ttgen', 'fsl',
                 t1w_image,
                 '5TT.mif']
                + (['-premasked'] \
                   if t1w_is_premasked \
                   else ['-mask', 'T1w_mask.mif']))
    if output_verbosity > 1:
        with open('5TT.json', 'w', encoding='utf-8') as out_5tt_json_file:
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
                         shared.freesurfer_subjects_dir / subdir,
                         subdir)

        # Run FreeSurfer pipeline on this subject's T1w image
        # If the pre-processed T1-weighted image is not brain-extracted,
        #   we'll use that here; but if it is, fingers crossed we have the
        #   raw T1-weighted image that was used to generate it...
        # TODO Improve on this
        if t1w_is_premasked:
            freesurfer_t1w_input = pathlib.Path('T1w_raw.nii')
            if not os.path.isfile('T1w_raw.nii'):
                raise MRtrixError(
                    'Cannot run FreeSurfer: '
                    'pre-processed T1-weighted image is skull-stripped')
        else:
            freesurfer_t1w_input = 'T1w.nii'
            run.command(f'mrconvert T1w.mif {freesurfer_t1w_input} -strides +1,+2,+3')
        run.command(['recon-all',
                     '-sd', app.SCRATCH_DIR,
                     '-subjid', 'freesurfer',
                     '-i', freesurfer_t1w_input])
        run.command(['recon-all',
                     '-sd', app.SCRATCH_DIR,
                     '-subjid', 'freesurfer',
                     '-all']
                     + shared.reconall_multithread_options)

        # Grab the relevant parcellation image and
        #   target lookup table for conversion
        parc_image_path = pathlib.Path(app.SCRATCH_DIR, 'freesurfer', 'mri')
        if shared.parcellation == 'desikan':
            parc_image_path = parc_image_path / 'aparc+aseg.mgz'
        elif shared.parcellation == 'destrieux':
            parc_image_path = parc_image_path / 'aparc.a2009s+aseg.mgz'
        else:
            # Non-standard parcellations are not applied as part of
            #   the recon-all command; need to explicitly map them to
            #   the subject
            # This requires SUBJECTS_DIR to be set;
            #   commands don't have a corresponding -sd option like recon-all
            env = run.shared.env
            env['SUBJECTS_DIR'] = str(app.SCRATCH_DIR)
            if shared.parcellation == 'brainnetome246fs':
                for index, hemi in enumerate(['l', 'r']):
                    run.command([
                        'mris_ca_label',
                        '-l',
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'label',
                                     f'{hemi}h.cortex.label'),
                        'freesurfer',
                        f'{hemi}h',
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'surf',
                                     f'{hemi}h.sphere.reg'),
                        shared.brainnetome_cortex_gcs_paths[index],
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'label',
                                     f'{hemi}h.BN_Atlas.annot')],
                        env=env)
                    run.command([
                        'mri_label2vol',
                        '--annot',
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'label',
                                     f'{hemi}h.BN_Atlas.annot'),
                        '--temp',
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'mri',
                                     'brain.mgz'),
                        '--o',
                        pathlib.Path(app.SCRATCH_DIR,
                                     'freesurfer',
                                     'mri',
                                     f'{hemi}h.BN_Atlas.mgz'),
                        '--subject', 'freesurfer',
                        '--hemi', f'{hemi}h',
                        '--identity',
                        '--proj', 'frac', '0', '1', '.1'],
                        env=env)
                run.command([
                    'mri_ca_label',
                    pathlib.Path(app.SCRATCH_DIR,
                                 'freesurfer',
                                 'mri',
                                 'brain.mgz'),
                    pathlib.Path(app.SCRATCH_DIR,
                                 'freesurfer',
                                 'mri',
                                 'transforms',
                                 'talairach.m3z'),
                    shared.brainnetome_sgm_gca_path,
                    pathlib.Path(app.SCRATCH_DIR,
                                 'freesurfer',
                                 'mri',
                                 'BN_Atlas_subcortex.mgz')],
                    env=env)
                parc_image_path = parc_image_path / 'aparc.BN_Atlas+aseg.mgz'
                # Need to deal with prospect of overlapping mask labels
                # - Any overlap between the two hemisphere ribbons
                #   = set to zero
                # - Any overlap between cortex and sub-cortical
                #   = retain cortex
                run.command(['mrcalc',
                            [pathlib.Path('freesurfer',
                                          'mri',
                                          f'{hemi}h.BN_Atlas.mgz')
                             for hemi in ['l', 'r']],
                            '-mult',
                            'cortex_overlap.mif',
                            '-datatype', 'bit'])
                run.command(['mrcalc',
                            [pathlib.Path('freesurfer',
                                          'mri',
                                          f'{hemi}h.BN_Atlas.mgz')
                             for hemi in ['l', 'r']],
                            '-add',
                            pathlib.Path('freesurfer',
                                         'mri',
                                         'BN_Atlas_subcortex.mgz'),
                            '-mult',
                            'sgm_overlap.mif',
                            '-datatype', 'bit'])
                run.command(['mrcalc',
                            [pathlib.Path('freesurfer',
                                          'mri',
                                          f'{hemi}h.BN_Atlas.mgz')
                             for hemi in ['l', 'r']],
                            '-add',
                            '1.0',
                            'cortex_overlap.mif',
                            '-sub',
                            '-mult',
                            pathlib.Path('freesurfer',
                                         'mri',
                                         'BN_Atlas_subcortex.mgz'),
                            '1.0',
                            'sgm_overlap.mif',
                            '-sub',
                            '-mult',
                            '-add',
                            parc_image_path])
                app.cleanup('cortex_overlap.mif')
                app.cleanup('sgm_overlap.mif')

            elif shared.parcellation == 'hcpmmp1':
                parc_image_path = parc_image_path / 'aparc.HCPMMP1+aseg.mgz'
                for index, hemi in enumerate(['l', 'r']):
                    run.command(['mri_surf2surf',
                                 '--srcsubject', 'fsaverage',
                                 '--trgsubject', 'freesurfer',
                                 '--hemi', f'{hemi}h',
                                 '--sval-annot',
                                 shared.hcpmmp1_annot_paths[index],
                                 '--tval',
                                 pathlib.Path(app.SCRATCH_DIR,
                                              'freesurfer',
                                              'label',
                                              f'{hemi}h.HCPMMP1.annot')],
                                env=env)
                run.command(['mri_aparc2aseg',
                             '--s', 'freesurfer',
                             '--old-ribbon',
                             '--annot', 'HCPMMP1',
                             '--o', parc_image_path],
                            env=env)
            elif shared.parcellation in ['yeo7fs', 'yeo17fs']:
                num = '7' if shared.parcellation == 'yeo7fs' else '17'
                parc_image_path = parc_image_path / f'aparc.Yeo{num}+aseg.mgz'
                for index, hemi in enumerate(['l', 'r']):
                    run.command(['mri_surf2surf',
                                 '--srcsubject', 'fsaverage5',
                                 '--trgsubject', 'freesurfer',
                                 '--hemi', f'{hemi}h',
                                 '--sval-annot',
                                 shared.yeo_annot_paths[index],
                                 '--tval',
                                 pathlib.Path(app.SCRATCH_DIR,
                                              'freesurfer',
                                              'label',
                                              f'{hemi}h.Yeo{num}.annot')],
                                env=env)
                run.command(['mri_aparc2aseg',
                             '--s', 'freesurfer',
                             '--old-ribbon',
                             '--annot', f'Yeo{num}',
                             '--o', parc_image_path],
                            env=env)
            else:
                assert False

        if shared.mrtrix_lut_file:
            # If necessary:
            # Perform the index conversion
            run.command(['labelconvert',
                         parc_image_path,
                         shared.parc_lut_file,
                         shared.mrtrix_lut_file,
                         'parc_init.mif'])
            # Substitute the sub-cortical grey matter parcellations
            #   with estimates from FSL FIRST
            run.command(['labelsgmfix',
                         'parc_init.mif',
                         freesurfer_t1w_input,
                         shared.mrtrix_lut_file,
                         'parc.mif'])
            app.cleanup('parc_init.mif')
        else:
            # Non-standard sub-cortical parcellation;
            #   labelsgmfix not applicable
            run.command(f'mrconvert {parc_image_path} parc.mif '
                        '-datatype uint32')
        app.cleanup('freesurfer')


    elif shared.do_mni:
        app.console('Registering to MNI template and transforming grey '
                    'matter parcellation back to subject space')

        # Use non-dilated brain masks for performing
        #   histogram matching & linear registration
        t1w_histmatched_path = pathlib.Path('T1w_histmatch.nii')
        histmatched_strides = '-1,+2,+3' \
                              if shared.template_registration_software == 'fsl' \
                              else '+1,+2,+3'
        run.command(f'mrhistmatch linear {t1w_image} {shared.template_image_path}'
                    ' -mask_input T1w_mask.mif'
                    f' -mask_target {shared.template_mask_path} - | '
                    f'mrconvert - {t1w_histmatched_path}'
                    f' -strides {histmatched_strides}')

        assert shared.template_registration_software
        if shared.template_registration_software == 'ants':

            # Use ANTs SyN for registration to template
            # From Klein and Avants, Frontiers in Neuroinformatics 2013:
            ants_prefix = 'template_to_t1_'
            run.command([
                'antsRegistration',
                '--dimensionality', '3',
                '--output', ants_prefix,
                '--use-histogram-matching', '1',
                '--initial-moving-transform',
                f'[{t1w_histmatched_path},{shared.template_image_path},1]',
                '--transform', 'Rigid[0.1]',
                '--metric',
                f'MI[{t1w_histmatched_path},{shared.template_image_path},1,32,Regular,0.25]',
                '--convergence', '1000x500x250x100',
                '--smoothing-sigmas', '3x2x1x0',
                '--shrink-factors', '8x4x2x1',
                '--transform', 'Affine[0.1]',
                '--metric',
                f'MI[{t1w_histmatched_path},{shared.template_image_path},1,32,Regular,0.25]',
                '--convergence', '1000x500x250x100',
                '--smoothing-sigmas', '3x2x1x0',
                '--shrink-factors', '8x4x2x1',
                '--transform', 'BSplineSyN[0.1,26,0,3]',
                '--metric',
                f'CC[{t1w_histmatched_path},{shared.template_image_path},1,4]',
                '--convergence', '100x70x50x20',
                '--smoothing-sigmas', '3x2x1x0',
                '--shrink-factors', '6x4x2x1'])
            transformed_atlas_path = 'atlas_transformed.nii'
            run.command([
                'antsApplyTransforms',
                '--dimensionality', '3',
                '--input', shared.parc_image_path,
                '--reference-image', t1w_histmatched_path,
                '--output', transformed_atlas_path,
                '--n', 'GenericLabel',
                '--transform', f'{ants_prefix}1Warp.nii.gz',
                '--transform', f'{ants_prefix}0GenericAffine.mat',
                '--default-value', '0'])
            app.cleanup(glob.glob(f'{ants_prefix}*'))

        elif shared.template_registration_software == 'fsl':

            # Subject T1w, brain masked; for flirt -in
            flirt_in_path = pathlib.Path(t1w_histmatched_path)
            if not t1w_is_premasked:
                while flirt_in_path.suffix:
                    flirt_in_path = flirt_in_path.with_suffix('')
                flirt_in_path = f'{flirt_in_path}_masked.nii'
                run.command(f'mrcalc {t1w_histmatched_path} T1w_mask.mif -mult {flirt_in_path}')
            # Template T1w, brain masked; for flirt -ref
            flirt_ref_path = pathlib.Path('template_masked.nii')
            run.command(f'mrcalc {shared.template_image_path} {shared.template_mask_path}'
                        ' -mult - | '
                        f'mrconvert - {flirt_ref_path} -strides -1,+2,+3')
            # Now have data required to run flirt
            run.command([shared.flirt_cmd,
                         '-ref', flirt_ref_path,
                         '-in', flirt_in_path,
                         '-omat', 'T1w_to_template.mat',
                         '-dof', '12',
                         '-cost', 'leastsq'])
            if not t1w_is_premasked:
                app.cleanup(flirt_in_path)
            app.cleanup(flirt_ref_path)

            # If possible, use dilated brain masks for non-linear
            #   registration to mitigate mask edge effects;
            #   if T1-weighted image is premasked, can't do this
            fnirt_in_path = t1w_histmatched_path
            fnirt_ref_path = shared.template_image_path
            if t1w_is_premasked:
                fnirt_in_mask_path = pathlib.Path('T1w_mask.nii')
                run.command(f'mrconvert T1w_mask.mif {fnirt_in_mask_path} -strides -1,+2,+3')
                fnirt_ref_mask_path = shared.template_mask_path
            else:
                fnirt_in_mask_path = pathlib.Path('T1w_mask_dilated.nii')
                run.command('maskfilter T1w_mask.mif dilate - -npass 3 | '
                            f'mrconvert - {fnirt_in_mask_path} -strides -1,+2,+3')
                fnirt_ref_mask_path = pathlib.Path('template_mask_dilated.nii')
                run.command(f'maskfilter {shared.template_mask_path} dilate {fnirt_ref_mask_path}'
                            ' -npass 3')

            run.command([
                shared.fnirt_cmd,
                f'--config={shared.fnirt_config_basename}',
                f'--ref={fnirt_ref_path}',
                f'--in={fnirt_in_path}',
                '--aff=T1w_to_template.mat',
                f'--refmask={fnirt_ref_mask_path}',
                f'--inmask={fnirt_in_mask_path}',
                '--cout=T1w_to_template_warpcoef.nii'])
            app.cleanup(fnirt_in_mask_path)
            if not t1w_is_premasked:
                app.cleanup(fnirt_ref_mask_path)
            app.cleanup('T1w_to_template.mat')
            fnirt_warp_subject2template_path = \
                fsl.find_image('T1w_to_template_warpcoef')

            # Use result of registration to transform atlas
            #   parcellation to subject space
            run.command([
                shared.invwarp_cmd,
                f'--ref={t1w_histmatched_path}',
                f'--warp={fnirt_warp_subject2template_path}',
                '--out=template_to_T1w_warpcoef.nii'])
            app.cleanup(fnirt_warp_subject2template_path)
            fnirt_warp_template2subject_path = \
                fsl.find_image('template_to_T1w_warpcoef')
            run.command([
                shared.applywarp_cmd,
                f'--ref={t1w_histmatched_path}',
                f'--in={shared.parc_image_path}',
                f'--warp={fnirt_warp_template2subject_path}',
                '--out=atlas_transformed.nii',
                '--interp=nn'])
            app.cleanup(fnirt_warp_template2subject_path)
            transformed_atlas_path = fsl.find_image('atlas_transformed')

        app.cleanup(t1w_histmatched_path)

        if shared.parc_lut_file and shared.mrtrix_lut_file:
            run.command(['labelconvert',
                         transformed_atlas_path,
                         shared.parc_lut_file,
                         shared.mrtrix_lut_file,
                         'parc.mif'])
        else:
            # Not all parcellations need to go through the labelconvert step;
            #   they may already be numbered incrementally from 1
            run.command(f'mrconvert {transformed_atlas_path} parc.mif')
        app.cleanup(transformed_atlas_path)


    if output_verbosity > 2 and shared.parcellation != 'none':
        if shared.mrtrix_lut_file:
            label2colour_lut_option = ['-lut', shared.mrtrix_lut_file]
        elif shared.parc_lut_file:
            label2colour_lut_option = ['-lut', shared.parc_lut_file]
        else:
            # Use random colouring if no LUT available, but
            #   still generate the image
            label2colour_lut_option = []
        run.command(['label2colour', 'parc.mif', 'parcRGB.mif']
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
        tractogram_filepath = pathlib.Path(f'tractogram_{num_streamlines}.tck')
        run.command(f'tckgen FOD_WM.mif {tractogram_filepath}'
                    ' -act 5TT.mif -backtrack -crop_at_gmwmi'
                    ' -maxlength 250'
                    ' -power 0.33'
                    f' -select {num_streamlines}'
                    ' -seed_dynamic FOD_WM.mif')

        # Step 6: Use SIFT2 to determine streamline weights
        app.console('Running the SIFT2 algorithm to assign '
                    'weights to individual streamlines')
        fd_scale_gm_option = []
        if not multishell:
            fd_scale_gm_option = ['-fd_scale_gm']
        # If SIFT2 fails, reduce number of streamlines and try again
        while num_streamlines:
            try:
                run.command(['tcksift2',
                             tractogram_filepath,
                             'FOD_WM.mif',
                             'weights.csv',
                             '-act', '5TT.mif',
                             '-out_mu', 'mu.txt']
                            + fd_scale_gm_option)
                break
            except run.MRtrixCmdError:
                app.warn('SIFT2 failed, likely due to running out of RAM; '
                         'reducing number of streamlines and trying again')
                num_streamlines = int(num_streamlines // 2)
                new_tractogram_filepath = \
                    pathlib.Path(f'tractogram_{num_streamlines}.tck')
                run.command(f'tckedit {tractogram_filepath} {new_tractogram_filepath}'
                            f' -number {num_streamlines}')
                app.cleanup(tractogram_filepath)
                tractogram_filepath = new_tractogram_filepath
        if not num_streamlines:
            raise MRtrixError('Unable to run SIFT2 algorithm for '
                              'any number of streamlines')
        run.function(shutil.move, tractogram_filepath, 'tractogram.tck')
        tractogram_filepath = pathlib.Path('tractogram.tck')


        if output_verbosity > 2:
            # Generate TDIs:
            # - A TDI at DWI native resolution, with SIFT mu scaling,
            #   and precise mapping
            #     (for comparison to WM ODF l=0 term, to
            #     verify that SIFT2 has worked correctly)
            app.console('Producing Track Density Images (TDIs)')
            # TODO This may have a header
            with open('mu.txt', 'r', encoding='utf-8') as f:
                mu = float(f.read())
            # In the space of the DWI image
            run.command(f'tckmap {tractogram_filepath} -'
                        ' -tck_weights_in weights.csv'
                        ' -template FOD_WM.mif'
                        ' -precise'
                        ' | '
                        f'mrcalc - {mu} -mult tdi_dwi.mif')
            # In the space of the T1-weighted image
            run.command(f'tckmap {tractogram_filepath} -'
                        ' -tck_weights_in weights.csv'
                        f' -template {t1w_image}'
                        ' -precise'
                        ' | '
                        f'mrcalc - {mu} -mult tdi_T1w.mif')
            # - Conventional TDI at super-resolution
            #   (mostly just because we can)
            run.command(f'tckmap {tractogram_filepath} tdi_hires.mif'
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
            ['-assignment_radial_search', '5'] \
            if shared.parcellation in ['yeo7mni', 'yeo17mni'] \
            else []
        run.command(['tck2connectome', tractogram_filepath, 'parc.mif', 'connectome.csv',
                     '-tck_weights_in', 'weights.csv',
                     '-out_assignments', 'assignments.csv']
                    + assignment_option)
        run.command(['tck2connectome', tractogram_filepath, 'parc.mif', 'meanlength.csv',
                     '-tck_weights_in', 'weights.csv',
                     '-scale_length',
                     '-stat_edge', 'mean']
                    + assignment_option)

        if output_verbosity > 2:
            # Produce additional data that can be used for
            #   visualisation within mrview's connectome toolbar
            app.console('Generating geometric data for '
                        'enhanced connectome visualisation')
            run.command(f'connectome2tck {tractogram_filepath} assignments.csv exemplars.tck'
                        ' -tck_weights_in weights.csv'
                        ' -exemplars parc.mif'
                        ' -files single')
            run.command('label2mesh parc.mif nodes.obj')
            run.command('meshfilter nodes.obj smooth nodes_smooth.obj')
            app.cleanup('nodes.obj')



    # Prepare output path for writing
    app.console(f'Processing for session "{session_label}" completed;'
                ' writing results to output directory')
    subdirs_to_make = ['anat', 'dwi', 'tractogram']
    if shared.parcellation != 'none':
        subdirs_to_make.insert(1, 'connectome')
    for subdir in subdirs_to_make:
        full_subdir_path = output_subdir / subdir
        if full_subdir_path.exists():
            run.function(shutil.rmtree, full_subdir_path)
        run.function(os.makedirs, full_subdir_path)

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
            output_analysis_level_path / \
            f'{parc_string[1:]}_lookup{lut_export_file.suffix}'
        try:
            shutil.copy(lut_export_file, lut_export_path)
        except OSError:
            pass

    # Copy / convert necessary files to output directory
    for scratch_file, output_item in output_items.items():
        if output_verbosity >= output_item.min_verbosity \
                and (multishell or not output_item.needs_multishell):
            full_output_path = output_subdir / output_item.path
            if output_item.is_image:
                run.command(f'mrconvert {scratch_file} {full_output_path}'
                            ' -clear_property comments'
                            + (' ' + output_item.options
                               if output_item.options
                               else ''),
                            force=full_output_path.exists())
            else:
                run.function(shutil.copyfile,
                             scratch_file,
                             full_output_path)

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
