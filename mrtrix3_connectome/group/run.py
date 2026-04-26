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
from .. import IS_CONTAINER
from ..sessions import get_sessions

GROUP_BRAINMASKS_DIR = pathlib.Path('brainmasks')
GROUP_BZEROS_DIR = pathlib.Path('bzeros')
GROUP_CONNECTOMES_DIR = pathlib.Path('connectomes')
GROUP_FA_DIR = pathlib.Path('fa')
GROUP_RESPONSES_DIR = pathlib.Path('responses')
GROUP_WMVOXELS_DIR = pathlib.Path('wmvoxels')
GROUP_WARPS_DIR = pathlib.Path('warps')

NORMALISATION_JSON_DATA = {
    'rf': {
        'Description': 'Multiplication term based on the difference in '
                       'magnitude between the white matter response '
                       'function used during independent participant-'
                       'level analysis, and the group average white '
                       'matter response function generated during group-'
                       'level analysis'
    },
    'mu': {
        'Description': 'Value of the "proportionality coefficient" '
                       'within the SIFT model',
        'Units': 'FD/mm'
    },
    'vol': {
        'Description': 'Volume of DWI voxels',
        'Units': 'mm^3'
    },
    'wmb0': {
        'Description': 'Multiplication term based on the median '
                       'intensity of the b=0 image within white matter, '
                       'compared to the mean of this value across '
                       'subjects'
    },
    'norm': {
        'Description': 'Normalisation factor applied to session '
                       'connectome data, calculated as the product of '
                       'the "rf", "mu", "vol" and "wmb0" terms'
    }
}

def run_group(bids_dir, output_verbosity, output_app_dir):

    preproc_dir = output_app_dir / 'MRtrix3_connectome-preproc'
    participant_dir = output_app_dir / 'MRtrix3_connectome-participant'
    group_dir = output_app_dir / 'MRtrix3_connectome-group'

    # Participant-level analysis no longer generates FA and mean b=0 images
    # These really should not be that expensive to compute in series,
    #   and will keep the output directory cleaner

    # Check presence of all required input files before proceeding
    # Pre-calculate paths of all files since many will be used in
    #   more than one location
    class SessionPaths:
        def __init__(self, session):
            session_label = '_'.join(session)
            preproc_root = pathlib.Path(preproc_dir, *session)
            participant_root = pathlib.Path(participant_dir, *session)
            group_root = pathlib.Path(group_dir, *session)
            # Get input DWI path here rather than in function
            in_dwi_image_list = (preproc_root / 'dwi').glob('*_dwi.nii*')
            if not in_dwi_image_list:
                raise MRtrixError(f'No DWI data found for session "{session_label}" output')
            if len(in_dwi_image_list) > 1:
                raise MRtrixError('More than one DWI mage found'
                                  f' in session "{session_label}" output')
            self.in_dwi = in_dwi_image_list[0]
            if not '_desc-preproc_' in self.in_dwi:
                raise MRtrixError(f'DWI image in output directory for session "{session_label}"'
                                  ' not flagged as pre-processed')
            in_dwi_prefix = pathlib.Path(self.in_dwi)
            while in_dwi_prefix.suffix:
                in_dwi_prefix = in_dwi_prefix.with_suffix('')
            self.in_bvec = in_dwi_prefix.with_suffix('.bvec')
            self.in_bval = in_dwi_prefix.with_suffix('.bval')
            self.in_rf = participant_root / 'dwi' / f'{session_label}_tissue-WM_response.txt'
            connectome_files = list((participant_root / 'connectome')
                                    .glob(f'{session_label}_desc-*_connectome.csv'))
            if not connectome_files:
                raise MRtrixError('No participant-level connectome file '
                                  f'found for session "{session_label}"')
            if len(connectome_files) > 1:
                raise MRtrixError('Connectomes from multiple parcellations'
                                  f' detected for session "{session_label}";'
                                  ' this is not yet supported')
            self.in_connectome = connectome_files[0]
            self.in_mu = participant_root / 'tractogram' / f'{session_label}_mu.txt'

            for entry in vars(self).values():
                if not entry.exists():
                    raise MRtrixError('Unable to find critical data '
                                      f'for session "{session_label}"'
                                      f'(expected location: {entry})')

            self.grad_import_option = ['-fslgrad', self.in_bvec, self.in_bval]

            self.bvalues = [float(value) for value in \
                    run.command(['mrinfo',
                                 self.in_dwi,
                                 '-shell_bvalues']
                                + self.grad_import_option).stdout.split()]

            self.parcellation = \
                re.findall('(?<=_desc-)[a-zA-Z0-9]*', self.in_connectome.name)[0]

            # Permissible for this to not exist at either location
            self.in_mask = preproc_root / 'dwi' / f'{session_label}_desc-brain_mask.nii.gz'
            if not self.in_mask.is_file():
                self.in_mask = participant_root / 'dwi' / f'{session_label}_desc-brain_mask.nii.gz'
                if not self.in_mask.exists():
                    self.in_mask = None

            # Not guaranteed to exist
            # Also needs to not just be a directory present, but also
            #   have the "eddy_quad" contents present (if EddyQC is not
            #   installed, that directory will still be constructed, it
            #   just will only contain contents from "eddy" itself)
            self.in_eddyqc_dir = preproc_root / 'dwi' / 'eddyqc'
            in_eddyqc_file = self.in_eddyqc_dir / 'qc.json'
            if not in_eddyqc_file.is_file():
                self.in_eddyqc_dir = None

            self.mu = matrix.load_vector(self.in_mu)[0]
            self.rf = matrix.load_matrix(self.in_rf)

            self.temp_mask = GROUP_BRAINMASKS_DIR / f'{session_label}.mif'
            self.temp_fa = GROUP_FA_DIR / f'{session_label}.mif'
            self.temp_bzero = GROUP_BZEROS_DIR / f'{session_label}.mif'
            self.temp_warp = GROUP_WARPS_DIR / f'{session_label}.mif'
            self.temp_voxels = GROUP_WMVOXELS_DIR / f'{session_label}.mif'
            self.temp_rf = GROUP_RESPONSES_DIR / f'{session_label}.txt'
            self.median_bzero = 0.0
            self.dwiintensitynorm_factor = 1.0
            self.rf_multiplier = 1.0
            self.volume_multiplier = 1.0
            self.global_multiplier = 1.0
            self.temp_connectome = GROUP_CONNECTOMES_DIR / f'{session_label}.csv'
            self.out_dir = group_root
            self.out_connectome_data = \
                group_root / 'connectome' / self.in_connectome.name
            self.out_connectome_json = \
                self.out_connectome_data.with_suffix('.json')

            self.session_label = session_label

    session_list = get_sessions(participant_dir)
    if not session_list:
        raise MRtrixError(
            'No processed session data found'
            f' in output directory "{participant_dir}"'
            ' for group analysis')
    if len(session_list) == 1:
        app.warn('Only one session present in participant directory; '
                 'some group-level analysis steps will be skipped')
    if os.path.exists(group_dir):
        app.warn('Output directory for group-level analysis already exists;'
                 ' all contents will be erased when this execution completes')


    bids_session_list = get_sessions(bids_dir)
    not_processed = [session for session in bids_session_list \
                     if session not in session_list]
    if not_processed:
        app.warn(f'{len(not_processed)} session{"s" if len(not_processed) > 1 else ""}'
                 ' present in BIDS directory'
                 f' {"have" if len(not_processed) > 1 else "has"}'
                 ' not yet undergone participant-level processing:'
                 f' {", ".join("_".join(session) for session in not_processed)}')

    sessions = []
    for session in session_list:
        sessions.append(SessionPaths(session))



    # Connectome-based calculations can only be performed if the
    #   parcellation is consistent across all sessions
    parcellation = sessions[0].parcellation
    consistent_parcellation = \
        all(s.parcellation == parcellation for s in sessions)
    out_connectome_path = group_dir / f'desc-{parcellation}_connectome.csv' \
                          if consistent_parcellation \
                          else None

    app.activate_scratch_dir()

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
            raise MRtrixError('Excessive deviation of b-values:'
                              f' mean across subjects b={shell_mean};'
                              f' range {min(shell)}-{max(shell)}')


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
        if s.in_mask is None:
            run.command(['dwi2mask',
                         'legacy',
                         s.in_dwi,
                         s.temp_mask]
                        + s.grad_import_option)
        else:
            run.command(['mrconvert',
                         s.in_mask,
                         s.temp_mask,
                         '-datatype', 'bit'])
        run.command(['dwiextract', s.in_dwi, '-bzero', '-']
                    + s.grad_import_option +
                    ['|',
                     'mrmath', '-', 'mean', s.temp_bzero,
                     '-axis', '3'])
        run.command(['dwi2tensor', s.in_dwi, '-',
                     '-mask', s.temp_mask]
                    + s.grad_import_option +
                    ['|',
                     'tensor2metric', '-',
                     '-fa', s.temp_fa,
                     '-mask', s.temp_mask])
        run.function(shutil.copy, s.in_rf, s.temp_rf)
        progress.increment()
    progress.done()

    # First group-level calculation:
    # Generate the population FA template
    if len(sessions) == 1:
        app.console('Duplicating single-subject FA image as '
                    'population template image')
        run.function(shutil.copyfile,
                     sessions[0].temp_fa,
                     'template.mif')
    else:
        app.console('Generating population template for '
                    'intensity normalisation WM mask derivation')
        run.command(['population_template',
                     GROUP_FA_DIR,
                     'template.mif',
                     '-mask_dir', GROUP_BRAINMASKS_DIR,
                     '-warp_dir', GROUP_WARPS_DIR,
                     '-type', 'rigid_affine_nonlinear',
                     '-rigid_scale', '0.25,0.5,0.8,1.0',
                     '-affine_scale', '0.7,0.8,1.0,1.0',
                     '-nl_scale', '0.5,0.75,1.0,1.0,1.0',
                     '-nl_niter', '5,5,5,5,5',
                     '-linear_no_pause'])
    app.cleanup(GROUP_FA_DIR)
    app.cleanup(GROUP_BRAINMASKS_DIR)

    # Generate the group average response function
    if len(sessions) == 1:
        app.console('Duplicating single-subject WM response function'
                    ' as group-average response function')
        run.function(shutil.copyfile,
                     sessions[0].temp_rf,
                     'response.txt')
    else:
        app.console('Calculating group-average WM response function')
        run.command(['responsemean',
                     [s.temp_rf for s in sessions],
                     'response.txt'])
    app.cleanup(GROUP_RESPONSES_DIR)
    mean_rf = matrix.load_matrix('response.txt')
    mean_rf_lzero = [line[0] for line in mean_rf]

    # Second pass through subject data in group analysis:
    #     - Warp template FA image back to subject space &
    #       threshold to define a WM mask in subject space
    #     - Calculate the median subject b=0 value within this mask
    #     - Store this in a file, and contribute to calculation of the
    #       mean of these values across subjects
    #     - Contribute to the group average response function
    if len(sessions) == 1:
        app.console('Calculating N=1 intensity normalisation factor')
        run.command('mrthreshold template.mif voxels.mif -abs 0.4')
        sessions[0].median_bzero = image.statistics(sessions[0].temp_bzero,
                                                    mask='voxels.mif').median
        app.cleanup(sessions[0].temp_bzero)
        sum_median_bzero = sessions[0].median_bzero
        app.cleanup('voxels.mif')
    else:
        progress = app.ProgressBar('Generating intensity normalisation factors',
                                   len(sessions))
        run.function(os.makedirs, GROUP_WMVOXELS_DIR)
        sum_median_bzero = 0.0
        for s in sessions:
            run.command(['mrtransform', 'template.mif', '-',
                         '-warp_full', s.temp_warp,
                         '-from', '2',
                         '-template', s.temp_bzero,
                         '|',
                         'mrthreshold', '-', s.temp_voxels,
                         '-abs', '0.4'])
            s.median_bzero = image.statistics(s.temp_bzero,
                                              mask=s.temp_voxels).median
            app.cleanup(s.temp_bzero)
            app.cleanup(s.temp_voxels)
            app.cleanup(s.temp_warp)
            sum_median_bzero += s.median_bzero
            progress.increment()
        progress.done()

    app.cleanup(GROUP_BZEROS_DIR)
    app.cleanup(GROUP_WMVOXELS_DIR)
    app.cleanup(GROUP_WARPS_DIR)
    app.cleanup('template.mif')


    # Second group-level calculation:
    # - Calculate the mean of median b=0 values
    mean_median_bzero = sum_median_bzero / len(sessions)

    # Third pass through session data in group analysis:
    # - Scaling factors for connectome strengths:
    #   - Multiply by SIFT proportionality coefficient mu
    #   - Multiply by (mean median b=0) / (subject median b=0)
    #   - Multiply by (subject RF size) / (mean RF size)
    #     (needs to account for multi-shell data)
    #   - Multiply by voxel volume
    progress = app.ProgressBar('Computing normalisation scaling factors'
                               ' for subject connectomes',
                               len(sessions))
    run.function(os.makedirs, GROUP_CONNECTOMES_DIR)
    # Determine, from the minimum connectivity value that can be represented
    #   in a streamlines-based representation, the maximum across sessions
    min_connectivity = 0.0
    for s in sessions:
        rf_lzero = [line[0] for line in s.rf]
        s.rf_multiplier = 1.0
        for (mean, subj) in zip(mean_rf_lzero, rf_lzero):
            s.rf_multiplier = s.rf_multiplier * subj / mean
        # Don't want to be scaling connectome independently for
        #   differences in RF l=0 terms across all shells;
        #   use the geometric mean of the per-shell scale factors
        s.rf_multiplier = math.pow(s.rf_multiplier, 1.0 / len(mean_rf_lzero))

        s.bzero_multiplier = mean_median_bzero / s.median_bzero

        # Calculate voxel volume
        for spacing in image.Header(s.in_dwi).spacing()[0:3]:
            s.volume_multiplier *= spacing

        s.global_multiplier = s.mu \
                              * s.bzero_multiplier \
                              * s.rf_multiplier \
                              * s.volume_multiplier

        # Minimum connectivity value that can be reasonably represented is
        #   1 streamline prior to scaling
        min_connectivity = max(min_connectivity, s.global_multiplier)

        progress.increment()
    progress.done()

    # Third group-level calculation:
    # Compute normalised connectomes, and generate the group mean connectome
    # Can only do this if the parcellation is identical across subjects;
    #     this needs to be explicitly checked
    # Use geometric mean for averaging across subjects, since variance
    #   across sessions is closer to multiplicative than additive
    if consistent_parcellation:
        progress = app.ProgressBar('Normalising subject connectomes, '
                                   'applying group-wise minimum connectivity, '
                                   'and calculating group mean connectome',
                                   len(sessions)+1)
        mean_connectome = []
        for s in sessions:
            connectome_prenorm = matrix.load_matrix(s.in_connectome)
            connectome_postnorm = [[max(v*s.global_multiplier,
                                        min_connectivity)
                                    for v in line]
                                   for line in connectome_prenorm]
            matrix.save_matrix(s.temp_connectome, connectome_postnorm)

            if mean_connectome:
                mean_connectome = [[c1+math.log(c2)
                                    for c1, c2 in zip(r1, r2)]
                                   for r1, r2 in zip(mean_connectome,
                                                     connectome_postnorm)]
            else:
                mean_connectome = [[math.log(v)
                                    for v in row]
                                   for row in connectome_postnorm]
            progress.increment()

        mean_connectome = [[math.exp(v/len(sessions))
                            for v in row]
                           for row in mean_connectome]
        progress.done()
    else:
        app.warn('Different parcellations across sessions, '
                 'cannot calculate a group mean connectome; '
                 'normalising and applying minimum connectivity '
                 'independently for each session')
        connectome_prenorm = matrix.load_matrix(s.in_connectome)
        connectome_postnorm = [[max(v, 1.0)*s.global_multiplier
                                for v in line]
                               for line in connectome_prenorm]
        matrix.save_matrix(s.temp_connectome, connectome_postnorm)


    # Run EddyQC group-level analysis if available
    # Do this LAST, as it writes back to the preproc EddyQC directories
    #   if successful
    do_squad = bool(shutil.which('eddy_squad'))
    if do_squad:
        quad_dirs = [s.in_eddyqc_dir for s in sessions if s.in_eddyqc_dir]
        missing_sessions = [s.session_label for s in sessions \
                            if not s.in_eddyqc_dir]
        if quad_dirs:
            if missing_sessions:
                app.warn('Some sessions do not contain EddyQC subject data, '
                         'and will be omitted from the group-level analysis: '
                         f'{missing_sessions}')
            run.command(['eddy_squad', quad_dirs])
        else:
            app.warn('No pre-processed sessions contain EddyQC data; '
                     '"eddy_squad" skipped')
            do_squad = False
    else:
        app.warn('EddyQC command "eddy_squad" not available; skipped')

    # Write results of interest back to the output directory;
    #     both per-subject and group information
    progress = app.ProgressBar('Writing results to output directory',
                               len(sessions)+3)
    if group_dir.exists():
        run.function(shutil.rmtree, group_dir)
    run.function(os.makedirs, group_dir)
    for s in sessions:
        run.function(os.makedirs,
                     s.out_dir / 'connectome')
        run.function(shutil.copyfile,
                     s.temp_connectome,
                     s.out_connectome_data)
        json_data = {'Contributions': {
                        'RFMagnitude': s.rf_multiplier,
                        'SIFTMu': s.mu,
                        'VoxelVolume': s.volume_multiplier,
                        'WMIntensity': s.dwiintensitynorm_factor},
                     'Multiplier': s.global_multiplier}
        with open(s.out_connectome_json, 'w', encoding='utf-8') as json_file:
            json.dump(json_data, json_file)
        progress.increment()
    app.cleanup(GROUP_CONNECTOMES_DIR)

    matrix.save_matrix(group_dir / 'tissue-WM_response.txt',
                       mean_rf,
                       force=IS_CONTAINER)
    progress.increment()
    if consistent_parcellation:
        matrix.save_matrix(out_connectome_path,
                           mean_connectome,
                           force=IS_CONTAINER)
        with open(group_dir / 'normalisation.tsv', 'w', encoding='utf-8') as tsv_file:
            tsv_file.write('session_id\trf\tmu\tvol\twmb0\tnorm\n')
            for s in sessions:
                tsv_file.write(f'{s.session_label}\t'
                               f'{s.rf_multiplier}\t'
                               f'{s.mu}\t'
                               f'{s.volume_multiplier}\t'
                               f'{s.dwiintensitynorm_factor}\t'
                               f'{s.global_multiplier}\n')
        with open(group_dir / 'normalisation.json',
                  'w',
                  encoding='utf-8') as json_file:
            json.dump(NORMALISATION_JSON_DATA, json_file)
    progress.increment()
    if do_squad:
        run.function(os.makedirs,
                     group_dir / 'eddyqc')
        for filename in ['group_db.json', 'group_qc.pdf']:
            run.function(shutil.copyfile,
                         filename,
                         group_dir / 'eddyqc' / filename)
    progress.done()

    # For group-level analysis, function is only executed once, so
    #   no need to bypass the default scratch cleanup
    # Only exception is if we want to capture the whole scratch directory
    #   in the output path
    if output_verbosity == 4:
        app.console('Copying scratch directory to output location')
        run.function(shutil.copytree,
                     app.SCRATCH_DIR,
                     group_dir / 'scratch')
