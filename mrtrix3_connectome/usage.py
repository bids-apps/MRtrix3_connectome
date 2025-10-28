from mrtrix3 import app

from . import IS_CONTAINER
from . import OPTION_PREFIX
from . import __version__

COPYRIGHT = '''Copyright (c) 2016-2025 The Florey Institute of Neuroscience
and Mental Health.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Covered Software is provided under this License on an "as is"
basis, without warranty of any kind, either expressed, implied, or
statutory, including, without limitation, warranties that the
Covered Software is free of defects, merchantable, fit for a
particular purpose or non-infringing.
See the Mozilla Public License v. 2.0 for more details.'''

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

    cmdline.set_copyright(COPYRIGHT)

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
            type=app.Parser.Int(0),
            metavar='number',
            dest='nthreads',
            help='Use this number of threads in MRtrix3 '
                 'multi-threaded applications '
                 '(0 disables multi-threading)')
        cmdline._action_groups[2].add_argument(
            '-scratch', '--scratch',
            dest='scratch',
            type=app.Parser.DirectoryOut(),
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
        'DWI: denoising (if applicable); '
        'Gibbs ringing removal (if applicable); '
        'motion, eddy current and EPI distortion correction '
        'and outlier detection & replacement; '
        'gradient non-linearity distortion correction '
        '(if available & necessary); '
        'brain masking, bias field correction and intensity normalisation; '
        'rigid-body registration & transformation to T1-weighted image (if available). '
        'T1-weighted image: '
        'gradient non-linearity distortion correction '
        '(if available & necessary); '
        'bias field correction; '
        'brain masking.')
    cmdline.add_description(
        '"participant": '
        'DWI: Response function estimation; FOD estimation. '
        f'T1-weighted image (if {OPTION_PREFIX}parcellation '
        'is not none): '
        'Tissue segmentation; grey matter parcellation. '
        f'Combined (if {OPTION_PREFIX}parcellation is not none, '
        f'or {OPTION_PREFIX}streamlines is provided): '
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
        f'{OPTION_PREFIX}participant_label and '
        f'{OPTION_PREFIX}session_label options '
        'correspond(s) to sub-<participant_label> and '
        'ses-<session_label> from the BIDS spec (so they do _not_ '
        'include "sub-" or "ses-"). Multiple participants / sessions '
        'can be specified with a space-separated list.')
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
        f'The {OPTION_PREFIX}output_verbosity option principally '
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
    cmdline.add_description(
        'As the production of individual connectomes is split into two '
        'stages, being preproc-level and participant-level analyses, '
        'there is scope for erroneous usage to lead to sub-optimal results. '
        'In particular, it is possible for registration between DWI and '
        'T1-weighted images to be either absent from the pre-processing '
        'step, or for different T1-weighted images to be utilised '
        'between the two stages of analysis, resulting in misalignment '
        'between DWI and anatomical information. It is up to the user '
        'to ensure that the expectation of the participant-level analyses '
        'of such alignment having already been achieved is not violated.'
    )

    cmdline.add_argument(
        'bids_dir',
        type=app.Parser.DirectoryIn(),
        help='The directory with the input dataset formatted '
             'according to the BIDS standard.')
    cmdline.add_argument(
        'output_dir',
        type=app.Parser.DirectoryIn(),
        help='The existing directory where BIDS Derivatives will be written')
    cmdline.add_argument(
        'analysis_level',
        help='Level of analysis that will be performed; '
             f'options are: {", ".join(ANALYSIS_CHOICES)}.',
        choices=ANALYSIS_CHOICES)

    cmdline.add_argument(
        f'{OPTION_PREFIX}output_verbosity',
        type=app.Parser.Int(1, 4),
        default=2,
        help='The verbosity of script output (number from 1 to 4).')

    batch_options = cmdline.add_argument_group(
        'Options specific to the batch processing of participant data')
    batch_options.add_argument(
        f'{OPTION_PREFIX}participant_label',
        nargs='+',
        help='The label(s) of the participant(s) that should be analyzed.')
    batch_options.add_argument(
        f'{OPTION_PREFIX}session_label',
        nargs='+',
        help='The session(s) within each participant that should be analyzed.')

    preproc_options = \
        cmdline.add_argument_group(
            'Options that are relevant to preproc-level analysis')
    preproc_options.add_argument(
        f'{OPTION_PREFIX}gdc',
        metavar='path',
        type=app.Parser.DirectoryIn(),
        help='Provide a directory containing pre-computed warps '
             'for gradient non-linearity distortion correction')
    preproc_options.add_argument(
        f'{OPTION_PREFIX}concat_denoise',
        choices=('before', 'after'),
        default='before',
        help='Specify whether, '
             'in the presence of multiple DWI series for a session, '
             'one prefers concatenation of those files into a single series '
             'before vs. after denoising '
             '(not guaranteed to be used depending on data)')

    preproc_participant_options = \
        cmdline.add_argument_group(
            'Options that are relevant to both preproc-level and '
            'participant-level analyses')
    preproc_participant_options.add_argument(
        f'{OPTION_PREFIX}t1w_preproc',
        metavar='path',
        type=app.Parser.DirectoryIn(),
        help='Provide a path by which pre-processed T1-weighted image data '
             'may be found for the processed participant(s) / session(s)')

    participant_options = \
        cmdline.add_argument_group(
            'Options that are relevant to participant-level analysis')
    if not IS_CONTAINER:
        participant_options.add_argument(
            f'{OPTION_PREFIX}atlas_path',
            metavar='path',
            type=app.Parser.DirectoryIn(),
            help='The filesystem path in which to search for atlas '
                 'parcellation files.')
    participant_options.add_argument(
        f'{OPTION_PREFIX}parcellation',
        help='The choice of connectome parcellation scheme '
             '(compulsory for participant-level analysis); '
             f'options are: {", ".join(PARCELLATION_CHOICES)}.',
        choices=PARCELLATION_CHOICES)
    participant_options.add_argument(
        f'{OPTION_PREFIX}streamlines',
        type=app.Parser.Int(0),
        default=0,
        help='The number of streamlines to generate for each subject '
             '(will be determined heuristically if not explicitly set).')
    participant_options.add_argument(
        f'{OPTION_PREFIX}template_reg',
        metavar='software',
        help='The choice of registration software for mapping subject to '
             'template space; '
             f'options are: {", ".join(REGISTRATION_CHOICES)}.',
        choices=REGISTRATION_CHOICES)
