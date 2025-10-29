import os
import pathlib
import shutil
from mrtrix3 import CONFIG
from mrtrix3 import MRtrixError
from mrtrix3 import app
from mrtrix3 import run
from mrtrix3 import utils
from . import IS_CONTAINER
from . import OPTION_PREFIX
from .preproc.shared import Shared as PreprocShared
from .preproc.run import run_preproc
from .participant.shared import Shared as ParticipantShared
from .participant.run import run_participant
from .group.run import run_group
from .sessions import get_sessions

def execute(): #pylint: disable=unused-variable

    # TODO Should be able to use pathlib.Path.resolve()
    app.ARGS.bids_dir = pathlib.Path(os.path.abspath(app.ARGS.bids_dir))
    app.ARGS.output_dir = pathlib.Path(os.path.abspath(app.ARGS.output_dir))
    app.ARGS.t1w_preproc = pathlib.Path(os.path.abspath(app.ARGS.t1w_preproc)) \
        if app.ARGS.t1w_preproc \
        else None

    # If running within a container, and the --debug option has been
    #     provided, modify the interlly-stored MRtrix3 configuration
    #     contents, so that any temporary directories will be constructed
    #     within the mounted output directory, and therefore temporary
    #     directory contents will not be lost upon container instance
    #     destruction if the script fails at any point.
    if IS_CONTAINER and app.ARGS.debug:
        app.DO_CLEANUP = False
        if 'ScriptScratchDir' not in CONFIG and not app.ARGS.scratch:
            CONFIG['ScriptScratchDir'] = str(app.ARGS.output_dir)

    if utils.is_windows():
        raise MRtrixError(
            'Script cannot be run on Windows due to FSL dependency')

    if app.ARGS.skipbidsvalidator:
        app.console('Skipping BIDS validation based on user request')
    elif shutil.which('bids-validator'):
        run.command(['bids-validator', app.ARGS.bids_dir])
    else:
        app.warn('BIDS validator script not installed; '
                 'proceeding without validation of input data')

    if app.ARGS.output_verbosity < 1 or app.ARGS.output_verbosity > 4:
        raise MRtrixError(f'Valid values for {OPTION_PREFIX}output_verbosity'
                          ' option are from 1 to 4')

    # At output verbosity level 4 we retain all data and move the
    #   scratch directory to the output
    if app.ARGS.output_verbosity == 4:
        app.DO_CLEANUP = False

    if app.ARGS.output_dir.name.lower() == f'mrtrix3_connectome-{app.ARGS.analysis_level}':
        output_app_path = app.ARGS.output_dir.parent
    else:
        output_app_path = app.ARGS.output_dir
        if output_app_path.is_file():
            raise MRtrixError('Output path cannot be an existing file')

    sessions_to_analyze = None
    if app.ARGS.analysis_level in ['preproc', 'participant']:
        sessions_to_analyze = get_sessions(
            app.ARGS.bids_dir,
            participant_label=app.ARGS.participant_label,
            session_label=app.ARGS.session_label)

    # TODO All three stages may involve an invocation of dwi2mask;
    #   need to integrate into Shared the ability to store which algorithm is used,
    #   and ideally also a command-line option to control what is used

    if app.ARGS.analysis_level == 'preproc':

        preproc_shared = PreprocShared(app.ARGS.gdc,
                                       app.ARGS.concat_denoise,
                                       app.ARGS.eddy_cubicflm,
                                       app.ARGS.eddy_mbs)

        for session_to_process in sessions_to_analyze:
            app.console(f'Commencing execution for session: {"_".join(session_to_process)}')
            run_preproc(app.ARGS.bids_dir,
                        session_to_process,
                        preproc_shared,
                        app.ARGS.t1w_preproc,
                        app.ARGS.output_verbosity,
                        output_app_path)

    if app.ARGS.analysis_level == 'participant':

        participant_shared = \
            ParticipantShared(getattr(app.ARGS, 'atlas_path', None),
                              app.ARGS.parcellation,
                              app.ARGS.streamlines,
                              app.ARGS.template_reg)

        for session_to_process in sessions_to_analyze:
            app.console(f'Commencing execution for session: {"_".join(session_to_process)}')
            run_participant(app.ARGS.bids_dir,
                            session_to_process,
                            participant_shared,
                            app.ARGS.t1w_preproc,
                            app.ARGS.output_verbosity,
                            output_app_path)

    elif app.ARGS.analysis_level == 'group':

        if app.ARGS.participant_label:
            raise MRtrixError(f'Cannot use {OPTION_PREFIX}participant_label option '
                              'when performing group-level analysis')
        if app.ARGS.session_label:
            raise MRtrixError(f'Cannot use {OPTION_PREFIX}session_label option '
                              + 'when performing group-level analysis')

        run_group(app.ARGS.bids_dir,
                  app.ARGS.output_verbosity,
                  output_app_path)
