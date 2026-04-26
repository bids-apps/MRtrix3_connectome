import os
from mrtrix3 import MRtrixError
from mrtrix3 import app

# Examine the contents of a directory (whether the raw BIDS dataset or a
#   derivatives directory), and return a list of sessions.
# This list may optionally be filtered based on the use of batch processing
#   command-line options; e.g. resticting the participant or session IDs.
def get_sessions(root_dir, **kwargs):

    participant_labels = kwargs.pop('participant_label', None)
    session_labels = kwargs.pop('session_label', None)
    if kwargs:
        raise TypeError(f'Unsupported keyword arguments passed to get_session(): {kwargs}')

    # Perform a recursive search through the BIDS dataset directory,
    #   looking for anything that resembles a BIDS session
    # For any sub-directory that itself contains directories "anat/" and "dwi/",
    #   store the list of sub-directories required to navigate to that point
    # This becomes the list of feasible processing targets for any level
    #   of analysis
    # From there:
    #   - "--participant_label" can be used to remove entries from the list
    #   - "--session_label" can be used to remove entries from the list
    all_sessions = []
    for dir_name, subdir_list, _ in os.walk(root_dir):
        subdir_list[:] = [entry \
                          for entry in subdir_list \
                          if entry.lower() != 'derivatives']
        if all(item in subdir_list for item in ('anat', 'dwi')):
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

    invalid_sessions = []
    for session in all_sessions:
        session = os.path.normpath(session).split(os.sep)
        if all(any(subdir.startswith(prefix) for prefix in ['sub-', 'ses-'])
               for subdir in session):
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
        else:
            invalid_sessions.append(os.sep.join(session))

    if invalid_sessions:
        app.warn(f'Entr{"ies" if len(invalid_sessions) > 1 else "y"}'
                 f' in "{root_dir}"'
                 ' found with valid anat/ and dwi/ sub-directories,'
                 f'but invalid directory name{"s" if len(invalid_sessions) > 1 else ""}:'
                 f' {invalid_sessions}')

    if not result:
        raise MRtrixError('No sessions were selected for processing')

    app.console(f'{len(all_sessions)} total sessions found'
                f' in directory "{root_dir}";'
                f' {"all" if len(result) == len(all_sessions) else str(len(result))}'
                ' will be processed')
    sub_not_found = [key for key, value in sub_found.items() if not value]
    if sub_not_found:
        app.warn(f'{len(sub_not_found)} nominated participant'
                 f' label{"s were" if len(sub_not_found) > 1 else " was"}'
                 ' not found in input dataset:'
                 f' {", ".join(sub_not_found)}')
    ses_not_found = [key for key, value in ses_found.items() if not value]
    if ses_not_found:
        app.warn(f'{len(ses_not_found)} nominated session'
                 f' label{"s were" if len(ses_not_found) > 1 else " was"}'
                 ' not found in input dataset:'
                 f' {", ".join(ses_not_found)}')

    app.debug(str(result))
    return result
