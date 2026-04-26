import os
IS_CONTAINER = os.path.exists('/version') \
               and os.path.exists('/mrtrix3_version')
OPTION_PREFIX = '--' if IS_CONTAINER else '-'

# pylint: disable=consider-using-with
__version__ = 'BIDS-App \'MRtrix3_connectome\'' \
              f'version {open("/version", "r", encoding="utf-8").read()}' \
              if IS_CONTAINER \
              else 'BIDS-App \'MRtrix3_connectome\' standalone'
