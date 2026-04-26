from .usage import usage as my_usage
from .execute import execute as my_execute

def usage(cmdline): #pylint: disable=unused-variable
    return my_usage(cmdline)

def execute(): #pylint: disable=unused-variable
    return my_execute()
