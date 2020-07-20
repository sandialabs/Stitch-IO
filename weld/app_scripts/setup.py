import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    from distutils.sysconfig import get_python_inc
    config=Configuration('app_scripts',parent_package,top_path)
    config.add_include_dirs(get_python_inc())

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

      
