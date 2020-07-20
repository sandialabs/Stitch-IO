import os
from os.path import join

def get_libraries_list(lapack_info):
    try:
        libraries=list(lapack_info['libraries'])
    except KeyError:
        libraries=list()
    return libraries

def get_library_dirs_list(lapack_info):
    try:
        library_dirs=list(lapack_info['library_dirs'])
    except KeyError:
        library_dirs=list()
    return library_dirs

def get_include_dirs_list(lapack_info):
    try:
        include_dirs=list(lapack_info['include_dirs'])
    except KeyError:
        include_dirs=list()
    return include_dirs

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info
    from distutils.sysconfig import get_python_inc
    from numpy.distutils.misc_util import Configuration
    config=Configuration('potts',parent_package,top_path)
    config.add_include_dirs(get_python_inc())

    lapack = get_info('lapack_opt')
    blas = get_info('blas_opt')
    include_dirs=get_include_dirs_list(lapack)
    library_dirs=get_library_dirs_list(lapack); 
    libraries=get_libraries_list(lapack); 

    # package
    weld_geometry_sources=['_weld_geometrymodule.cxx']

    # package
    config.add_extension('_weld_geometry',
                         sources=[join('',x) for x in weld_geometry_sources],
                         include_dirs=include_dirs,
                         libraries=libraries,
                         library_dirs=library_dirs,
                         extra_compile_args=['-std=c++11'],
                         extra_link_args=None)
    # package
    teardrop_sources=['_teardropmodule.cxx']

    # package
    config.add_extension('_teardrop',
                         sources=[join('',x) for x in teardrop_sources],
                         include_dirs=include_dirs,
                         libraries=libraries,
                         library_dirs=library_dirs,
                         extra_compile_args=['-std=c++11'],
                         extra_link_args=None)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
