#!/usr/bin/env python
import warnings
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import (get_info, AtlasNotFoundError,
                                         BlasNotFoundError)

def configuration(parent_package='', top_path=None):
    blas_info = get_info('atlas')
    if not blas_info:
        warnings.warn(AtlasNotFoundError.__doc__)
        blas_info = get_info('blas')
        if not blas_info:
            # Blas is required
            print "\nError:\n%s\n" % BlasNotFoundError.__doc__
            raise SystemExit(1)

    config = Configuration('bvp1lg', parent_package, top_path)
    config.add_extension('_colnew',
                         sources=['../../lib/colnew.pyf',
                                  '../../lib/colnew.f',
                                  '../../lib/dgesl.f',
                                  '../../lib/dgefa.f'],
                         libraries=blas_info['libraries'],
                         library_dirs=blas_info['library_dirs'])
    config.add_extension('_mus',
                         sources=['../../lib/mus.pyf',
                                  '../../lib/mus1.f',
                                  '../../lib/mus2.f',
                                  '../../lib/mus3.f',
                                  '../../lib/dgesl.f',
                                  '../../lib/dgefa.f'],
                         libraries=blas_info['libraries'],
                         library_dirs=blas_info['library_dirs'])

    config.add_data_dir('tests')

    return config
