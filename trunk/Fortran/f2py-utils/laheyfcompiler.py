import os
import sys

from cpuinfo import cpu
from fcompiler import FCompiler

class LaheyFCompiler(FCompiler):

    compiler_type = 'lahey'
    version_pattern =  r'Lahey/Fujitsu Fortran 95 Compiler Release (?P<version>[^\s*]*)'

    executables = {
        'version_cmd'  : ["lf95", "-version"],
        'compiler_f77' : ["lf95", "-fix"],
        'compiler_fix' : ["lf95", "-fix"],
        'compiler_f90' : ["lf95"],
        'linker_so'    : ["lf95"],
        'linker_exe'   : ["lf95"],
        'archiver'     : ["ar", "-cr"],
        'ranlib'       : ["ranlib"]
        }

    module_dir_switch = '-mod '  #XXX Fix me
    module_include_switch = '' #XXX Fix me
    obj_extension = ".obj"
    shared_lib_extension = ".obj"
    static_lib_extension = ".lib"
    def get_flags(self):
        return ['-ZERO']
    def get_flags_opt(self):
        return [
                '-O1',
                '-nap',
                '-nchk',
                '-nchkglobal',
                '-ncover',
                '-ng',
                '-npca',
                '-nsav',
                '-nstchk',
                '-ntrace'
                ]
    def get_flags_debug(self):
        return [
                '-g',
                '-chk(a,e,s,u,x)',
                '-chkglobal',
                '-pca',
                '-stchk',
                '-trace',
                '-w',
                '-info'
                ]
    def get_library_dirs(self):
        opt = []
        d = os.environ.get('LAHEY')
        if d:
            opt.append(os.path.join(d,'lib'))
        return opt
    def get_libraries(self):
        opt = []
        #CF -- Changed libs for Lahey95
        opt.extend([
            'f90gl',
            'f90glu',
            'f90glut',
            'f90SQL',
            'fj90dcw',
            'fj90f',
            'fj90g',
            'fj90i',
            'fj90vi',
            'fj90vsw',
            'fjfdcw',
            'fjfi',
            'fjfvi',
            'fjfvsw',
            'lf95w',
            'ncos',
            'SSL2',
            'winter',
            'user32',
            'libcmodified',
            ])
        #

        return opt

if __name__ == '__main__':
    from distutils import log
    log.set_verbosity(2)
    from fcompiler import new_fcompiler
    #compiler = new_fcompiler(compiler='lahey')
    compiler = LAHEYFCompiler()
    compiler.customize()
    print compiler.get_version()
