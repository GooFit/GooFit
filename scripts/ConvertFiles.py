#!/usr/bin/env python

from __future__ import print_function
import re
try:
    from plumbum import local, cli
except ImportError:
    print("This file uses the plumbum library. Install with pip or conda (user directory or virtual environment OK).")
    raise

FILE = local.path(__file__).dirname

base = FILE / '../include'
files = base // '*/*.h'
files += base // 'goofit/*/*.h'
rel = [f - base for f in files]

def fix_text(contents):
    """
    >>> text = '''
    ... #include "Variable.h"
    ... #include "GlobalCudaDefines.h"
    ... #include "goofit/FitControl.h"
    ... #include <set>
    ... #include <goofit/BinnedDataSet.h>
    ... #include <UnbinnedDataSet.h>
    ... '''

    >>> corrected = fix_text(text)
      Converting GlobalCudaDefines.h -> goofit/GlobalCudaDefines.h
      Converting UnbinnedDataSet.h -> goofit/UnbinnedDataSet.h
      Converting Variable.h -> goofit/Variable.h
    >>> print(corrected.strip())
    #include "goofit/Variable.h"
    #include "goofit/GlobalCudaDefines.h"
    #include "goofit/FitControl.h"
    #include <set>
    #include <goofit/BinnedDataSet.h>
    #include "goofit/UnbinnedDataSet.h"
    """

    for r in rel:
        before = r.parts[-1]
        after = str(r)
        re_before = re.compile(r'["<]'+before+'h?[">]')
        if re_before.search(contents):
            print('  Converting {0} -> {1}'.format(before, after))
            contents = re_before.sub('"'+after+'"',contents)
    return contents

def fix_files(src):
    for name in src:
        print('Converting: {0}'.format(name))
        with name.open('r') as f:
            contents = f.read()
        contents = fix_text(contents)
        with name.open('w') as f:
            f.write(contents)



class ConvertFiles(cli.Application):

    @cli.positional(cli.ExistingFile)
    def main(self, *src):
        fix_files(src)


if __name__ == '__main__':
    ConvertFiles()

