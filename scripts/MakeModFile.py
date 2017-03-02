#!/usr/bin/env python

from __future__ import print_function

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

class MakeModFile(cli.Application):

    def main(self):
        print("conversion = [")
        for r in rel:
            before = r.parts[-1]
            after = '"' + str(r) + '"'
            re_before = (r'["<]'+before+'h?[">]')
            print(repr((re_before, after)), ',', sep='')
        print("]")


if __name__ == '__main__':
    MakeModFile.run()
