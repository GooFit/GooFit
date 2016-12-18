#!/usr/bin/env python

# This is a python script to test the output of the files.
# It requires plumbum (pip install plumbum or conda install -c conda-forge plumbum)

from __future__ import print_function
from plumbum import local, cli, TEE, colors

LOCAL_DIR = local.path(__file__).dirname

def test(filename):
    command = local[LOCAL_DIR / filename / filename]
    colors.info.print('Running', command)
    with local.cwd(LOCAL_DIR / filename):
        code, stdout, stderr = command & TEE(retcode=None)
        if code==0:
            colors.success.print(filename, 'Successful')
        else:
            colors.fatal.print(filename, 'Failed with status code:', code)
    return code


class TestAll(cli.Application):

    def main(self):
        filelist = [f.name for f in LOCAL_DIR // '*' if f.is_dir()]
        codes = [test(fname) for fname in filelist]
        failed_names = [name for code, name in zip(codes, filelist) if code !=0]
        if failed_names:
            colors.fatal.print("Failed:", failed_names)
        else:
            colors.success.print("All programs completed.")

if __name__ == '__main__':
    TestAll()
