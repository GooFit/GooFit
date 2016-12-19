#!/usr/bin/env python

# This is a python script to test the output of the files.
# It requires plumbum (pip install plumbum or conda install -c conda-forge plumbum)

from __future__ import print_function
from plumbum import local, cli, TEE, colors
import time

# Simple timer
class Timer:
    def __enter__(self):
        self.start = time.time()
        return self
    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

LOCAL_DIR = local.path(__file__).dirname

def test(filename):
    command = local[LOCAL_DIR / filename / filename]
    colors.info.print('Running', command)
    with local.cwd(LOCAL_DIR / filename):
        with Timer() as t:
            code, stdout, stderr = command & TEE(retcode=None)
        if code==0:
            colors.success.print(filename, 'Successful')
        else:
            colors.fatal.print(filename, 'Failed with status code:', code)
    return dict(name=filename, code=code, time=t.interval, stdout=stdout, stderr=stderr)


class TestAll(cli.Application):

    def main(self):
        filelist = [f.name for f in LOCAL_DIR // '*' if f.is_dir()]
        results = [test(fname) for fname in filelist]
        failed = [result for result in results if result['code'] !=0]
        if failed:
            colors.fatal.print("Failed:")
            for result in failed:
                colors.fatal.print(result['name'], result['code'])
        else:
            colors.success.print("All programs completed.")
        print()
        colors.info.print('Resulting times:')
        for result in results:
            print((colors.success if result['code'] == 0 else colors.warn) | '{0[name]}: {0[time]}'.format(result))

if __name__ == '__main__':
    TestAll()
