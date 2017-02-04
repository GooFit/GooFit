#!/usr/bin/env python

from __future__ import print_function, division

from plumbum import local, cli, FG
from plumbum.cmd import curl

FILES = [
    'https://raw.githubusercontent.com/henryiii/CLI11/master/include/CLI.hpp',
]

DIR = local.path(__file__).dirname

def download_file(path):
    name = path.split('/')[-1]
    (curl[path] > name) & FG

class UpdateCLI(cli.Application):
    def main(self):
        with local.cwd(DIR / '../include/goofit/detail'):
            for f in FILES:
                download_file(f)

if __name__ == "__main__":
    UpdateCLI()
