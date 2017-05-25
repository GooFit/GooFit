#!/usr/bin/env python

from __future__ import print_function

from plumbum import local, cli
import re

expr = re.compile(r"#include")

def add_namespace(lines, ns):
    n = 0;
    for v,line in enumerate(lines):
        # Quit if already contains namespace
        if 'namespace ' + ns in line:
            return None

        # Add namespace after last include
        if ('#include' in line
         or '#pragma once' in line):
            n = v
    lines.insert(n+1, '\nnamespace '+ns+' {\n\n')
    lines.append(r"} // namespace " + ns + "\n\n")
    return lines

class AddNamespace(cli.Application):
    namespace = cli.SwitchAttr('--name', str, default='GooFit', help='namespace to add')

    @cli.positional(cli.ExistingFile)
    def main(self, *filenames):
        for fi in filenames:
            with open(fi) as f:
                lines = f.readlines()
            lines = add_namespace(lines, self.namespace)
            if lines is None:
                print("Skipping", fi)
                continue
            with open(fi, 'w') as f:
                f.write(''.join(lines))
            print("Added namespace to", fi)

if __name__ == '__main__':
    AddNamespace()
