#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 17:12:53 2016

@author: thorsten
"""

import os
import glob

DOC_DIR  = os.path.dirname(__file__)
ROOT_DIR = os.path.join(DOC_DIR, '..')

def main():
    in_path = os.path.join(ROOT_DIR, 'Tutorials')

    fns = glob.glob(os.path.join(in_path, '*.py'))

    for fn in fns:
        bn = os.path.basename(fn)
        out_fn = os.path.join(DOC_DIR, 'Tutorials', '__' + bn.replace('.py', '.txt'))

        in_code_block   = False
        in_ignore_block = False
        out_fh = open(out_fn, 'w')
        for line in open(fn, 'r'):
            if in_ignore_block==False and line.startswith('"""'):
                in_ignore_block = True
                in_code_block   = False
                continue
            elif in_ignore_block==True and line.startswith('"""'):
                in_ignore_block = False
                in_code_block   = False
                continue
            elif in_ignore_block==True:
                in_code_block   = False
                continue
            elif line.startswith('# -*-'):
                continue
            elif not line.startswith('##'):
                if not in_code_block:
                    if len(line.strip())==0:
                        continue
                    out_fh.write('\n.. code-block:: python\n\n')
                    in_code_block = True
                out_fh.write('    ' + line)
            elif line.startswith('###'):
                if in_code_block:
                    out_fh.write('\n')
                in_code_block   = False
                line = line.replace('#','').strip()
                out_fh.write('**' + line + '**\n\n')
#                out_fh.write('"'*len(line) + '\n')
            elif line.startswith('##'):
                if in_code_block:
                    out_fh.write('\n')
                in_code_block   = False
                out_fh.write(line.replace('#','').strip() + '\n')
        out_fh.close()

if __name__ == '__main__':
    main()
