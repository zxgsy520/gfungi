#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys

cov_file = sys.argv[1]
kat_mx_file = sys.argv[2]
name = sys.argv[3]
kat_path = sys.argv[4]


ngs_depth = 200

for line in open(cov_file, 'r'):
    if line.startswith('Coverage'):
        ngs_depth = line.strip().split()[-1]
        if int(float(ngs_depth)) > 200:
            ngs_depth = int(float(ngs_depth)) + 20
        else:
            ngs_depth = 200

os.system('%s plot spectra-cn -x %s -o %s_kat %s' % (kat_path, ngs_depth, name, kat_mx_file))
