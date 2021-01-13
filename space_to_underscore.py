# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 12:33:11 2021

@author: David
"""

import glob, os

for filename in glob.glob('*.csv'):
    new_name = filename.replace(' ', '_')
    os.rename(filename, new_name)