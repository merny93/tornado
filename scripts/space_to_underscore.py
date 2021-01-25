# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 12:33:11 2021

@author: David
"""

import glob, os

# White spaces to underscore for one file
# for filename in glob.glob('../first_week_raw_data/data/Angles/55/00 Other Sources/Ba-133 000.Chn'):
#     new_name = filename.replace(' ', '_')
#     os.rename(filename, new_name)
    
def replace(parent):
    '''
    Convert white spaces to underscores in all the file names contained in the path provided.
    '''
    for path, folders, files in os.walk(parent):
        for f in files:
            os.rename(os.path.join(path, f), os.path.join(path, f.replace(' ', '_')))
        for i in range(len(folders)):
            new_name = folders[i].replace(' ', '_')
            os.rename(os.path.join(path, folders[i]), os.path.join(path, new_name))
            folders[i] = new_name
            
replace('../data')