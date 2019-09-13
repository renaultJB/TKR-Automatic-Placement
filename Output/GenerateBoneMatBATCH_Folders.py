# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 08:40:04 2017

@author: Jean-Baptiste RENAULT

Generate the Bonemat Batch file for all the inp files in the
current Directory subfolders
"""

import os

cwd = os.getcwd()


INPs_WM = cwd + '\\INPs_WM'

FolderList = [ fldr for fldr in next(os.walk('.'))[1] if fldr != 'INPs_WM' and 'DES_R' in fldr]

if not os.path.exists(INPs_WM):
    os.makedirs(INPs_WM)


fout = open('BoneMatBatch_ALL.txt','w')
for Folder in FolderList:
    # =============================================================================
    # list files in folder
    # =============================================================================
    FileList = os.listdir(Folder)
    
    # =============================================================================
    # Identify the inp files in the directory
    # =============================================================================
    OutputFileID = [f for f in FileList if f[0:6] in ['Tibia_','TIBIA_','tibia_']  and f[-3:] in ['inp','INP'] and '_cut_' not in f]
    OutputFileID.sort()
    
    # =============================================================================
    # Get the subject code (without and with the side [L or R])
    # =============================================================================
    SubjectCode = OutputFileID[0].split('_')[1]
    SubjectCodeLR = '_'.join(OutputFileID[0].split('_')[1:3])
    
    # =============================================================================
    # Create the bonemat batch output file
    # =============================================================================
    
    for fileID in OutputFileID :
        
        SubjectCodeLRAlpha = 'Tib_' + SubjectCodeLR + '_' + '.'.join(fileID.split('_')[-1].split('.')[0:2])
        #SubjectCodeLRAlpha = '_'.join(SubjectCodeLRAlpha.split('.'))
        
        fout.write(cwd+'\\'+ Folder +'\\'+fileID+'\n')
        fout.write(cwd+'\\'+ Folder +'\\'+ SubjectCode +'\\'+ SubjectCode +'.vtk\n')
        fout.write(cwd+'\\'+ 'JBMAFfile.xml\n')
        fout.write(INPs_WM +'\\'+SubjectCodeLRAlpha+'_WM.inp'+'\n\n')

fout.close()
