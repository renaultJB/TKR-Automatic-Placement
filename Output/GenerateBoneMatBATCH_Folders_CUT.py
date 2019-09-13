# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 10:00:17 2017

@author: Jean-Baptiste

Generate the Bonemat Batch file for all the inp files in the
current Directory subfolders
"""


import os

cwd = os.getcwd()


INPs_WM = cwd + '\\INPs_WM'
#if os.path.exists(INPs_WM):
#    os.removedirs(INPs_WM)

FolderList = [ fldr for fldr in next(os.walk('.'))[1] if fldr != 'INPs_WM' and 'GUI_R' in fldr]

if not os.path.exists(INPs_WM):
    os.makedirs(INPs_WM)


fout = open('BoneMatBatch_ALL.txt','w')
for Folder in FolderList:
    # =========================================================================
    # list files in folder
    # =========================================================================
    FileList = os.listdir(Folder)
    
    # =========================================================================
    # Identify the inp files in the directory
    # =========================================================================
    OutputFileID = [f for f in FileList if f[0:6] in ['Tibia_','TIBIA_','tibia_']  and f[-3:] in ['inp','INP'] and '_cut_' in f]
    OutputFileID.sort()
    
    # =========================================================================
    # Get the subject code (without and with the side [L or R])
    # =========================================================================
    SubjectCode = OutputFileID[0].split('_')[1]
    SubjectCodeLR = '_'.join(OutputFileID[0].split('_')[1:3])
    
    # =========================================================================
    # Create the bonemat batch output file
    # =========================================================================
    
    for fileID in OutputFileID :
        
        SubjectCodeLRAlpha = 'Tib_' + SubjectCodeLR + '_' + '.'.join(fileID.split('_')[-1].split('.')[0:2])
        SubjectCodeLRAlpha = ''.join(SubjectCodeLRAlpha.split('.'))
        
        fout.write(cwd+'\\'+ Folder +'\\'+fileID+'\n')
        fout.write(cwd+'\\'+ Folder +'\\'+ SubjectCode +'\\'+ SubjectCode +'.vtk\n')
        fout.write(cwd+'\\'+ 'JBMAFfile.xml\n')
        fout.write(cwd+'\\'+ Folder +'\\'+SubjectCodeLRAlpha+'_cut_WM.inp'+'\n\n')

fout.close()
