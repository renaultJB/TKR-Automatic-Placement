# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 08:40:04 2017

@author: Jean-Baptiste RENAULT

Generate the Bonemat Batch file for all the inp files in the
current Directory subfolders
"""
import time
import os


HU_E_relationship = ['Carter1977', 'Morgan2003_Snyder1991']

print('Select HU to Rho to Elastic modulus relationship:')
time.sleep(0.5)
for i, relationship in enumerate(HU_E_relationship) :
    print(str(i+1) + ' --> ' + relationship)
    time.sleep(0.5)
    
i = input("Enter the number corresponding to the selected relationship: ")
selected_HU_E_relationship = HU_E_relationship[i-1]
print(selected_HU_E_relationship + ' selected...')
          
cwd = os.getcwd()

INPs_WM = cwd + '\\INPs_WM'

FolderList = [ fldr for fldr in next(os.walk('.'))[1] if fldr != 'INPs_WM' ] #and 'DES_R' in fldr

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
    OutputFileID = [f for f in FileList if f[0:6] in ['Tibia_','TIBIA_','tibia_']  and f[-3:] in ['inp','INP'] ] #and '_cut_' not in f
    OutputFileID.sort()

    
    # =============================================================================
    # Get the subject code (without and with the side [L or R])
    # =============================================================================
    if OutputFileID :
        SubjectCode = OutputFileID[0].split('_')[1]
        SubjectCodeLR = '_'.join(OutputFileID[0].split('_')[1:3])
        
    
    # =============================================================================
    # Create the bonemat batch output file
    # =============================================================================
    
    for fileID in OutputFileID :
        Cut = '_cut' if '_cut_' in fileID else ''
        
        SubjectCodeLRAlpha = 'Tib_' + SubjectCodeLR + '_' + '.'.join(fileID.split('_')[-1].split('.')[0:2])
        #SubjectCodeLRAlpha = '_'.join(SubjectCodeLRAlpha.split('.'))
        
        fout.write(cwd+'\\'+ Folder +'\\'+fileID+'\n')
        fout.write(cwd+'\\'+ Folder +'\\'+ SubjectCode +'\\'+ SubjectCode +'.vtk\n')
        fout.write(cwd+'\\'+ 'MAFfile_'+selected_HU_E_relationship+'.xml\n')
        fout.write(INPs_WM +'\\'+SubjectCodeLRAlpha+Cut+'_WM.inp'+'\n\n')


fout.close()
time.sleep(0.5)
print('Batch file created... closing ...')
time.sleep(2.5)
