# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 08:56:47 2017

@author: Jean-Baptiste
"""

import os
import subprocess

# =============================================================================
# list files in folder
# =============================================================================

#wd = 'C:\Users\Jean-Baptiste\Documents\These\Outils\Briques Codes\Generate Mesh from Step'
#DirList = os.listdir(wd)
DirList = os.listdir(os.getcwd())

# =============================================================================
# Identify correct files
# =============================================================================
CorrectFileID = [f for f in DirList if f[0:6] in ['Tibia_','TIBIA_','tibia_'] and f[-3:] in ['tep','stp'] and '_cut_' in f] #
CorrectFileID.sort()

# =============================================================================
# Get the files containing the altitude of the the implant
# =============================================================================
OutputFileID = [f for f in DirList if f[0:7] in ['Output_','OUTPUT_','output_'] and f[-3:] in ['txt','Txt']]
OutputFileID.sort()

# =============================================================================
# Generate abaqus inp files
# =============================================================================
i=0
for fileID in CorrectFileID:
    shortID = '.'.join(fileID.split('.')[0:-1])+'_CUT
    
    f = open('TibiaCutted_MESH_3D.geo','r')
    filedata = f.read()
    f.close()
    
    f = open(OutputFileID[i],'r')
    line = f.readline()
    while 'ZALT' not in line :
        line = f.readline()
    line = f.readline()
    Z_Alt = float(line)
    f.close()
    i+=1
    
    newdata = filedata.replace("Tibia.step" , fileID)
    newdata = newdata.replace("ZALT" , str(Z_Alt))
    finaldata = newdata.replace("Tibia.inp" , shortID + ".inp")
    
    
    f = open(shortID + '_inp.geo','w')
    f.write(finaldata)
    f.close()
    
    subprocess.call(["gmsh", shortID + '_inp.geo'], shell=True)
    os.remove(shortID + '_inp.geo')
    print("Mesh generated")

# =============================================================================
# print("Script that tries to delete 1D and 2D elements of current folder Abaqus input files")
# print("Only works with Tetrahedral Elements")
# =============================================================================
#DELETE Non Volumetric elements of Abaqus Inp File
# !!! Works only for tetrahedral 3D Elements

cwdf = os.listdir(os.getcwd())
inp_files = []
for files in cwdf :
    if files[len(files)-4:len(files)]== '.inp':
        print("Found .inp files : "+files)
        output=''
        with open(files, "r") as f:
            line = f.readline()
            while line[0:4] not in ['*ELE','*ELS','*SUR','*VOL']:
                output = output + line
                line = f.readline()
                #print(line)
            #output = output + line  
            while line[0:18] != '*ELEMENT, type=C3D':
                line = f.readline()
            output = output + line
            line = f.readline()
            output = output + line
            while line[0:4] not in ['*ELE','*ELS','*SUR','*VOL'] and line:
                line = f.readline()
                output = output + line
        with open(files, "w") as f:
            print(files+" ---> done reading.")
            f.write(output)
            print(files+" ---> done writing.")
