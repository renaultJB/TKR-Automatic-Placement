# -*- coding: utf-8 -*-
"""
Created on Sun Oct 6 12:29:49 2019

@author: Jean-Baptiste Renault
"""
# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import os

os.chdir(r'C:\Users\Doctorant\Documents\GitHub\TKR-Automatic-Placement\Output')


import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import optimization
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import odbAccess
import numpy as np
import AbqsMdlFuns as mdl_funs

#Find the subject folders
ParentFolder = os.getcwd()
FolderList = os.listdir(ParentFolder)
SubjectFolder = [fldr for fldr in FolderList if len(fldr)==5  and '_' in fldr and '.' not in fldr]

#Find the odb output folders
OutputInpsFolder = os.listdir(ParentFolder + '\\OutputInps')
ODBFilesAll = [f for f in OutputInpsFolder if f.endswith('.odb') ]

ODBFilesOK = [f for f in ODBFilesAll if f+'_f' not in OutputInpsFolder]
ODBFilesNotOK = ['_'.join(f.split('_')[0:3]) for f in ODBFilesAll if f+'_f' in OutputInpsFolder]

#Write a file with all the mis-executed models
NotOkODBFile = open("NotOkODBFile.txt",'w')
NotOkODBFile.write("\n".join(ODBFilesNotOK))
NotOkODBFile.close() 


# =============================================================================
# Open output file and write first line (Header)
# =============================================================================
with open('OutputLoadCentrality.txt','w') as fCScore :
    fCScore.write('{},{},{},{},{},{},'.format('subjectCode','LegSide','ImplantSize','Alignement','alpha','beta'))
    fCScore.write('{},{},{},{},{},{},'.format('CV','Min.Max','Min.Mean','mMPTA','TibialSlopeDia','TibialSlopeMech'))
    fCScore.write('{},{},{},{},{},{},{},{},'.format('Lsup','Lflank','Ltip','Ltot','Lshare','LTipCv','LTipStd','LTipMax'))
    fCScore.write('{},{},{},{},{},{},{},'.format('Mean.Distance','Wdth.Implt.ML','Wdth.Implt.AP','Wdth.Tib.ML','Wdth.Tib.AP','Ratio.AP','Ratio.ML'))
    fCScore.write('{},{}\n'.format('LoadStep','Epoch'))

    # =============================================================================
    # Read odb file and look for the 
    # =============================================================================
    for ODBf in ODBFilesOK:
        print(ODBf + ' is currently being processed...')
        subjectCodeLR = '_'.join(ODBf.split('_')[0:2])
        LegSide = subjectCodeLR.split('_')[1]
        shortMdlName = '_'.join(ODBf.split('_')[0:3])+'.'+ODBf.split('_')[3]
        mdlFolder = ParentFolder + '\\' + subjectCodeLR + '\\'
        if 'NoOp' in ODBf :
            NoOp = True
            epoch = -1 # Boolean to check if its a pre operative ODB file as OFB structure varies a little
        else :
            epoch = int(ODBf.split('_')[-1].split('.')[0])
            NoOp = False
        
        dictFile = mdlFolder + 'Dict_' + shortMdlName +'.txt'
        Data = mdl_funs.readData(dictFile)
        # Normal of the implant tibial plateau equivalent to normal of the tibial plateau cut plan
        nxp = Data['Nxp']
        
        # Get Centrality
        with open(mdlFolder + 'Output_' + shortMdlName + '.txt','r') as centralityf:
            GoodLine = False
            Infos = dict()
            ScoreName = ''
            for line in centralityf:
                if GoodLine :
                    Infos[ScoreName] = float(line)
                    GoodLine = False
                if 'Centrality ' in line:
                    GoodLine = True
                    ScoreName = 'Ctrly.' + line.split(' ')[2]
                if '(pente tibiale)' in line:
                    GoodLine = True
                    ScoreName = 'TibSlopeDia'
                if 'Coverage ' in line:
                    GoodLine = True
                    ScoreName = 'Coverage'
                if 'Rotation ' in line:
                    GoodLine = True
                    ScoreName = 'MalRot'

        
        with open(mdlFolder + 'Centrality_'+shortMdlName+'.txt','r') as centralityf :

            line = centralityf.readline()
            line = centralityf.readline()
            Cntrlty = {key: None for key in ['Xring','Yring','Zring','Xtib','Ytib','Ztib','Distance']}
            
            ringPts = []
            tibPts = []
            Distance = []
            
            while len(line.split(',')) > 5 :
                ringPts.append([float(i) for i in line.split(',')[0:3]])
                tibPts.append([float(i) for i in line.split(',')[3:6]])
                Distance.append(float(line.split(',')[6]))
                line = centralityf.readline()
            
        ringPts = np.array(ringPts)
        tibPts = np.array(tibPts)
        Distance = np.array(Distance)
        
        Width_Tib_ML = np.dot(tibPts,Data['Ymech'])
        Width_Tib_AP = np.dot(tibPts,Data['Xmech'])
        
        Width_Implt_ML = np.dot(ringPts,Data['Ymech'])
        Width_Implt_AP = np.dot(ringPts,Data['Xmech'])
        
        Infos['Mean.Distance'] = np.mean(Distance)
        Infos['Width.Implt.ML'] = np.max(Width_Implt_ML) - np.min(Width_Implt_ML)
        Infos['Width.Implt.AP'] = np.max(Width_Implt_AP) - np.min(Width_Implt_AP)
        Infos['Width.Tib.ML'] = np.max(Width_Tib_ML) - np.min(Width_Tib_ML)
        Infos['Width.Tib.AP'] = np.max(Width_Tib_AP) - np.min(Width_Tib_AP)

        Ratio_AP = [abs(np.min(Width_Tib_AP) - np.min(Width_Implt_AP)),abs(np.max(Width_Tib_AP) - np.max(Width_Implt_AP))]
        Infos['Ratio.AP'] = min(Ratio_AP) / max(Ratio_AP)

        Ratio_ML = [abs(np.min(Width_Tib_ML) - np.min(Width_Implt_ML)) , abs(np.max(Width_Tib_ML) - np.max(Width_Implt_ML))]
        Infos['Ratio.ML'] = min(Ratio_ML) / max(Ratio_ML)
            
    # =============================================================================
    #     Calculate the load share the proportion of load transiting through the stem tip
    # =============================================================================
        MyOdb = session.openOdb(name=ParentFolder + '\\OutputInps\\' + ODBf)

        lCtctNSTibia = ['NS-XP','NS-STEM-FLANK','NS-STEM-TIP','NS-CUT']
        for stepName in MyOdb.steps.keys() :
            LoadStepFrame = MyOdb.steps[stepName].frames[-1]
            
            LoadShare = dict()
            if NoOp : ################################## A verifier si le ctct n'est pas iverser pour noOP ###################################
                CnormF = LoadStepFrame.fieldOutputs['CNORMF   ASSEMBLY_TIBIA_CUT-1_SURF-CUT/ASSEMBLY_TIBIA-1_SURF-CUT']
                CshearF = LoadStepFrame.fieldOutputs['CSHEARF  ASSEMBLY_TIBIA_CUT-1_SURF-CUT/ASSEMBLY_TIBIA-1_SURF-CUT']
            else :
                CnormF = LoadStepFrame.fieldOutputs['CNORMF   ASSEMBLY_TIBIA-1_SURF-CUT/ASSEMBLY_IMPLT-1_SURF-CTCT']
                CshearF = LoadStepFrame.fieldOutputs['CSHEARF  ASSEMBLY_TIBIA-1_SURF-CUT/ASSEMBLY_IMPLT-1_SURF-CTCT']

            # Get the contact force Field
            CF = CnormF + CshearF
            AxForceShareStep = []
            for ns in lCtctNSTibia :
                TNS = MyOdb.rootAssembly.instances['TIBIA-1'].nodeSets[ns]
                CF_TNS = CF.getSubset(region=TNS)
                AxForce = [nxp[0]*val.data[0]+nxp[1]*val.data[1]+nxp[2]*val.data[2] for val in CF_TNS.values] #Axial Force are carried by the Z axis
                LoadShare[ns] = np.sum(np.array(AxForce))
                if 'NS-STEM-TIP' in ns :
                    Forces = [val.magnitude for val in CF_TNS.values]
                    LoadShare['STD_Forces_Tip'] = np.std(Forces) 
                    ImaxForces = Forces.index(np.max(Forces))
                    LoadShare['Max_Tip'] = np.max(Forces)
                    LoadShare['CV_Forces_Tip'] = LoadShare['STD_Forces_Tip']/np.mean(Forces)

            if not NoOp :
                TibBrdrNS = MyOdb.rootAssembly.instances['TIBIA-1'].nodeSets['NS-CUT-BORDER']
                ImpltBrdrNS = MyOdb.rootAssembly.instances['IMPLT-1'].nodeSets['NS-PLATE-BORDER']
            

            
        #    Write too output file
            fCScore.write('{},{},{},{},{},{},'.format(subjectCodeLR,LegSide,Data['ImpltType'],Data['Alignment'],Data['alpha'],Data['beta']))
            fCScore.write('{},{},{},{},{},{},'.format(Infos['Ctrly.CV'],Infos['Ctrly.Min_Max'],Infos['Ctrly.Min_Mean'],Data['mMPTA'],Infos['TibSlopeDia'],Data['TPslope']))
            for ns in lCtctNSTibia:
                fCScore.write('{},'.format(LoadShare[ns]))
            fCScore.write('{},'.format(LoadShare['NS-STEM-TIP']/LoadShare['NS-CUT']))
            fCScore.write('{},{},{},'.format(LoadShare['CV_Forces_Tip'],LoadShare['STD_Forces_Tip'],LoadShare['Max_Tip']))
            fCScore.write('{},{},{},{},{},{},{},'.format(Infos['Mean.Distance'],Infos['Width.Implt.ML'],Infos['Width.Implt.AP'],Infos['Width.Tib.ML'],Infos['Width.Tib.AP'],Infos['Ratio.AP'],Infos['Ratio.ML']))
            fCScore.write('{},{}\n'.format(stepName,epoch))
    #    Close ODB and go to next file
        MyOdb.close()


