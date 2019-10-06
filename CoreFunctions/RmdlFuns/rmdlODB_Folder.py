# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 19:38:15 2019

@author: Jean-Baptiste
"""
# -*- coding: mbcs -*-
from abaqus import *
from abaqusConstants import *
import __main__
import odbAccess
import os
import numpy as np
import time
import subprocess
import math
import re


os.chdir(r'C:\Users\Doctorant\Documents\JB\Rmdl_Macro\TRO_L')
# cwd = os.getcwd() + '\\'
#-----------------------------------------------------------------
# Parameters
matKeyword = 'MAT_'
elSetKeyword = 'SET_'
partName = 'TIBIA-1'
nCPUs = 8
rmdl_T = 0.40
dt = 1./2.
law = 'Morgan2003' # Rho to E law used in the models
#-----------------------------------------------------------------
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import odbAccess
import rmdlFun as rmdl_funs
import creaPostOpFun
import creaPreOpFun



# Get the models that have been assigned material
Folder = os.getcwd()
FileList = os.listdir(Folder)
mdlNames_WM = [f for f in FileList if f.startswith('Tib_') and f.endswith('_WM.inp') and 'cut' not in f ]

for mdlName in mdlNames_WM :
    # Read the original model before implantation
    mdlName = mdlName[0:-4]
    #mdlName.replace('.', '_')
    mdlName_short = '_'.join(mdlName.split('_')[1:4])
    mdlName_NoOp = mdlName_short.replace('.', '_') + '_NoOp'

    # Check if PreOp analysis has already been executed, launched it otherwise
    if not os.path.isfile(mdlName_NoOp+'.odb'):
        if not os.path.isfile(mdlName_NoOp+'.inp'):
            creaPreOpFun.creaPreOpMdl(mdlName,nCPUs)
            print('PreOp Inp file of' + mdlName + 'has been created.')
        rmdl_funs.launch_inp(mdlName_NoOp, nCPUs)
        print('PreOp analysis of' + mdlName + 'has been launched.')
        # Wait for job to be completed -----
        rmdl_funs.check_analysis_completed(mdlName_NoOp)
    print('PreOp Simulation of ' + mdlName_NoOp + ' done...')

    preOp = session.openOdb(name= mdlName_NoOp + '.odb', readOnly=FALSE)
    TibRA_preOp = preOp.rootAssembly.instances['TIBIA-1']
    TibRA_preOp_ES_All = TibRA_preOp.elementSets['ES-2RMDL']

    # Read inp of initial file 
    inpText_NoOp = rmdl_funs.get_inp_text(mdlName_NoOp)

    #Find Mat Ppties with regex
    Dict_Mat_E = rmdl_funs.get_Mat_From_inp(inpText_NoOp,matKeyword)

    #FindElset with regex
    dict_Elset = rmdl_funs.get_ES_From_inp(inpText_NoOp,elSetKeyword)
    ElementList = sorted([ el for k in dict_Elset.keys() for el in dict_Elset[k]])

    # FindSection : Assocition between ElSet and Materials
    dict_ELset_E, dict_Mat_Elset = rmdl_funs.get_ES2MAT_From_inp(inpText_NoOp,matKeyword,elSetKeyword,Dict_Mat_E)
    dict_E_Elset = { dict_ELset_E[elset] : elset for elset in dict_ELset_E} # Invert dict

    # Get elements elastic modulus correspondance : 
    dict_El_E = dict()
    for k in dict_ELset_E.keys():
        for el in dict_Elset[k] :
            dict_El_E[el] = dict_ELset_E[k]

            
    # Get referentiel strain energy 'massic' density
    Dict_Sref = {el.label : [] for el in TibRA_preOp_ES_All.elements}
    Data_Sref = []
    elmtData= []
    Data_rho = []
    for stepName in preOp.steps.keys() :
        preOp_SED = preOp.steps[stepName].frames[-1].fieldOutputs['ESEDEN']
        SED_preOp =  preOp_SED.getSubset(region=TibRA_preOp_ES_All)
        if 'Sref' not in preOp.steps[stepName].frames[-1].fieldOutputs.keys() :
            SEMD_preOp = preOp.steps[stepName].frames[-1].FieldOutput(name='Sref',
                        description='Strain Energy Massic Density', type=SCALAR)
        else :
            SEMD_preOp = preOp.steps[stepName].frames[-1].fieldOutputs['Sref']

        if 'rho' not in preOp.steps[stepName].frames[-1].fieldOutputs.keys() :
            RHO_preOp = preOp.steps[stepName].frames[-1].FieldOutput(name='rho',
                        description='Apparent Density of TB', type=SCALAR)
        else :
            RHO_preOp = preOp.steps[stepName].frames[-1].fieldOutputs['rho']
        
        for i, val in enumerate(SED_preOp.values) :
            E = dict_El_E[val.elementLabel]
            rho = rmdl_funs.rho_from_E(E,law)
            S = val.data/rho
            Data_Sref.append((S,))
            elmtData.append(val.elementLabel)
            Data_rho.append((rho,))
            Dict_Sref[val.elementLabel].append(S)
            
        SEMD_preOp.addData(position=WHOLE_ELEMENT, instance=TibRA_preOp,
                                        labels=elmtData, data=Data_Sref)
        RHO_preOp.addData(position=WHOLE_ELEMENT, instance=TibRA_preOp,
                                        labels=elmtData, data=Data_rho)

    epoch = 0

    Sref_Equi = { el : np.mean(s) for el, s in Dict_Sref.items() }

    while epoch < 1 :
        name_curr = mdlName_short.replace('.', '_') +'_Op_'+ str(epoch)
        # Check if current epoch post Op analysis has already been executed, launched it otherwise
        if not os.path.isfile(name_curr+'.odb'):
            if not os.path.isfile(name_curr+'.inp') and epoch == 0:
                creaPostOpFun.creaPostOpMdl(mdlName,nCPUs)
                print('PostOp Inp file of' + name_curr + 'has been created.')
            rmdl_funs.launch_inp(name_curr, nCPUs)
            print('PostOp analysis of' + name_curr + ' has been launched.')
            # Wait for job to be completed -----
            rmdl_funs.check_analysis_completed(name_curr)
        print('Simulation of ' + mdlName + ' epoch ' + str(epoch) + ' done...')

        #-----------------------------------------------------------------
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #               Remodelling of the tibia bulk volume
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #-----------------------------------------------------------------

        # -----------------------------------------------------------
        # Dictionnary of elements and associated elastic modulus
        #------------------------------------------------------------
        # READ original INP file
        inpText = rmdl_funs.get_inp_text(name_curr)

        #Find Mat Ppties with regex
        Dict_Mat_E = rmdl_funs.get_Mat_From_inp(inpText,matKeyword)
        
        #FindElset with regex
        dict_Elset = rmdl_funs.get_ES_From_inp(inpText,elSetKeyword)
        ElementList = sorted([ el for k in dict_Elset.keys() for el in dict_Elset[k]])

        # FindSection : Assocition between ElSet and Materials
        dict_ELset_E, dict_Mat_Elset = rmdl_funs.get_ES2MAT_From_inp(inpText,matKeyword,elSetKeyword,Dict_Mat_E)
        dict_E_Elset = { dict_ELset_E[elset] : elset for elset in dict_ELset_E} # Invert dict

        # Get elements elastic modulus correspondance : 
        dict_El_E = dict()
        for k in dict_ELset_E.keys():
            for el in dict_Elset[k] :
                dict_El_E[el] = dict_ELset_E[k]


        # -----------------------------------------------------------
        # Get the values of the signal S
        #------------------------------------------------------------
        postOp = session.openOdb(name= name_curr + '.odb', readOnly=FALSE)
        TibRA = postOp.rootAssembly.instances['TIBIA-1']
        TibRA_ES_ALL = TibRA.elementSets['ES-2RMDL']
        TibRA_ES_TBCMT = TibRA.elementSets['ES-LAYER-ALL']


        # Get SEMD for each steps
        Dict_S = {el.label : [] for el in TibRA_ES_ALL.elements}
        for stepName in postOp.steps.keys() :
            preOp_SED = preOp.steps[stepName].frames[-1].fieldOutputs['ESEDEN']
            SED_preOp =  preOp_SED.getSubset(region=TibRA_preOp_ES_All)
            preOp_Sref = preOp.steps[stepName].frames[-1].fieldOutputs['Sref']
            SREF =  preOp_Sref.getSubset(region=TibRA_preOp_ES_All)
            #
            postOp_SED = postOp.steps[stepName].frames[-1].fieldOutputs['ESEDEN']
            SED =  postOp_SED.getSubset(region=TibRA_ES_ALL)
            #
            #postOp_Sref = postOp.steps[stepName].frames[-1].fieldOutputs['Sref']
            #Sref =  postOp_Sref.getSubset(region=TibRA_ES_ALL)
            
            # Compute Strain Energy Density Shielding field for the current model relative to the préop situation
            SEMD = postOp.steps[stepName].frames[-1].FieldOutput(name='SEMD',
                                description='Strain Energy Massic Density', type=SCALAR)
            SEDSField = postOp.steps[stepName].frames[-1].FieldOutput(name='SEDSh',
                                description='SED Shielding', type=SCALAR)
            SEMDsh = postOp.steps[stepName].frames[-1].FieldOutput(name='SEMDsh',
                                description='Strain Energy Massic Density shielding', type=SCALAR)
            RHO = postOp.steps[stepName].frames[-1].FieldOutput(name='rho',
                                description='Apparent Density of TB', type=SCALAR)
            E_MOD = postOp.steps[stepName].frames[-1].FieldOutput(name='E_Mod',
                                description='Elastic Modulus of material', type=SCALAR)
            Data_S = []
            Data_deltaS = []
            Data_SEDsh = []
            Data_Rho = []
            Data_E_Mod = []
            elmtData = []
            for i, val in enumerate(SED.values) :
                E = dict_El_E[val.elementLabel]
                rho = rmdl_funs.rho_from_E(E,law)
                Data_Rho.append((rho,))
                Data_E_Mod.append((E,))
                
                S = val.data/rho
                Data_S.append((S,))
                Dict_S[val.elementLabel].append(S)
                
                sedSh = [val.data - SED_preOp.values[i].data]
                Data_SEDsh.append(tuple(sedSh))
                
                deltaS = [S - SREF.values[i].data]
                Data_deltaS.append(tuple(deltaS))
                       
                elmtData.append(val.elementLabel)
                
            # Write new fields
            SEMD.addData(position=WHOLE_ELEMENT, instance=TibRA,
                labels=elmtData, data=Data_S)
            SEDSField.addData(position=WHOLE_ELEMENT, instance=TibRA,
                labels=elmtData, data=Data_SEDsh)
            SEMDsh.addData(position=WHOLE_ELEMENT, instance=TibRA,
                labels=elmtData, data=Data_deltaS)
            RHO.addData(position=WHOLE_ELEMENT, instance=TibRA, labels=elmtData, data=Data_Rho)
            E_MOD.addData(position=WHOLE_ELEMENT, instance=TibRA, labels=elmtData, data=Data_E_Mod)
        # Get the equivalent strain enery density 1st strategy : Mean diff
        S_Equi = { el : np.mean(s) for el, s in Dict_S.items() }
        

        # -----------------------------------------------------------
        # Get the new value of the elastic modulus
        #------------------------------------------------------------
        dict_El_E_Rmdl = dict()
        dict_Elset_Rmdl = {k: [] for k in dict_Elset.keys()}

        for el in dict_El_E:
            dict_El_E_Rmdl[el] = rmdl_funs.bone_remodeling(dict_El_E[el], S_Equi[el], Sref_Equi[el], dt, law)
            E_closest = rmdl_funs.find_nearest_E_group(sorted(dict_E_Elset.keys()),dict_El_E_Rmdl[el])
            dict_Elset_Rmdl[dict_E_Elset[E_closest]].append(el)
                
        #-----------------------------------------------------------------
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Remodelling of the Trabecular bone interdifitated in the cement
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #-----------------------------------------------------------------
        # Get SvM for each steps
        TBCMT_El_SvM = {el.label : [] for el in TibRA_ES_TBCMT.elements}
        for stepName in postOp.steps.keys() :
            postOp_SvM = postOp.steps[stepName].frames[-1].fieldOutputs['MISESONLY']
            SvM =  postOp_SvM.getSubset(region=TibRA_ES_TBCMT)
            for i, val in enumerate(SvM.values) :
                TBCMT_El_SvM[val.elementLabel].append(val.data)
        # Get the equivalent von Mises Stress 1st strategy : Mean diff
        TBCMT_El_MaxSvM = { el : np.mean(SvMList) for el, SvMList in TBCMT_El_SvM.items() }
        TBCMT_El_SvM = TBCMT_El_MaxSvM
        
        # -----------------------------------------------------------
        # Dictionnary of elements and associated elastic modulus
        #-----------------------------------------------------------------
        # READ original INP file
        inpText = rmdl_funs.get_inp_text(name_curr)
        
        #Find Mat Ppties with regex
        TBCMT_Mat_E = rmdl_funs.get_Mat_From_inp(inpText,'TB-PMMA')

        #FindElset with regex
        TBCMT_Elset = rmdl_funs.get_ES_From_inp(inpText,'SECT_TB-PMMA')

        # FindSection : Assocition between ElSet and Materials
        TBCMT_ELset_E, TBCMT_Mat_Elset = rmdl_funs.get_ES2MAT_From_inp(inpText,'TB-PMMA','SECT_TB-PMMA',TBCMT_Mat_E)
        TBCMT_E_Elset = { TBCMT_ELset_E[elset] : elset for elset in TBCMT_ELset_E} # Invert dict

        # Get elements elastic modulus correspondance :
        TBCMT_El_E = dict()
        for k in TBCMT_ELset_E.keys():
            for el in TBCMT_Elset[k] :
                TBCMT_El_E[el] = TBCMT_ELset_E[k]

        # Get knew value of elastic modulus
        TBCMT_El_E_Rmdl = dict()
        TBCMT_Elset_Rmdl = {k: [] for k in TBCMT_Elset.keys()}

        
        for el in TBCMT_El_E:
            TBCMT_El_E_Rmdl[el] = rmdl_funs.tbcmt_remodeling(TBCMT_El_E[el],TBCMT_El_SvM[el],dt)
            E_closest = rmdl_funs.find_nearest_E_group(sorted(TBCMT_E_Elset.keys()),TBCMT_El_E_Rmdl[el])
            TBCMT_Elset_Rmdl[TBCMT_E_Elset[E_closest]].append(el)



        for stepName in postOp.steps.keys() :
            Data_E_Mod = []
            elmtData = []
            if 'E_Mod' not in postOp.steps[stepName].frames[-1].fieldOutputs.keys() :
                E_Mod = postOp.steps[stepName].frames[-1].FieldOutput(name='E_Mod',
                                description='Elastic Modulus of material', type=SCALAR)
            else :
                E_Mod = postOp.steps[stepName].frames[-1].fieldOutputs['E_Mod']
                
            for el, E in TBCMT_El_E.iteritems(): 
                Data_E_Mod.append((E,))
                elmtData.append(el)
                
            E_MOD.addData(position=WHOLE_ELEMENT, instance=TibRA, labels=elmtData, data=Data_E_Mod)
        
        #-----------------------------------------------------------------
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #               Write updated Inp file with new elset
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #-----------------------------------------------------------------
        pElset2 = re.compile('elset=('+elSetKeyword+'[0-9]+)',re.IGNORECASE)
        pElset3 = re.compile('elset=('+'SECT_TB-PMMA'+'[0-9]+)',re.IGNORECASE)
        name_rmdl = mdlName_short.replace('.', '_') +'_Op_'+ str(epoch+1)                
        fout = open(name_rmdl + '.inp','w')
        with open(name_curr + '.inp', 'r') as f:
            writeLine = True
            for line in f.readlines() :
                if '*' in line :
                    writeLine = True
                if line.startswith('*Elset'):
                    elset = pElset2.findall(line,re.I)
                    elset3 = pElset3.findall(line,re.I)
                    fout.write(line)
                    if elset :
                        writeLine = False
                        elmts = dict_Elset_Rmdl[elset[0].upper()]
                        if elmts :
                            #break elmts list
                            elmts = [str(el) for el in elmts]
                            elmts_chunk = [elmts[i:i + 10] for i in xrange(0, len(elmts), 10)]
                            for c in elmts_chunk:
                                fout.write(', '.join(c)+'\n')
                        else :
                            fout.write('  \n')
                    elif elset3 :
                        writeLine = False
                        elmts = TBCMT_Elset_Rmdl[elset3[0].upper()]
                        if elmts :
                            #break elmts list
                            elmts = [str(el) for el in elmts]
                            elmts_chunk = [elmts[i:i + 10] for i in xrange(0, len(elmts), 10)]
                            for c in elmts_chunk:
                                fout.write(', '.join(c)+'\n')
                        else :
                            fout.write('  \n')
                elif writeLine :
                    fout.write(line)
                
                    
        fout.close()
        postOp.save()
        postOp.close()
        epoch += 1

    preOp.save()
    preOp.close()
