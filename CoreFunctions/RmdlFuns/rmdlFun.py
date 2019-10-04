# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 19:38:15 2019

@author: Jean-Baptiste
"""
from abaqus import *
from abaqusConstants import *
import time
import numpy as np

def check_analysis_completed(jobName):
    checkFile = jobName+'.log'
    print('looking for: ' + checkFile)
    while not os.path.exists(checkFile):
        time.sleep(10)
    print('found file: ' + checkFile)
    print('found file: ' + checkFile + ' waiting for analysis completion')
    test = -1
    while test < 0 :
        f1 = open(checkFile,'r')
        txt = f1.read()
        f1.close()
        time.sleep(10)
        test = txt.find(jobName+' COMPLETED\n')
        if 'error' in txt:
            print('Abaqus/Analysis exited with errors')
            break

    print('Job ' + jobName + ' terminated')

def find_nearest_E_group(array,value):
    array = sorted(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def exterior_elSet_from_ndSet(part, ndSetName):
    NS_ext = part.sets[ndSetName]
    els = part.elements

    if 'NS' in ndSetName.upper():
        elSetName = 'ES' + ndSetName.split('NS')[-1]
    else :
        elSetName = 'ES' + ndSetName

    list_ext_layer = [el for nd in NS_ext.nodes for el in nd.getElements()]
    list_ext_layer_label = [el.label for el in list_ext_layer]
    el_ext_layer_seq = els.sequenceFromLabels(tuple(set(list_ext_layer_label)))
    ES_ext_layer = part.Set(elements=el_ext_layer_seq, name=elSetName)
    return list_ext_layer_label
            
def write_job(jobName, mdlName , nCPUS=8):
    mdb.Job(name=jobName, model=mdlName, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCPUS,
        numDomains=nCPUS, numGPUs=0)

def write_job_launch_it(jobName, mdlName , nCPUS=8):
    write_job(jobName, mdlName , nCPUS=8)
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    print('Job: ' + jobName + ' launched')

def launch_inp(jobName, nCPUS=8):
    import os
    cmd  = 'abq2019 job='
    cmd += jobName
    cmd += ' cpus='+str(nCPUS)
    cmd += ' ask_delete=OFF'
    os.system(cmd)

    
def read_inp_write_job_launch_it(jobName , nCPUS=8):
    # Import inp file, write corresponding job, and run it
    mdb.ModelFromInputFile(name=jobName, 
        inputFileName=jobName + '.inp')
    write_job_launch_it(jobName, jobName , nCPUS=8)

def get_inp_text(mdlName):
    # READ original INP file
    with open(mdlName + '.inp', 'r') as inpFile :
        inpText = inpFile.read()
    return inpText

def get_Mat_From_inp(inpText,matKeyword):
    # Read inp file to get the material properties
    import re
    pattern = re.compile("Material, name=(.*)\n[*](?:Density.*\n.*\n.*|Elast.*)\n[ ]*([0-9]*\.[0-9]*).*\n",re.IGNORECASE)
    pElasMod = re.compile("Material, name="+matKeyword+".*\n[*](?:Density.*\n.*\n.*|Elast.*)\n[ ]*([0-9]*\.[0-9]*).*\n",re.IGNORECASE)
    elasMod = [float(E) for E in pElasMod.findall(inpText,re.I)]
    pMat = re.compile("Material, name=("+matKeyword+".*)\n[*](?:Density.*\n.*\n.*|Elast.*)\n[ ]*([0-9]*\.[0-9]*).*\n",re.IGNORECASE)
    Dict_Mat_E = dict(pMat.findall(inpText,re.I))
    # convert elastic modulus to float
    for k in Dict_Mat_E.keys(): 
        Dict_Mat_E[k] = float(Dict_Mat_E[k])
    return Dict_Mat_E


def get_ES_From_inp(inpText,elSetKeyword):
    #FindElset with regex
    import re
    pElset = re.compile('[*]*Elset, elset=('+elSetKeyword+'[0-9]+)\n((?:(?:[ ]*[0-9]+,)*[ ]*[0-9]+(?:\n)*)*)',re.IGNORECASE)
    elset = pElset.findall(inpText,re.I)
    pEl = re.compile('([0-9]+)')
    dict_Elset = dict(elset)
    # convert elements to integer
    for k in dict_Elset.keys():
        dict_Elset[k] = [int(el) for el in pEl.findall(dict_Elset[k])]
    return dict_Elset

def get_ES2MAT_From_inp(inpText,matKeyword,elSetKeyword,Dict_Mat_E):
    # Find the association between the element set and the material properties
    import re
    # FindSection : Assocition between ElSet and Materials
    pSection = re.compile('[*][ ]*Solid Section, elset=('+elSetKeyword+'[0-9]+), material=('+matKeyword+'[0-9]+)',re.IGNORECASE)
    section = pSection.findall(inpText,re.I)
    dict_Mat_Elset = dict(section)
    dict_ELset_E = dict()
    for k in dict_Mat_Elset.keys():
        dict_ELset_E[k] = Dict_Mat_E[dict_Mat_Elset[k]]
    return dict_ELset_E, dict_Mat_Elset

def a_rho(rho,law='Adams2014'):
    # Specific surface for the remodelling equation (see Huiskes et Al. 1992 p 127)
    if law == 'Adams2014':
        # normalized so that max(a_rho) = 1
        p_i = [-0.0693, -0.1595, 1.005, -2.138, 2.357]
    elif law == 'Martin' :
        p_i = []
    # Get the value of a_rho from the fifth order polynomial fit
    return np.sum( [p*rho**(5-i) for i, p in enumerate(p_i)] )


def rho_from_E(E,law='Carter77') :
    if law == 'Carter77' :
        # E = 3790.rho^3
        rho = (E/3790.0)**(1.0/3.0)
    elif law == 'Morgan2003' :
        if E < 5280.4 :
            rho = (E/15520.0)**(1.0/1.93)
        elif E > 12261.5 :
            rho = (E/10714.0)**(1.0/0.74)
        else :
            rho = (E/9965.0)**(1.0/1.137)

    # Ensure that rho is always between 0.01 and 1.85 g/cm^3
    rho = max(rho,0.01)
    rho = min(rho,1.85)
    return rho
    
def rho_to_E(rho,law='Carter77') :
    if law == 'Carter77' :
        E = 3790.0*rho*3.0
    elif law == 'Morgan2003' :
        if rho < 0.572 :
            E = 15520.0*rho**1.93
        elif rho > 1.2 :
            E = 10714.0*E**0.74
        else :
            E = 9965.0*E**1.137
    # Ensure that E is always positive
    E = max(0.5,E)
    return E

def bone_remodeling(E, S, Sref, dt, law='Carter77') :
    # Function that control the remodelling of the bulk bone
    # i : current state
    # ip1 : i+1 -> next step of forward euler scheme
    lz = 0.75 # Lazy zone coefficient -> Huiskes 1992
    tau = 115 # remodelling rate
    rho_i = rho_from_E(E,law)
    a_rho = a_rho(rho_i, 'Adams2014') # Specific surface
    
    if S < (1-lz)*Sref :
        delta_rho =  tau * a_rho * ( S - (1-lz)*Sref )
    elif  S > (1+lz)*Sref :
        delta_rho =  tau * a_rho * ( S - (1+lz)*Sref )
    else :
        delta_rho = 0

    rho_ip1 = rho_i + delta_rho
    # Ensure that rho is always between 0.01 and 1.85 g/cm^3
    rho_ip1 = max(rho_ip1,0.01)
    rho_ip1 = min(rho_ip1,1.85)
    E_ip1 = rho_to_E(rho_ip1,law)
    return E_ip1

def tbcmt_remodeling(E,SvM,dt) :
    # Function that control the remodelling of the bulk bone
    from math import exp
    # Fitted parameters
    x = [3.0910 , 0.4648 , 5.5258 , 0.9260 , 2.2579]
    BVTV = 8.3436e-10*E**2.6133
    V_rmdl_max = 0.28
    SvM0 = x[0]*BVTV**x[1]
    dSvM = SvM - SvM0;
    lmbda = x[2]*BVTV**x[3]
    diff_E_fit = -(V_rmdl_max*dt*E*exp(-dSvM/lmbda)/((1+exp(-dSvM/lmbda))))*x[4]
    newE = E + diff_E_fit
    return newE
