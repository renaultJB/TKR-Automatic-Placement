import numpy as np
import os

from abaqus import *
from abaqusConstants import *
import __main__
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

def creaPreOpMdl(mdlName,nCPUs):
    import AbqsMdlFuns as MdlFuns
    cwd = os.getcwd()
    #=============================================================================#
    #-----------------------------------------------------------------------------#
    #                   Start    of    code    main   part
    #-----------------------------------------------------------------------------#
    #=============================================================================#
    dictFile = 'Dict_' + '_'.join(mdlName.split('_')[1:4])+'.txt'
    mdlName_r = '_'.join(mdlName.split('.')) #removeDots from mdlName


    tibiaCutInpFileName = '_'.join(mdlName.split('_')[0:4]) + '_cut_WM.inp'

    if not os.path.exists(cwd+'/CAEs/'):
        os.makedirs(cwd+'/CAEs/')

    # Read data
    Data = MdlFuns.readData(dictFile)
    SelectedForces, SelectedMoments = MdlFuns.GetLoadBCsFromGaitPCA(0.0)
        
    # Get Body Weight of current subject to normalise the forces
    BW = Data['BW'][mdlName.split('_')[1]]

    #====================================================
    # Place the implant by changing the node coordinates 
    #====================================================
    impltInpFileName = MdlFuns.transformInpNodesCoor(Data['ImpltType'] +'_surf.inp',Data['T'])

    #==========================================
    # Import models and parts and rename them
    #==========================================
    mdl1 = mdb.ModelFromInputFile(name = mdlName_r, inputFileName = mdlName + '.inp')
    RA = mdl1.rootAssembly
    del RA.features['PART-1-1']

    mdl1.parts.changeKey(fromName='PART-1', toName='TIBIA')
    tib = mdl1.parts['TIBIA']
    #----------------------------
    # Create dictionnary of element - section - material correspondance of tibia
    ES_Name_Labels, Sect_ES_Name, Mat_Sect, EMod_Mat = MdlFuns.getElsetSectMat(mdl1,tib)

    # Import the tibial implants -> Useful to create surfaces and sets at the interface
    mdl1.PartFromInputFile( inputFileName = impltInpFileName)
    mdl1.parts.changeKey(fromName='PART-1', toName='IMPLT')
    implt = mdl1.parts['IMPLT']

    # Import the volume of bone that will be cut and replaced by the implant
        # The associated material must be imported too
    ##mdl1.PartFromInputFile( inputFileName = tibiaCutInpFileName)
    ##mdl1.parts.changeKey(fromName='PART-1', toName='TIBIA_CUT')
    ##tibCut = mdl1.parts['TIBIA_CUT']
    mdlCut = mdb.ModelFromInputFile(name = 'mdlCut', inputFileName = tibiaCutInpFileName)

    dic_SectName_SAId = { sa.sectionName : i for i,sa in enumerate(mdlCut.parts['PART-1'].sectionAssignments)}
    dic_Mat_SectName = {mdlCut.sections[s].material : s for s in mdlCut.sections.keys()}
    for mat in mdlCut.materials.keys():
        mdlCut.materials.changeKey(fromName=mat,toName='MATCUT'+mat[3:])
        sectName = dic_Mat_SectName[mat]
        setName = sectName.split('-')[-1]
        setNameNew = setName.split('_')[0]+'CUT_'+setName.split('_')[1]
        mdlCut.parts['PART-1'].sets.changeKey(fromName=setName, toName=setNameNew)
        setRegion = mdlCut.parts['PART-1'].sets[setNameNew]
        mdlCut.sections[sectName].setValues(material='MATCUT'+mat[3:], thickness=None)
        mdlCut.sections.changeKey(fromName=sectName, toName='SectCUT'+sectName[7:])
        SAId = dic_SectName_SAId[sectName]
        mdlCut.parts['PART-1'].sectionAssignments[SAId].setValues(sectionName='SectCUT'+sectName[7:],
                                                                  offsetType=MIDDLE_SURFACE,
                                                                  offsetField='', offset=0.0, region=setRegion)

    mdl1.Part('TIBIA_CUT', mdlCut.parts['PART-1'])
    mdl1.copyMaterials(sourceModel=mdlCut)
    mdl1.copySections(sourceModel=mdlCut)
    tibCut = mdl1.parts['TIBIA_CUT']

    del mdb.models['mdlCut']
    for key in mdb.models.keys() :
        if key.startswith('Model-') :
            del mdb.models[key]

    # Delete transformed implant Inp file
    os.remove(impltInpFileName)
    #==========================================
    # Get surface and nodes Set of Implants
    #==========================================
    implt_SFC = MdlFuns.getAllExtSFC(implt)
    implt_NS_CTCT = MdlFuns.getNSClose2PlanFromSFC(
                        implt, implt_SFC, 'NS-CTCT', Data['Nxp'],
                        Data['Pt_xp'], 0.025, -1000.)
    implt_NS_PROX_FACE = MdlFuns.getNSClose2PlanFromSFC(implt, implt_SFC, 'NS-PROX-FACE',
                                                        Data['Nxp'],Data['Pt_xp'],4.025,3.975)
    implt_NS_XP = MdlFuns.getNSClose2PlanFromSFC(
                        implt, implt_SFC, 'NS-XP', Data['Nxp'],
                        Data['Pt_xp'], 0.025, -0.025)
    implt_NS_STEM_TIP = MdlFuns.getNSClose2PlanFromSFC(
                        implt, implt_SFC, 'NS-STEM-TIP', Data['Nst'],
                        Data['Pt_StemTip'],1.85, -5.)
    implt_NS_STEM_FLANK = implt.SetByBoolean(name='NS-STEM-FLANK', operation=DIFFERENCE,
                                             sets=(implt_NS_CTCT, implt_NS_XP,implt_NS_STEM_TIP, ))

    implt_NS_STEM_ALL = implt.SetByBoolean(name='NS-STEM-ALL', operation=UNION,
                                             sets=(implt_NS_STEM_TIP, implt_NS_STEM_FLANK, ))

    implt_SFC_PROX_FACE = MdlFuns.getSubSurfFromNS(implt,'NS-PROX-FACE',implt_SFC)
    implt_SFC_CTCT = MdlFuns.getSubSurfFromNS(implt,'NS-CTCT',implt_SFC)
    implt_SFC_XP = MdlFuns.getSubSurfFromNS(implt,'NS-XP',implt_SFC)
    implt_SFC_STEM_TIP = MdlFuns.getSubSurfFromNS(implt,'NS-STEM-TIP',implt_SFC)
    implt_SFC_STEM_FLANK = implt.SurfaceByBoolean(name='SURF-STEM-FLANK', operation=DIFFERENCE,
                                                  surfaces=(implt_SFC_CTCT, implt_SFC_XP,
                                                            implt_SFC_STEM_TIP, ))


    #==========================================
    #      Create surface and Node Sets 
    #             of the tibia 
    #==========================================
    # Get an element Set of the whole volume
    tib_ES_ALL = tib.SetFromElementLabels(name = 'ES-ALL' ,
                               elementLabels = tuple([el.label for el in tib.elements])) 


    # Get the whole exterior surface of the tibia and associated node Set
    tib_SFC = MdlFuns.getAllExtSFC(tib)

    tib_NS_EXT = tib.SetFromNodeLabels(name = 'NS-EXT-ALL' ,
                               nodeLabels = tuple([nd.label for nd in tib_SFC.nodes]))
    #----------------------------
    # Get distal face for future encastrement
    tib_NS_Dist_Face = MdlFuns.getDistalFaceNS(tib,tib_SFC,'NS-DIST-FACE')

    #----------------------------
    # Get node sets of the cut plan on the tiba
    tib_NS_XP = MdlFuns.getNSClose2PlanFromSFC(
                        tib, tib_SFC, 'NS-XP', Data['Nxp'],
                        Data['Pt_xp'], 0.025, -0.025)
    tib_SFC_XP = MdlFuns.getSubSurfFromNS(tib,'NS-XP',tib_SFC)

    #----------------------------
    # Get node sets of ctct with implant
    tib_NS_CTCT = MdlFuns.getNSclose2NS(tib,tib_NS_EXT,implt_NS_CTCT,'NS-CTCT',1.5)
    tib_NS_CTCTonXP = tib.SetByBoolean(name='NS-CTCT-ON-XP', operation=INTERSECTION,
                                sets=(tib_NS_XP, tib_NS_CTCT, )) 
    # Get node sets of Stem cut in tibia
    tib_NS_STEM_0 = MdlFuns.getNSclose2NS(tib,tib_NS_EXT,implt_NS_STEM_ALL,'NS-STEM-0',1.5)

    #----------------------------
    # Get node set and surface of the surfaces cut for implantation
    tib_NS_CUT = tib.SetByBoolean(name='NS-CUT', operation=UNION,
                                             sets=(tib_NS_XP, tib_NS_STEM_0, ))
    tib_SFC_CUT = MdlFuns.getSubSurfFromNS(tib,'NS-CUT',tib_SFC)

    # Get Node Set of Stem
    tib_NS_STEM = tib.SetByBoolean(name='NS-STEM', operation=DIFFERENCE,
                                             sets=(tib_NS_STEM_0, tib_NS_XP,  ))
    tib_SFC_STEM = tib.SurfaceByBoolean(name='SURF-STEM', operation=DIFFERENCE,
                                                  surfaces=(tib_SFC_CUT, tib_SFC_XP, ))
    del tib.sets['NS-STEM-0']

    #----------------------------
    # Get node set of Stem Tip
    tib_NS_STEM_TIP = MdlFuns.getNSClose2PlanFromSFC(
                        tib, tib_SFC_CUT, 'NS-STEM-TIP', Data['Nst'],
                        Data['Pt_StemTip'], 1.85, -5.)
    tib_SFC_STEM_TIP = MdlFuns.getSubSurfFromNS(tib,'NS-STEM-TIP',tib_SFC_CUT)

    #----------------------------
    # Get Stem Flank Node Set and Surface
    tib_NS_STEM_FLANK = tib.SetByBoolean(name='NS-STEM-FLANK', operation=DIFFERENCE,
                                sets=(tib_NS_CUT, tib_NS_XP, tib_NS_STEM_TIP, ))
    tib_SFC_STEM_FLANK = tib.SurfaceByBoolean(name='SURF-STEM-FLANK', operation=DIFFERENCE,
                                                  surfaces=(tib_SFC_CUT, tib_SFC_STEM_TIP,
                                                            tib_SFC_XP, ))
    #----------------------------
    # Exterior Surface of the bone
    tib_SFC_CORT = tib.SurfaceByBoolean(name='SURF-CORT', operation=DIFFERENCE,
                                                  surfaces=(tib_SFC,tib_SFC_CUT, ))
    #==========================================
    #      Create elSet of the furure   
    #           cement-bone layer
    #==========================================
    #----------------------------
    # Cement-Bone Global
    tib_ES_LAYER_ALL = MdlFuns.ESlayerFromNS(tib,tib_NS_CTCT,'ES-LAYER-ALL')

    tib_ES_TO_BE_RMDL = tib.SetByBoolean(name='ES-2RMDL', operation=DIFFERENCE,
                                             sets=(tib_ES_ALL, tib_ES_LAYER_ALL,  ))


    #==========================================
    #      Create surface and Node Sets 
    #             of the tibia cut
    #==========================================
    # Get the whole exterior surface of the tibia and associated node Set
    tibCut_SFC = MdlFuns.getAllExtSFC(tibCut)

    tibCut_NS_EXT = tibCut.SetFromNodeLabels(name = 'NS-EXT-ALL' ,
                               nodeLabels = tuple([nd.label for nd in tibCut_SFC.nodes]))
    #----------------------------
    # Get distal face for future encastrement
    tibCut_NS_Dist_Face = MdlFuns.getDistalFaceNS(tibCut,tibCut_SFC,'NS-DIST-FACE')

    #----------------------------
    # Get node sets of the cut plan on the tiba
    tibCut_NS_XP = MdlFuns.getNSClose2PlanFromSFC(
                        tibCut, tibCut_SFC, 'NS-XP', Data['Nxp'],
                        Data['Pt_xp'], 0.025, -0.025)
    tibCut_SFC_XP = MdlFuns.getSubSurfFromNS(tibCut,'NS-XP',tibCut_SFC)

    #----------------------------
    # Get node sets of ctct with implant
    tibCut_NS_CTCT = MdlFuns.getNSclose2NS(tibCut,tibCut_NS_EXT,implt_NS_CTCT,'NS-CTCT',1.5)
    tibCut_NS_CTCTonXP = tibCut.SetByBoolean(name='NS-CTCT-ON-XP', operation=INTERSECTION,
                                sets=(tibCut_NS_XP, tibCut_NS_CTCT, )) 
    # Get node sets of Stem cut in tibia
    tibCut_NS_STEM_0 = MdlFuns.getNSclose2NS(tibCut,tibCut_NS_EXT,implt_NS_STEM_ALL,'NS-STEM-0',1.5)

    #----------------------------
    # Get node set and surface of the surfaces cut for implantation
    tibCut_NS_CUT = tibCut.SetByBoolean(name='NS-CUT', operation=UNION,
                                             sets=(tibCut_NS_XP, tibCut_NS_STEM_0, ))
    tibCut_SFC_CUT = MdlFuns.getSubSurfFromNS(tibCut,'NS-CUT',tibCut_SFC)

    # Get Node Set of Stem
    tibCut_NS_STEM = tibCut.SetByBoolean(name='NS-STEM', operation=DIFFERENCE,
                                             sets=(tibCut_NS_STEM_0, tibCut_NS_XP,  ))
    tibCut_SFC_STEM = tibCut.SurfaceByBoolean(name='SURF-STEM', operation=DIFFERENCE,
                                                  surfaces=(tibCut_SFC_CUT, tibCut_SFC_XP, ))
    del tibCut.sets['NS-STEM-0']

    #----------------------------
    # Get node set of Stem Tip
    tibCut_NS_STEM_TIP = MdlFuns.getNSClose2PlanFromSFC(
                        tibCut, tibCut_SFC_CUT, 'NS-STEM-TIP', Data['Nst'],
                        Data['Pt_StemTip'], 1.85, -5.)
    tibCut_SFC_STEM_TIP = MdlFuns.getSubSurfFromNS(tibCut,'NS-STEM-TIP',tibCut_SFC_CUT)

    #----------------------------
    # Get Stem Flank Node Set and Surface
    tibCut_NS_STEM_FLANK = tibCut.SetByBoolean(name='NS-STEM-FLANK', operation=DIFFERENCE,
                                sets=(tibCut_NS_CUT, tibCut_NS_XP, tibCut_NS_STEM_TIP, ))
    tibCut_SFC_STEM_FLANK = tibCut.SurfaceByBoolean(name='SURF-STEM-FLANK', operation=DIFFERENCE,
                                                  surfaces=(tibCut_SFC_CUT, tibCut_SFC_STEM_TIP,
                                                            tibCut_SFC_XP, ))
    #----------------------------
    # Exterior Surface of the bone
    tibCut_SFC_CORT = tibCut.SurfaceByBoolean(name='SURF-CORT', operation=DIFFERENCE,
                                                  surfaces=(tibCut_SFC,tibCut_SFC_CUT, ))
    #---------------------------
    # Get articular surface of bone
    tibCut_NS_CONDYLES = tibCut_NS_XP = MdlFuns.getNSClose2PlanFromSFC(
                        tibCut, tibCut_SFC, 'NS-CONDYLES', Data['Ztp'],
                        Data['O_knee'], 2.0, -1.5)
    tibCut_SFC_CONDYLES = MdlFuns.getSubSurfFromNS(tibCut,'NS-CONDYLES',tibCut_SFC)


    ## Delete the implant that is not useful anymore
    del mdl1.parts['IMPLT']
    #=============================================================================
    # End of operation on parts
    #=============================================================================
    #=============================================================================
    # Building assembly
    #=============================================================================
    #-+-+-+-+-+-+-+-+-+-+-+-+
    # Add parts to assembly
    #-+-+-+-+-+-+-+-+-+-+-+-+
    RA.Instance(name='TIBIA-1', part=tib, dependent=ON)
    RA.Instance(name='TIBIA_CUT-1', part=tibCut, dependent=ON)
    RA.regenerate()

    tibRA = RA.instances['TIBIA-1']
    tibCutRA = RA.instances['TIBIA_CUT-1']

    LSign = 1.0 if Data['LegSide'] == 'R' else -1.0


    #==========================================
    #    Interaction 
    #==========================================
    # Create Csys of implant and tibia
    CsysTib = MdlFuns.createCSysTibia(RA,Data)
    R_Tib = RA.datums[CsysTib.id]

    CsysTibialPlateau = MdlFuns.createCSysTibialPlateau(RA,Data)
    R_TP = RA.datums[CsysTibialPlateau.id]

    # Implant to bone interaction
    Surf_tibCutCUT = tibCutRA.surfaces['SURF-CUT']
    Surf_tibCUT = tibRA.surfaces['SURF-CUT']

    # Add tie Constraint
    ##mdl1.Tie(name='TIE_CTCT', master=Surf_tibCutCUT, 
    ##        slave=Surf_tibCUT, positionToleranceMethod=COMPUTED, adjust=ON, 
    ##        tieRotations=ON, thickness=ON)

    # Add Cohesive Contact
    mdl1.ContactProperty('Coh_LC_5000')
    mdl1.interactionProperties['Coh_LC_5000'].NormalBehavior(
            pressureOverclosure=LINEAR, contactStiffness=5000.0, 
            constraintEnforcementMethod=DEFAULT)
    mdl1.interactionProperties['Coh_LC_5000'].CohesiveBehavior(
            defaultPenalties=OFF, table=((5000.0, 5000.0, 5000.0), ))
    mdl1.SurfaceToSurfaceContactStd(name='STS_Cohesive', 
            createStepName='Initial', master=Surf_tibCUT, slave=Surf_tibCutCUT, 
            sliding=SMALL, thickness=ON, interactionProperty='Coh_LC_5000', 
            adjustMethod=OVERCLOSED, initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None, tied=OFF)

    # - - - - - - - - - - - - - - - - - - - - 
    # Points to Tibia Cut interactions
    # - - - - - - - - - - - - - - - - - - - -
        # Define the reference points on the implants and associated Geometric Sets
    RA_Rpts = RA.referencePoints
    RPtOriKnee = RA.ReferencePoint(point=tuple(Data['O_knee']))
    RA.features.changeKey(
            fromName=RPtOriKnee.name , toName='O_Knee')
    GS_RPtOriKnee = RA.Set(referencePoints=(RA_Rpts[RPtOriKnee.id],), name='GS_O_Knee')
                                        #-----#
    RPtMedCond = RA.ReferencePoint(point=tuple(Data['Pt_MCC']))
    RA.features.changeKey(
            fromName=RPtMedCond.name , toName='RPt_Med_Cond')
    GS_RPtMedCond = RA.Set(referencePoints=(RA_Rpts[RPtMedCond.id],), name='GS_RPtMedCond')
                                        #-----#
    RPtLatCond = RA.ReferencePoint(point=tuple(Data['Pt_LCC']))
    RA.features.changeKey(
            fromName=RPtLatCond.name , toName='RPt_Lat_Cond')
    GS_RPtLatCond = RA.Set(referencePoints=(RA_Rpts[RPtLatCond.id],), name='GS_RPtLatCond')

    GS_RPtMedLatConds = RA.Set(referencePoints=(RA_Rpts[RPtMedCond.id], RA_Rpts[RPtLatCond.id],),
                               name='GS_RPtMedLatCond')

        # Get the proximal surface of the implant

    SFC_CONDYLES = tibCutRA.surfaces['SURF-CONDYLES']

        # Define the coupling of the med and lat points with the surface
        # Radius determined from normal surface contact LIU2010 with 20% diffusion
    mdl1.Coupling(name='PtMedCond', 
        controlPoint= GS_RPtMedCond , surface= SFC_CONDYLES , influenceRadius=12.5, 
        couplingType=DISTRIBUTING, weightingMethod=CUBIC, localCsys=R_TP, 
        u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF) 

    mdl1.Coupling(name='PtLatCond', 
        controlPoint= GS_RPtLatCond , surface= SFC_CONDYLES , influenceRadius=12.0, 
        couplingType=DISTRIBUTING, weightingMethod=CUBIC, localCsys=R_TP, 
        u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)   

        # Define the coupling of the Implant Origin with the med and lat points with the surface
    mdl1.Coupling(name='PtOrigin2MedLatConds', 
        controlPoint= GS_RPtOriKnee , surface =  GS_RPtMedLatConds,  
        couplingType=DISTRIBUTING, influenceRadius=WHOLE_SURFACE,
        localCsys=R_TP, 
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)  

    # =============================================================================
    # Boundary conditions
    # =============================================================================



    # Encastrement
    NS_ENCASTRE_RA = tibRA.sets['NS-DIST-FACE']
    mdl1.EncastreBC(name='ENCASTRE_TIB_DISTAL', 
            createStepName='Initial', region=NS_ENCASTRE_RA, localCsys=None)

    # --------------------------
    # LoadSteps
    steps = ['GC_15','GC_50','CU','SU1']
    i=0
    allSteps = ['Initial'] + steps
    for s in allSteps[1:] :
        mdl1.StaticStep(name=s, previous=allSteps[i], nlgeom=OFF)
        i += 1

    # Apply Implant Forces
    i=0
    for stepName in steps :
        F_x, F_y, F_z =  list(SelectedForces[stepName]*BW)
        M_x, M_y, M_z =  list(SelectedMoments[stepName]*BW*1000.) #Convert to N.mm

        if 'GC' in stepName :
            #Adapt M_y to HKA:
            U = R_TP.axis2.direction
            M_y = MdlFuns.M_y_from_MFR(F_z,'PreOp',Data,Data['Pt_LCC'],Data['Pt_MCC'],U)
        
        F_x, F_y, F_z = [round(F_x,1), round(F_y,1), round(F_z,1)]
        M_x, M_y, M_z = [round(M_x,1), round(M_y,1), round(M_z,1)]
        # M_x was removed from forces because it is affected by the point of expression of the torseur

        if i==0:
            mdl1.ConcentratedForce(name='TKRForces', 
                createStepName=stepName, region=GS_RPtOriKnee, cf1=F_x, cf2=LSign*F_y, cf3=F_z, 
                distributionType=UNIFORM, field='', localCsys=R_TP)

            mdl1.Moment(name='TKRMoments', 
                createStepName=stepName, region=GS_RPtOriKnee, cm1=LSign*0.1, cm2=M_y, cm3=LSign*0.1, 
                distributionType=UNIFORM, field='', localCsys=R_TP)
        else :
            mdl1.loads['TKRForces'].setValuesInStep(stepName=stepName, cf1=F_x, cf2=LSign*F_y, cf3=F_z)           
            mdl1.loads['TKRMoments'].setValuesInStep(stepName=stepName, cm1=LSign*0.1, cm2=M_y, cm3=LSign*0.1)

        i += 1
            

    # Apply MuscleForces
    forces = MdlFuns.importMusclesForces(Data)
    sfcExtTib = tibRA.surfaces['SURF-CORT']
    dict_MuscleRpts = MdlFuns.generateMusclesLoads(mdl1,steps,sfcExtTib,Data,forces,BW,CsysTib,0.5)

    # Biceps femoris is applied on fibula -> proportion that transfers to the tibia is unknown
    mdl1.loads['fBF'].suppress() # We decided to suppress it
                
    # =============================================================================
    # Add fieldoutputs
    # =============================================================================
    mdl1.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'MISES', 'MISESONLY', 'E', 'LE', 'U', 'RF', 'CF', 'P', 
        'CSTRESS', 'CDISP', 'CFORCE', 'ELEDEN'))

    mdl1.historyOutputRequests['H-Output-1'].setValues(
        variables=('CFNM', 'CFN1', 'CFN2', 'CFN3', 'ALLAE', 'ALLCD', 'ALLDMD', 
        'ALLEE', 'ALLFD', 'ALLIE', 'ALLJD', 'ALLKE', 'ALLKL', 'ALLPD', 'ALLQB', 
        'ALLSE', 'ALLSD', 'ALLVD', 'ALLWK', 'ETOTAL'))

    # =============================================================================
    # Create job and write input and save models as .cae
    # =============================================================================
    jobName = mdl1.name[4:-3]+'_NoOp'
    mdb.Job(name=jobName, model=mdl1.name, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCPUs, 
        numDomains=nCPUs, numGPUs=0)

    mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
    mdb.saveAs(pathName=cwd+'/CAEs/'+jobName)
    
