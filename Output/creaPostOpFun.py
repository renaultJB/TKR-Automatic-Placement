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

def creaPostOpMdl(mdlName,nCPUs) :
    import AbqsMdlFuns as MdlFuns
    cwd = os.getcwd()
    #=============================================================================#
    #-----------------------------------------------------------------------------#
    #                   Start    of    code    main   part
    #-----------------------------------------------------------------------------#
    #=============================================================================#
    dictFile = 'Dict_' + '_'.join(mdlName.split('_')[1:4])+'.txt'
    mdlName_r = '_'.join(mdlName.split('.')) #removeDots from mdlName


    # Read data
    Data = MdlFuns.readData(dictFile)
    SelectedForces, SelectedMoments = MdlFuns.GetLoadBCsFromGaitPCA(0.0)
        
    # Get Body Weight of current subject to normalise the forces
    BW = Data['BW'][mdlName.split('_')[1]]

    #====================================================
    # Place the implant by changing the node coordinates 
    #====================================================
    impltInpFileName = MdlFuns.transformInpNodesCoor(Data['ImpltType'] +'_surf.inp',Data['T'])
    ELabels_HMWPE_implt = MdlFuns.getESLabelFromInp(Data['ImpltType'] +'_surf.inp','ELSET=Volume2')

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

    mdl1.PartFromInputFile( inputFileName = impltInpFileName)
    mdl1.parts.changeKey(fromName='PART-1', toName='IMPLT')
    implt = mdl1.parts['IMPLT']

    # Delete transformed implant Inp file
    os.remove(impltInpFileName)
    #==========================================
    # Add parts to assembly
    #==========================================
    RA.Instance(name='TIBIA-1', part=tib, dependent=ON)
    RA.Instance(name='IMPLT-1', part=implt, dependent=ON)
    RA.regenerate()

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
    implt_SFC_NO_CONTACT = implt.SurfaceByBoolean(name='SURF-NO-CTCT-BONE', operation=DIFFERENCE,
                                                  surfaces=(implt_SFC, implt_SFC_CTCT, ))

    #----------------------------
    # Border of the cut surface on the bone
    implt_NS_PLATE_BORDER = MdlFuns.getNSSfcIntersect(implt,implt_SFC_CTCT,implt_SFC_NO_CONTACT,'NS-PLATE-BORDER')

    #==========================================
    #      Assign section to Implant
    #          Get UHMPWE section
    #==========================================
    # Get elements sets
    implt_ES_ALL = implt.SetFromElementLabels(name = 'ES-ALL' ,
                               elementLabels = tuple([el.label for el in implt.elements]))

    ##implt_ES_HMWPE = implt.SetFromElementLabels(name = 'ES-HMWPE' ,
    ##                           elementLabels = tuple(ELabels_HMWPE_implt))

    implt_ES_HMWPE = MdlFuns.getESinCyl(implt,'ES-HMWPE',Data['Pt_Cyl'],Data['rCyl'],
                                        Data['Ucyl'],Data['heightCyl'])

    ##implt_ES_TA6V = implt.SetByBoolean( name='ES-TA6V', operation=DIFFERENCE,
    ##                                              sets=(implt_ES_ALL, implt_ES_HMWPE,) )

    implt_ES_CRCO = implt.SetByBoolean( name='ES-CRCO', operation=DIFFERENCE,
                                                  sets=(implt_ES_ALL, implt_ES_HMWPE,) )

    # Create material and section
    dict_E_Sect_TBCMT =  MdlFuns.CreateMaterialsSections(mdl1) # also get a dict of TB-PMMA mix section elastic modulus

    ##implt.SectionAssignment(region=implt_ES_TA6V, sectionName='SECT_TA6V', offset=0.0, 
    ##        offsetType=MIDDLE_SURFACE, offsetField='', 
    ##        thicknessAssignment=FROM_SECTION)

    implt.SectionAssignment(region=implt_ES_CRCO, sectionName='SECT_CRCO', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)

    implt.SectionAssignment(region=implt_ES_HMWPE, sectionName='SECT_UHMWPE', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)


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

    #----------------------------
    # Border of the cut surface on the bone
    tib_NS_CUT_BORDER = MdlFuns.getNSSfcIntersect(tib,tib_SFC_CORT,tib_SFC_CUT,'NS-CUT-BORDER')

    #==========================================
    #      Create the cement-bone layer  
    #             element set
    #==========================================
    #----------------------------
    # Cement-Bone Global
    ES_LAYER_ALL = MdlFuns.ESlayerFromNS(tib,tib_NS_CTCT,'ES-LAYER-ALL')
    #----------------------------
    # Cement-Bone Layer at XP
    ES_LAYER_XP = MdlFuns.ESlayerFromNS(tib, tib_NS_CTCTonXP, 'ES-LAYER-XP')
    #----------------------------
    # Cement-Bone on Stem Flank
    ES_LAYER_STEM_FLANK = MdlFuns.ESlayerFromNS(tib, tib_NS_STEM_FLANK,
                                                'ES-LAYER-STEM-FLANK-0')
    ES_LAYER_STEM_FLANK = tib.SetByBoolean(name='ES-LAYER-STEM-FLANK', operation=DIFFERENCE,
                                                  sets=(ES_LAYER_STEM_FLANK, ES_LAYER_XP, ))
    del tib.sets['ES-LAYER-STEM-FLANK-0']
    #----------------------------
    # Cement-Bone Layer at Stem Tip
    ES_LAYER_STEM_TIP = tib.SetByBoolean(name='ES-LAYER-STEM-TIP', operation=DIFFERENCE, sets=(ES_LAYER_ALL, ES_LAYER_STEM_FLANK, ES_LAYER_XP, ))

    #----------------------------
    # Get the vol of bone that is not mixed with bone cement and that will be remodelled
    tib_ES_TO_BE_RMDL = tib.SetByBoolean(name='ES-2RMDL', operation=DIFFERENCE,
                                             sets=(tib_ES_ALL, ES_LAYER_ALL,  ))

    #==========================================
    #    Find the density of the elements  
    #       within the cmt-TB layer
    #   and assign material to the layer
    #==========================================
    #----------------------------
    # Substract ES layer from dict of Elset:
    ElmtLayerSet = set([el.label for el in ES_LAYER_ALL.elements])
    ES_Name_Labels_updated = {es : list(set(elLabels)-ElmtLayerSet) for es, elLabels in ES_Name_Labels.iteritems()}

    # Update SET by removing elements xithin TB-CMT layer
    MdlFuns.updateCreateElset(tib,ES_Name_Labels_updated)

    # Build element -> elastic moduli dictionnary
    ES_Name_Sect = {v : k for k,v in Sect_ES_Name.iteritems()}
    Sect_Mat = { v : k for k,v in Mat_Sect.iteritems()}
    Mat_EMod = { v : k for k,v in EMod_Mat.iteritems()}
    EMod_Labels = { Mat_EMod[Sect_Mat[ES_Name_Sect[k]]] : v for k,v in ES_Name_Labels.iteritems()}
    elmt_EMod = { el : E for E,Labels in EMod_Labels.iteritems() for el in Labels}

    # assign material to the layer elements
    dict_E_Elmts_CMTTB = MdlFuns.assignMat2TBCMTLayer(tib,ES_LAYER_ALL,dict_E_Sect_TBCMT,elmt_EMod)


    #=============================================================================
    # End of operation on parts
    #=============================================================================
    #=============================================================================
    # Building assembly
    #=============================================================================
    RA.regenerate()

    tibRA = RA.instances['TIBIA-1']
    impltRA = RA.instances['IMPLT-1']
    LSign = 1.0 if Data['LegSide'] == 'R' else -1.0


    #==========================================
    #    Interaction 
    #==========================================
    # Create Csys of implant and tibia
    CsysImplt = MdlFuns.createCSysImplant(RA,Data)
    R_Implt = RA.datums[CsysImplt.id]

    CsysTib = MdlFuns.createCSysTibia(RA,Data)
    R_Tib = RA.datums[CsysTib.id]

    # Implant to bone interaction
    Surf_impltCTCT = impltRA.surfaces['SURF-CTCT']
    Surf_tibCUT = tibRA.surfaces['SURF-CUT']

    LCF_name = MdlFuns.createLCFproperty(mdl1,1000.0,0.3)

    STS_InterA = mdl1.SurfaceToSurfaceContactStd(name='STS-'+LCF_name, 
        createStepName='Initial', master=Surf_impltCTCT, slave=Surf_tibCUT, sliding=FINITE, 
        thickness=ON, interactionProperty=LCF_name, adjustMethod=OVERCLOSED, 
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF)
    # - - - - - - - - - - - - - - - - - - - - 
    # Points to implants interactions
    # - - - - - - - - - - - - - - - - - - - -
        # Define the reference points on the implants and associated Geometric Sets
    RA_Rpts = RA.referencePoints
    RPtOriImplt = RA.ReferencePoint(point=tuple(Data['O_Implt']))
    RA.features.changeKey(
            fromName=RPtOriImplt.name , toName='O_Implt')
    GS_RPtOriImplt = RA.Set(referencePoints=(RA_Rpts[RPtOriImplt.id],), name='GS_O_Implt')
                                        #-----#
    RPtMedImplt = RA.ReferencePoint(point=tuple(Data['Pt_Med_Implt']))
    RA.features.changeKey(
            fromName=RPtMedImplt.name , toName='RPt_Med_Implt')
    GS_RPtMedImplt = RA.Set(referencePoints=(RA_Rpts[RPtMedImplt.id],), name='GS_RPtMedImplt')
                                        #-----#
    RPtLatImplt = RA.ReferencePoint(point=tuple(Data['Pt_Lat_Implt']))
    RA.features.changeKey(
            fromName=RPtLatImplt.name , toName='RPt_Lat_Implt')
    GS_RPtLatImplt = RA.Set(referencePoints=(RA_Rpts[RPtLatImplt.id],), name='GS_RPtLatImplt')

    GS_RPtMedLatImplt = RA.Set(referencePoints=(RA_Rpts[RPtMedImplt.id], RA_Rpts[RPtLatImplt.id],),
                               name='GS_RPtMedLatImplt')

        # Get the proximal surface of the implant
    SFC_Prox_Face_implt = impltRA.surfaces['SURF-PROX-FACE']

        # Define the coupling of the med and lat points with the surface
    mdl1.Coupling(name='PtMedImplt', 
        controlPoint= GS_RPtMedImplt , surface= SFC_Prox_Face_implt , influenceRadius=12.0, 
        couplingType=DISTRIBUTING, weightingMethod=CUBIC, localCsys=RA.datums[CsysImplt.id], 
        u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)    

    mdl1.Coupling(name='PtLatImplt', 
        controlPoint= GS_RPtLatImplt , surface= SFC_Prox_Face_implt , influenceRadius=10.0, 
        couplingType=DISTRIBUTING, weightingMethod=CUBIC, localCsys=RA.datums[CsysImplt.id], 
        u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)   

        # Define the coupling of the Implant Origin with the med and lat points with the surface
    mdl1.Coupling(name='PtOrigin2MedLatImplt', 
        controlPoint= GS_RPtOriImplt , surface =  GS_RPtMedLatImplt,  
        couplingType=DISTRIBUTING, influenceRadius=WHOLE_SURFACE,
        localCsys=RA.datums[CsysImplt.id], 
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
            U = R_Implt.axis1.direction
            M_y = MdlFuns.M_y_from_MFR(F_z,Data['Alignment'],Data,Data['Pt_Med_Implt'],Data['Pt_Lat_Implt'],U)

        F_x, F_y, F_z = [round(F_x,1), round(F_y,1), round(F_z,1)]
        M_x, M_y, M_z = [round(M_x,1), round(M_y,1), round(M_z,1)]
        # M_x was removed from forces because it is affected by the point of expression of the torseur
        
        if i==0:
            mdl1.ConcentratedForce(name='TKRForces', 
                createStepName=stepName, region=GS_RPtOriImplt, cf1=F_x, cf2=LSign*F_y, cf3=F_z, 
                distributionType=UNIFORM, field='', localCsys= R_Implt)

            mdl1.Moment(name='TKRMoments', 
                createStepName=stepName, region=GS_RPtOriImplt, cm1=LSign*0.1, cm2=M_y, cm3=LSign*M_z, 
                distributionType=UNIFORM, field='', localCsys= R_Implt)
        else :
            mdl1.loads['TKRForces'].setValuesInStep(stepName=stepName, cf1=F_x, cf2=LSign*F_y, cf3=F_z)           
            mdl1.loads['TKRMoments'].setValuesInStep(stepName=stepName, cm1=LSign*0.1, cm2=M_y, cm3=LSign*M_z)

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
    jobName = mdl1.name[4:-3]+'_Op_0'
    mdb.Job(name=jobName, model=mdl1.name, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCPUs, 
        numDomains=nCPUs, numGPUs=0)

    mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
    mdb.saveAs(pathName=cwd+'/CAEs/'+jobName)

