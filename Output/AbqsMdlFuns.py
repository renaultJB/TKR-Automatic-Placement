# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:29:55 2017
Updated on  Thu Oct 3 16:30:30 2019

@author: Jean-Baptiste RENAULT
"""
import numpy as np

def GetLoadBCsFromGaitPCA(scaleFactorSD):
    # Get force
    # Force Unit N/Kg
    # Moment Unit N/Kg.m
    # Data are for a right leg
    # Read the PCA results file
    with open('PCA_8GaitLoads.txt','r') as f:
        Latent = [float(l) for l in f.readline().split("\t")[0:-1]]
        mu = [float(l) for l in f.readline().split("\t")[0:-1]]
        eigenVctrs = list()
        test = False
        while not test:
            eigenVctrs.append([float(l) for l in f.readline().split("\t")[0:-1]])
            test = not eigenVctrs[-1]
        
        eigenVctrs.pop()
        
    
    # convert to numpy arrays
    Latent = np.array(Latent)
    mu = np.array(mu)
    eigenVctrs = np.array(eigenVctrs)
    
    # Generate a random Gait Cycle Load from PCA decomposition first ncomp PCs
    ncomp = 4
    
    SampledValueSD = scaleFactorSD*np.random.randn(1,ncomp)
    randMatrix = SampledValueSD*np.sqrt(Latent[0:ncomp])
    
    Curve = np.transpose(mu) + np.dot(randMatrix,eigenVctrs[0:ncomp])

    #Curve = Curve[0]
    
    # Axes Directions X, Y, Z:	Lateral, Anterior, Superior

    # Extract the Loads at key points of Gait Cycle
    #    0% 15% 30% 50% 65% 95%
    
    # Compute resultant force and moment
    Forces = np.concatenate((Curve[:,0::6],Curve[:,1::6],Curve[:,2::6]), axis=0)
    Moments = np.concatenate((Curve[:,3::6],Curve[:,4::6],np.concatenate((Curve[:,5::6],Curve[:,-1:]),axis=1)), axis=0)
    
    FRes = np.sqrt(np.sum(Forces**2,axis=0))
    MRes = np.sqrt(np.sum(Moments**2,axis=0))
    
    #matplotlib.pyplot.plot(FRes) 
    
    SelectedForces = dict()
    SelectedMoments = dict()
    
    SelectedForces['GC_00'] = Forces[:,0]
    SelectedMoments['GC_00'] = Moments[:,0]
    
    for key in [15, 50] :
        #find peak load
        Imax = np.argmax(FRes[int((key/100.-0.05)*len(FRes)):int((key/100.+0.05)*len(FRes))])
        Imax += int((key/100.-0.05)*len(FRes))
        
        SelectedForces['GC_'+str(key)] = Forces[:,Imax]
        SelectedMoments['GC_'+str(key)] = Moments[:,Imax]
         

    for key in [30, 65, 95] :
        SelectedForces['GC_'+str(key)] = Forces[:,int((key/100.)*len(FRes))]
        SelectedMoments['GC_'+str(key)] = Moments[:,int((key/100.)*len(FRes))]
    
    
    SelectedForces['Forces'] = Forces
    SelectedForces['ForceRes'] = FRes
    
    SelectedMoments['Moments'] = Moments
    SelectedMoments['MomentRes'] = MRes
    
    SelectedForces['PCA_Sample'] = str(SampledValueSD[0])
    SelectedMoments['PCA_Sample'] = str(SampledValueSD[0])

    # Forces for stair Up peak 1 and 2 and Chair Up activity in N/Kg
    # Mass in Kg of each subject is required
    SelectedForces['SU1'] = np.array([1.14, -1.81, -30.37])
    SelectedForces['SU2'] = np.array([-0.15, -2.29, -29.85])
    SelectedForces['CU'] = np.array([0.322, -0.192, -27.48])

    # Moments in N/Kg.Meters
    SelectedMoments['SU1'] = np.array([0.147, -0.23, -0.0073])
    SelectedMoments['SU2'] = np.array([0.182, -0.193, -0.072])
    SelectedMoments['CU'] = np.array([0.12, 0.099, -0.054])
    
    return SelectedForces, SelectedMoments

def M_y_from_MFR(Fz,case,data,Pt1,Pt2,U):
    # Compute the moment in N.mm induced by the HKA
    # Assume that all the deformation is localised at the tibia level
    Pt1Pt2 = np.array(Pt1)-np.array(Pt2)
    d = np.dot(U,Pt1Pt2)
    d = abs(d)
    Fz = abs(Fz)
    #d2 = ((Pt1[0]-Pt2[0])**2+(Pt1[1]-Pt2[1])**2+(Pt1[2]-Pt2[2])**2)
    #d = d2**0.5
    # Get HKA preOp or Post Op depending on alignement
    if case == 'Mech' :
        HKA = data['alpha']
    elif case == 'Kine' :
        HKA = data['alpha'] + 2 # Physio MMPTA = 88°
    elif case == 'PreOp' : 
        HKA = 88 - data['mMPTA'] # Physio MMPTA = 88°
    # Linear regression obtained from gait cycles  
    M_y = (0.053*HKA-0.374)*Fz*d/2.
    return M_y

def CreateMaterialsSections(Mdl):
    Mdl.Material(name='TA6V')
    Mdl.materials['TA6V'].Elastic(table=((110000.0, 0.34), ))
    Mdl.materials['TA6V'].Density(table=((4.4, ), ))
    Mdl.HomogeneousSolidSection(name='SECT_TA6V',material='TA6V', thickness=None)
    # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mtp641

    Mdl.Material(name='CRCO')
    Mdl.materials['CRCO'].Elastic(table=((240000.0, 0.3), ))
    Mdl.materials['CRCO'].Density(table=((8.3, ), ))
    Mdl.HomogeneousSolidSection(name='SECT_CRCO',material='CRCO', thickness=None)
    # https://www.hainaut-prothese.com/assets/doc/zirkonzahn/CrCoData.pdf

    Mdl.Material(name='UHMWPE')
    Mdl.materials['UHMWPE'].Elastic(table=((862.0, 0.45), ))
    Mdl.materials['UHMWPE'].Density(table=((0.95, ), ))
    Mdl.HomogeneousSolidSection(name='SECT_UHMWPE',material='UHMWPE', thickness=None)
    
    Mdl.Material(name='PMMA')
    Mdl.materials['PMMA'].Elastic(table=((1837, 0.23), ))
    Mdl.materials['PMMA'].Density(table=((1.18, ), ))
    Mdl.HomogeneousSolidSection(name='SECT_PMMA',material='PMMA', thickness=None)
#   Properties of an injectable low modulus PMMA bone cement for osteoporotic bone
#   Boger, Andreas Bohner, Marc Heini, Paul Verrier, Sophie Schneider, Erich
#   Journal of Biomedical Materials Research  2008
#   DOI: 10.1002/jbm.b.31044

    dict_E_Sect = dict()
    TB_PMMA_E = list(np.linspace(1.0, 330., num=223)**1.5)
    TB_PMMA_E += list(np.linspace(6500., 20000., num=28))
    for i,e in enumerate(TB_PMMA_E):
        Mdl.Material(name='TB-PMMA'+str(i))
        Mdl.materials['TB-PMMA'+str(i)].Elastic(table=((e, 0.3), ))
        Mdl.HomogeneousSolidSection(name='SECT_TB-PMMA'+str(i),material='TB-PMMA'+str(i), thickness=None)
        dict_E_Sect[e] = 'SECT_TB-PMMA'+str(i)
        
    
    return dict_E_Sect
    

def getAllExtSFC(p,SET0=[]) :
    from abaqusConstants import S1,S2,S3,S4
    from collections import Counter
    import numpy
    
    elements = p.elements

    npFace   = 6 # Number of nodes per face. If using C3D4 elems, this will be 3
    FacesPerEl = 4 # Number of faces per element. If using C3DX elems, this will be 4
    
    # Get element connectivity (elements connect to nodes)
    elemConn = {}
    for elem in elements:
        e     = elem.label
        nconn = elem.connectivity
        for n in nconn:
            if not elemConn.has_key(n):
                elemConn[n] = []
            elemConn[n].append(e)
    
    # Now get list of all neighbouring elements for each element
    # Use the dictionary value to tally how many nodes are shared. If this
    # number is equal to the number of nodes per face, then faces must also
    # be shared
    
    elemNbr = {}
    for elem in elements:
        e     = elem.label
        nconn = elem.connectivity
        if not elemNbr.has_key(e):
            elemNbr[e] = {}
        for n in nconn:
            for ec in elemConn[n]:
                if not elemNbr[e].has_key(ec):
                    elemNbr[e][ec] =  1
                else: 
                    elemNbr[e][ec] += 1
    
    #elemNbr : Element Neighbour, its a dictionnay that lists in another dictionnary
    #    all the neighbour element and the number of common node found
    
    # Delete element label from its own list of neighbours. Also delete key if
    # number of nodes is not equal to npFace
    for k1 in elemNbr.keys():
        del elemNbr[k1][k1]
        for k2,v2 in elemNbr[k1].items():
            if v2!=npFace: 
                del elemNbr[k1][k2]
        elemNbr[k1] = elemNbr[k1].keys()
    
    
    # find elements that have one or several free face, lying on a surface
    Face1El = []
    Face2El = []
    Face3El = []
    Face4El = []
    
    elemSurf = []
    if not SET0 :
        elemSurf = [el for el in elemNbr.keys() if len(elemNbr[el])<FacesPerEl]
    
    else :
        SET0Label = [el.label for el in SET0.elements ]
        IntersectSet0 = numpy.intersect1d(numpy.array(elemNbr.keys()), numpy.array(SET0Label), assume_unique=False).tolist()
        elemSurf = []
        for el in IntersectSet0:
            if len(elemNbr[el])<FacesPerEl :
                elemSurf.append(el)
    
    
    # Identify the face identifier of the element with a free face
    for elsurf in elemSurf:
        nodesID = []
        elsurfNd = elements.getFromLabel(elsurf).connectivity[0:4]

        # Get the list of all nodes of the neighbour elements of the current elements
        # Keep only the nodes that are shared with the current elements
        # If a node is present 3 time in the list it means that this node is not on the ext surface
        # One element can have more than one exterior faces...
        for elnbr in elemNbr[elsurf]:
            nodesID.append(elements.getFromLabel(elnbr).connectivity[0:4])
        nodesID = [ v for conNdList in nodesID for v in conNdList if v in elsurfNd]
        
    #   Only one Face  
        if 3 in Counter(nodesID).values() :
            NdIDNotSurface = Counter(nodesID).most_common(1)[0][0]
            OppositeNodePlace = elsurfNd.index(NdIDNotSurface)
            
            if OppositeNodePlace == 0:
                Face3El.append(elsurf)
            elif OppositeNodePlace == 1:
                Face4El.append(elsurf)  
            elif OppositeNodePlace == 2:
                Face2El.append(elsurf)
            elif OppositeNodePlace == 3:
                Face1El.append(elsurf)
            else :
                print('error face out of boundary')
             
    #    Dual Faces
        elif 3 not in Counter(nodesID).values() and len(Counter(nodesID))!=3:
            for ndIDMC in Counter(nodesID).most_common(2) :
                OppositeNodePlace = elsurfNd.index(ndIDMC[0])
                if OppositeNodePlace == 0:
                    Face3El.append(elsurf)
                elif OppositeNodePlace == 1:
                    Face4El.append(elsurf)  
                elif OppositeNodePlace == 2:
                    Face2El.append(elsurf)
                elif OppositeNodePlace == 3:
                    Face1El.append(elsurf)
                else :
                    print('error face out of boundary')
                    
    #    Tri Faces
        elif len(Counter(nodesID))==3:
            for ndIDMC in Counter(nodesID).keys() :
                OppositeNodePlace = elsurfNd.index(ndIDMC)
                if OppositeNodePlace == 0:
                    Face3El.append(elsurf)
                elif OppositeNodePlace == 1:
                    Face4El.append(elsurf)  
                elif OppositeNodePlace == 2:
                    Face2El.append(elsurf)
                elif OppositeNodePlace == 3:
                    Face1El.append(elsurf)
                else :
                    print('error face out of boundary')
            
            print('Tri faces')
        
        else :
            print('error in element face identification')
    
    Face1ElSet = p.SetFromElementLabels(name = '_ELSETF1', elementLabels = tuple(Face1El))
    Face2ElSet = p.SetFromElementLabels(name = '_ELSETF2', elementLabels = tuple(Face2El))
    Face3ElSet = p.SetFromElementLabels(name = '_ELSETF3', elementLabels = tuple(Face3El))
    Face4ElSet = p.SetFromElementLabels(name = '_ELSETF4', elementLabels = tuple(Face4El))
    
    SetAllSurf = p.MeshSurfaceFromElsets(name='SURF-EXT-ALL',elementSetSeq=((Face1ElSet,S1),
                                (Face2ElSet,S2),(Face3ElSet,S3),(Face4ElSet,S4)))

    del Face1ElSet
    del Face2ElSet
    del Face3ElSet
    del Face4ElSet

    del p.sets['_ELSETF1']
    del p.sets['_ELSETF2']
    del p.sets['_ELSETF3']
    del p.sets['_ELSETF4']

    
    return SetAllSurf


def getNSClose2PlanFromSFC(p, SFC, NS_out_Name, nrml, ptOnPlan, distPos, distNeg):
    # Convert to numpy arrays
    nrml = np.array(nrml)
    ptOnPlan = np.array(ptOnPlan)
                  
    #Get exterior surface nodes
    # - - - - - - - - - - - - - - - - - -
##    if 'NS-EXT-SFC' not in p.sets.keys() :
##        extNodesLabel = np.array([nd.label for nd in SetAllSFC.nodes])
##        NS_Ext_Sfc = p.SetFromNodeLabels(name = 'NS-EXT-SFC' , nodeLabels = tuple(extNodesLabel))
##    else :
##        NS_Ext_Sfc = p.sets['NS-EXT-SFC']
##        extNodesLabel = np.array([nd.label for nd in NS_Ext_Sfc.nodes])

    # Get node labels from surfaces
    # - - - - - - - - - - - - - - - - - -
    extNodesLabel = np.array([nd.label for nd in SFC.nodes])

    #Get cut surface nodes
    # - - - - - - - - - - - - - - - - - -
    extNodesCoor = np.array([nd.coordinates for nd in SFC.nodes])

    # Get nodes distance to the resected surface plan
    distNd2CutSurf = np.dot( (extNodesCoor - ptOnPlan) , nrml )

    # Get label of elements on resected surface
    goodNodesLogical = np.logical_and( distNd2CutSurf < distPos  ,  distNd2CutSurf > distNeg )
    goodNodesLabel = extNodesLabel[goodNodesLogical]

    # Generate Node Set
    NS_Close2Plan_from_SFC = p.SetFromNodeLabels(name = NS_out_Name , nodeLabels = tuple(goodNodesLabel))

    return NS_Close2Plan_from_SFC


def getNSclose2NS(p, NS0, NS_target, NS_out_Name, distMax = 1.5):
    from abaqusConstants import INTERSECTION,DIFFERENCE

    # Get Stem nodes
    # - - - - - - - - - - - - - - - - - -
    # Get NS_target Bounding box (To reduce compuatution time
    NSCoor = np.array([nd.coordinates for nd in NS_target.nodes])
##    NSBB_corners = [c for c in np.min(NSCoor,axis=0)]
##    NSBB_corners += [c for c in np.max(NSCoor,axis=0)]
    NSBB_corners = [c - distMax*0.5 for c in np.min(NSCoor,axis=0)]
    NSBB_corners += [c + distMax*0.5 for c in np.max(NSCoor,axis=0)]

    # Get a NodeSet of part p nodes within bounding box -> NS_BB
    Nds_in_BB = p.nodes.getByBoundingBox(NSBB_corners[0],NSBB_corners[1],NSBB_corners[2],
                                            NSBB_corners[3],NSBB_corners[4],NSBB_corners[5])
    Nds_in_BB_label = np.array([nd.label for nd in Nds_in_BB])
    NS_BB = p.SetFromNodeLabels(name = 'NS-BB-TMP-ALL' , nodeLabels = tuple(Nds_in_BB_label))

    # Get NS0 nodes within bounding box
    NS0_tmp = p.SetByBoolean(name='NS-Close2NS-SFC-0', operation=INTERSECTION, sets=(NS_BB, NS0, ))

    NS0_tmp_Lab = np.array([nd.label for nd in NS0_tmp.nodes])
    NS0_tmp_Coor = np.array([nd.coordinates for nd in NS0_tmp.nodes])

    # Dist between
    dist2target = np.array([np.sqrt(np.min(np.sum((NSCoor - x)**2, axis=1))) for x in NS0_tmp_Coor])
    NS_Ok_Lab = NS0_tmp_Lab[dist2target<distMax]
    
    NS_Close_2_Target_NS = p.SetFromNodeLabels(name = NS_out_Name , nodeLabels = tuple(NS_Ok_Lab))
    
    #NS_Stem_Sfc = p.SetByBoolean(name='NS-STEM-SFC', operation=DIFFERENCE, sets=(NS_Stem_Sfc_0, NS_Cut_Sfc, ))

    del p.sets['NS-BB-TMP-ALL']
    del p.sets['NS-Close2NS-SFC-0']

    return NS_Close_2_Target_NS    
    

###########################################################################
def getSubSurfFromNS(p,NS_SubSfc_name,SetAllSurf):
###########################################################################
    from abaqusConstants import FACE1,FACE2,FACE3,FACE4
    from abaqusConstants import S1,S2,S3,S4

    NS_SubSfc = p.sets[NS_SubSfc_name]

    #### /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ ####
    #### /!\ /!\ /!\ nd.label-1 used might not work if element label is different from index + 1
    Set_Nds_SubSfc = set([nd.label-1 for nd in NS_SubSfc.nodes])
    #### /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ ####
    #### /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ ####


    Face1El = []
    Face2El = []
    Face3El = []
    Face4El = []


    #### /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ ####
    ####  Warning only works if 'FACESi' are sorted in SetAllSurf.sides  ####
    #### /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ ####
    NbFace1 = 0 + len([1 for side in SetAllSurf.sides if side == FACE1])
    NbFace2 = NbFace1 + len([1 for side in SetAllSurf.sides if side == FACE2])
    NbFace3 = NbFace2 + len([1 for side in SetAllSurf.sides if side == FACE3])
    NbFace4 = NbFace3 + len([1 for side in SetAllSurf.sides if side == FACE4])

    for el in SetAllSurf.elements[0:NbFace1]:
        elConnFace1 = set((el.connectivity[0],el.connectivity[1],el.connectivity[2]))
        ndNotonSurface = elConnFace1 - Set_Nds_SubSfc
        # If there is no node that are not on the surface this a correct face
        if not ndNotonSurface :
            Face1El.append(el.label)

    for el in SetAllSurf.elements[NbFace1:NbFace2]:
        elConnFace2 = set((el.connectivity[0],el.connectivity[1],el.connectivity[3]))
        ndNotonSurface = elConnFace2 - Set_Nds_SubSfc
        # If there is no node that are not on the surface this a correct face
        if not ndNotonSurface :
            Face2El.append(el.label)

    for el in SetAllSurf.elements[NbFace2:NbFace3]:
        elConnFace3 = set((el.connectivity[1],el.connectivity[2],el.connectivity[3]))
        ndNotonSurface = elConnFace3 - Set_Nds_SubSfc
        # If there is no node that are not on the surface this a correct face
        if not ndNotonSurface :
            Face3El.append(el.label)

    for el in SetAllSurf.elements[NbFace3:NbFace4]:
        elConnFace4 = set((el.connectivity[0],el.connectivity[2],el.connectivity[3]))
        ndNotonSurface = elConnFace4 - Set_Nds_SubSfc
        # If there is no node that are not on the surface this a correct face
        if not ndNotonSurface :
            Face4El.append(el.label)
    
    Face1ElSet = p.SetFromElementLabels(name = '_ELSETF1', elementLabels = tuple(Face1El))
    Face2ElSet = p.SetFromElementLabels(name = '_ELSETF2', elementLabels = tuple(Face2El))
    Face3ElSet = p.SetFromElementLabels(name = '_ELSETF3', elementLabels = tuple(Face3El))
    Face4ElSet = p.SetFromElementLabels(name = '_ELSETF4', elementLabels = tuple(Face4El))

    SetSubSurf = p.MeshSurfaceFromElsets(name='SURF-'+NS_SubSfc_name[3:],elementSetSeq=((Face1ElSet,S1),
                                (Face2ElSet,S2),(Face3ElSet,S3),(Face4ElSet,S4)))

    del Face1ElSet
    del Face2ElSet
    del Face3ElSet
    del Face4ElSet


    del p.sets['_ELSETF1']
    del p.sets['_ELSETF2']
    del p.sets['_ELSETF3']
    del p.sets['_ELSETF4']
    
    return SetSubSurf

def getNSSfcIntersect(p,sfc1,sfc2,NS_out_Name) :
    # Get the Node Set of the intersection of the nodes of two surfaces sfc1 and sfc2
    nd_lab1 = set([nd.label for nd in sfc1.nodes])
    nd_lab2 = set([nd.label for nd in sfc2.nodes])
    nd_lab_itrsct = nd_lab1.intersection(nd_lab2)
    NS_itrsct = p.SetFromNodeLabels(name = NS_out_Name , nodeLabels = tuple(nd_lab_itrsct))
    return NS_itrsct

    


def getESinCyl(p,ES_name,centrPt,radius,u,height):
    #             ^  Axis : u      
    #             |                
    #       ------+------               
    #             |
    #             |
    #             | Height
    #             |
    #             |
    #       ------+------  radius
    #           centrPt
    #
    #
    #
    #

    centrPt=np.array(centrPt)

    # Get the centroid of the elements of the part
    #dict_ElLab_Centroid = {el.label : np.mean([p.nodes[nd].coordinates for nd in el.connectivity], axis = 0) for el in p.elements}
    Elmts_Label = np.array([el.label for el in p.elements])
    Elmts_Centroid = np.array([ np.mean([p.nodes[nd].coordinates for nd in el.connectivity[0:4]], axis = 0) for el in p.elements ])


    # Get orthogonal projection matrix
    u = np.array([u,])
    u = np.transpose(u)/np.linalg.norm(u)
    u_T = np.transpose(u)

    P = np.dot(u,np.dot(np.linalg.inv(np.dot(u_T,u)),u_T))

    # Project elements on axis
    Centr_Center = np.transpose(Elmts_Centroid-centrPt)
    Elmts_Centroid_Proj = np.transpose( np.dot(P, Centr_Center ) ) + centrPt

    # Get Distance to axis and distance to origin point 
    dist2Axis = np.sqrt(np.sum((Elmts_Centroid_Proj-Elmts_Centroid)**2, axis = 1))
    dist2CenterPt = np.dot(u_T,np.transpose(Elmts_Centroid_Proj-centrPt))[0]

    # Identify element within the cylinder
    # if element form cylinder border its centroid must be inside
    inside_Cyl_El_Lab_Logical = np.logical_and( dist2Axis < radius  ,  dist2CenterPt < height )

    El_Lab_in_Cylinder = Elmts_Label[inside_Cyl_El_Lab_Logical]

    ES_in_Cylinder = p.SetFromElementLabels(name = ES_name , elementLabels = tuple(El_Lab_in_Cylinder))

    return ES_in_Cylinder    

def getDistalFaceNS(p,SURF_Ext,NSname='NS-DIST-FACE'):
    BB = p.elements.getBoundingBox()
    BB = tuple(np.array(BB['low']+BB['high'])+
                  np.array([-3,-3,-10,3,3,
                10-float(BB['high'][2])+float(BB['low'][2])]))
    NodesInBB = SURF_Ext.nodes.getByBoundingBox(*BB)
    pts = np.array([ nd.coordinates for nd in NodesInBB ])
    pts_0 = pts
    for delta in [10.,8.5,7.,5.,2.5,2.,1.5,1.,0.75,0.5,0.25,0.1,0.1] :
        center = np.mean(pts,axis=0)
        # Find least square plan
        covMat = np.cov(pts.T)
        w,v = np.linalg.eig(covMat)
        nrml = v.T[2]
        if nrml[2] < 0 :
            nrml = - nrml
        # Get nodes distance to the resected fitted plan
        distNd2FittedPlan = np.dot( (pts - center) , nrml )
        # Get label of elements on resected surface
        pts = pts[distNd2FittedPlan<=float(delta)]

    NS_Dist_Face = getNSClose2PlanFromSFC(p, SURF_Ext, NSname, nrml, center, 0.05, -5)

    if len(NS_Dist_Face.nodes)<150 :
        print('Low number of nodes on distal face '+ NSname )
        print('Method was changed to Bounding Box Method')
        center = np.mean(pts_0,axis=0)
        nrml = [0,0,1]
        NS_Dist_Face = getNSClose2PlanFromSFC(p, SURF_Ext, NSname, nrml, center, 1., -5.)

    return NS_Dist_Face

def transformInpNodesCoor(inpFile,T=[[1.,0,0,0],[0,1.,0.,0],[0.,0.,1.,0.],[0.,0,0,1]]):
    # Change the node coordinates of the mesh within the inp file by applying
    # the affine transformatio matrix T (Roto-Translation Matrix)
    
    #inpFile = 'S4_surf_4GMSH.inp'
    inpOutputFile = inpFile[:-4]+'_T.inp'
    
    output=''
    
    # Read inp file and transform the coordinates
    # write updated lines in output
    with open(inpFile,'r') as f:
        line = f.readline()
        while '*NODE' not in line.upper():
            output += line
            line = f.readline()
        output +=  line
        line = f.readline()
        while not line.startswith('*'):
            lineSplit = line.split(', ')
            coor = [float(c) for c in lineSplit[1:4]]
            coor += [1]
            coorNp = np.transpose([coor,])
            coorNp_trans = np.transpose(np.dot(T,coorNp))
            updatedCoorList = [str(c) for c in coorNp_trans[0][0:3]]  
            updatedLine = ', '.join([lineSplit[0]]+updatedCoorList)+'\n'
            output += updatedLine
            line = f.readline()
        
        output += line
        for line in f.readlines():
            output += line
        
    with open(inpOutputFile, 'w') as f:
        f.write(output)
        print('Inp File Transformed by applying Transformation Matrix T')
        
    return inpOutputFile


def getESLabelFromInp(inpFile, ESname = 'ELSET=Volume2'):
    # Only works with Elements Connectivity - elset codefinition
    ES_labels = []
    #ES_Cyl = '*ELEMENT, type=C3D10, ELSET=Volume2'.upper()
    with open(inpFile,'r') as f:
        line = f.readline()
        while ESname.upper() not in line.upper():
            line = f.readline()
        line = f.readline()
        while not line.startswith('*'):
            lineSplit = line.split(', ')
            ES_labels.append(int(lineSplit[0]))
            line = f.readline()
    return ES_labels

def ESlayerFromNS(p,NS,ES_name):
    #Get the layer of elements that got at least a node within the node set NS 
    els = p.elements
    #ES_label = [el.label for i, nd in enumerate(NS.nodes) for el in NS.nodes[i].getElements()]
    elmts_label = [el.label for nd in NS.nodes for el in nd.getElements()]
    elmts_Seq = els.sequenceFromLabels(tuple(set(elmts_label)))
    ES_LayerFromNS = p.Set(elements=elmts_Seq, name=ES_name)
    
    return ES_LayerFromNS


def getElsetSectMat(mdl,p,prefix = 'SET_'):
    
    # Dict Elset Name -> Label of elmts within set
    ES_Name_Labels = dict()
    for es in p.sets.keys():
        if es.startswith(prefix):
            ES_Name_Labels[es]=[el.label for el in p.sets[es].elements]
            
    # Dict Section -> Elset Name
    Sect_ES_Name = dict()
    for sa in p.sectionAssignments:
        Sect_ES_Name[sa.sectionName] = sa.region[0]
    
    # Dict Mat -> Section
    Mat_Sect = dict()
    for s in mdl.sections.keys():
        Mat_Sect[mdl.sections[s].material] = s
        
    # Dict Elastic Modulus -> Mat
    EMod_Mat = dict()
    for m in  mdl.materials.keys():
        ElasticModulus = mdl.materials[m].elastic.table[0][0]
        EMod_Mat[ElasticModulus] = m
            
    return [ES_Name_Labels, Sect_ES_Name, Mat_Sect, EMod_Mat ]
    
    
    
def rhoFromE(E,law='Carter77'):
    # Get volumic mass from elastic modulus
    if law == 'Morgan2003' :
        paramE = dict()
        paramE[1] = {'a' : 0 , 'b' : 15520 , 'c' : 1.93}
        paramE[2] = {'a' : 0 , 'b' : 9965 , 'c' : 1.137}
        paramE[3] = {'a' : 0 , 'b' : 10714 , 'c' : 0.74}
        paramE['t1'] = 5280
        paramE['t2'] = 12260
        if E < paramE['t1'] :
            i = 1
        elif E > paramE['t2']:
            i = 3
        else :
            i = 2
        rho = ( (E-paramE[i]['a'])/paramE[i]['b'] )**( 1./paramE[i]['c'] )
    elif law == 'Carter77' :
        rho = (E/3790.0)**(1.0/3.0)
    return rho

def EFromRhoTBCMT(rho):
    # Get elastic modulus from volumic mass
    # Parameters were fitted from micro-FE models
    b = 6639.
    c = 0.7074
    BVTV = rho/1.9
    E = b*BVTV**c   
    return E
    
def updateCreateElset(p,dict_Elset):
    #dict_Elset a dictionnary with Elset names as keys and a list of associated
    #   element labels as values
    for es in dict_Elset.keys():
        if es in p.sets.keys():
            del p.sets[es]
        p.SetFromElementLabels(name = es , elementLabels = tuple(dict_Elset[es]))
        

def assignMat2TBCMTLayer(p,ES_Layer,dict_E_Sect,elmt_EMod,law='Carter77'):
    from abaqusConstants import MIDDLE_SURFACE,FROM_SECTION
    # dict_E_Sect # Dictionnary Elastic moudilus -> SectionName of CTMTB
    dict_E_Elmts = { E:[] for E in dict_E_Sect.keys()}
    E_CMTTB = sorted(dict_E_Sect.keys())
    for el in ES_Layer.elements:
        rho = rhoFromE( elmt_EMod[el.label], law)
        E = EFromRhoTBCMT(rho)
        Ebin = find_nearest(E_CMTTB,E)
        dict_E_Elmts[Ebin].append(el.label) 
    for E, sect in dict_E_Sect.iteritems():
        if dict_E_Elmts[E] :
            ES_curr = p.SetFromElementLabels(name = sect , elementLabels = tuple(dict_E_Elmts[E]))
            p.SectionAssignment(region=ES_curr, sectionName=sect,
                            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                            thicknessAssignment=FROM_SECTION)
    print('Assign section to TBCMT --> done.')
    return dict_E_Elmts
        

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]
        
def readData(fileName):
    data = dict()
    with open(fileName) as f:
        line = f.readline()
        ls = line.split(', ')
        data[ls[0]] = str(ls[1].split('"')[1]) #Prosthesis Name
        print(data[ls[0]])
        line = f.readline()
        ls = line.split(', ')
        data[ls[0]] = str(ls[1].split('"')[1]) #LegSide
        print(data[ls[0]])
        line = f.readline()
        ls = line.split(', ')
        data[ls[0]] = str(ls[1].split('"')[1]) #Alignment
        print(data[ls[0]])
        for line in f.readlines():
            ls = line.split(', ')
            data[ls[0]] = [float(nbr) for nbr in ls[1:]]
            if len(data[ls[0]])==1:
                data[ls[0]] = data[ls[0]][0]    

    T = data['T']
    data['T'] = [T[0:4],T[4:8],T[8:12],T[12:16]]

    LegSide = 1 if data['LegSide']=='R' else -1
    print(LegSide)

    SubjBW = {'BEN':122. , 'GIL':97., 'MOR':105., 'DES':76., 'GUI':68.,'IPO':102.,
          'ISN':63.,'LAU':76.,'LEC':80.,'LIM':91.,'MIL':70.,'NOR':83.,'ALO':99.,
          'TAR':84.,'THI':80.,'THI':80.,'TRO':78.,'VIT':68.,'MEL':80.}
    data['BW']=SubjBW
    if data['ImpltType'] == 'S4':
        data['Pt_Cyl'] = [0.002,29.421,41.404]
        data['Pt_Med_Implt'] = [LegSide*0.26*66.25, 22.0, 0.0]
        data['Pt_Lat_Implt'] = [-LegSide*0.26*66.25, 22.0, 0.0]
        data['O_Implt'] = [0.0 , 26.04, -6.975]
        print('Match with implant S4')
    elif data['ImpltType'] == 'S5':
        data['Pt_Cyl'] = [0.002,29.647,43.298]
        data['Pt_Med_Implt'] = [LegSide*0.26*75. , 22.325, 0.0]
        data['Pt_Lat_Implt'] = [-LegSide*0.26*75., 22.325, 0.0]
        data['O_Implt'] = [0.0 , 26.32, -7.05]
        print('Match with implant S5')
    elif data['ImpltType'] == 'S6' :
        data['Pt_Cyl'] = [0.002,29.647,43.298]
        data['Pt_Med_Implt'] = [LegSide*0.26*75. , 23.75, 0.0]
        data['Pt_Lat_Implt'] = [-LegSide*0.26*75. , 23.75, 0.0]
        data['O_Implt'] = [0.0 , 28.0, -7.5]
        print('Match with implant S6')
    else :
        print('No match with implant codes')
    data['Pt_Cyl'] = transformPt(data['Pt_Cyl'],data['T'])
    data['Pt_Med_Implt'] = transformPt(data['Pt_Med_Implt'],data['T'])
    data['Pt_Lat_Implt'] = transformPt(data['Pt_Lat_Implt'],data['T'])
    data['O_Implt'] = transformPt(data['O_Implt'],data['T'])
    data['rCyl'] = 4.5
    data['Ucyl'] = data['Nst']
    data['heightCyl'] = 22.0
    return data

def transformPt(pt,T):
    pt = list(pt) + [1.0,]
    coorNp = np.transpose([pt,])
    coorNp_trans = np.transpose(np.dot(T,coorNp))
    pt_transformed = [c for c in coorNp_trans[0][0:3]]
    return pt_transformed
       
def createCsysOnPartRA(p,O,X,Y,nameCsys):
    pt_origin = p.DatumPointByCoordinate(coords=tuple(O))
    pt_Xpos = p.DatumPointByCoordinate(coords=tuple(X))
    pt_Ypos = p.DatumPointByCoordinate(coords=tuple(Y))

    # Create Points
    PtYPos = p.DatumPointByMidPoint(point1=
        ImplantPart.vertices[41], point2=
        ImplantPart.vertices[19])
    
    PtYNeg = p.DatumPointByMidPoint(point1=
        ImplantPart.vertices[42], point2=
        ImplantPart.vertices[31])
    
    PtMiddleImplant = p.DatumPointByMidPoint(point1=
        ImplantPart.datums[PtYPos.id], point2=
        ImplantPart.datums[PtYNeg.id])
    
    
    # Create Csys
    CSys = p.DatumCsysByThreePoints(coordSysType=CARTESIAN
        , name=nameCsys, origin=p.datums[pt_origin.id], 
        point1=p.datums[PtYPos.id], point2=p.datums[PtYNeg.id])
    return CSys

def createLCFproperty(mdl,k,f):
    from abaqusConstants import FRACTION,OFF,PENALTY,ISOTROPIC,LINEAR,DEFAULT
    Iname = 'LC'+str(int(k))+'_f'+str(int(10*f))
    mdl.ContactProperty(Iname)
    mdl.interactionProperties[Iname].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((float(f), ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdl.interactionProperties[Iname].NormalBehavior(
        pressureOverclosure=LINEAR, contactStiffness=float(k), 
        constraintEnforcementMethod=DEFAULT)
    return Iname

def createCSysImplant(p,data):
    from abaqusConstants import CARTESIAN
    LegSide = 1 if data['LegSide'] == 'R' else -1
    O = np.array(data['O_Implt'])
    T = np.array(data['T'])
    Nxp = np.array(data['Nxp'])
    R = T[0:3,0:3];
    U_X = -R[:,0]*LegSide
    U_Y = np.cross(Nxp,U_X)
    Pt_X = O + 10.0*U_X
    Pt_Y = O + 10.0*U_Y
    CSysImplant = p.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='CSysImplant',
                                           origin=tuple(O), point1=tuple(Pt_X), point2=tuple(Pt_Y))
    return CSysImplant

def createCSysTibia(p,data):
    from abaqusConstants import CARTESIAN
    O = np.array(data['O_knee'])
    U_X = np.array(data['Xmech'])
    U_Y = np.array(data['Ymech'])
    Pt_X = O + 10.0*U_X
    Pt_Y = O + 10.0*U_Y
    CSysTibia = p.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='CSysTibia',
                                           origin=tuple(O), point1=tuple(Pt_X), point2=tuple(Pt_Y))
    return CSysTibia

def createCSysTibialPlateau(p,data):
    from abaqusConstants import CARTESIAN
    import numpy as np
    O = np.array(data['O_knee'])
    Ztp = np.array(data['Ztp'])
    U_Y = np.array(data['Xmech'])
    U_X = np.cross(U_Y,Ztp)
    U_X = U_X/np.linalg.norm(U_X)
    U_Y = np.cross(Ztp,U_X)
    Pt_X = O + 10.0*U_X
    Pt_Y = O + 10.0*U_Y
    CSysTibialPlateau = p.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='CSysTibialPlateau',
                                           origin=tuple(O), point1=tuple(Pt_X), point2=tuple(Pt_Y))
    return CSysTibialPlateau

def rodriguesRot(v,k,t):
    # See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    # v : Vector to be rotated around vector u by an angle t(theta) in degree
    v = np.array(v)
    k = np.array(k)
    t = np.radians(t)
    t1 = v*np.cos(t)
    t2 = np.cross(k,v)*np.sin(t)
    t3 = np.dot(k,np.dot(k,v))*(1-np.cos(t))
    vrot = t1 + t2 + t3
    return vrot

def importMusclesForces(data):
    LegSide = 1 if data['LegSide'] == 'R' else -1
    forceList = ['fPL','fSM','fBF','fST','fSL','fTA','fTP']
    forces = { f:dict() for f in forceList}
    steps = ['GC_00','GC_15','GC_30','GC_50','GC_65','GC_95','CU','SU1','SU2']
    stepsGC = [0.0,15.0,30.0,50.0,65.0,95,-1000.0,-1000.0,-1000.0]
    stepsAngle = {'CU' : 90.0, 'SU1' : 52.5 , 'SU2' : 20.0}
    #
    # Forces values from Winby 2009 and Trinker 2019 and PL (Patellar Ligament) Force Calculated with Defrate 2007
    # Forces values for FTP from Masouros 2010 and Ellis 1984 [9.,14.]
    GC_percent = list(np.arange(0.,100.01,2.5))
    KneeFlex_GC = [10.14, 13.10, 15.38, 18.24, 20.34, 21.60, 21.88, 21.40, 20.35,
                19.04, 17.91, 16.79, 15.26, 14.20, 13.31, 12.90, 12.87, 13.28,
                13.46, 13.89, 14.89, 16.75, 19.92, 25.13, 32.32, 40.90, 49.10,
                55.68, 60.04, 62.00, 61.74, 59.45, 55.39, 49.63, 42.44, 34.07,
                24.90, 16.49, 10.08, 7.20, 8.53]
    # Force between the
    # van Eijden 1987 :
    KneeFlex_QPL = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.]
    PLF_QF = [1.07, 1.05, 1, 0.95, 0.87, 0.8, 0.7, 0.6, 0.53, 0.56, 0.58, 0.61, 0.65]
    sagAnglePL = [23., 13., 6.5, -3.5, -9.2]  # -------------
    corAnglePL = [17.7, 7.7, 4.2, 1.5, 1.8] # Defrate et Al. 2007
    anglePL = [0.0,30.0,60.0,90.0,110.0]    # -------------
    #
    # Angle in degree with Z axis of the tibia, in the XZ plan and YZ plan
    forceDir = { 'fPL' : [2,[]],
                 'fSM' : [1,[400,-10,-50]],
                 'fBF' : [1,[300,60,0]],
                 'fST' : [1,[420,30,-70]],
                 'fSL' : [0,[-300,0,0]],
                 'fTA' : [0,[-300,10,5]],
                 'fTP' : [0,[-300,15,10]] }
    # Forces intensity in terme of BW
    forceInt = [0.11, 0.07, 0.04, 0.07, 0.00, 0.14, 0.00, #---
                0.73, 0.25, 0.22, 0.11, 0.48, 0.21, 0.29, # |
                0.06, 0.09, 0.14, 0.05, 1.18, 0.12, 0.27, # | Tinler
                0.77, 0.00, 0.16, 0.00, 2.95, 0.00, 0.34, # |   2019
                0.32, 0.00, 0.04, 0.00, 0.00, 0.09, 0.00, # |
                0.00, 0.35, 0.28, 0.18, 0.00, 0.22, 0.00, #---
                3.42, 0.55, 0.41, 0.14, 0.00, 0.00, 0.00, # Navacchia 2019
                1.62, 0.45, 0.83, 0.15, 0.00, 0.00, 0.00, # Rasnick 
                1.32, 0.53, 0.62, 0.14, 0.00, 0.00, 0.00] #     2016

    # Dir in tibial Coordinate Syst
    Zmech = np.array([0.,0.,1.])
    Ymech = np.array([0.,1.,0.])
    Xmech = np.array([1.,0.,0.])
    i=0
    for j,s in enumerate(steps):
        # Get Flexion angle
        if stepsGC[j] < -999. :
            flexA = stepsAngle[s]
        else :
            flexA = np.interp(stepsGC[j],GC_percent,KneeFlex_GC)
        #
        for f in forceList :
            # Get force direction for tibial muscles
            if forceDir[f][0] == 0:
                coefs = forceDir[f][1]
                #
                u_01 = coefs[0]*Zmech
                u_02 = coefs[1]*Ymech
                u_03 = coefs[2]*LegSide*Xmech
                u_0 = u_01 + u_02 + u_03 
                #
                u = u_0/np.linalg.norm(u_0)
                #
                force = forceInt[i]*u
                forces[f].update({s:force})
                i += 1
                #
            # Get force direction for hamstring muscles
            elif forceDir[f][0] == 1:
                coefs = forceDir[f][1]
                u_01 = coefs[0]*Zmech
                u_02 = coefs[1]*Ymech
                u_03 = coefs[2]*LegSide*Xmech
                u_0 = u_01 + u_02 + u_03 
                #
                rotAxisDir = -LegSide*Ymech
                u_rot = rodriguesRot(u_0,rotAxisDir,flexA)
                u = u_rot/np.linalg.norm(u_rot)
                #
                force = forceInt[i]*u
                forces[f].update({s:force})
                i += 1
                #
            # Get force direction for patellar ligament
            elif forceDir[f][0] == 2:
                ratioPLQ = np.interp(flexA,KneeFlex_QPL,PLF_QF)
                #print(ratioPLQ)
                alpha = np.radians( np.interp(flexA,anglePL,sagAnglePL) )
                #print(np.degrees(alpha))
                beta = np.radians( np.interp(flexA,anglePL,corAnglePL) )
                #print(np.degrees(beta))
                u = np.array( [LegSide*np.sin(alpha), np.sin(beta), 0] )
                u[2] = np.sqrt(1.000-np.sum(u**2))
                #
                #print(flexA)
                #print(u)
                force = ratioPLQ*forceInt[i]*u
                forces[f].update({s:force})
                i += 1
                #
    #
    return forces

def generateMusclesLoads(mdl,stepList,sfc,data,forces,BW,CSys,frr):
    # frr : Force Reduction Ratio
    import regionToolset
    from abaqusConstants import DISTRIBUTING,UNIFORM,CUBIC,QUADRATIC,ON,OFF,WHOLE_SURFACE,STRUCTURAL
    RA = mdl.rootAssembly
    CSysDatum = RA.datums[CSys.id]
    Pt_TT = np.array(data['Pt_TT'])
    Pt_SM = np.array(data['Pt_SM'])
    Pt_SM0 = Pt_SM - 7.*np.array(data['Ymech'])
    Pt_SM1 = Pt_SM + 7.*np.array(data['Ymech'])
    Pt_ST = np.array(data['Pt_ST'])
    Pt_Fib = np.array(data['Pt_Fib'])
    Pt_Sol = np.array(data['Pt_Sol'])
    Pt_TA = np.array(data['Pt_TA'])
    Pt_TP = np.array(data['Pt_TP'])
    forcesAPts = {'fPL':Pt_TT ,'fSM':[Pt_SM0,Pt_SM1,Pt_SM,Pt_SM]
                  ,'fBF':Pt_Fib, 'fST':Pt_ST, 'fSL':Pt_Sol,
                  'fTA':Pt_TA, 'fTP':Pt_TP}
    forcesRadius = {'fPL':9,'fSM':8,'fBF':8,'fST':8,'fSL':20,'fTA':20,'fTP':20}
    rPts = dict()
    reg_RPts = dict()
    for f,v in forcesAPts.iteritems() :
        if isinstance(v,list):
            subPt = []
            for i,pt in enumerate(v[0:-1]):
                RP = RA.ReferencePoint(point=tuple(pt))
                RA.features.changeKey(
                    fromName=RP.name , toName=f+str(i+1))
                reg_RP = regionToolset.Region(referencePoints=
                                        (RA.referencePoints[RP.id], ))
                #
                mdl.Coupling(name='Pt_'+f+str(i), controlPoint = reg_RP ,
                              surface= sfc, influenceRadius= forcesRadius[f],
                              couplingType=DISTRIBUTING, weightingMethod=QUADRATIC, localCsys=None,
                              u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF, adjust=True)
                #
                subPt.append(RA.referencePoints[RP.id])
            rPts[f] = RA.ReferencePoint(point=tuple(v[-1]))
            RA.features.changeKey(
                fromName=rPts[f].name , toName=f)
            reg_RPts[f] = regionToolset.Region(referencePoints = (RA.referencePoints[rPts[f].id], ))
            subRefPoints= tuple(subPt)
            print(subRefPoints)
            region2=RA.Set(referencePoints=subRefPoints, name='Set_'+f)
            mdl.Coupling(name='Pt_'+f, controlPoint=reg_RPts[f],
                          surface=region2, influenceRadius=WHOLE_SURFACE, 
                          couplingType=STRUCTURAL, localCsys=CSysDatum,
                          u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
            RA.regenerate()
        else :
            rPts[f] = RA.ReferencePoint(point=tuple(v))
            RA.features.changeKey(
                fromName=rPts[f].name , toName=f)
            reg_RPts[f] = regionToolset.Region(referencePoints=(RA.referencePoints[rPts[f].id], ))
            mdl.Coupling(name='Pt_'+f, controlPoint= reg_RPts[f] ,
                surface= sfc, influenceRadius=forcesRadius[f], 
                couplingType=DISTRIBUTING, weightingMethod=CUBIC, localCsys=CSysDatum, 
                u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
            RA.regenerate()
    #
    for s in stepList :
        # Create loads
        for f, v in forces.iteritems() :
            F = frr * 9.81 * BW * forces[f][s] + 0.01 #Forces cannot be 0.0, frr : Force Reduction Ratio
            if s == stepList[0] : # Create forces at first step
                mdl.ConcentratedForce(name=f,
                                       createStepName=s, region = reg_RPts[f], 
                                       cf1 = F[0], cf2 = F[1], cf3= F[2],
                                       distributionType = UNIFORM, 
                                       field='', localCsys=CSysDatum)
            else : #Update Force for other steps
                mdl.loads[f].setValuesInStep(stepName=s,
                      cf1 = F[0], cf2 = F[1], cf3= F[2])   
        #
        RA.regenerate()
    #
    return rPts



