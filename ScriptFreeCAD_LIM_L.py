# -*- coding: utf-8 -*-
import FreeCAD
import Part
import Import
import Mesh

App.ActiveDocument=None

doc = FreeCAD.newDocument("Noname") 
App.setActiveDocument("Noname")
App.ActiveDocument = App.getDocument("Noname")




TibiaStp = Part.read("C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Tibia_LIM_L.stp")
CementStp = Part.read("C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\data\\Prosthesis1\\Cement\\C_S5_15mm.stp")
CutProthStp = Part.read("C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\data\\Prosthesis1\\Cut\\Cut_S5_15mm.stp")
ProthStp = Part.read("C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\data\\Prosthesis1\\Implant\\Implant1_S5.stp")
CondylesMsh = Mesh.read("C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Condyles4Mimics.stl")

doc.recompute()

Tibia = doc.addObject("Part::Feature","Tibia")
Tibia.Shape = TibiaStp

Tibia2 = doc.addObject("Part::Feature","Tibia2")
Tibia2.Shape = TibiaStp

Cement = doc.addObject("Part::Feature","Cement")
Cement.Shape = CementStp

CutProth = doc.addObject("Part::Feature","CutProth")
CutProth.Shape = CutProthStp

CutProth2 = doc.addObject("Part::Feature","CutProth2")
CutProth2.Shape = CutProthStp

Proth = doc.addObject("Part::Feature","Proth")
Proth.Shape = ProthStp

Condyles = doc.addObject("Mesh::Feature","Condyles")
Condyles.Mesh = CondylesMsh
matTsfrm=FreeCAD.Matrix()
matTsfrm.scale(82.1289385714,59.3507637344,1.0)
CdlMesh=doc.Condyles.Mesh.copy()
CdlMesh.transform(matTsfrm)
doc.recompute() 

newplace=FreeCAD.Matrix(0.89050932,0.35461302,-0.28503116,0.00000000,0.37375576,-0.92742328,0.01388135,0.00000000,-0.25942202,-0.11889351,-0.95841773,0.00000000,65.52826636,55.28115614,-571.24396378,1.00000000)


newplaceAnat=FreeCAD.Matrix(0.95083972,0.28503721,-0.12106866,70.93020786,0.28871915,-0.95731625,0.01366896,37.42834364,-0.11200483,-0.04795183,-0.99255002,-561.87311883,0.00000000,0.00000000,0.00000000,1.00000000)


doc.recompute() 

CutProth.Placement = newplace
CutProth2.Placement = newplace
Proth.Placement = newplace
Cement.Placement = newplace
CdlMesh.Placement = FreeCAD.Placement(newplaceAnat)
doc.recompute()
Mesh.show(CdlMesh)
doc.recompute()

App.activeDocument().addObject("Part::Cut","TibiaCutted")
App.activeDocument().TibiaCutted.Base = Tibia
App.activeDocument().TibiaCutted.Tool = CutProth


App.activeDocument().addObject("Part::MultiCommon","TibialCut")
App.activeDocument().TibialCut.Shapes = [Tibia2, CutProth2,]


doc.recompute()

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("TibiaCutted"))
Import.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\Tibia_LIM_L_cutted_a-10.0.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("TibialCut"))
Import.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\Tibia_LIM_L_cut_a-10.0.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("Proth"))
Import.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\Implant_S5_LIM_L_a-10.0.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("Cement"))
Import.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\C_S5_LIM_L_a-10.0.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("CutProth"))
Import.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\CutProth_S5_LIM_L_a-10.0.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("Mesh"))
import Mesh
Mesh.export(__objs__,u"C:\\Users\\Jean-Baptiste\\Documents\\These\\Methodes\\AACS_Knee\\TKR-Automatic-Placement\\Output\\LIM_L\\Condyle_LIM_L_a-10.0.stl")
del __objs__

quit()


