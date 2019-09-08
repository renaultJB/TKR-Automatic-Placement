# -*- coding: utf-8 -*-
import FreeCAD
import Part
import Import
import Mesh

App.ActiveDocument=None

doc = FreeCAD.newDocument("Noname") 
App.setActiveDocument("Noname")
App.ActiveDocument = App.getDocument("Noname")




TibiaStp = Part.read("CWD\Tibia_NAME.stp")
CutProthStp = Part.read("PATH2CUT")
ProthStp = Part.read("PATH2IMPLANT")
CondylesMsh = Mesh.read("CWD\Condyles4Mimics.stl")

doc.recompute()

Tibia = doc.addObject("Part::Feature","Tibia")
Tibia.Shape = TibiaStp

Tibia2 = doc.addObject("Part::Feature","Tibia2")
Tibia2.Shape = TibiaStp

CutProth = doc.addObject("Part::Feature","CutProth")
CutProth.Shape = CutProthStp

CutProth2 = doc.addObject("Part::Feature","CutProth2")
CutProth2.Shape = CutProthStp

Proth = doc.addObject("Part::Feature","Proth")
Proth.Shape = ProthStp

Condyles = doc.addObject("Mesh::Feature","Condyles")
Condyles.Mesh = CondylesMsh
matTsfrm=FreeCAD.Matrix()
matTsfrm.scale(ML_WIDTH,AP_WIDTH,1.0)
CdlMesh=doc.Condyles.Mesh.copy()
CdlMesh.transform(matTsfrm)
doc.recompute() 

NEW_PLACEMENT_MATRIX

NEW_PLACEMENT_ANAT_MATRIX

doc.recompute() 

CutProth.Placement = newplace
CutProth2.Placement = newplace
Proth.Placement = newplace
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
Import.export(__objs__,u"CWD\Output\NAME\Tibia_NAME_cutted_aALPHA.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("TibialCut"))
Import.export(__objs__,u"CWD\Output\NAME\Tibia_NAME_cut_aALPHA.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("Proth"))
Import.export(__objs__,u"CWD\Output\NAME\Implant_SIZE_WLS_NAME_aALPHA.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("CutProth"))
Import.export(__objs__,u"CWD\Output\NAME\CutProth_SIZE_WLS_NAME_aALPHA.step")
del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Noname").getObject("Mesh"))
import Mesh
Mesh.export(__objs__,u"CWD\Output\NAME\Condyle_NAME_aALPHA.stl")
del __objs__

quit()


