import os
import subprocess
from script_checks import matlab_module_exists,find_ProsthFiles
matlab_module_exists('matlab.engine')
import matlab.engine


#### Paramaters ####
PatientCode = str.upper(raw_input("Enter Patient Code (eg : AM_R) : "))
if '_L' in PatientCode :
    RightKnee = 0
elif '_R' in PatientCode :
    RightKnee = 1
else :
    LeftKnee = input("What knee side ? (1 for Right, 0 for Left) : ")
    RightKnee = 0
    
LongStem = input("Use of long Stem ? (1 for yes) : ")
#Size = str.upper(raw_input("Enter Component Size (eg : S5) : "))
MeshExist = input("Meshes already created ? (1 for Yes): ")
PosFiles = input("Generate GMSH .pos files ? (1 for Yes): ")
if PosFiles!=1 :
    PosFiles=0

if MeshExist!=1 :
    ElmtSize2D = input("Enter Element Size in for Surface Mesh (in mm) : ")
    ElmtSize3D = input("Enter Element Size in for Volume Mesh (in mm) : ")
TypeProth = input("Prosthesis Type ? (1 for Old, 2 for recent): ")

print("\n###########################################\n######### END OF USER INPUT ###############\n###########################################\n")

cwd = os.getcwd()
subprocess.call(["cd", cwd], shell=True)
if not os.path.exists(cwd+"\\Output\\"+PatientCode):
    os.makedirs(cwd+"\\Output\\"+PatientCode)
    print("Patient folder created.")
else :
    print("Patient folder already exists.")

if MeshExist==1 :

    print("Mesh already exist, won't be generated again")
    os.rename(cwd+"\\Output\\"+PatientCode+"\\Tibia_" + PatientCode + "_2D.msh",cwd+"\\"+"Tibia_" + PatientCode + "_2D.msh")
    os.rename(cwd+"\\Output\\"+PatientCode+"\\Tibia_" + PatientCode + "_3D.msh",cwd+"\\"+"Tibia_" + PatientCode + "_3D.msh")


else :
        
    for dim in ["Initial_MESH_2D","Initial_MESH_3D"]:
               
        
        f = open(dim+'.geo','r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace("Esize2D = ;" , "Esize2D = " + str(ElmtSize2D) + ";")
        newdata = newdata.replace("Esize3D = ;" , "Esize3D = " + str(ElmtSize3D) + ";")
        newdata = newdata.replace("Tibia.stp" , "Tibia_" + PatientCode + ".stp")
        newdata = newdata.replace("Tibia_2D.msh" , "Tibia_" + PatientCode + "_2D.msh")
        finaldata = newdata.replace("Tibia_3D.msh" , "Tibia_" + PatientCode + "_3D.msh")


        f = open(dim + '_' + PatientCode +'.geo','w')
        f.write(finaldata)
        f.close()

        subprocess.call(["F:\\Implantation\\gmsh\\gmsh.exe", dim + "_"+PatientCode +".geo"], shell=True)
        os.remove(dim + "_"+ PatientCode +'.geo')
        print("Mesh generated")

eng = matlab.engine.start_matlab()

# [-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0]
# [-4.5,-1.5,1.5,4.5]
for alpha in [0.0] :
    if TypeProth == 2:
        T, Tanat, ML_Width , AP_Width, ProstName = eng.PositionProth2(PatientCode,alpha,RightKnee,LongStem,PosFiles,nargout=5)
    elif TypeProth == 1 :
        T, Tanat, ML_Width , AP_Width, ProstName = eng.PositionProth1(PatientCode,alpha,RightKnee,LongStem,PosFiles,nargout=5)
    print(ML_Width)    

    f = open('ScriptFreeCAD.py','r')
    filedata = f.read()
    f.close()
    
    paths2STPfiles = find_ProsthFiles(os.getcwd(),ProstName,TypeProth)

    newdata = filedata.replace("NAME",PatientCode)
    newdata = newdata.replace("PATH2CEMENT",str(paths2STPfiles[0]))
    newdata = newdata.replace("PATH2CUT",str(paths2STPfiles[1]))
    newdata = newdata.replace("PATH2IMPLANT",str(paths2STPfiles[2]))
    newdata = newdata.replace("CWD",cwd)
    newdata = newdata.replace("ALPHA",str(alpha))
    newdata = newdata.replace("AP_WIDTH",str(AP_Width))
    newdata = newdata.replace("ML_WIDTH",str(ML_Width))
    newdata = newdata.replace("NEW_PLACEMENT_MATRIX",T)
    newdata = newdata.replace("NEW_PLACEMENT_ANAT_MATRIX",Tanat)
    if LongStem != 1:
        newdata = newdata.replace("_WLS","")
    finaldata = newdata.replace("SIZE",ProstName)

    f = open('ScriptFreeCAD_'+PatientCode +'.py','w')
    f.write(finaldata)
    f.close()
    subprocess.call(["FreeCADcmd", 'ScriptFreeCAD_'+PatientCode +'.py'], shell=True)

eng.quit()
##os.remove('ScriptFreeCAD_'+PatientCode +'.py')
os.rename(cwd+"\\"+"Tibia_" + PatientCode + "_2D.msh",cwd+"\\Output\\"+PatientCode+"\\Tibia_" + PatientCode + "_2D.msh")
os.rename(cwd+"\\"+"Tibia_" + PatientCode + "_3D.msh",cwd+"\\Output\\"+PatientCode+"\\Tibia_" + PatientCode + "_3D.msh")


quit()
