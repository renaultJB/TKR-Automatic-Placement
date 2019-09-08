import os
import subprocess
from script_checks import matlab_module_exists,find_ProsthFiles
matlab_module_exists('matlab.engine')
import matlab.engine


#### Paramaters ####
SubjectCode = str.upper(raw_input("Enter Subject Code (eg : AM_R) : "))
if '_L' in SubjectCode :

    RightKnee = 0
elif '_R' in SubjectCode :
    RightKnee = 1
else :
    LeftKnee = input("What knee side ? (1 for Right, 0 for Left) : ")
    RightKnee = 0
    
# LongStem = input("Use of long Stem ? (1 for yes) : ")
LongStem = 0

MeshExist = input("Meshes already created ? (1 for Yes): ")
##PosFiles = input("Generate GMSH .pos files ? (1 for Yes): ")
##if PosFiles!=1 :
##    PosFiles=0
PosFiles=0

if MeshExist!=1 :
    ElmtSize2D = 0.5
TypeProth = input("Prosthesis Type ? (1 for Sym, 2 for alt Sym, 3 for ASym): ")

print("\n###########################################\n######### END OF USER INPUT ###############\n###########################################\n")

cwd = os.getcwd()
subprocess.call(["cd", cwd], shell=True)
if not os.path.exists(cwd+"\\Output\\"+SubjectCode):
    os.makedirs(cwd+"\\Output\\"+SubjectCode)
    print("Subject folder created.")
else :
    print("Subject folder already exists.")

if MeshExist==1 :

    print("Mesh already exist, won't be generated again")
    os.rename(cwd+"\\Output\\"+SubjectCode+"\\Tibia_" + SubjectCode + ".msh",cwd+"\\"+"Tibia_" + SubjectCode + ".msh")
    os.rename(cwd+"\\Output\\"+SubjectCode+"\\DistTibia_" + SubjectCode + ".msh",cwd+"\\"+"DistTibia_" + SubjectCode + ".msh")


else :
        
    for TibiaPart in ["Tibia_","DistTibia_"]:
        
        f = open('Initial_MESH.geo','r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace("Esize2D = ;" , "Esize2D = " + str(ElmtSize2D) + ";")
        newdata = newdata.replace("Tibia.stp" , TibiaPart + SubjectCode + ".stp")
        finaldata = newdata.replace("Tibia.msh" , TibiaPart + SubjectCode + ".msh")

        fName = 'Initial_MESH_' + TibiaPart + SubjectCode +'.geo'
        f = open(fName,'w')
        f.write(finaldata)
        f.close()

        subprocess.call(["gmsh", fName], shell=True)
        os.remove(fName)
        print("Mesh generated")
        
os.chdir(cwd+'\\CoreFunctions')

eng = matlab.engine.start_matlab()

os.chdir(cwd)
# [-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0]
# [-4.5,-1.5,1.5,4.5]
for alpha in [0.0 , 100.0] :
    os.chdir(cwd+'\\CoreFunctions')
    beta = 6.0
    T, Tanat, ML_Width , AP_Width, ProstName = eng.PlacementTI(SubjectCode,alpha,TypeProth,LongStem,beta,nargout=5)
    print(ML_Width) 
    
    os.chdir(cwd)

    f = open('ScriptFreeCAD_NoCMT.py','r')
    filedata = f.read()
    f.close()

    paths2STPfiles = find_ProsthFiles(os.getcwd(),ProstName,TypeProth,True)

    newdata = filedata.replace("NAME",SubjectCode)
    newdata = newdata.replace("PATH2CUT",str(paths2STPfiles[0]))
    newdata = newdata.replace("PATH2IMPLANT",str(paths2STPfiles[1]))
    newdata = newdata.replace("CWD",cwd)
    newdata = newdata.replace("ALPHA",str(alpha))
    newdata = newdata.replace("AP_WIDTH",str(AP_Width))
    newdata = newdata.replace("ML_WIDTH",str(ML_Width))
    newdata = newdata.replace("NEW_PLACEMENT_MATRIX",T)
    newdata = newdata.replace("NEW_PLACEMENT_ANAT_MATRIX",Tanat)
    if LongStem != 1:
        newdata = newdata.replace("_WLS","")
    newdata = newdata.replace("SIZE",ProstName)
    finaldata = newdata.replace("\\","\\\\")

    f = open('ScriptFreeCAD_NoCMT_'+SubjectCode +'.py','w')
    f.write(finaldata)
    f.close()
    subprocess.call(["FreeCADcmd", 'ScriptFreeCAD_NoCMT_'+SubjectCode +'.py'], shell=True)

    if alpha%1 == 0 :
        os.rename(cwd+"\\"+"Output_" + SubjectCode + "_alpha" + str(int(alpha))+ ".txt",cwd+"\\Output\\"+SubjectCode+"\\Output_" + SubjectCode + "_alpha" + str(int(alpha))+ ".txt")
        os.rename(cwd+"\\"+"Centrality_" + SubjectCode + "_alpha" + str(int(alpha))+ ".txt",cwd+"\\Output\\"+SubjectCode+"\\Centrality_" + SubjectCode + "_alpha" + str(int(alpha))+ ".txt")
    else :
        os.rename(cwd+"\\"+"Output_" + SubjectCode + "_alpha" + str(alpha)+ ".txt",cwd+"\\Output\\"+SubjectCode+"\\Output_" + SubjectCode + "_alpha" + str(alpha)+ ".txt")
        os.rename(cwd+"\\"+"Centrality_" + SubjectCode + "_alpha" + str(alpha)+ ".txt",cwd+"\\Output\\"+SubjectCode+"\\Centrality_" + SubjectCode + "_alpha" + str(alpha)+ ".txt")
        
            

eng.quit()


os.remove('ScriptFreeCAD_'+SubjectCode +'.py')
os.rename(cwd+"\\"+"Tibia_" + SubjectCode + ".msh",cwd+"\\Output\\"+SubjectCode+"\\Tibia_" + SubjectCode + ".msh")
os.rename(cwd+"\\"+"DistTibia_" + SubjectCode + ".msh",cwd+"\\Output\\"+SubjectCode+"\\DistTibia_" + SubjectCode + ".msh")


quit()
