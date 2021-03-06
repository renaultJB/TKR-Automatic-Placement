def matlab_module_exists(module_name):
    import subprocess
    
    try:
        __import__(module_name)
    except ImportError:
        tmp = subprocess.check_output(["where", "matlab"])
        tmp = tmp.decode("utf-8")
        matlabRoot = tmp.split('bin\\matlab')[0]
        matlabEngineRoot = matlabRoot+'extern\\engines\\python'
        matlabEngineSetup = matlabEngineRoot + '\\setup.py'
        subprocess.call(['python', 'setup.py' , 'install'], cwd = matlabEngineRoot)
        try:
            __import__(module_name)
        except ImportError:
            print('Matlab Engine Module is not installed')
            print('A file to explain how to install it will be created')
            f = open('MatlabEngineInstall.txt','w')
            lines = ['open a cmd prompt as admin and execute the following commands :']
            lines.append('cd '+ matlabEngineRoot)
            lines.append('python setup.py install')
            lines.append('_______________________')
            lines.append('https://fr.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html')
            f.write('\n'.join(lines))
            f.close()
            print(' The \'MatlabEngineInstall.txt\' file was created')
            wait = raw_input("PRESS ENTER TO CONTINUE.")
        else:
            print('Matlab Engine Module for python has been installed and imported')
    else:
        print('Matlab Engine Module for python already installed and imported')


def find_ProsthFiles(directory,Pname,Ptype,CmtOnly=False):
    import os
    Ptype = int(Ptype)
    os.chdir(directory)
    files_list = []
    prosthType = 'Prosthesis'+str(Ptype)
    if CmtOnly :
        filesPref = ['Cut_','Implant'+str(Ptype)]
        suffix = '.stp' 
    else :
        filesPref = ['C_','Cut_','Implant'+str(Ptype)]
        suffix = '.stp'
    for path, subdirs, files in os.walk(directory):
        for name in files:
            files_list.append(os.path.join(path, name))
    filesList = []
    for filePref in filesPref :
        A = [fdir for fdir in files_list  if prosthType in fdir and Pname in fdir and filePref in fdir and suffix in fdir]
        filesList.append(str(A[0]))
    return filesList
