import os 

cwd = os.getcwd()

inpFiles = [f for f in os.listdir(cwd) if f.endswith('.inp')]

with open('Batch_HeND.bat', 'w') as outfile:
    outfile.write('(\n')
    outfile.write('cd ' + cwd + '\n')
    for inpf in inpFiles :
         outfile.write('abq6142 job={} cpus=8 interactive \n'.format(inpf[:-4]))
    outfile.write(')')
