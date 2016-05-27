from qutip import *
import os

folder = os.path.join('twirl','test_sim','data')
fileprefix = 'run' # or 'iteration'
num_of_files = 4
filenames = [os.path.join(fileprefix+str(i+1)) for i in range(num_of_files)]
fullstates_filename = 'final_states'
header_filename = 'header.txt'

#result_states = qload(os.path.join('twirl','test_sim','data','run1'))

def prepare_bloch():
    db=Bloch()
    colors = ["g"]#,"r","g","#CC6600"]#,"r","g","#CC6600"]
    db.point_color = "r"
    db.point_marker = ['o']
    return db

def plot_bloch(db,result_states): 
    for t in range(0,len(result_states),10):
        db.add_states(result_states[t].ptrace(0), kind='point')
        db.add_states(result_states[t].ptrace(1), kind='point')    


for file in filenames:
    result_states = qload(os.path.join(folder,file))
    db = prepare_bloch()
    plot_bloch(db,result_states)
    #db.show()
    db.save(os.path.join(folder,file+'.pdf'))
    
print('Bloch spheres saved')
final_states = qload(os.path.join(folder,fullstates_filename))
with open(os.path.join(folder,header_filename)) as f:
    header = dict()
    header.update(yaml.load(f))

Plotter(final_states, header,filetype='.pdf',display=False)
print('Plotter run')