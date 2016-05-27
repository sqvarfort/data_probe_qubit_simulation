from qutip import *
import os
import yaml
from Plotter import *

info = os.path.join('twirl','twirl_config.yml')
lind_args = dict()
if type(info) == str:
    with open(info) as cfile:
        lind_args.update(yaml.load(cfile))
elif type(config) == dict:
    lind_args.update(info)
else:
    print "info is wrong type"

#path = os.path.join('twirl','displaced')
#num_of_files = 4
path = os.path.join(lind_args.get('folder'),lind_args.get('subfolder'))
num_of_files = lind_args.get('runs')
data_directory = 'data'
save_directory = 'plots'
data_path = os.path.join(path,data_directory)
save_path = os.path.join(path,save_directory)
fileprefix = 'run' # or 'iteration'

filenames = [os.path.join(fileprefix+str(i+1)) for i in range(num_of_files)]
fullstates_filename = 'final_states'
header_filename = 'header.txt'

# Make sure save directory exists
if not os.path.exists(save_path):
    os.mkdir(save_path)

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
    result_states = qload(os.path.join(data_path,file))
    db = prepare_bloch()
    plot_bloch(db,result_states)
    #db.show()
    db.save(os.path.join(save_path,file+'.pdf'))
    
print('Bloch spheres saved')
final_states = qload(os.path.join(data_path,fullstates_filename))
with open(os.path.join(data_path,header_filename)) as f:
    header = dict()
    header.update(yaml.load(f))

# correct header to save to right directory
header['folder'] = path
header['subfolder'] = save_directory
Plotter(final_states, header,filetype='.pdf',display=False)
print('Plotter run')