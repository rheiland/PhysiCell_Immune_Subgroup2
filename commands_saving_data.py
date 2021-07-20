import sys
import os
import glob
import numpy as np
from pyMCDS_cells import pyMCDS_cells
from pyMCDS import pyMCDS
import matplotlib.pyplot as plt

"""""
timetotal = 72;
timetotal = 3;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

for tcount = 1:timetotal
   % clf
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    deltax = abs(MCDS.mesh.X(1,1,k)-MCDS.mesh.X(1,2,k));
    deltay = abs(MCDS.mesh.Y(1,1,k)-MCDS.mesh.Y(2,1,k));

    
    virion(tcount) = sum(sum(MCDS.continuum_variables(7).data(:,:,k)))*deltax*deltay*20;%virion
    IFN(tcount) = sum(sum(MCDS.continuum_variables(2).data(:,:,k)))*deltax*deltay*20;%debris
    cytokine(tcount) = sum(sum(MCDS.continuum_variables(3).data(:,:,k)))*deltax*deltay*20;%proinflamcytokine virion
    chemokine(tcount) = sum(sum(MCDS.continuum_variables(4).data(:,:,k)))*deltax*deltay*20;%chemokine
    debris(tcount) = sum(sum(MCDS.continuum_variables(5).data(:,:,k)))*deltax*deltay*20;%debris
    ROS(tcount) = sum(sum(MCDS.continuum_variables(6).data(:,:,k)))*deltax*deltay*20;%debris
    
    total_cells = length(MCDS.discrete_cells.live_cells);
    
    CD8(tcount) = length(find( MCDS.discrete_cells.metadata.type == 3)); %CD8    
    neutrophils(tcount) = length(find( MCDS.discrete_cells.metadata.type == 5)); %neutrophils
    DC(tcount) = length(find( MCDS.discrete_cells.metadata.type == 6)); %DCs
    CD4(tcount) = length(find( MCDS.discrete_cells.metadata.type == 7)); %CD4
    
    inactivated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==0);
    activated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==1);
    
    macrophagesinactive(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),inactivated_immune)); %macs inactive
    macrophagesactive(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),activated_immune)); %macs active

    dead(tcount) = length(intersect(MCDS.discrete_cells.dead_cells,find(MCDS.discrete_cells.metadata.type == 1)));
    those_withvirus = find(MCDS.discrete_cells.custom.Vnuc>=1/8000);
    those_notdead_withvirus = those_withvirus;%intersect(MCDS.discrete_cells.live_cells,those_withvirus);
    those_notantiviral = find(MCDS.discrete_cells.custom.antiviral_state<0.5);
    infected(tcount) = length(intersect(those_notdead_withvirus,those_notantiviral));
    uninfected(tcount) = length(find( MCDS.discrete_cells.metadata.type == 1))-dead(tcount)-infected(tcount);
    antiviral_cells(tcount) = length(find(MCDS.discrete_cells.custom.antiviral_state>0.5));
          
 end
 
 save('simulated_data.mat','virion','IFN','cytokine','chemokine','debris','ROS','CD8','neutrophils','DC','CD4','macrophagesinactive','macrophagesactive','dead','infected','uninfected','antiviral_cells')
"""""

argc = len(sys.argv)-1
print("# args=",argc)

#data_dir = 'output'
if (argc < 1):
#  data_dir = int(sys.argv[kdx])
  print("Usage: provide output subdir")
  sys.exit(-1)

kdx = 1
data_dir = sys.argv[kdx]
print('data_dir = ',data_dir)
os.chdir(data_dir)
xml_files = glob.glob('output*.xml')
os.chdir('..')
xml_files.sort()
print('xml_files = ',xml_files)

ds_count = len(xml_files)
print("----- ds_count = ",ds_count)
#mcds0 = pyMCDS_cells(xml_files[0], data_dir)
mcds0 = pyMCDS(xml_files[0], data_dir)
print("substrates: ", mcds0.get_substrate_names())

mcds = [pyMCDS(xml_files[i], data_dir) for i in range(ds_count)]

tval = np.linspace(0, mcds[-1].get_time(), ds_count)
print('tval= ',tval)

# DC(tcount) = length(find( MCDS.discrete_cells.metadata.type == 6)); %DCs

# # mac,neut,cd8,DC,cd4,Fib
# Macs
yval4 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 4) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )
print("Macs: ",yval4)

# count Neuts
yval5 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 5) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count CD8
yval6 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 3) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count DC
yval7 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 6) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )
print("DC: ",yval7)

# count CD4
yval8 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 7) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count Fibroblasts
yval9 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 8) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

#f = np.array([ (mcds[idx].get_concentrations('virion')).sum()  for idx in range(ds_count)] )

print("VTEST for 0th: ",mcds[0].get_concentrations('VTEST'))
