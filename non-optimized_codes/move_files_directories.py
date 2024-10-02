import pandas as pd
import os
import glob
import shutil

#bad files
name_bad_files=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_files_elec/test',header=None)
cut_name_bad_files=name_bad_files[0].apply(lambda x: x[-40:])

#All files 6000 
directory_all_files='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_10PeV_Fe_900m'

name_all_files= pd.DataFrame({'files':glob.glob(os.path.join(directory_all_files,'*'))})
cut_name_all_files=name_all_files['files'].apply(lambda x: x[-40:])

#bad files in 6000 
directory_bad_in_all=[]
for i in range(0,len(cut_name_bad_files)):
    directory_bad_in_all.append(name_all_files['files'].iloc[cut_name_all_files[cut_name_all_files.str.contains(cut_name_bad_files.iloc[i])].index[0]])

for i in range(0,len(directory_bad_in_all)):
    shutil.copy(directory_bad_in_all[i], '/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_files_elec/bad_files')
