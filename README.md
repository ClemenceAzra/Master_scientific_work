# Master's scientific work of Clemence Azra, 2024

The work has the title "Estimation of the EAS axis position and the arrival direction based on the Cherenkov light characteristic in the SPHERE experiment".

All of the used codes are in the repository /optimized_codes/

  main.py : dashboard of the algorithm. Here you have to input the directory names, characteristics of EAS events (telescope, altitude, energy, nuclei). 
            Each dirname of events need to have the same structure: /Data/{telescope}/{type_analyse}/{type_analyse}_{En}PeV_{nuclei}_{H}m/
            A loop browses the directory files. If the file is not an event file, it is skipped.
            You must load the py files : Initial_values, Openfiles, Functions_algorithm, Conditions.
            Conditions skip the events based on event selection criteria.
            By using fnc_algorithm.name_of_function (see functions_algorithm.py), load the results (axis, angles...).
            
  open_files.py : 
  conditions.py
  functions_algorithm.py
  initial_values.py
  open_files.py
  results.py
