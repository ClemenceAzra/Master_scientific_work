# Master's scientific work of Clemence Azra, 2024

The work has the title "Estimation of the EAS axis position and the arrival direction based on the Cherenkov light characteristic in the SPHERE experiment".

## Getting EAS_angles_sphere

The latest release version is available and can be installed using

```
pip install EAS_angles_sphere
```

https://github.com/ClemenceAzra/Master_scientific_work/tree/main/EAS_angles_sphere
The latest source code is available from the github repository at [(https://github.com/ClemenceAzra/Master_scientific_work/tree/main/EAS_angles_sphere
)]([https://github.com/spacepy/spacepy](https://github.com/ClemenceAzra/Master_scientific_work/tree/main/EAS_angles_sphere
).

## Algorithm 

All of the used codes are in the repository EAS_angles_sphere/

- main.py : dashboard of the algorithm
  - You have to input the directory names, characteristics of EAS events (telescope, altitude, energy, nuclei). 
  - Each dirname of events need to have the same structure: /Data/{telescope}/{type_analyse}/{type_analyse}_{En}PeV_{nuclei}_{H}m/
  - A loop browses the directory files. If the file is not an event file, it is skipped.
  - You must load the py files : Initial_values, Openfiles, Functions_algorithm, Conditions.
  - Conditions skip the events based on event selection criteria.
  - By using fnc_algorithm.name_of_function (see functions_algorithm.py), load the results (axis, angles...).
            
- initial_values.py: all initial values (telescope geometry, thresholds...)
  - Some are constant and not changeable (telescope geometry) : #Unchangeable
  - Others are can be changed (thresholds, response...) : #Changeable

- open_files.py : create a dataframe of size: (N bins of response * N PMT)
  
- functions_algorithm.py: all of the created functions to finally calculate the arrival angles
  - sum_of_impulse: return the summed impulse over all channels in the response
  - intensity: return the maximum amplitude of the summed impulse (N photons), and the intensity (multiplicity max amplitude / background) 
  - diapason: return the new temporal diapason of the impulse in the response: list of N photons, list of time (ns), list of background values. len = number of bins
  - translation_snow_mos: return the translated front (xi, yi, ti) from mosaic (mm) to snow (m)
  - axis: return the EAS axis on snow (xi, yi) and the distance from the telescope axis in mm
  - length: return the duration of the impulse in ns
  - amplification: return the amplified "diapason" signal by the sliding window method
  - noise_pmt: return the background in "piedestal" for each channel
  - front_up_noise: return the channels with the maximum impulse above the background
  - DFS: return the channels with EAS impulse
  - neighbors: return the channels with at least N neighbors
  - angles: return both zenith (theta) and azimthal angles (phi) in rad, and front parameters a0, a1, a2 

- conditions.py : all the event selection criteria

- results.py: calculate the error values, plot figures, save data
  - delta: return angles errors in deg
  - save_results: save values (e.g. values of angles, axis...)
  - fig_sum_impulse: plot the summed impulse in the response
  - front: plot the front (xi, yi, ti) on the snow

## Examples of algorithm operation

Open the repository Examples/

- 10 files of EAS events are loaded in the directory: data:SPHERE-3:only_sig:only_sig_10PeV_P_500m

- the file main_and_results.py contains the dashboard opening these files and printed results for selected events


