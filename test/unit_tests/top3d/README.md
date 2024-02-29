### Notes
1. The timer values being displayed at the end when we run the top3d.jl file is not accurate since various File IO operations are happening which are increasing the overall time taken by a function.
2. The values of all the MATLAB variables are being stored in different files in ```matlab_values``` directory.
3. The values of all the Julia variables are being stored in different files in ```julia_values``` directory.
4. Corresponding variable values of MATLAB and Julia are being compared in the Julia file ```top3d.jl``` and if they are not matching then their names are printed in the file ```not_matching.txt```. And if they are matching then their names are printed in the file ```matching.txt```.
5. As of now, the tests start breaking after we enter into the optimization loop so many variables used within the optimization loop are present in the file ```not_matching.txt```.
6. The MATLAB code runs for 146 iterations while the Julia code is running for 200 iterations so after 146 iterations in Julia, the variables stop getting printed in the ```not_matching.txt``` file since there is no corresponding MATLAB value to compare to after that iteration.
7. We have also displayed the results of MATLAB iterations at the end for the ease of comparision.
