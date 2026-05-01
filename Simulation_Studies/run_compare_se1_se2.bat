@echo off
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::install('C:/Users/stat-user/iCloudDrive/Desktop/EBMR/EBMRalgorithmFast4', quiet=TRUE, upgrade='never')" > "C:\Users\stat-user\iCloudDrive\Desktop\EBMR\Simulation_Studies\compare_se1_se2_install.txt" 2>&1
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "C:\Users\stat-user\iCloudDrive\Desktop\EBMR\Simulation_Studies\compare_se1_se2.R" > "C:\Users\stat-user\iCloudDrive\Desktop\EBMR\Simulation_Studies\compare_se1_se2_output.txt" 2>&1
