@echo off
cd /d "C:\Users\stat-user\iCloudDrive\Desktop\EBMR"
echo === Step 1: Installing EBMRalgorithmFast4 ===
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::install('EBMRalgorithmFast4', quiet=TRUE, upgrade='never')"
echo === Step 1 done, exit code: %ERRORLEVEL% ===
echo === Step 2: Running compare_se1_se2.R ===
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "C:\Users\stat-user\iCloudDrive\Desktop\EBMR\Simulation_Studies\compare_se1_se2.R"
echo === Step 2 done, exit code: %ERRORLEVEL% ===
