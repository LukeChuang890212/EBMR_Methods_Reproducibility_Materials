Set objShell = CreateObject("WScript.Shell")
Set objFSO = CreateObject("Scripting.FileSystemObject")

Dim outFile
outFile = "C:\Users\stat-user\iCloudDrive\Desktop\EBMR\compare_se1_se2_output.txt"

' Step 1: Install package
Dim installCmd
installCmd = """C:\Program Files\R\R-4.5.2\bin\Rscript.exe"" -e ""devtools::install('C:/Users/stat-user/iCloudDrive/Desktop/EBMR/EBMRalgorithmFast4', quiet=TRUE, upgrade='never')"""
objShell.Run "cmd.exe /c " & installCmd & " >> """ & outFile & """ 2>&1", 0, True

' Step 2: Run comparison script
Dim runCmd
runCmd = """C:\Program Files\R\R-4.5.2\bin\Rscript.exe"" ""C:\Users\stat-user\iCloudDrive\Desktop\EBMR\Simulation_Studies\compare_se1_se2.R"""
objShell.Run "cmd.exe /c " & runCmd & " >> """ & outFile & """ 2>&1", 0, True

MsgBox "Done! Output saved to: " & outFile
