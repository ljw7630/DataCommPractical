DataCommPractical
=================

A visual studio project for data communication practical

1. compile the project as Release version
2. go to Realease folder, put hadsund.dhm, mhad435.10, receiver_preprocess.py run.bat and plot_result.py into the folder
3. in command window, run "python batch_run.py" to generate experiment results.

Details of file information in DataCommPractical/DataCommunicationPractical/
1. batch_run.py: Generate five outputs for each number of receiver pairs. Then number of receiver pairs is vary from [8, 10, 12, .., 24].
2. Calculator.h: Calculate all results using the error output from localization.cpp. The error evaluating methods are provided in the paper.
3. FilePrinter.h: Print the result to files for further processing.
4. hadsund.dhm: The terrain profile.
5. Localization.cpp: The origin program that given by Eamonn that implement localization for single transmitter.
6. mhad435.10m: The measurements for the whole terrain.
7. plot_result.py: A python script to generate experiment result diagrams.
8. recevier_preprocess.py: A python script to randomly generate arbitrary number of receiver pairs.
9. ReceiverReader.h: Read the output of receiver_preprocess.py, which is a number of receiver pairs, load it as C++ vector.
10. run.bat: run the whole program for once.
11. DataCommunicationPractical.sln: The visual studio solution file
12. README.md: This readme file.
13. VisualStudio.gitignore: ignore file list for git source control.
