This R script allows to calculate closed-system 230Th-U ages.
Simply, set the working directory “path”, name of the sample to solve “sample_name” (can solve several samples if they share the same character string as defined by the user), and run the code.
Data need to be in a tab-separated .txt file “IoliteExport_All_Integrations.txt” (although the file name can be changed), with sample names, measured (230Th/238U), (234U/238U) activity ratios  and their 2s errors in columns named respectively “X”, “Th230_U238_CORR”, “U234_U238_CORR”, “Th230_U238_CORR_Int2SE” and “U234_U238_CORR_Int2SE”.
You can also set the number of optimisation for each sample, “nbit”, and the lower and upper bounds, respectively “lowerbound” and “upperbound”, for log10 of the age (in log10(yr)) and initial (234U/238U).

The code returns ages (in kyr), initial (234U/238U) and their 2s errors in a comma-separated .csv file with “sample_name” as file name.
