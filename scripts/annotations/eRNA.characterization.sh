## log of the commands used to analyze the TNE characterization step 


# eRNA.characterize.xls include both plus and minus strand TNEs (unfiltered)

awk '$24==3{print $0}' ./input_files/characterization/eRNA.characterize.xls | wc -l 
# class 3 = 64,423

awk '$24==2{print $0}' ./input_files/characterization/eRNA.characterize.xls | wc -l 
# class 2 = 19,139

awk '$24==1{print $0}' ./input_files/characterization/eRNA.characterize.xls | wc -l 
# class 1 = 26,035

# total number of TNEs : 109,597

# eRNA.plus.characterize.xls 

awk '$24==3{print $0}' ./input_files/characterization/eRNA.plus.characterize.xls | wc -l 
# class 3 = 42008

awk '$24==2{print $0}' ./input_files/characterization/eRNA.plus.characterize.xls | wc -l 
# class 2 = 10644

awk '$24==1{print $0}' ./input_files/characterization/eRNA.plus.characterize.xls | wc -l 
# class 1 = 14130

# total number of pus TNEs : 66,782

# eRNA.minus.characterize.xls 

awk '$24==3{print $0}' ./input_files/characterization/eRNA.minus.characterize.xls | wc -l 
# class 3 = 22415

awk '$24==2{print $0}' ./input_files/characterization/eRNA.minus.characterize.xls | wc -l 
# class 2 = 8495

awk '$24==1{print $0}' ./input_files/characterization/eRNA.minus.characterize.xls | wc -l 
# class 1 = 11905

# total number of pus TNEs : 42,815