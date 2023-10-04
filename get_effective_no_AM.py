import sys,fileinput
import numpy as np

counter = 0

sample_size = sys.argv[2]

fin = open(sys.argv[1],"r")
fin2 = fin.readlines()
fin3 = [row.rstrip().split() for row in fin2]
data_matrix = np.array(fin3)
data_matrix2 = data_matrix.transpose()

for idx,line in enumerate(data_matrix2): # inferences file
        if (idx+1)%2 != 0:
                hap1 = line
        else:
                counter +=1
                hap2 = line
                prev1 = hap1[0] ; prev2 = hap2[0]
                for i in range(1,len(hap1)):
                     if (hap1[i] == prev1) and (hap2[i] == prev2):
                          pass
                     elif (hap1[i] == prev2) and (hap2[i] == prev1):
                          pass
                     else:
                          counter +=1
                     prev1 = hap1[i] ; prev2 = hap2[i]
no_of_test = round(counter/float(sample_size),0)
print("Effective no. of tests: ",no_of_test)


