from random import uniform
from random import randint
import argparse

file = open("text.vcf", "w")
total_induced_mutations = 0
total_lines = 0
sample_number = -1

file_path = input("FILE PATH HERE: ")

with open(file_path) as search:
    for line in search:
        if "#CHROM" == line[:6]:
            file.write(line)
            headers = line.split()
            samples = headers[headers.index('FORMAT')+1:]
            sample_number = len(samples)
            print("\n\n" +str(sample_number)+" SAMPLES DETECTED: "+ str(samples) + "\n\n")
        elif "##" != line[:2]:
            #elements = line.split()[-1*sample_number:]
            #elements contains all the 0|0 1|0 0|1 1|1 stuff
            elements = line.split()
            vt = elements[7].split(";")
            b = randint(0,2)
            if b==1:
                vt[len(vt)-1] = "VT=INDEL"
            elements[7] = ";".join(vt)

            for x in range(0,sample_number-1):
                if uniform(0,9) < 2: #induce a mutation
                    mutation_type = randint(0, 2)
                    mut_type_str = ""
                    
                    if mutation_type == 0:
                        mut_type_str = "1|0"
                    elif mutation_type == 1:
                        mut_type_str = "0|1"
                    else:
                        mut_type_str = "1|1"
                    
                    #line = line[:len(line)-((x+1)*3+x)] + mut_type_str + line[(len(line)-(x+1)*3+x)+3:]
                    elements[len(elements)-x-1] = mut_type_str

                    total_induced_mutations = total_induced_mutations + 1
                    total_lines = total_lines + 1
                    file.write("\t".join(elements) + "\n")
                else:
                    total_lines = total_lines + 1
                    file.write(line)
            
        else:
            file.write(line)

print("TOTAL INDUCED MUTATIONS: "+ str(total_induced_mutations) + " TOTAL LOCATIONS POSSIBLE: "+str(total_lines) + "\n\n")


