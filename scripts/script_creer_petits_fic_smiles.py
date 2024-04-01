#!/usr/bin/env python
import csv
import sys
import os

file_data = sys.argv[1]

with open("scripts/script_smi_to_mol.sh", "w") as f_script :
	with open(file_data, 'r') as f :
		csvreader = csv.reader(f)
		i = 0
		for row in csvreader :
			if row[1] != "name" :
				petit_smile = row[6]
				name = "_".join(row[1].split(" "))
				if name+".mol" not in os.listdir("data/smi_files_reduit/") :
					with open("data/smi_files_reduit/%s.smi"%name, 'w') as f :
						f.write("%s %s\n"%(petit_smile, name))
					i += 1
				
					f_script.write("obgen data/smi_files_reduit/%s.smi > data/smi_files_reduit/%s.mol\n"%(name, name))

