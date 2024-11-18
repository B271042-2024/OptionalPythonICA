#!/usr/bin/python3

#--------------------IMPORT MODULES

import subprocess
import os
import re
import shutil
import sys

#--------------------STATING PURPOSE

print("\n- - This pipeline involves multiple analyses - -\n\n1	clustalO [clustalo -i {input_file} -o {inp_fil}.aln --force --threads=30]\n2	plotcon [plotcon -sequences {inp_fil}.aln -winsize 15 -graph png]\n3	patmatmotifs [patmatmotifs -sequence {input_file} -outfile {inp_fil}.motif]\n")

#--------------------GET PATH/TO/FASTA WITH INPUT()

print("Please insert your FASTA file that containing a set of protein sequences.")
input_file = input("path/to/input_fasta_file: ")	#input: file.fa
qc1_1 = input("Filter out sequences tagged 'PREDICTED' y/n: ")
qc1_2 = input("Filter out sequences tagged 'LOW QUALITY PROTEIN' y/n: ")

#--------------------GET FILE NAME (WITHOUT PATH/EXTENSION)

inp_fil_1 = input_file.split('.')[0]	#grab the text before '.'
inp_fil_ext = input_file.split('.')[1]	#grab extension
inp_fil_0 = inp_fil_1.split('/')
inp_fil = inp_fil_0[-1]	#grab all text after the last '/'
inp_fil_path = '/'.join(inp_fil_0[:-1])	#path

#--------------------QC0: if file exist

if not os.path.isfile(input_file):
	print('File does not exist\n')
	sys.exit()

#--------------------MOVE TO WORKING DIRECTORY & CREATE NEW DIRECTORY

if len(inp_fil_path) > 0:
	os.chdir(inp_fil_path)

if os.path.isdir(f"{inp_fil}_protein/supplementary"):
	shutil.rmtree(f"{inp_fil}_protein")
os.makedirs(f"{inp_fil}_protein/supplementary")

#copy fasta file
shutil.copy(input_file, f"{inp_fil}_protein/{inp_fil}_copy.{inp_fil_ext}")

#--------------------QC1: USER OPTION TO REMOVE BASED ON THE NAME OF EACH SEQUENCE

#create a new edited file
edi_file = f"{inp_fil}_protein/{inp_fil}_edited.{inp_fil_ext}"
os.system(f"touch {edi_file}")
fil_file = f"{inp_fil}_protein/{inp_fil}_filtered.{inp_fil_ext}"

#merge all lines under each id
with open(input_file, "r") as orifile, open(edi_file, "w") as newfile:
	seq = ""
	for i in orifile:
		i = i.strip()
		if re.search(r"^>", i):
			if seq:
				newfile.write(seq + "\n\n")
			newfile.write(i + "\n")
			seq = ""
		else:
			seq+=i
	newfile.write(seq)

#if y for both predicted and lower quality
if qc1_1.lower() == "y" and qc1_2.lower() == "y":
	print("\n1. Filtering out sequences tagged 'PREDICTED' and 'LOW QUALITY PROTEIN'")
	os.system(f"touch {fil_file}")
	with open(edi_file, "r") as edifile, open(fil_file, "w") as filtfile:
		lines = iter(edifile)
		for i in lines:
			if not re.search(r"\bPREDICTED\b", i) and not re.search(r"\bLOW QUALITY PROTEIN\b", i):
				if re.search(r"^>", i):
					filtfile.write(i)
					nextline = next(lines)
					filtfile.write(nextline + "\n")

	print(f"[FILE]	{input_file}")
	print(f"[OUTFILE]	{edi_file}")
	print(f"[OUTFILE]	{fil_file}")

if qc1_1.lower() == "y" and not qc1_2.lower() == "y":
	print("\n1. Filtering out sequences tagged 'PREDICTED'")
	os.system(f"touch {fil_file}")
	with open(edi_file, "r") as edifile, open(fil_file, "w") as filtfile:
		lines = iter(edifile)
		for i in lines:
			if not re.search(r"\bPREDICTED\b", i):
				if re.search(r"^>", i):
					filtfile.write(i)
					nextline = next(lines)
					filtfile.write(nextline + "\n")
	print(f"[FILE]  {input_file}")
	print(f"[OUTFILE]        {edi_file}")
	print(f"[OUTFILE]        {fil_file}")

if not qc1_1.lower() == "y" and qc1_2.lower() == "y":
	print("\n1. Filtering out sequences tagged 'LOW QUALITY PROTEIN'")
	os.system(f"touch {fil_file}")
	with open(edi_file, "r") as edifile, open(fil_file, "w") as filtfile:
		lines = iter(edifile)
		for i in lines:
			if not re.search(r"\bLOW QUALITY PROTEIN\b", i):
				if re.search(r"^>", i):
					filtfile.write(i)
					nextline = next(lines)
					filtfile.write(nextline + "\n")
	print(f"[FILE]  {input_file}")
	print(f"[OUTFILE]        {edi_file}")
	print(f"[OUTFILE]        {fil_file}")

if not qc1_1.lower() == "y" and not qc1_2.lower() == "y":
	print("\n1. DEFAULT: No Filter is SET")

#--------------------ALIGNMENT WITH CLUSTALO [file.fasta]

print("\n2. Performing alignment [software: ClustalO]")

if not os.path.isfile(fil_file):
	print(f"[FILE] {input_file}")
	subprocess.call(f"clustalo -i {input_file} -o {inp_fil}_protein/supplementary/{inp_fil}.aln --force --threads=30", shell=True)	#align file.fa via clustalo
else:
	print(f"[FILE] {fil_file}")
	subprocess.call(f"clustalo -i {fil_file} -o {inp_fil}_protein/supplementary/{inp_fil}.aln --force --threads=30", shell=True)
print(f"[OUTFILE]	{inp_fil}_protein/supplementary/{inp_fil}.aln")

#--------------------RE-FORMAT ALIGNED SEQUENCES WITH EMBOSS:SHOWALIGN [file.aln]

#original
print("\n3. Aligning sequences in parallel [software: showalign, EBLOSUM62]: original")
subprocess.call(f"showalign -sequence {inp_fil}_protein/supplementary/{inp_fil}.aln -show=a -outfile {inp_fil}_protein/supplementary/{inp_fil}_showalign_ori.aln", shell=True)
print(f"[FILE]	{inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]	{inp_fil}_protein/supplementary/{inp_fil}_original.aln")

#identity: show characters that align
print("\n4. Aligning sequences in parallel [software: showalign, EBLOSUM62]: identity")
subprocess.call(f"showalign -sequence {inp_fil}_protein/supplementary/{inp_fil}.aln -show=i -outfile {inp_fil}_protein/supplementary/{inp_fil}_identity.aln", shell=True)
print(f"[FILE]  {inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]       {inp_fil}_protein/supplementary/{inp_fil}_identity.aln")

#non-identities: show characters that don't align
print("\n5. Aligning sequences in parallel [software: showalign, EBLOSUM62]: non-identity")
subprocess.call(f"showalign -sequence {inp_fil}_protein/supplementary/{inp_fil}.aln -show=n -outfile {inp_fil}_protein/supplementary/{inp_fil}_nonidentity.aln", shell=True)
print(f"[FILE]  {inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]       {inp_fil}_protein/supplementary/{inp_fil}_nonidentity.aln")

#similarity
print("\n6. Aligning sequences in parallel [software: showalign, EBLOSUM62]: similarities (based on function/biochemical properties)")
subprocess.call(f"showalign -sequence {inp_fil}_protein/supplementary/{inp_fil}.aln -show=s -outfile {inp_fil}_protein/supplementary/{inp_fil}_similarity.aln", shell=True)
print(f"[FILE]  {inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]       {inp_fil}_protein/supplementary/{inp_fil}_similarity.aln")

#dissimilarities/less strict: show characters that differ in function, biochemical properties [based on EBLOSUM62]
print("\n7. Aligning sequences in parallel [software: showalign, EBLOSUM62]: dissimilarities (based on function/biochemical properties)")
subprocess.call(f"showalign -sequence {inp_fil}_protein/supplementary/{inp_fil}.aln -show=d -outfile {inp_fil}_protein/supplementary/{inp_fil}_dissimilarity.aln", shell=True)
print(f"[FILE]  {inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]       {inp_fil}_protein/supplementary/{inp_fil}_dissimilarity.aln")

#--------------------CONSERVATION PLOT WITH EMBOSS:PLOTCON [file.aln]

print("\n8. Plotting conservation graph [software: plotcon]")
subprocess.call(f"plotcon -sequences {inp_fil}_protein/supplementary/{inp_fil}.aln -winsize 10 -graph png", shell=True)
subprocess.call(f"plotcon -sequences {inp_fil}_protein/supplementary/{inp_fil}.aln -winsize 10 -graph data", shell=True)
os.system(f"mv plotcon.1.png {inp_fil}_protein/supplementary")
os.system(f"mv plotcon1.dat {inp_fil}_protein/supplementary")
os.rename(f"{inp_fil}_protein/supplementary/plotcon.1.png", f"{inp_fil}_protein/supplementary/{inp_fil}_plotcon.png")
os.rename(f"{inp_fil}_protein/supplementary/plotcon1.dat", f"{inp_fil}_protein/supplementary/{inp_fil}_plotcon.data")
print(f"Renamed plotcon.1.png to {inp_fil}_plotcon.png")
print(f"Renamed plotcon1.dat to {inp_fil}_plotcon.data")
print(f"[FILE]	{inp_fil}_protein/supplementary/{inp_fil}.aln")
print(f"[OUTFILE]	{inp_fil}_plotcon.png")
print(f"[OUTFILE]	{inp_fil}_plotcon.data")

#--------------------GET MOTIFS WITH EMBOSS:PROSITE [file.fasta]

print("\n9. Searching for motifs")
if not os.path.isfile(fil_file):
	subprocess.call(f"patmatmotifs -sequence {input_file} -outfile {inp_fil}_protein/supplementary/{inp_fil}.motif", shell=True)
	print(f"[FILE]	{input_file}")
else:
	subprocess.call(f"patmatmotifs -sequence {fil_file} -outfile {inp_fil}_protein/supplementary/{inp_fil}.motif", shell=True)
	print(f"[FILE]	{fil_file}")
print(f"[OUTFILE]	{inp_fil}_protein/supplementary/{inp_fil}.motif")

#--------------------BUILD A REPORT

#create file
report = f"{inp_fil}_protein/Report_{inp_fil}.txt"
subprocess.call(f"touch {report}",shell=True)
print("\nAnalysis done")
with open(report, "w") as openreport, open(f"{inp_fil}_protein/supplementary/{inp_fil}_plotcon.data", "r") as opendata, open(f"{inp_fil}_protein/supplementary/{inp_fil}.motif") as openmotif:
	openreport.write("REPORT SUMMARY\n\n")
	openreport.write("File Paths:\n")
	openreport.write(f"[FOLDER] {inp_fil}_protein\n")
	openreport.write(f" |\n")
	openreport.write(f"  ----- [FILE] {inp_fil}_copy.fasta\n")
	openreport.write(f"  ----- [FILE] {inp_fil}_edited.fasta\n")
	openreport.write(f"  ----- [FILE] {inp_fil}_filtered.fasta\n")
	openreport.write(f"  ----- [FILE] Report_{inp_fil}.txt\n")
	openreport.write(f"  ----- [FOLDER] supplementary\n")
	openreport.write(f"         |\n")
	openreport.write(f"          ----- [FILE] {inp_fil}.aln\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_original.aln\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_identity.aln\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_nonidentity.aln\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_dissimilarity.aln\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_plotcon.png\n")
	openreport.write(f"          ----- [FILE] {inp_fil}_plotcon.data\n")
	openreport.write(f"          ----- [FILE] {inp_fil}.motif\n\n\n")
	openreport.write("DESCRIPTION:\n")
	openreport.write(f'''Analysis focuses on protein conservation in a certain taxonomic group, comparing a group of organisms listed in the input fasta file.
1	A copy of fasta file was made.	[output: {inp_fil}_copy.fasta]
2	Sequences were cleaned to produce 2 lines per sequence and then, filtered. (optional: with tags 'PREDICTED' and/or 'LOW QUALITY PROTEIN')
	[output: {inp_fil}_edited.fasta (cleaning/re-formatted: 1 line=title, 1 line=sequences) & {inp_fil}_filtered.fasta (filtered sequences)]
3	Sequences were aligned with ClustalO.	[output: {inp_fil}.aln]	(MODE: --force --threads=30)
4	Aligned sequences were re-formatted using EMBOSS showalign to produce files below:
		- Re-formats.	[output: {inp_fil}_original.aln] (MODE: -show=a)
		- Shows identical sequences. (based on sequence comparison)	[output: {inp_fil}_identity.aln] (MODE: -show=i)
		- Generates non-identical report. (based on sequence comparison) 	[output: {inp_fil}_nonidentity.aln] (MODE: -show=n)
		- Generates similarity report. (based on function/biochemical)	[output: {inp_fil}_similarity.aln] (MODE: -show=s)
		- Generates dissimilarity report. (based on function/biochemical)	[output: {inp_fil}_dissimilarity.aln] (MODE: -show=d)
5	Protein conservation plot was generated using EMBOSS plotcon.	[output: {inp_fil}_plotcon.png & {inp_fil}_plotcon.data]
6	Motif search was done against PROSITE database using EMBOSS patmatmotif.	[output: {inp_fil}.motif]\n\n\n''')

	openreport.write("DATA OVERVIEW\n\n")
	openreport.write("Evolutionary conservation details:\n")
	var_data = ['Points', 'XminA', 'Xmin', 'ScaleXmin']
	for line in opendata:
		for item in var_data:
			if item in line:
				openreport.write(line)
				break

	#openreport.write("\n\nMotifs\n")
	#var_motif = ['Sequence','Motif']
	#for line in openmotif:
	#	for item in var_motif:
	#		if item in line:
	#			openreport.write(line)
	#			break

sys.exit()
