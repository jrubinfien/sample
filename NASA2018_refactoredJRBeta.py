#!/usr/bin/env python3


#to-do: 
	#QC at end -- output chart
	#format all outputs (krona, stacked bar chart, qc) into single HTML
	#color grade krona by q-score
	#
	#eventually, stop using subprocess.run (unnecessary forking, leaking so much memory):
		#see if i can write a bash script to work simultaneously, in sync w/ python

import pickle


import os, subprocess, sys, time, re, testPlots
from tqdm import tqdm
import pandas as pd

#this is the directory that I move all my seq files to from the minIT
path = r"/Users/Julian/Desktop/NASA/testP/"

#pull list of sub directories
folder = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x))]

work_run = "BLANK"

def choose_run():
	print("These are your runs: \n \n")
	for s in folder:
		print(s)
	print("\n \n")
	
	#allow user to select which run folder they want to work out of (forms the base of everything else)	
	
	global work_run
	work_run = str(input("Please enter the name of the run to be analyzed (or enter \"exit\" to exit): ")).strip()
	
	global pathToRun
	pathToRun = path + work_run	
	
	#allow user to specify a sub-directory in their run folder for pulling fastqs
		#e.g. if you basecall manually
	global fastqOrigin
	fastqOrigin = []
	fqOriginRequest = str(input(" \nIf you would like to specify a folder of fastq files, please enter the absolute path with slashes now (or press enter to skip)"))
	
	#allow user to specify location for the sequencing summary, if it's somewhere weird
	global seq_sum
	seq_sum = str(input(" \nIf you would like to specify a sequencing summary file for filtering, please enter the absolute path with slashes now (or press enter to skip)"))

	#look for sequencing summary file in main directory or in 
	#a sequencing_summary folder, or exit
	if len(seq_sum) == 0:
		if os.path.exists(pathToRun + "/sequencing_summary"):
			for x in os.listdir(pathToRun + "/sequencing_summary/"):
				if x.endswith('sequencing_summary.txt'):
					seq_sum = pathToRun + "/sequencing_summary/" + x
					break
		for item in os.listdir(pathToRun):
			if item.endswith('sequencing_summary.txt'):
				seq_sum = pathToRun + "/" + item
				break
		if len(seq_sum) == 0:
			print("\n \n CAN'T FIND SEQUENCING SUMMARY. PLEASE ENTER A PATH NEXT TIME.")
			sys.exit()

		print("\n Using " + seq_sum + " as the sequencing summary file. \n")

	#list all directories user specifies for pulling fastqs from
	if len(fqOriginRequest) > 0:
		fastqOrigin = [x.strip() for x in re.findall(r"[\w']+",fqOriginRequest)]

	#if user doesn't specify, assume default locations
	if len(fqOriginRequest) == 0:
		fastqOrigin.append(pathToRun + "/fastq_fail/")
		fastqOrigin.append(pathToRun + "/fastq_pass/")

	#if the user wants to nanofilter demultiplexed reads, they can input specs
	global nFiltInfo
	nFiltInfo = []
	nFiltInfoReq = str(input(" \nIf you would like to Nanofilter your reads, please specify your minimum and maximum desired read lengths (you must enter 2 numbers) separated by a comma (or press enter to skip)"))
	if len(nFiltInfoReq) > 0:
		nFiltInfo = [x.strip() for x in re.findall(r"[\w']+",nFiltInfoReq)]

	#if the main directory exists, go ahead, otherwise try again or leave
	if any(work_run == s for s in folder):
		print("\n \n Got it. Merging fastq files...")

	elif work_run == "exit":
		sys.stdout.write('\n \n Program terminated. \n \n')
		sys.exit(0)		
	else:
		print("\n You entered: ")
		print(work_run)
		print("\n")
		return choose_run()




def start():
	choose_run()

	#create a log file in the main directory of all operatiions
	log = open(pathToRun + "/log.txt", mode='a+')

	print("\n \n \n")
	# ##### starts below

	#make the subdirectory "working_files"
	global MakeWkDirCmd
	MakeWkDirCmd = "mkdir " + pathToRun + "/working_files"
	subprocess.run([MakeWkDirCmd], shell=True)

	#merges all fastq files from the default fastq origin directories (fail and pass) 
	#or from directories specified by user and then sends the output to working_files
	CatFilesCmdBeginning = "cat "
	for i in fastqOrigin:	
		CatFilesCmdBeginning = CatFilesCmdBeginning + i + "*.fastq "
	CatFilesCmd = CatFilesCmdBeginning + "> " + pathToRun + "/working_files/merged_all.fastq"
	log.write("[" + str(time.time()) + "]" + "merge fastq files:   " + CatFilesCmd + "\n \n")
	subprocess.run([CatFilesCmd], shell=True)

	#counts all reads in merged fastq file
	print("\n \n Merged all fastq files. Counting total reads...")

	countTotalReads = "awk '{s++}END{print s/4}' " + pathToRun + "/working_files/merged_all.fastq"
	totalReads = os.popen(countTotalReads).read()



	print(totalReads.strip("\n") + " total read(s) exist")
	log.write("[" + str(time.time()) + "]" + totalReads.strip("\n") + " total read(s) exist" + "\n \n")


	#this makes the demux subdirectory and starts demuxing
	MakeDemuxDir = "mkdir " + pathToRun + "/working_files/demux"
	subprocess.run([MakeDemuxDir], shell=True)

	print("\n \n Starting to demultiplex. This may take a while. \n \n ")
	log.write("\n \n Starting to demultiplex. This may take a while. \n \n ")

	DemuxCmd = "nohup qcat -f " + pathToRun + "/working_files/merged_all.fastq -b " + pathToRun + "/working_files/demux --kit RAB204 --trim | tee -a " + pathToRun + "/log.txt"
	subprocess.run([DemuxCmd], shell=True)
	print("\n \n Qcat demultiplexing is complete. Starting to count barcoded reads...\n \n ")


	#now will recursively count and print reads qcat assigned to each barcode
	demuxFolder = pathToRun + "/working_files/demux/"
	mybcodes = [x for x in os.listdir(demuxFolder) if os.path.isfile(os.path.join(demuxFolder,x))]
	
	for i in tqdm(mybcodes):
		cmd = "awk '{s++}END{print s/4}' " + pathToRun + "/working_files/demux/" + i
		numb = os.popen(cmd).read()
		print(i[:-6] + " has " + numb.strip("\n") + " reads \n \n")
		log.write("[" + str(time.time()) + "]" + i[:-6] + " has " + numb.strip("\n") + " reads \n \n")
	

	#####NANOFILTERING
	MakeFilteredDir = "mkdir " + pathToRun + "/working_files/demux/filtered_data"
	
	subprocess.run([MakeFilteredDir], shell=True)
	print("MADE FILTERED DIR")
	log.write("[" + str(time.time()) + "]" + "\n MADE FILTERED DIR\n ")

	filtered_barcodes = []

	#if the user specified nanofiltering parameters, slice them in and call nfilt
	if len(nFiltInfo) > 0:
		print("\n \n Counting complete. Filtering out reads <" + str(nFiltInfo[0]) + "bp and >"+ str(nFiltInfo[1]) +"bp. This may take a while. \n \n")	
		log.write("[" + str(time.time()) + "]" + "\n \n Counting complete. Filtering out reads <" + str(nFiltInfo[0]) + "bp and >"+ str(nFiltInfo[1]) +"bp. This may take a while. \n \n")	
		
		for i in tqdm(mybcodes):
			filtered_barcodes.append(work_run + "_" + i[:-6] + "_nf_qcat_" + str(nFiltInfo[0]) + "_" + str(nFiltInfo[1]) + ".fastq")
			cmd = "cat '" + pathToRun + "/working_files/demux/" + i + "' | NanoFilt --length "+ str(nFiltInfo[0]) + " -s " + seq_sum + " --maxlength " + str(nFiltInfo[1]) + " --readtype 1D > " + pathToRun + "/working_files/demux/filtered_data/" + work_run + "_" + i[:-6] + "_nf_qcat_" + str(nFiltInfo[0]) + "_" + str(nFiltInfo[1]) + ".fastq"
			subprocess.run([cmd], shell=True)
		print("Nanofiltering finished. Counting files that passed.")
		log.write("[" + str(time.time()) + "]" + "Nanofiltering finished. Counting files that passed.")

		#count the result
		counter = 0
		for i in tqdm(filtered_barcodes):	
			cmd = "awk '{s++}END{print s/4}' " + pathToRun + "/working_files/demux/filtered_data/" + i
			numb = os.popen(cmd).read()
			print(mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that passed filtering \n \n")
			log.write("[" + str(time.time()) + "]" + mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that passed filtering \n \n")
			counter += 1
		print("\n \n Counting complete.")
		log.write("[" + str(time.time()) + "]" + "\n \n Counting complete.")


	#if user didn't ask for nfilt, still move files to filt folder to 
	#simplify downstream processing
	if len(nFiltInfo) == 0:
		print("Counting complete. Not Nanofiltering. Moving your barcode files around.")
		log.write("Counting complete. Not Nanofiltering. Moving your barcode files around.")
		for i in tqdm(mybcodes):
			filtered_barcodes.append(work_run + "_" + i[:-6] + "_unFilt.fastq")
			cmd = "cp " + pathToRun + "/working_files/demux/" + i + " " + pathToRun + "/working_files/demux/filtered_data/" + work_run + "_" + i[:-6] + "_unFilt.fastq"
			subprocess.run([cmd], shell=True)
		print("Moving complete. \n \n")
		log.write("[" + str(time.time()) + "]" + "Moving complete. \n \n")



		
	print("\n \nStarting Minimap2 alignment to 16S sequence. This will take a while. \n \n")
	log.write("[" + str(time.time()) + "]" + "\n \n Starting Minimap2 alignment to 16S sequence. This will take a while. \n \n")


	time.sleep(5)
	
	#Make analysis_mapping directory
	MakeMappingDir = "mkdir " + pathToRun + "/working_files/analysis_mapping"
	subprocess.run([MakeMappingDir], shell=True)
	log.write("[" + str(time.time()) + "]" + "MADE MAPPING DIR \n \n")
#####

	#make a list to track which barcodes I'm actually continuing with,
	#so i don't needlessly process barcodes w/o data
	sxtnEssAlgn = []

	#IF there is actually data in the filtered files, then continue, otherwise DONT MAP IT
	#realized problem is that if you put an empty file into alignment, it won't tell you...
	# then you get samtools trunc. file error
	# so i get around this by filtering out files with no reads (via file size)

	#however the 573311 filter only works if the header size stays the same, which it should for
	#these 16S reads. If they grow longer, then this won't work anymore
	#I could attempt to use try/except or assert, but not sure if i can easily
	#get a samtools error to raise a real error without parsing stdout and looking for it
	#hopefully nothing will cause empty filtered fastq files to hold data.

	for i in tqdm(filtered_barcodes):
		if os.stat(pathToRun + "/working_files/demux/filtered_data/" + i).st_size > 0:
			sxtnEssAlgn.append(i[:-6] + ".sam")

			print(sxtnEssAlgn)

			cmd = "minimap2 -ax map-ont 16SMicrobial.fa " + pathToRun + "/working_files/demux/filtered_data/" + i + " > " + pathToRun + "/working_files/analysis_mapping/" + i[:-6] + ".sam | tee -a " + pathToRun + "/log.txt"
			subprocess.run([cmd], shell=True)
		else:
			continue

	log.write("[" + str(time.time()) + "]" + str(sxtnEssAlgn))

	###count all mapping reads
	
	print("\n \n Mapping complete. Counting reads per barcode that mapped to 16S. \n \n")
	log.write("[" + str(time.time()) + "]" + "\n \n Mapping complete. Counting reads per barcode that mapped to 16S. \n \n")

	counter = 0
	for i in tqdm(sxtnEssAlgn):
		cmd = "samtools view -S -c " + pathToRun + "/working_files/analysis_mapping/" + i
		numb = os.popen(cmd).read()
		print(mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that map to 16S \n")
		log.write("[" + str(time.time()) + "]" + mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that map to 16S \n")		
		counter += 1

	###alignment to honeybee
	
	print("\n \n Counting complete. Starting Minimap2 alignment to positive control sequence. \n \n")
	log.write("[" + str(time.time()) + "]" + "\n \n Counting complete. Starting Minimap2 alignment to positive control sequence. \n \n")
	

	time.sleep(5)
	
	POSalignment_files = []
	
	#map non-filtered, straight-from-demux reads to positive control, and count
	for i in tqdm(mybcodes):
		if os.stat(pathToRun + "/working_files/demux/" + i).st_size > 0:
			POSalignment_files.append(i[:-6] + "_POShoneybee.sam")
			cmd = "minimap2 -ax map-ont apismellifera.fa " + pathToRun + "/working_files/demux/" + i + " > " + pathToRun + "/working_files/analysis_mapping/" + i[:-6] + "_POShoneybee.sam | tee -a " + pathToRun + "/log.txt"
			subprocess.run([cmd], shell=True)
		else:
			continue


	print("Mapping complete. Counting reads per barcode that mapped to positive control. \n \n")
	log.write("[" + str(time.time()) + "]" + "Mapping complete. Counting reads per barcode that mapped to positive control. \n \n")

	counter = 0
	for i in tqdm(POSalignment_files):
		cmd = "samtools view -S -c " + pathToRun + "/working_files/analysis_mapping/" + i
		numb = os.popen(cmd).read()
		print(mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that map to the positive control \n")
		log.write("[" + str(time.time()) + "]" + mybcodes[counter][:-6] + " has " + numb.strip("\n") + " reads that map to the positive control \n")
		counter += 1



	print("\n \n Counting complete. Extracting primary/unique mapped reads.")
	log.write("[" + str(time.time()) + "]" + "\n \n Counting complete. Extracting primary/unique mapped reads.")


	###extract primary mapped reads, send to new sam file
	##USING 2308 PER THE OLD VERSION OF SAMTOOLS. IF SAMTOOLS IS UPGRADED TO VERSION >1,
	##THEN THE EXTRACTION MAY BE CHAINED
	
	primaryMappedReads = []
	
	for i in tqdm(sxtnEssAlgn):
		if os.stat(pathToRun + "/working_files/analysis_mapping/" + i).st_size > 537311:
			primaryMappedReads.append(i[:-4] + "_allaligned.sam")		
			cmd = "samtools view -S -h -F 2308 " + pathToRun + "/working_files/analysis_mapping/" + i + " > " + pathToRun + "/working_files/analysis_mapping/" + i[:-4] + "_allaligned.sam"
			subprocess.run([cmd], shell=True)
		else:
			continue
		
	print("\n \n Primary extracton complete. Extracting mapped reads that meet Q1, Q10, Q20, and Q30 quality.")
	log.write("[" + str(time.time()) + "]" + "\n \n Primary extracton complete. Extracting mapped reads that meet Q1, Q10, Q20, and Q30 quality.")

#########
	print("\n \n Calculating MarginStats on primary aligned reads.")
	log.write("[" + str(time.time()) + "]" + "\n \n Calculating MarginStats on primary aligned reads.")

	for i in tqdm(primaryMappedReads):
		if os.stat(pathToRun + "/working_files/analysis_mapping/" + i).st_size > 537311:
			filt_fastq = [x for x in filtered_barcodes if x[:-6] == i[:-15]]
			print(i)
			print([x[:-6] for x in filtered_barcodes])
			print([j[:-15] for j in primaryMappedReads])
			print(filt_fastq[0] + "\n \n ")
			cmd = "./marginAlign/marginStats " + pathToRun + "/working_files/analysis_mapping/" + i + " " + pathToRun + "/working_files/demux/filtered_data/" + filt_fastq[0] + " " + "16SMicrobial.fa --readIdentity --alignmentIdentity --readCoverage --mismatchesPerAlignedBase --deletionsPerReadBase --insertionsPerReadBase --localAlignment --readLength --printValuePerReadAlignment > " + pathToRun + "/working_files/analysis_mapping/" + i[:-4] + "_margStats.txt"
			print(cmd + "\n \n ")			
			subprocess.run([cmd], shell=True)
			#now plot
			cmd = "Rscript plots.R " + pathToRun + "/working_files/analysis_mapping/" + i[:-4] + "_margStats.txt " + pathToRun + "/working_files/analysis_mapping/" + i[:-4] + "_margStats_plot.pdf" 
			subprocess.run([cmd], shell=True)
		else:
			continue
#########



	#just make [1,10, 20, 30] for Q1, Q10, Q20, and Q30. And make separate subdirectories for each
	listmapQs = [1,10,20,30]
	for i in listmapQs:
		MkQDirCmd = "mkdir " + pathToRun + "/working_files/analysis_mapping/Q" + str(i)
		subprocess.run([MkQDirCmd], shell=True)


	#new list to continue to track barcodes that progress
	Q_reads = []

	#if primary mapped reads exist, extract by Q1, Q10, Q20, Q30 scores from each, and
	#put into corresponding dir
	for i in tqdm(primaryMappedReads):
		if os.stat(pathToRun + "/working_files/analysis_mapping/" + i).st_size > 537311:
			Q_reads.append(i[:-4] + "_Q_mapped.sam")
			for x in listmapQs:
				cmd = "samtools view -S -h -q " + str(x) + " " + pathToRun + "/working_files/analysis_mapping/" + i + " > " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i[:-4] + "_Q_mapped.sam"
				print(cmd)
				subprocess.run([cmd], shell=True)
		else:
			continue
		

	print("\n \n Q-score extracton complete. Now counting Q-scored reads.")
	log.write("[" + str(time.time()) + "]" + "\n \n Q-score extracton complete. Now counting Q-scored reads.")



	#Count each barcode three times for each Q-score w/ nested loop
	counter = 0
	for i in tqdm(Q_reads):
		for x in listmapQs:
			if os.stat(pathToRun + "/working_files/analysis_mapping/Q"+ str(x) + "/" + i).st_size > 537311:
				cmd = "samtools view -S -c " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i
				numb = os.popen(cmd).read()

				index = i.find("barcode")
				jump = 9
				if index == -1:
					index = i.find("none")
					jump = 4

				print(i[index:index+jump] + " has " + numb.strip("\n") + " reads that pass Q" + str(x) + " quality \n \n")
				log.write("[" + str(time.time()) + "]" + i[index:index+jump] + " has " + numb.strip("\n") + " reads that pass Q" + str(x) + " quality \n \n")
				counter += 1

			else:
				index = i.find("barcode")
				jump = 9
				if index == -1:
					index = i.find("none")
					jump = 4

				print(i[index:index+jump] + " has no reads that pass Q" + str(x) + " \n \n ")
				log.write(i[index:index+jump] + " has no reads that pass Q" + str(x) + " \n \n ")
				continue


	print("\n \n Counting complete. Now extracting top 50 aligned reads per barcode.")
	log.write("[" + str(time.time()) + "]" + "\n \n Counting complete. Now extracting top 50 aligned reads per barcode.")


	all_top50 = []

	#extract top 50 refseq IDs by count, move to txt file
	for i in tqdm(Q_reads):
		for x in listmapQs:
			if os.stat(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i).st_size > 537311:
				all_top50.append(i[:-4] + "_top50_mapped.txt")			
				cmd = "samtools view -S " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i + " | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 50 > " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i[:-4] + "_top50_mapped.txt"
				subprocess.run([cmd], shell=True)
			else:
				continue
		
	print("\n \n Top 50 extracton complete. Merging with organism reference names.")
	log.write("[" + str(time.time()) + "]" + "\n \n Top 50 extracton complete. Merging with organism reference names.")


	merged_with_refs = []

	#merge top 50 refseq IDs w/ latin names
	#should just import the actual python script -- no need for subprocess
	for i in tqdm(all_top50):
		for x in listmapQs:
			if os.path.exists(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i):
				merged_with_refs.append(i[:-4] + "_refNames.txt")			
				cmd = "python intersect_sam_id_ref_id.py --in " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i + " --ref 16SMicrobial.fa.headers.txt --out " + pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i[:-4] + "_refNames.txt"
				subprocess.run([cmd], shell=True)

	print("\n \n Top 50 reads per barcode merged with reference names.")
	log.write("[" + str(time.time()) + "]" + "\n \n Top 50 reads per barcode merged with reference names.")

	#plot stacked bar plot of IDs in each output ref ID file


	#remove duplicates from merged_with_refs so I don't calculate redundantly
	merged_with_refs = list(dict.fromkeys(merged_with_refs))



	#import plotter, call directly
	print("\n \n Making stacked barplots \n \n ")
	for i in tqdm(merged_with_refs):
		for x in listmapQs:
			if os.path.exists(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i):
				testPlots.plotter(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i)

	#import kronaplotter, though it calls with subprocess
	print("Making Krona plots \n \n ")
	for i in tqdm(merged_with_refs):
		for x in listmapQs:
			if os.path.exists(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i):
				testPlots.kronaplotter(pathToRun + "/working_files/analysis_mapping/Q" + str(x) + "/" + i)



start()

