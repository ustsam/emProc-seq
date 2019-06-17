import sys,gzip

##
#Known bug, it wont return the alignment for last group of alignment in the file.
##
workdir = sys.argv[1]

infiles = ["/9_alignment.sam.gz","/10_alignment.sam.gz"]
outfile = gzip.open(workdir + "/11_alignment.sam.gz","wb")
outfile1 = gzip.open(workdir + "/13_Enderror.sam.gz","wb")
outfile2 = gzip.open(workdir + "/14_hardclipped.sam.gz","wb")

for file in infiles:
	infile = gzip.open(workdir + file,"rb")
	
##running lists
	Alignment = []
	AlignmentScore = []
	nM = []	#For if condition to screen Spliced read with difference AS
##running lists	
#current read
	SequenceID = ""
	
	for line in infile:
		RawLine = line
		line = line.split()
#		print RawLine
		if line[0] == SequenceID:
			k = 11
			AS = []
			MuTag = []
			while k < len(line):
				if line[k].split(":")[0] == "AS":
					AS = line[k].split(":")
					k+=1
				elif line[k].split(":")[0] == "NM":
					MuTag = line[k].split(":")
					k+=1
				elif len(AS) != 0 and len(MuTag) != 0:
					break
				else:
					k+=1
			AlignmentScore.append(int(AS[2]))
#			print "ASscore"
#			print AlignmentScore
			numbers = ""
			letters = []
			for character in line[5]:       #Chop according to the number field 6
				if character == "M" or character == "S" or character == "N" or character == "I" or character == "D" or character == "H" or character == "P":
					letters.append(character)
					numbers += "X"
				else:
					numbers += character
			numbers = numbers.split("X")
#			print "MuTag[2]"
#			print MuTag[2]
#			print "numbers[]"
#			print numbers
			i = 0
			while i < len(letters):
				if letters[i] == "S" or letters[i] == "H":
					MuTag[2] = int(MuTag[2]) + int(numbers[i])    #!int! is needed since what append in above is string type.
				i+=1
			nM.append(int(MuTag[2]))
			Alignment.append(RawLine)
			
		else:

			if len(Alignment) > 0 :
				OriginalIndex = []
				i = 0
				while i < len(AlignmentScore):
				        OriginalIndex.append(i)
				        i+=1
				ASsorted = sorted(AlignmentScore,reverse=True)
##This chunk is used to sorted list before, it would have bug if The list contain the same value, it would sort the list by the value in the list. So i change the sorting method as below, sort Index first then nM according to index sorted in IndexSorted
#				nMsortedbyAS = sorted(nM, key=lambda x: AlignmentScore.index(x))
#				IndexSorted = sorted(OrginalIndex, key=lambda x: AlignmentScore.index(x))
#				nMsortedbyAS = [x for _,x in sorted(zip(nM, AlignmentScore),reverse=True)]
#				IndexSorted = [x for _,x in sorted(zip(OriginalIndex, AlignmentScore),reverse=True)]
#				nMsortedbyAS = [a for b,a in sorted(zip(AlignmentScore,nM),reverse=True)]
#				IndexSorted = [a for b,a in sorted(zip(AlignmentScore,OriginalIndex),reverse=True)]
##
				IndexSorted = sorted(range(len(AlignmentScore)), key=lambda k: AlignmentScore[k], reverse=True)
				pix = 0
				nMsortedbyAS = []
				while pix < len(IndexSorted):
					nMsortedbyAS.append(nM[IndexSorted[pix]])
					pix+=1
##Debug field
#				print AlignmentScore
#				print nM
#				print "---------------------------------------------------"
#				print ASsorted
#				print nMsortedbyAS
#				print IndexSorted
#				print "+++++++++++++++++++++++++++++++++++++++++++++++++++"
##
				BestIndex = []
				g = 0
#				print "Min(nM)"+str(min(nM))
				while g < len(ASsorted):
#					print "nMsortedbyAS[g]:"+str(nMsortedbyAS[g])
#					print "min(nM):"+str(min(nM))
					if str(nMsortedbyAS[g]) == str(min(nM)):
						BestIndex.append(IndexSorted[g])
						break
#					print g
					g+=1
#				print len(BestIndex)
#				print BestIndex

				if len(BestIndex) > 0 :
					Splitchecker = Alignment[BestIndex[0]]
					Splitchecker = Splitchecker.split()
					if Splitchecker[5].count("S") == 0 and Splitchecker[5].count("H") == 0:
						outfile.write(Alignment[BestIndex[0]])
					elif Splitchecker[5].count("S") > 0 and Splitchecker[5].count("H") == 0:
						outfile1.write(Alignment[BestIndex[0]])
					elif Splitchecker[5].count("H") > 0:
						outfile2.write(Alignment[BestIndex[0]])

			SequenceID = line[0]
#			print "SequenceID"
#			print SequenceID

			Alignment = []
			AlignmentScore = []
			nM = []

			k = 11
			AS = []
			MuTag = []
			while k < len(line):
				if line[k].split(":")[0] == "AS":
					AS = line[k].split(":")
					k+=1
				elif line[k].split(":")[0] == "NM":
					MuTag = line[k].split(":")
					k+=1
				elif len(AS) != 0 and len(MuTag) != 0:
					break
				else:
					k+=1

			AlignmentScore.append(int(AS[2]))
			numbers = ""
			letters = []
			for character in line[5]:       #Chop according to the number field 6
				if character == "M" or character == "S" or character == "N" or character == "I" or character == "D" or character == "H" or character == "P":
					letters.append(character)
					numbers += "X"
				else:
					numbers += character
			numbers = numbers.split("X")

			i = 0
			while i < len(letters):
				if letters[i] == "S" or letters[i] == "H":
					MuTag[2] = int(MuTag[2]) + int(numbers[i])    #!int! is needed since what append in above is string type.
				i+=1
			nM.append(int(MuTag[2]))
			Alignment.append(RawLine)


