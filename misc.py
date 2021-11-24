#######################################################################################################
'Sequence Mutation - GAP - others'
'Single mutation of sequence '
#Functions: mutate, tangorun  ..  creates single fasta file with all mutants
def mutate(infile, start = 1, end = 0, mutAA = 'A'):
	#Input Parameters
	#start = 1 # Default 1
	#end = 0   # 0 not defined: auto set to max (total length)
	#mutAA = 'A'

	#Verifying the argument
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	pro = temp[0] + "_"

	#Reading Input file 
	F = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = F.read().splitlines()

	i = 1
	seq = ""
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:			
		if i == 1 and line[0] == '>':
			title = line
		elif i != 1 and line[0] == '>':
			print("Only one fasta sequence is expected. Plz check the file content for multiple sequence or blank lines !!")
			break
		else:
			seq = seq + line
			print(i-1, "  ", line)
		i = i + 1
		F.close()
	print("Input sequence file: ", IFILE, " read successfully.\nThe sequence titled as ", title)
	print("The sequence (", len(seq), "AA ) is ", seq)

	i = 1	
	#Creating folder and Defining Prefix for output file names
	temp = pro + 'mutants_' + mutAA + ".fasta" 
	while os.path.exists(temp):
		temp = pro + 'mutants_'+  mutAA + str(i) + ".fasta" 
		i = i + 1
	OFILE = temp

	mutno = 0
	#Creating Mutant Sequences
	i = start-1;
	if end == 0:
		end = len(seq)


	tmp = "X"
	F = open(OFILE, 'w')
	F.write(">mutant_" + str(mutno) + "_res" + str(0) + tmp + "\n")
	F.write(seq+ "\n")	
	
	while i < end:
#		print i;	
		#print "started"
		temp = list(seq)
		tmp = 'X'
		if (temp[i] != mutAA) & (i > -1):
		#Mutating the original sequence	using Alanine at position i			
			tmp = temp[i]
			temp[i] = mutAA
			mutseq = "".join(temp)
			mutno = mutno + 1
		#Writing mutated sequence to file
			#temp = title + " Mutant_" + str(i+1) + "\n"
			F.write(">mutant_" + str(mutno) + "_res" + str(i+1) + tmp + "\n")
			F.write(mutseq+ "\n")
		i = i + 1 
	F.close()
	print("job completed successfully !!", mutno,"Mutant sequences have been created in file: ", OFILE, 'by replacing residues with amino acid "'+mutAA+'"'			)
	return OFILE	

def tangorun(infile, ct="N",nt="N",ph=7.4,te=303,io=0.1):
	#Input Parameters
	curdir = os.getcwd();
	epath = 'C:/Users/rpkaran/Documents/GoogleDrive/Project/mylib/mutation_executables/Tango.exe'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
#	param = "N N 7 298 0.1"
	app = "tango"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			OUT.write(title + " " + param + " ")
			OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	OUT.write(title + " " + param + " ")
	OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	print("Input file for",app,"created successfully: ", OFILE, "In folder ", temp)

	cmd = epath + " -inputfile=" + OFILE + " > tango.log"
	print("executing ",cmd)
	status = subprocess.call(cmd, shell=True)
	os.chdir(curdir)
	if status == 0: 
		print(app,"successfully executed !!")
		#os.chdir(curdir)  #os.chdir('../')
		#print curdir 	print os.getcwd()
	else:
		print("Error in execution of", app)
	return 0	

def GAPrun(infile):
	#Input Parameters
	#epath = '/home/user1/Documents/Project/mylib/mutation_executables/tango_x86_64_release'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
#	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
	param = ""
#	param = "N N 7 298 0.1"
	app = "GAP"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File ##and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	#OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	#OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			#OUT.write(title + " " + param + " ")
			#OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	#OUT.write(title + " " + param + " ")
	#OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	#OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	#print "Input file for",app,"created successfully: ", OFILE, "In folder ", temp

	#cmd = epath + " -inputfile=" + OFILE + " > GAP.log"
#	print "executing ",cmd
	#status = subprocess.call(cmd, shell=True)
	status = accessGAP(infile, titlelist)
	if status == 0:
		print(app,"successfully executed !!")
		os.chdir('../')
		return titlelist
	else:
		print("Error in execution of", app)
		return -1
		
		
def accessGAP(infile, titlelist):
	outfile1 = "GAPout.html"
	outfile2 = "GAPout.txt"
	#os.mkdir("test")
	F = open("../"+infile, 'r')
	Seq = F.read()
	F.close()

	#Seq = '>sp|P37840|SYUA_HUMAN Alpha-synuclein OS=Homo sapiens GN=SNCA PE=1 SV=1' + '\n' + 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
	params = urllib.urlencode({'sequence': Seq, 'predict': 'Predict'})
	#, '@type': 'submit'
	headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/html","Accept-Encoding":"gzip, deflate", "Connection":"keep-alive","Host":"www.iitm.ac.in","Origin":"http://www.iitm.ac.in","Referer":"http://www.iitm.ac.in/bioinfo/GAP/"}
	conn = httplib.HTTPConnection("10.93.249.92:80")
	#"Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"

	#print(Seq)
	#print(params)
	#print(headers)

	conn.request("POST", 'http://www.iitm.ac.in/bioinfo/cgi-bin/GAP/cal_fasta.pl', params, headers)
	response = conn.getresponse()
	print(response.status, response.reason)
	#if response.headers().get('Content-Encoding') == 'gzip':
	buf = StringIO( response.read())
	f = gzip.GzipFile(fileobj=buf)
	conn.close()

	html = f.read()
	#html = html.decode('utf8')

	#html = open("foobar.html").read()
	#data = html2text.html2text(html)
	#print(data)
	#soup = BeautifulSoup(html)

	F = open(outfile1, 'w')
	F.write(html)
	F.close()

#	F = open(outfile2, 'w')
#	F.write(data)
#	F.close()
	
	html2dat(outfile1, titlelist)
	return 0
	
def html2dat(infile, titlelist):
	
	#infile = "out.html"
	outfile1 = "outtable1.txt"
	outfile2 = "outtable2.txt"
	outfile3 = "outtable3.txt"

	PATnumSeq = re.compile('(?<=Number of sequences: )\d*',re.IGNORECASE)
	PATSeq = re.compile('(?<=Input sequence</u>: )\w*',re.IGNORECASE)
	PATout = re.compile('(?<=<br>The given peptide is Amyloid \()\d*.\d*',re.IGNORECASE)

	PATtable = re.compile('<table.*?</table>',re.IGNORECASE)
	PATtr = re.compile('<tr>.*?</tr>',re.IGNORECASE)
	PATtd = re.compile('<td.*?</td>',re.IGNORECASE)
	PATtag = re.compile('<.*?>',re.IGNORECASE)

	theader = ("No","Peptide","Amyloid probability","Non-amyloid probability")
	gap = 6 - 1
	txt = ""
	ttxt = "\t" + " No " + "\t" + "Peptide" + "\t" + "Amyloid probability" + "\t" + "Non-amyloid probability" + '\t'
	t = list()
	data = list()

	F = open(infile, 'r')
	Seq = F.read()
	F.close()

	NumSeqs = int(PATnumSeq.findall(Seq)[0])
	Seqs = PATSeq.findall(Seq)
	outs = PATout.findall(Seq)

	print(NumSeqs)
	#print(Seqs)
	#print(outs)

	for i in range(0,NumSeqs):
		#print i
		txt = txt + titlelist[i] + '\t' + str(i) + '\t' + outs[i] + '\t' + Seqs[i] + '\n'
		#print txt
	
	with open(outfile1, "w") as output:
		output.write(txt)

	txt = ""
	seltxt = ""
	
	#tables = pat1.match(Seq)
	tables = PATtable.findall(Seq)
	if tables:
	#	print tables
	#	print tables.groups()
		i = -1
		for table in tables:		
			i = i + 1
			mutation = titlelist[i]
			#mutant_0_res0X
			sel = int(re.findall('(?<=res)\d*',mutation)[0])
	#		print i,table,'\n'
			trows = PATtr.findall(table)
			if trows:
				j = -1
				for trow in trows:
					j = j + 1
	#				print i,j,PATtag.sub("\t",trow) ,i
					temp = PATtag.sub("\t",trow) 
					temp = re.sub("\t\t","\t",temp)
					if (temp != ttxt):
						txt = txt + mutation + temp + '\n'
						if (j in range(sel - gap,sel+1)) or (sel == 0):
							seltxt = seltxt + mutation + temp + '\n'				
			
	if (NumSeqs != i+1):
		print("Warn: Out of ",NumSeq,", ",i+1," outputs were retrieved !!!. #########")

	#print(data)
	#print(type(txt))		
	#F.write(str(data))	
	#with open(outfile2, "w") as output:
	#	writer = csv.writer(output, lineterminator='\n')
	#	writer.writerows(data)			
	with open(outfile2, "w") as output:
		output.write(txt)
	with open(outfile3, "w") as output:
		output.write(seltxt)
	return(0)	
	
def PASTA2run(infile):
	perlf="PastaPairs.pl"
	pot="pot_pasta.dat"
	abspath = os.path.dirname(os.path.abspath(__file__)) + "/PASTA2II/"
	path= abspath 
	top=10
	flag= 2
	typ = "self"
	cutoff = -5.5
	base = ""
	#Input Parameters
	epath = '/home/user1/Documents/Project/mylib/pasta_exe/'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	app = "PASTA2"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	path = path + temp + '/'
	print("Created output folder: ", temp)
	

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	status = 0
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':			
			base = title
			os.mkdir(base)
			os.chdir(base)
			OFILE = "mutants_inputs" + str(seqno) + ".fasta"
			OUT = open(OFILE, 'w')			
			OUT.write(">" + title + "\n" + seq )
			OUT.close()			
			cmd = "(cd " + epath + " && perl " + perlf + " " + pot + " " + path + base + '/ ' + str(top) + " " + str(flag) + " " + str(typ) + " " + str(cutoff) + " " + base + ")" + ' > PASTA2.log'
			#print cmd
			print("Executing ", app,": mutant seqno -> ", seqno)
			status = subprocess.call(cmd, shell=True)
			if status != 0:				
				print("Error in execution of", app, ": sequenceno",seqno)
				break
			os.chdir(path)
			seq = ""
			seqno = seqno + 1
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
			i = i + 1
		IN.close()
	if status == 0:		
		base = title
		os.mkdir(base)
		os.chdir(base)
		OFILE = "mutants_inputs" + str(seqno) + ".fasta"
		OUT = open(OFILE, 'w')
		OUT.write(">" + title + "\n" + seq )
		OUT.close()
		cmd = "(cd " + epath + " && perl " + perlf + " " + pot + " " + path + base + '/ ' + str(top) + " " + str(flag) + " " + str(typ) + " " + str(cutoff) + " " + base + ")" + ' > PASTA2.log'
		status = subprocess.call(cmd, shell=True)
		if status != 0:				
			print("Error in execution of", app, ": sequenceno",seqno)
		seqno = seqno + 1
		titlelist.append(title)	
	
	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	print("Input file for",app,"created successfully in folder ", temp)
	os.chdir(abspath)
	return status	


def waltzrun(infile):
	#Input Parameters
	#epath = '/home/user1/Documents/Project/mylib/mutation_executables/tango_x86_64_release'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
#	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
	param = ""
#	param = "N N 7 298 0.1"
	app = "waltz"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File ##and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	#OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	#OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			#OUT.write(title + " " + param + " ")
			#OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	#OUT.write(title + " " + param + " ")
	#OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	#OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	#print "Input file for",app,"created successfully: ", OFILE, "In folder ", temp

	#cmd = epath + " -inputfile=" + OFILE + " > GAP.log"
#	print "executing ",cmd
	#status = subprocess.call(cmd, shell=True)
	status = accesswaltz(infile, titlelist)
	if status == 0:
		print(app,"successfully executed !!")
		os.chdir('../')
		return titlelist
	else:
		print("Error in execution of", app)
		return -1
		
		
def accesswaltz(infile, titlelist):
	outfile1 = "waltzout.html"
	outfile2 = "waltzout.txt"
	#os.mkdir("test")
	F = open("../"+infile, 'r')
	Seq = F.read()
	F.close()

	#Seq = '>sp|P37840|SYUA_HUMAN Alpha-synuclein OS=Homo sapiens GN=SNCA PE=1 SV=1' + '\n' + 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
	params = urllib.urlencode({'sequence': Seq, 'threshold': '92', 'custom_value' : 0, 'ph': 7.0, 'output': 'text_short', 'Submit': 'Submit sequences'})
	#, '@type': 'submit'
	headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/html","Accept-Encoding":"gzip, deflate", "Connection":"keep-alive","Host":"waltz.switchlab.org","Origin":"http://waltz.switchlab.org","Referer":"http://waltz.switchlab.org/results.cgi"}
	conn = httplib.HTTPConnection("134.58.50.6:80")
	#"Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"

	#print(Seq)
	#print(params)
	#print(headers)

	conn.request("POST", 'http://waltz.switchlab.org/results.cgi', params, headers)
	response = conn.getresponse()
	print(response.status, response.reason)
	#if response.headers().get('Content-Encoding') == 'gzip':
	buf = StringIO( response.read())
	f = gzip.GzipFile(fileobj=buf)
	conn.close()

	html = f.read()
	#html = html.decode('utf8')
	PATURL = re.compile('<a href="\.\.(.*?)">here</a>',re.IGNORECASE)
	URL = PATURL.findall(html)[0]		

	URL = 'http://waltz.switchlab.org/output' + URL
	print(URL)
	conn = httplib.HTTPConnection("134.58.50.6:80")
	conn.request("GET", URL)
	response = conn.getresponse()	
	print(response.status, response.reason)
#	if response.headers().get('Content-Encoding') == 'gzip':
#	buf = StringIO( response.read())
#	f = gzip.GzipFile(fileobj=buf)
#	html = f.read()
#	else:
	html = response.read()	
	conn.close()	

	
	F = open(outfile1, 'w')
	F.write(html)
	F.close()
	
	
	html2datwaltz(outfile1, titlelist)
	return 0
	
def html2datwaltz(infile, titlelist):	
	#infile = "out.html"
	outfile1 = "outtable1.txt"
	outfile2 = "outtable2.txt"
	outfile3 = "outtable3.txt"
	
	PATdiv = re.compile('<div id="bodyText">.*?</div>', re.DOTALL)
	PATtable = re.compile('<table.*?</table>',re.DOTALL)
	PATtr = re.compile('<tr>.*?</tr>', re.DOTALL)
	PATtd = re.compile('<td.*?</td>',re.IGNORECASE)
	PATtag = re.compile('<.*?>',re.IGNORECASE)
		
	txt = ""
	F = open(infile, 'r')
	Seq = F.read()
	F.close()
	theader = ("Job","Regions","Description")
	#print(Seq) 
	#print("######## seq ")
	div = PATdiv.findall(Seq)[0]
	#print div
	if div:
		table = PATtable.findall(div)[0]
	#print table				
	if table:
		trows = PATtr.findall(table)
	if trows:
		j = -1
		for trow in trows:
			j = j + 1			
			tcells = PATtd.findall(trow)
			if tcells:
				i = 0
				for tcell in tcells:				
					tcell = PATtag.sub("",tcell)					
					txt = txt + tcell + "\t"
					
				txt = txt + "\n"
							
	print("######", j ," outputs were retrieved !!!. #########")

	#print(data)
	#print(type(txt))		
	#F.write(str(data))	
	#with open(outfile2, "w") as output:
	#	writer = csv.writer(output, lineterminator='\n')
	#	writer.writerows(data)			
	with open(outfile1, "w") as output:
		output.write(txt)
	return(0)

####################################################################
def read_amylpred(filename):
    with open(filename, 'r') as ifile:
        tdict = {}
        line = ifile.readline()
        while line:
            if line[:16] == 'Sequence length:' :
                l = int(line.replace(' ', '').split(':')[1])
            if line[0:4] == 'HITS':
                line = ifile.readline()
                line = ifile.readline()
                temp = line.strip().split(':')
                while len(temp) == 2:                
                    #print(temp)
                    tdict[temp[0]] = temp[1]
                    line = ifile.readline()   
                    temp = line.strip().strip('>--> ').replace(' ', '').split(':')
            line = ifile.readline()
    tdict['seqlen'] = l
    tdict['filename'] = filename
    return tdict

def read_amylpred2(ifilename):
    pattern = re.compile("(\w+\s*\w*): ([0-9-,:# ]+)")
    flag = False
    dictlist = []
    dict1 = {}
    with open(ifilename, 'r') as ifile:
        line = ifile.readline()
        while line :
            match = pattern.findall(line)
            if len(match) > 0:
                #print( match )
                match = match[0]
                if match[0] == "Job ID":
                    if flag:
                        dictlist.append(dict1)
                        dict1 = {}
                    else:
                        flag = True
               # elif match[0] == "Sequence length":
                #    print(match)
                dict1[match[0]] = match[1]                
            line = ifile.readline()
        dictlist.append(dict1)
    return dictlist