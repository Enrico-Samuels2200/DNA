#   This program allows the user to find and view abnormalities within a DNA sequence.

DNA = "DNAFile.txt"
Normal = "normalDNA.txt"
Mutate = "mutateDNA.txt"

#   The keys in this dictionary are protiens and their values are the protiens that the protiens produce.
amino_acid = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	}

#   This function allows us to pass a file name into it and utilizes it.
def readSeq(file):
  with open(file, 'r') as line:
    codon = line.read()
    sequence = codon.replace("\n" , "").strip() #   Once we read the file we replace all line breaks with "". 
    return sequence

#   This function is responsible for writing two files. The one file contain the normal DNA sequence and
#   the other contain the mutated DNA sequence. This function is also responsible for pointing out where
#   the abnormality occured.
def mutate(data):
    find_char = data.find("a") #    This allows us to identify where the abnormality was found.
    
    with open("normalDNA.txt", "w") as normal:
        normal.write(data.replace("a", "A"))
    normal.close()
    
    with open("mutateDNA.txt", "w") as mutated:
        mutated.write(data.replace("a", "T"))
    mutated.close
    
    return find_char + 1

#   This function is responsible for reading the data passed into it and returning the protiens based off the
#   DNA sequence passed through it.
def txtTranslate(data):
    length = len(data)
    protiens = ''
    base_index = 0
    index = 3
    
    while len(protiens) < length: #   This set up loop and a factor of control
        #   I've utilized a try and catch method because if the length is equal to an odd number it will
        #   result in an error. Therefor if the error arises it'll fall to the except statement and by
        #   default the data that will be passed will be 'X'.
        
        try: 
            read_data = data[base_index : index]
            
            if len(read_data) % 3 == 0: #   This makes sure the length of data being passed is divisble by 3
                protiens += amino_acid[read_data] # amino_acid[read_data] is calling the key value
                
            #   This iterations allows us to move to the next three values of data passesd into the function
            base_index += 3
            index += 3
        except:
            protiens += 'X'
    return protiens

mutation_occurence = mutate(readSeq(DNA))

normal_DNA = txtTranslate(readSeq(Normal))
mutated_DNA = txtTranslate(readSeq(Mutate))

#   This prints the data in a specific order for the user to view.
print(f'''
      DNA Scan Results:
      -----------------
      
      [!] Abnormality occur at nucleotide {mutation_occurence}
      
      [*] Normal DNA Sequence:
          ====================
          {normal_DNA}
          
          
      [!] Mutated DNA Sequence:
          =====================
          {mutated_DNA}
      
     ''')