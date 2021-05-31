#Path_mips

import os
import argparse 
import re



# give indexes of regions which lack gaps in multiple sequence alignment
def remove_gaps(file_path):

    from Bio import AlignIO
    alignment = AlignIO.read(file_path, "fasta")

    nogap_alignment = []
    for i in range(len(alignment[0])):
        if '-' not in alignment[:,i]:
            nogap_alignment.append(i)

    return nogap_alignment



def select_index(nogap_alignment):
    from itertools import groupby
    from operator import itemgetter

    indexes =[]
    for k,g in groupby(enumerate(nogap_alignment),lambda x:x[0]-x[1]): 
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        indexes.append((group[0],group[-1]))

    final_indexes=[]
    for (i,j) in indexes:
        if j-i >= 100:
            final_indexes.append((i,j)) 

    return final_indexes

# Select regions lacking gaps from the sequence 
def load_seq(final_indexes,seq):
    seque = list(seq[i:j] for i,j in final_indexes)
    return seque

# open fasta file containing sequence
def open_seq(seq):
    from Bio import SeqIO
    import os
   
    seq = SeqIO.read(open(seq),'fasta')
    return str(seq.seq)


def count_repeat(seq, num_repeats=6):
    count_repeat_A = "A" * num_repeats
    count_repeat_T = "T" * num_repeats
    count_repeat_G = "G" * num_repeats
    count_repeat_C = "C" * num_repeats
    if count_repeat_A in seq or count_repeat_T in seq or count_repeat_G in seq or count_repeat_C in seq:
        return True
	#if re.search(r"A{6,}|T{6,}|G{6,}|C{6,}", seq):
 #   if re.search(r"G{num_repeats,}", seq) or re.search(r"A{num_repeats,}", seq) or re.search(r"C{num_repeats,}", seq) or re.search(r"T{num_repeats,}", seq):
#	return True
def count_repeats(seq, num_repeats):
    if re.search(r"A{%s,}|T{%s,}|G{%s,}|C{%s,}"%(num_repeats,num_repeats,num_repeats,num_repeats), seq):
        return True

def reverseCompl(s):
    complement = {"A":"T", "G":"C", "T":"A", "C":"G"}
    st = ""
    for i in s.upper():
        st = complement[i] + st # here complement bases are added in reverse order
    
    return st


def check_probes(data):
    import Bio.SeqUtils
    from Bio.SeqUtils import MeltingTemp as tm
    import pandas as pd

    tmext = []
    tmlig = []
    gcext = []
    gclig = []
    
    ind1 = []
    ind2 = []
    ind3 = []
    ind4 = []
    
    # check for GC content of all_probes (expected range: 40-70%)
    for ext, lig in zip(data['Extension probe'], data['Ligation probe']):
        if Bio.SeqUtils.GC(ext) < 40 or Bio.SeqUtils.GC(ext) > 70 or Bio.SeqUtils.GC(lig) < 40 or Bio.SeqUtils.GC(lig) > 70:
            # remove those which follows the above criteria
            ind1.append(int(data[data['Extension probe'] == ext].index.values))
            ind2.append(int(data[data['Ligation probe'] == lig].index.values))
            

    for ext, lig in zip(data['Extension probe'], data['Ligation probe']):
        if tm.Tm_NN(ext, nn_table=tm.DNA_NN2) < 60 or tm.Tm_NN(lig, nn_table=tm.DNA_NN2) < 60:
            # remove those which follows the above criteria
            ind3.append(int(data[data['Extension probe']== ext].index.values))
            ind4.append(int(data[data['Ligation probe']== lig].index.values))

    inde = set(ind1+ind2+ind3+ind4)       
    for i in inde:
        data.drop([i], inplace=True)
    
    for ext, lig in zip(data['Extension probe'], data['Ligation probe']):
        tmext.append(tm.Tm_NN(ext, nn_table=tm.DNA_NN2))
        tmlig.append(tm.Tm_NN(lig, nn_table=tm.DNA_NN2))
        gcext.append(Bio.SeqUtils.GC(ext)) 
        gclig.append(Bio.SeqUtils.GC(lig))
    
    data['tm_ex'] = tmext
    data['tm_li'] = tmlig
    data['GC_ex'] = gcext
    data['GC_li'] = gclig

    return data


def probes(seq, extension_probe_length, ligation_probe_length, insert_size):
    
    import pandas as pd
    # plus is 5'-3' and minus is 3'-5' strand
    # extension_probe_length = 18nts (ideally: 16-20nts)
    # ligation_probe_length = 20nts (ideally: 20-24nts)
    # insert_size = 112 nts (max size: 112nts)
    extension_probe_plus = []
    ligation_probe_plus = []
    extension_probe_minus = []
    ligation_probe_minus = []
    target_region_plus = []
    target_region_minus = []
    
    #all positions are zero-based
    extension_probe_pos_plus =[]
    ligation_probe_pos_plus = []
    
    extension_probe_pos_minus =[]
    ligation_probe_pos_minus = []
    
    
  #  steps = int(extension_probe_length/2) what to keep?
    steps = 2
    plus = 0
    minus = 0
    
    # selecting probes of defined length iteratively from 5'-3' of the input sequence 
    for i in range(0, len(seq), steps):
        for j in range((i + 1), (len(seq) + 1)):
            if len(seq[i:j]) == extension_probe_length and len(seq[j+insert_size:j+insert_size+ligation_probe_length]) == ligation_probe_length:
                # discard sequences with consecutive polyNs (=repeats)
                if count_repeat(seq[i:j]) == True or count_repeat(seq[j+insert_size:j+insert_size+ligation_probe_length]) == True:
                    break

                # if sequences do not contain consecutive polyNs, proceed
                else:
                    if minus <= plus:
                        # extract extension_probe of required length and store in a list
                        # note start and end indexes of extension_probe 
                        extension_probe_minus.append(seq[i:j])
                        extension_probe_pos_minus.append(f"start={i} end={j}")
                        
                        # skip 112 nts (=insert size) and extract following ligation_probe of required length

                        ligation_probe_minus.append(seq[j+insert_size:j+insert_size+ligation_probe_length])

                        # note start and end indexes of extension_probe and store in a list
                        ligation_probe_pos_minus.append(f"start={j+insert_size} end={j+insert_size+ligation_probe_length}")
                        target_region_minus.append(reverseCompl(seq[j:j+insert_size])) # direction of minus strand 5-3; written as reverse complement of plus strand
                        minus += 1
                        break

                    if plus < minus:
                        ligation_probe_plus.append(reverseCompl(seq[i:j]))
                        ligation_probe_pos_plus.append(f"start={i} end={j}")
                        extension_probe_plus.append(reverseCompl(seq[j+insert_size:j+insert_size+ligation_probe_length]))
                        extension_probe_pos_plus.append(f"start={j+insert_size} end={j+insert_size+ligation_probe_length}")
                        target_region_plus.append(seq[j:j+insert_size])
                        plus += 1
                        
                        break

                
                  
    plus_probes=pd.DataFrame({"Extension probe":pd.Series(extension_probe_plus), "extension_probe_pos_plus":pd.Series(extension_probe_pos_plus),"Target region Plus_strand":pd.Series(target_region_plus),"Ligation probe":pd.Series(ligation_probe_plus),"ligation_probe_pos_plus":pd.Series(ligation_probe_pos_plus)})
    minus_probes=pd.DataFrame({"Extension probe":pd.Series(extension_probe_minus), "extension_probe_pos_minus":pd.Series(extension_probe_pos_minus),"Target region Minus_strand":pd.Series(target_region_minus),"Ligation probe":pd.Series(ligation_probe_minus),"ligation_probe_pos_minus":pd.Series(ligation_probe_pos_minus)})

    return plus_probes, minus_probes


## index of probes needs to be changed!

def find_mips(seque, extension_probe_length, ligation_probe_length, insert_size):
    
    for seq in seque:
        name=seque.index(seq)
        plus_probes, minus_probes = probes(seq, extension_probe_length, ligation_probe_length, insert_size)
        final_plus_probes=check_probes(plus_probes)
        final_minus_probes=check_probes(minus_probes)
        
        final_plus_probes["MIPs plus strand"] = final_plus_probes["Ligation probe"]+"--linker--"+final_plus_probes["Extension probe"]
        final_minus_probes["MIPs minus strand"] = final_minus_probes["Ligation probe"]+"--linker--"+final_minus_probes["Extension probe"]
        final_plus_probes.to_csv(f"final_plus_probes_{name}.csv")
        final_minus_probes.to_csv(f"final_minus_probes_{name}.csv")

    

if __name__ == "__main__":
	# initialize your parser
    parser = argparse.ArgumentParser(description = "This script design molecular inversion probes")   

	# parse the arguments
    parser.add_argument("--path", type=str, help="Path of the directory where alignment files are stored")
 #   parser.add_argument("--align_file", type=str, help="Name of the alignment fasta file")
    parser.add_argument("--infolder", type=str, help="Path of the directory where fasta files containing gene sequence are stored")
    parser.add_argument("--seq",type=str, help="Name of the fasta file containing gene sequence to be analysed")
    parser.add_argument("--extension_probe_length",type=int, help="Size of extension_probe_length: 16-20nts")
    parser.add_argument("--ligation_probe_length",type=int, help="Size of ligation_probe_length: 20-24nts")
    parser.add_argument("--insert_size",type=int, help="Size of the target sequence: maximum = 112 nts")

    args = parser.parse_args()
    path = args.path
   # align_file = path.align_file
    infile = args.infolder
    seq = args.seq
    extension_probe_length = args.extension_probe_length
    ligation_probe_length = args.ligation_probe_length
    insert_size = args.insert_size

    os.chdir(infile)
    seq=open_seq(seq)
    
    os.chdir(path) 

	# iterate through all files in the directory
    for file in os.listdir(): 
	    # Check whether file is in correct format or not 
        if file.endswith(".fna") or file.endswith(".fasta"):
            file_path = f"{path}/{file}"
            new_path = f'{file_path}.MIPsout'
            os.makedirs(new_path)
            os.chdir(new_path)
  			
	    	# call functions
            nogap_alignment=remove_gaps(file_path) 
            final_indexes=select_index(nogap_alignment)
            seque=load_seq(final_indexes,seq)
            find_mips(seque,extension_probe_length, ligation_probe_length, insert_size)


            
