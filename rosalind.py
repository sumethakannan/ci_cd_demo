import end as end

s="CCAGAAAGGAAGTAACGAGCGAACTTGAGACGCATGAAGGTTATCATAGAGTGACGTTTTGCCTCTGACATTATGAGGTGCGTGGCCGGGCCGACATTCGACCTCTGTAACTATGCGGGATTCCGTAGTATAGCCAATCTGCGGACCTGTGGTGCGGCGGGTATTAGTGGTATACCAGTTCTACGTTATTGTATAACCTCCTTCCATCAGGGGGGATGGTAGTATAAGTTTTGCCGGCTTGATAGGACCTTGGCCTAATCATACAAAAGAACGAGGCGTTTGGGCATGACGAGACCGCGAATCAGATTAGACGCACCCGCGGAGCCACCTGATGGTCACCGGTTGGATGCCTAGGGAGCTGCTATGATCTCGAAGGCGAGGTTAGATGAGTATGTTGGGGAACGGGTGCCCGGAGGAGTGGAGTTGGTGCTTGCGTGAGCACGATCTCATAGAGGAATCTCGCGGCTCGATAACTGATTTAAAACGCACGGCCCTCCATCTCCTCTCCGGATGCTTAAAGAAATAGAAGCGTCGAGTGCGACAGGCATAAGACCATGGGGACCCTAAGTAGCCCCCTCGCCTCAAATATCTTCCTACTATGCCGGAAAAGATACGCCCTTTGCAAAAGCCATGTTGCATCCCAGTTTATAACAGTCCGCTAAACAGACCCTTACTTGGAAAGGCCGTCACCAATGACGCCCTTGGTTGATTGCTGGGTGCGCCGAATGAATACCCCATGGGTGTAATAGGATCGCTGAGCTCGTGAGTGAGATCGGTCTTCGCCGTCCGTCCCCGAGTGTACCATCAGGAACCAACTAGCAACATAGTCTCCATGGCTCGATCAAACGGTTCTATATCTGGGGATATGATGTCAGGAATTAGCGCTTCATATAACTACTTTGGCGTAGGAGGGTAACACCTAGGGATGCA"
a_count=0
g_count=0
t_count=0
c_count=0
for i in range(0,len(s)):
    if s[i] == 'A':
        print(s[i] == 'A')
        a_count +=1
    elif s[i] == 'T':
        t_count +=1
    elif s[i] == 'G':
        g_count +=1
    if s[i] == 'C':
        c_count += 1
print str(a_count)+" "+str(c_count)+" "+str(g_count)+" "+str(t_count)



dna_string="GGAATCGGGCTAACTACAAAACATGACCAGAATCATCCCGAGGATCGTACCAAACTGGGCGTCGTTCGCTGGGAAGGGAGTAATATCCCATCTATAACCGTCGGCCCATTATGCACAAGGCACCATCCGATAGACTCGGGCGTCGACCCGCTTCGAGGCCACAGATCTCAGGTGACCGGCCGACGATAGCGTTGCGTTGCGCTGTCTCTGTGGGTCCTCCTCTCTTGCTATGGGATCGAAAGCGATCAAGATTGCAAGAACATCCAAGCTGACCATTGCACGTATGCCTGAGAGCGCCGAAGTCTGGGATGAGACGCGGACGAATCCTGAGCATTGTTCCCTCTGGTCGAGAGGATCACCGTAGTTGCTGACACGAGCTGAAAGGTACCCTGAGGTAGGGAGTCGGCCGGGAACTAATTTTCTATTTTAGAACAAATCGGCGGACGGACGCCGTTTGCCCCGGCTAGCCGATTTACACCAACGCTAAGGTGGGTTGGTCGTTTGCATGGATCCATGTCAGACACCAAGTAGATTGGTAGCTTCTTTTCGCCAGCCAAACGCTGTTAGGAGACTATTGTGGGACTTGGGGGCTTTAGCTGGATGGATGGCCCTTCTAATGCCCCTCTGCCCCCGGGGAGGGTTAAAGCTCGCAGTGTAACCAGGGCGGGTTGTGCAGATGGTACGGCGAGGCTGTGGCGGCCATACTGAACAGGCACTGCGCCCGTGCTAGGAAATCTCTATTTAGGCTTGCCCCCCCTGAAACAACACGATCGGAGTGCTAGTCGAGACCATGCTACAGTCCATATTGCAATGATCTAACATCGTGCCCGGAGCCATTACAGTCTCCCCCGGAATTTTGGTCAAAAATGGTAAGATCGAGTGAATTGAACGCTTCAAGAGTAAT"
final_str=""
for str in dna_string:
    if str=="T":
        replace_str="U"
    else:
        replace_str=str
    final_str+=replace_str
print(final_str)


dna_string="CATTATGCGACAGGATTGCCAACATCCGTGTCACCAGCCATATCTGCACGACGCAGGGCCAACATCCTGACTAGGGCGGACGACCGGCAGGGCCCTGAACATCTTGCTCGAGAGCCAGGCATTAGCTTCGTATTCTGTCCGATACTCAATTGCAACTAAGACAGCTTTCGCCGTAGGTCTAGCCAACTAATGGATACCGTTTTGTCGAGAGGATTGGAGGCCAAGATGGCTCTGCAATATATAGCGGTATCAACAACCCACCTATCGTGTGTGCATCTGAGCAATGTTGAGCTACTACCACTGGACATTGCAAGGATGTGAGAATCTGGATCGTGTCTTTTAACGTGGTCTGATTGTAGCTCGCATTTAACCTAGCTTCTTGCGATGGTCAAGGATATACAGTTCAGACAACCGGCGATCAAACATGTGGCCTCTATGACCGATACCATTGCTGGGCATCTATCCCGGAAAGATTAGCCCCCCCGTTGATGGATGTGCCCGAAGCGGAAGGTGTGCATCGGGTCCGACCGATTTCAACTCAATTAGTGTTAAGAACCCACTAATTTCGATGGACGCACGTCCTTTTGGGGTCGGTATAGTGCGCCGAGATCCCGCCTGTTAAGAGCCCATCACTGTTTTCGTAACGAGTCTCTTCTGACCTCAGTCAAACAACTCTACTAGAACACTCGCTTTGCGACGATCCACCACTGAGCTCCGCGCGCAGAATACGACGACTGCTCGATTGTTTTAGCTGCCGCATTTCCTGGTAATCGTCATTTCAAGTACCAGCAGTTTACATCAACTCGATTAGTAGCGATATGCTGCCCGCAGAGTGCGCATATTATCGTAAGGGGGCCCCGCTGGTTTTTATCTGTCGTACCCAGGTTAG"
reverse_string= dna_string[::-1]
complement_dict={
    "A":"T",
    "T": "A",
    "G" :"C",
    "C" : "G"
}
reverse_complement=""
for char_str in reverse_string:
    reverse_complement_chr= complement_dict[char_str]
    reverse_complement+=reverse_complement_chr
reverse_complement


print "fobonocci-series"
def fibinocci(n,start):
    if not start:
        start=0
    # print start,"start"
    if start==0:
        elem=[start,start+1]
    else:
        elem=[1,1]
    for series in range(1,n-1):
        previous_elem=elem[-2]*start
        s= previous_elem+elem[-1]
        print(elem[-1],elem[-2],elem)
        elem.append(s)
    return elem
series=fibinocci(29,2)
print series

def fib(n, k):
    a, b = 1, 1
    for i in range(2, n):
        a, b = b, k*a + b
    return b
fib(29,2)


def fib(n, k):
    a, b = 1, 1
    for i in range(2, n):
        if i==k:
            a, b,c = b, (a+b)-1
        if i> k:
            a,b= b,a+b-1
        else:
            a,b=b,a+b
        print a,b
    return b
fib(85,18)
n = 94
m = 20
bunnies = [1, 1]
months = 2
while months < n:
    if months < m:
        bunnies.append(bunnies[-2] + bunnies[-1])
    elif months == m :
        bunnies.append(bunnies[-2] + bunnies[-1] - 1)
    else:
        # print(bunnies,m,bunnies[-(m + 1)])
        bunnies.append(bunnies[-2] + bunnies[-1] - bunnies[-(
            m + 1)])
    months += 1

print(bunnies[-1])


# gc content
gc_content_dict={}
dna_fasta_file=open("dna_fasta.fa","r")
line_str=""
fasta_seq_dict={}
for line in dna_fasta_file:
    line_str=""
    if line.startswith(">"):
        dna_id=line.strip(">").rstrip()
    else:
        line_str=line.strip()
    if dna_id not in fasta_seq_dict.keys():
        fasta_seq_dict[dna_id]=line_str
    else:
        fasta_seq_dict[dna_id]+=line_str
for dna_id, dna_str in fasta_seq_dict.items():
    length_str=len(dna_str)
    g_count=dna_str.count("G")
    c_count=dna_str.count("C")
    gc_content=g_count+c_count
    gc_content_percentage=(float(gc_content)/float(length_str))*100
    gc_content_dict[dna_id]=gc_content_percentage
print(max(gc_content_dict.items(), key=lambda k: k[1]))



#find hamming distance
print [ a!=b for (a, b) in zip(s1, s2)].count(True)
string1="CCCTCGCTATCGGACTCAGTAAATTCGGAGAGCTTGGTTACGTCCTATTGTAAGCTTATATGATAATATGTGCTGGTTATCGTCGCCACGTGATCCACGAGTACTACCAGATGACCTATCCTCCAAACGGTCTATGTGGCTACATAGCGTGCATTCAATCCCATTCTGAATGGACCAACTGAAAATGGTGAACATATAACGGGTCATCCTGGGTAGTATCCGGTGATCCAATACAATTTAGTACTACGGCTGCCAGCGTCCTCCGTCAAGCTGGACTGACGACTACGCACGCTCAGTTGCCGAAAGATGACCGACCGGCAGGAAAAGGCTATGGCAATTATTGCGGACGTACGCTGATTAGGGTCCCGCTAAGTCTCCCGCTGAAGGGCCACTACACGCCCCAACTGCCGTTTAAGTAAAAAGTACACAAAATGTATGGTTGGTTAACCTGAATCTAAGTCTTCCTCAGAGTCTAAGAAGTTAGTCATGTATTACTCCGCCGCGCAGCGCCTACAACTCCGCAATGTTCGAATTCTTCACTGCAACCCACGGTTGGACAGCGAGCTAATGACGCGCCTAGTTGGGCAAGCTGCCACGAATTAAGCGCGCCTCGATTGCGAGCAGCCCATTCCGATCCCGCGCACCGCGCGAGGATCCCCCCTCCAGTTCTCCTAGGGCCCATAGAGGCCTGAACTGGTCAGCAACGATAATGGGATTCCCAACGTGGCTTTTACCAGGGCCCTGCGGAGACCCAGTAATGCGCAGCCCATTGGAATTTACGGTATATTTCTCGTTAGATCTACCGATCGTCGACGCTTGGTAGGCATGAGTTTTGCTGCTAGAGCACTTTCATCATCCTGCACTAGCTCTACAGAATCGTGCCAAGTTACACCGAATCCTCACCATTAGAATCCCCATCGATCCCGGCTCGAGGCGA"
string2="CTCTTGATGTGGGTTCCCTGCAAGCCGGATCGCTCCGGCATGTCCTAGTTTAACTAGGTATCTCAGAAGGACCGTGGGGTCCTCGTCAGGCCTTCTACCAAAATTACCAGAAGAACTATTAGTCGAAATTTCGGAGTAACTACTATTCGTTGATTACTTACTTCTCTACACGTGCTTAATCAATACAGTGAACAAAAGACGTAGTATCCGTTTTTTTAGACGGTTATCCATTGCAAAGTCGTAGGATTGGCGCCTGTGTCCTGCGAAACTCAGTACGGCCACGTACGAACGTCCAGTTTTCGGAGTATTGCCAAAGGAGGTGCAAGGGCTATTGCCATTATAGAGTATCCCACTCTATCTTATCACCGCTAATGCCCGGGATAGGGCTCCATTGTCCGCCTGACCTCCAGTGGTAGTATTTCTTCCCCTCACGGTAAACTGGCAGACTCTCTTTACCGATTTTGATCGGTGACTATGAATTTAAGCAGTAATTACAGTGACTGGCCACGCTTTCCGTTTCTAAAAGCCCGACGTCCTTAATGGGAACAGCCGTTTGATTATGCTCCGAAGACGTGCTTAGTTGGTTTGGACGCCACGGATGGAGTGCACGACGCTTCCCAGCATCATCCCTCGGATCGGCGCACTGATTTGTTGTCTCTAGTACAATAGATCTTCGCCCCGTATGGGCCTGGATTGTAAACCCAGAACAGCGAAATACCTAACGTAACTTGACCACCCAACCTAGGGAGGCCTTCCGTAGTTAAGGCCATTTCACTCGGCCGTAAACATATCGGTTCTAATACCGAGGCTCTTCCGCTAGACGGCATTAGGTTCTTTCCTGGCATACTTTCCTTATGTCGTCAGAGCTAGGTCGATCCGTGCCAAACTAGAGCGAACCCTAAAATGTAGGATAAGCATACGCCCTAGAACGAATCCA"
mismatches=0
for i, (char1, char2) in enumerate(zip(string1, string2)):
    if char1!=char2:
        mismatches+=1
print mismatches

#mendelian inheritence
k=15
m=20
n=23
a = k + m + n
p_recessive = ((0.25*m*(m-1)) + (0.5*m*n) + (0.5*m*n) + (n*(n-1)))/(a*(a-1))
p_wanted = 1 - p_recessive
p_wanted = round(p_wanted, 5)
print(p_wanted)


#find the motif

sting1="ATAGGCTGGTAACTTGTAACTTGTTAACTTGCCATAACTTGTTTAACTTGGTAACTTGCTAACTTGGGGTAACTTGCTGAGTTCTAACTTGCTTTGTAAGTAACTTGTTAACTTGTTAACTTGTAACTTGCCGTAACTTGATTAACTTGGGTAACTTGCTAATAACTTGGCGCGACGGTAACTTGCTAACTTGTTAACTTGTAACTTGGTAACTTGGGGAAGATAACTTGGAGTAACTTGATTAACTTGTGTACATAACTTGTCCGTGCATAACTTGACACTTAACTTGCTAACTTGGGATATAACTTGTAACTTGGGTGTGTCTAACTTGTGATTAACTTGTAACTTGTAGTGGTAACTTGGTTAACTTGAACTAACTTGGTAACTTGCTAACTTGCCATACCGTTGCGTAACTTGTTTGCGAGGTAACTTGTTAACTTGGGTAACTTGAGTGGTAACTTGAAGAGTAACTTGTAACTTGGTAACTTGGTCTAACTTGATAACTTGTAACTTGGTTCTAACTTGTAACTTGTTAACTTGTAACTTGCGCCCGGTAACTTGAATAACTTGCTAACTTGCTAACTTGAATCGTAACTTGGTAACTTGATAACTTGACACATAACTTGTAACTTGTAACTTGCCTAACTTGATTAACTTGTTAACTTGTTCCCTAACTTGGCCTGTAACTTGGGTAACTTGATTCGAATTAACTTGACTAACTTGGCAGTAACTTGACCTGTTGATTAACTTGTAACTTGTCTAACTTGCCGTTAACTTGCGTTAACTTGGCTAACTTGCGTCATTTATCCGTATTAACTTGTTATAACTTGTAACTTGTAGAGGTTAACTTGGAATAACTTGGGGTAACTTGAACTAACTTGGATAACTTGCAGTTAACTTGCTTAACTTGACTAACTTGGTTAACTTGGGTAACTTGACGGGAGCGTAACTTGGTAACTTGTAACTTGTATAGCGTAAGTAACTTGTTTAACTTG"
motif="TAACTTGTA"

for char_idx in range(0,len(sting1)+1):
    if motif ==sting1[char_idx : char_idx +len(motif)]:
        print(char_idx+1 )


dna_seq=open("dna_Seq.fa",'r')

def fasta_to_list(dna_seq):
    fasta_seq_dict = {}
    for line in dna_seq:
        line_str=""
        if line.startswith(">"):
            dna_id=line.strip(">").rstrip()
        else:
            line_str=line.strip()
        if dna_id not in fasta_seq_dict.keys():
            fasta_seq_dict[dna_id]=line_str
        else:
            fasta_seq_dict[dna_id]+=line_str
    seq_list=fasta_seq_dict.values()
    return seq_list
seq_list=fasta_to_list(dna_seq)
result = [] # list to save nucleotide with max count from each "row"
n = len(seq_list[0]) # length of each sequence

profile_matrix = {
  'A': [0]*n,
  'C': [0]*n,
  'G': [0]*n,
  'T': [0]*n
}

for dna in seq_list:
    print(dna,enumerate(dna))
    for position, nucleotide in enumerate(dna):
        print(position, nucleotide)
        profile_matrix[nucleotide][position] += 1
profile_matrix
for position in range(n):
    max_count = 0
    max_nucleotide = None
    for nucleotide in ['A', 'C', 'G', 'T']:
        count = profile_matrix[nucleotide][position]
        if count > max_count:
            max_count = count
            max_nucleotide = nucleotide
    result.append(max_nucleotide)

consensus = ''.join(result)
consensus

def read_fasta_file(fa_file):
    dna_seq = open(fa_file, 'r')
    seq_list = []
    fasta_seq_dict = {}
    for line in dna_seq:
        line_str = ""
        if line.startswith(">"):
            dna_id = line.strip(">").rstrip()
        else:
            line_str = line.strip()
        if dna_id not in fasta_seq_dict.keys():
            fasta_seq_dict[dna_id] = line_str
        else:
            fasta_seq_dict[dna_id] += line_str
    return fasta_seq_dict

def is_k_overlap(s1, s2, k):
    return s1[-k:] == s2[:k]

import itertools
def overall_graph(dna_Seq_file,k):
    fasta_seq_dict=read_fasta_file(dna_Seq_file)
    edges=[]
    for u, v in itertools.combinations(fasta_seq_dict, 2):
        u_dna, v_dna = fasta_seq_dict[u], fasta_seq_dict[v]
        if is_k_overlap(u_dna, v_dna, k):
            edges.append((u, v))

        if is_k_overlap(v_dna, u_dna, k):
            edges.append((v, u))
    return edges



dna_Seq_file="overlap_graph.fa"
overlap_graph=overall_graph(dna_Seq_file,3)
for seq in overlap_graph:
    print(seq[0] + " "+ seq[1])



filepath="rosalind_iev.txt"
displayDominant = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
offspring = 2

with open(filepath) as file:
    parentCounts = [int(x) for x in file.read().split()]

print sum([offspring * x[0] * x[1] for x in zip(displayDominant, parentCounts)])


from collections import Counter
common_string_file=open("common_substring.fa","r")
common_string_list=fasta_to_list(common_string_file)
l = sorted(common_string_list, key=len)[-1]
string_length=(len(common_string_list[0]))
def common_substring(string_list):
    res= ""
    for i in range(string_length):
        for j in range(i + 1, string_length + 1):

            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = common_string_list[0][i:j]
            if len(stem)>2:
                stem=stem[-2:]
            k = 1
            for k in range(1, len(common_string_list)):

                # Check if the generated stem is
                # common to all words
                if stem not in common_string_list[k] or len(stem)==1:
                    # print(stem,common_string_list[k])
                    break
                elif len(res)<len(stem):
                    res=stem


            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if (k + 1 == len(common_string_list) and len(res) < len(stem)):
                print(stem, common_string_list[k])
                res = stem
    return res
common_substring(common_string_list)



def make_motifs(dna_string):
    motifs = []
    for i in range(0, len(dna_string)+1):
        for j in range(i+1, len(dna_string)+1):
            motifs.append(dna_string[i:j])
    return motifs

seqs = fastr.read("file.txt")
seqs=common_string_list
motifs = make_motifs(common_string_list[0])
saved_motif = ""
cont = 0

for motif in motifs:
    for i in range(1, len(seqs)):
        if seqs[i].find(motif) != -1:
            cont += 1
    if cont == len(seqs)-1 and len(saved_motif) < len(motif):
        saved_motif = motif
    cont = 0

print saved_motif



#independent varaibles

import math
k = 7
N = 34
P = 2**k
probability = 0
for i in range(N, P + 1):
    prob = (math.factorial(P) /
            (math.factorial(i) * math.factorial(P - i))) * (0.25**i) * (0.75**(P - i))
    probability += prob
print(probability)

rna_table = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
"UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
"UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
"UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
"CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
"CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
"CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
"CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
"AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
"ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
"AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
"AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
"GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
"GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
"GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
"GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

with open("rosalind_mrna.txt") as myfile:

  data = myfile.readlines()

charData = list(data[0].strip())


frequency_list = []
for i in range(0,len(charData),3):
    if i <=(len(charData)-3):
        print(charData[i:i+3],i,i+3)
        codon=rna_table["".join(charData[i:i+3])]
        if codon== 'STOP':
            break
        frequency_list.append(codon)
"".join(frequency_list)



# caculate permutation of a number

n=6
final_result=1
for i in range(1,n+1,1):
    final_result=final_result*i

final_result

import numpy as np
import itertools

n = 7
order = np.arange(1, n+1, 1)
count = 0

perm = list(itertools.permutations(order))
print(len(perm))

for i in perm:
    count += 1
    ordered_list= " "
    for x in range(n):
        if x < n - 1:
            print str(i[x]),
        else:
            print(str(i[x]))


monoisotopic_mass_table={
"A" :  71.03711,
"C" :  103.00919,
"D" :  115.02694,
"E" :  129.04259,
"F" : 147.06841,
"G" :  57.02146,
"H" :  137.05891,
"I" :  113.08406,
"K" :  128.09496,
"L" :  113.08406,
"M" :  131.04049,
"N" :  114.04293,
"P" :  97.05276,
"Q" :  128.05858,
"R" :  156.10111,
"S" :  87.03203,
"T" :  101.04768,
"V" :  99.06841,
"W" :  186.07931,
"Y" :  163.06333,
}
mass=0
sample="KRGWYAWDHQKEYQMMMKDMKKPDDTTRSIMGMMPCACIVKPDVRMAKIKHRVRCCVDDDRYQPREAFWDWRGYAVNKMIGEHRCTEVKIWPKGAPIMMHNITDAGITGDGQFQRLLGQRTFMVVDDHLCQLCFCQDGQIPANTRHTQWVLIIFQTIIGAGWMLTECPRQNGFLSVCRVYDIDYEELLMYFRVNTETHVTLYCEQGRMEGWVRSAHPQKVTGTQFLKTNVDKAETRWEWEFFSPNCRKAVSVHSHDKMSKGMNYAGACIKQPCHDDSKALGHYTCEMVVQMWKSFMQCAKLDFQAISQEVGREQTTYQIWVLRSLLDQEKWIMPFKPNIETNGPRLWVMFQHCVESGPGHHFKRCKEVTASNSLMEQQRNTWQAECALFPLHMWNVHWTDTEAQKRCDSLMLEYFGDAMIEPMGMHTSFLVDMRRYARGAPTIKIVYVVIIAAFRIKHVQENLELRWAIPTFDCNYMHLQAHFEPTKGYWTFDVTCTFSGLEFDFDDEEPGMKYCSHPHKYPHKPHQEHWQQGTIPWVPNQHSCLFRDEEGEASWSGLGGECWALGNMGKRPTPMNWVETSRHIVHTNEVGDWFVADFIYLKWEKTIQITDHSPAWVKGWEEHMQDQANFIRHEHVCLSLMVFQPNPTEFHTEEKMNQTEYLGYWDNWLNERDKKYYPGDIMNAVVAIVASCNTPLWYPLTMWWAWAGMAPLGVMKYTGGHARTNYLYLNFCLVGAGGLLYSAYTSRTPDCKAGMQIKRVKPSNWPQVQKFYISIHKWKKGEVMFYYIWESEFSGLEDMVMWQNHIHQTRIISHNISDYKRQCFVLGMICIPGFHCFIVYNENYYNSTDPPIEATRHRWDQPTRTANITFLHPHTSTFCNRFYSRQVVASAKGQASR"
for char_aa in sample:
    if monoisotopic_mass_table[char_aa]:
        mass+=monoisotopic_mass_table[char_aa]
round(mass,3)

dna_seq=open("rev_complement.fa",'r')

seq_list=fasta_to_list(dna_seq)
dna_string=seq_list[0]
Result=[]
for nuc in range(len(dna_string)-4):
    length = 3
    while length != 14:
        length += 1
        if (length + nuc) > len(dna_string):
            continue
        reverse = dna_string[nuc:nuc + length].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c')
        reverse = reverse.upper()
        relist = list(reverse)
        relist.reverse()
        reverse = ''.join(relist)
        if dna_string[nuc:nuc + length] == reverse:
            print(reverse, dna_string[nuc:nuc + length])
            Result.append([nuc + 1, length])
for i in range(len(Result)):
    seq = str(Result[i][0]) + ' ' + str(Result[i][1])
    print(seq)


dna_seq=open("exon_seq.fa",'r')
seq_list=fasta_to_list(dna_seq)
sorted_seq_list=sorted(seq_list, key=len)
# seq_list=[]
# charData=""
# for line in dna_seq:
#     if not line.startswith('>'):
#         seq_list.append(line.strip())
introns=sorted_seq_list[0:-1]
# introns.sort(reverse=True)

exon_data = sorted_seq_list[-1]
frequency_list = []
for intron in introns:
    if intron in exon_data:
        print(intron)
        exon_data=exon_data.replace(intron,'')
        print(exon_data)
exon_data=exon_data.replace('T','U')

for i in range(0,len(exon_data),3):
    if i <= (len(exon_data)-3):
        codon=rna_table[exon_data[i:i+3]]
        if codon== 'STOP':
            break
        print(codon,exon_data,i,exon_data[i:i+3])
        frequency_list.append(codon)

"".join(frequency_list)


#enumerate k-mers

sample_string = "A B C D E F G H I"
k_mer=3
from itertools import permutations
sample_string=sample_string.replace(" ","")*k_mer
permutations(sample_string, k_mer)
kmer_combinations=[]
for kmer_str in permutations(sample_string, k_mer):
    kmer_combinations.append("".join(kmer_str))
sorted_unique_kmer_combinations=sorted(list(set(kmer_combinations)))
for combination_kmer in sorted_unique_kmer_combinations:
    print(combination_kmer)

# Using flask to make an api
# import necessary libraries and functions
from flask import Flask, jsonify, request

# creating a Flask app
app = Flask(__name__)


# on the terminal type: curl http://127.0.0.1:5000/
# returns hello world when we use GET.
# returns the data that we send when we use POST.
@app.route('/', methods=['GET', 'POST'])
def home():
    if (request.method == 'GET'):
        data = "hello world"
        return jsonify({'data': data})