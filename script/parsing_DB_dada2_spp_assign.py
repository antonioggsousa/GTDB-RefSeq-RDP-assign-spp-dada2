import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


parser = argparse.ArgumentParser(description='Parsing assignTaxonomy() dada2 \
											formatted databases for addSpecies() \
											assignment')
											
parser.add_argument('-i', action="store", dest="input", type=str, help='Input database file path to parse')
parser.add_argument('-o', action="store", dest="output", type=str, help='Output database file path parsed')
parser.add_argument('-f', action="store_true", default=False, help='Trim/Filter database')

args = parser.parse_args()
#print('i     =', args.input)
#print('o   =', args.output)


DB = list(SeqIO.parse(args.input, "fasta"))

fasta_headers = []
fasta_seqs = []

for seq in DB:
	fasta_head = seq.id
	fasta_seq = seq.seq
	fasta_seqs.append(fasta_seq)
	fasta_headers.append(fasta_head)

fasta_id = []
species_name = []
count = 0 

for header in fasta_headers:
	count += 1
	if len(header.split(';')) < 7:
		spp = 'unclassified species'
		seq_id = 'UNKNOW_'+str(count)		
	else: 
		sp = header.split(';')[6]
		sp_seq_id = re.split('\(|\)', sp)
		if sp_seq_id[0] is '':
			sp_seq_id[0] = 'unclassified_species'
		spp = sp_seq_id[0].replace('_', ' ', 2 )
		spp = spp.replace('"', "") # replace the character '"' by nothing for RefSeq-RDP DB
		spp = spp.replace("'", "")
		seq_id = sp_seq_id[1]
	fasta_id.append(seq_id)
	species_name.append(spp)

species_fasta = []

for spp_fasta in range(len(species_name)): 
	species_fasta.append(SeqRecord(fasta_seqs[spp_fasta], id=fasta_id[spp_fasta], description=species_name[spp_fasta]))
	

unclassified_spp = []
classified_spp_cld = []

for sp in species_fasta:
	if sp.description == 'unclassified species':
		unclassified_spp.append(sp)
	else: 
		classified_spp_cld.append(sp)


upper_case = list(range(65,91))
lower_case = list(range(97, 123))	

unclassified_cld = []
classified_spp = []


for sp in classified_spp_cld: 
	if (ord(sp.description[0]) in upper_case) and (ord(sp.description[1]) in lower_case):
		classified_spp.append(sp)
	else:
		sp.description = sp.description + str(count)
		unclassified_cld.append(sp)


SeqIO.write(species_fasta, args.output + '_ALL.fa', "fasta")
SeqIO.write(classified_spp_cld, args.output + '_CLASS_SPP_CLD.fa', "fasta")
SeqIO.write(unclassified_spp, args.output + '_UNCLASS_SPP.fa', "fasta")
SeqIO.write(classified_spp, args.output + '_CLASS_SPP.fa', "fasta")
SeqIO.write(unclassified_cld, args.output + '_UNCLASS_CLD.fa', "fasta")


if args.f is True:
	
	seq_len = []
	
	for len_seq in classified_spp:
		seq_len.append(len(len_seq.seq))

	trim_seq_len = sum(seq_len)/len(seq_len) * 1/2

	classified_spp_trimmed = []

	for len_seq in classified_spp:
		if (len(len_seq.seq) > trim_seq_len) and (len_seq.seq.count('N') == 0): 
			classified_spp_trimmed.append(len_seq)
			
	SeqIO.write(classified_spp_trimmed, args.output + '_TRIM_CLASS_SPP.fa', "fasta")

