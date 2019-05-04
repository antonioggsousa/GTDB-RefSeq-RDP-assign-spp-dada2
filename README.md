## GTDB and RefSeq-RDP DBs parsed for the *addSpecies()* dada2 function

<br>

This tutorial describes how to parse dada2 formatted databases, such as 16S rRNA gene **GTDB** and **RefSeq-RDP** databases, 
for [*addSpecies*](https://rdrr.io/bioc/dada2/man/addSpecies.html) [**dada2**](https://benjjneb.github.io/dada2/tutorial.html) function.

The focus will be **GTDB** and **RefSeq-RDP** databases, but it can be generalised to any other dada2 formatted database. 
 
[dada2 formatted databases](https://benjjneb.github.io/dada2/training.html) consist of a multi-fasta file where each fasta sequence is identified by 
the full taxonomy path (with each taxonomic rank separated by **;**) in the fasta header, like this: 


    >Kingdom;Phylum;Class;Order;Family;Genus;Species(seq_id)
    ATGCATCGGCTAGACTGATGACGTGCTGCTAGTTAGCTGCCTAGTCGTAGCTAGCTAGCTGTAGCATGCATCGTGCTCGATCGTAGCT 

If any database formatted for dada2 contains references sequences for which the taxonomic path is shorter, i.e, the reference sequence 
was only confidently assigned up to Class, Order or Family level, any environmental sequence matching these references will 
be only classified up to that level, e.g., Class, Order or Family level. The **seq_id** (between brackets) stands for the accession number of the 
reference sequence in the database. 

This format works well with the [*assignTaxonomy*](https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/assignTaxonomy) **dada2** function that makes use of **RDP classifier** to classify environmental sequences. 
If the reference sequences contain the full taxonomic path up to Species level, the environmental sequences matching this will be classified up to 
Species level, otherwise they will be classified up to the taxonomic level available (e.g., Class, Order or Family level). 

However, if **dada2** users prefer to use the *addSpecies* function that assigns taxonomy only at Species level based on 
100% exact matches, as far as I know, there is only available the [**16S rRNA gene Silva v132**](https://zenodo.org/record/1172783#.XMte3SOZPOQ). 
Therefore, this tutorial describes a **python** script to parse dada2 formatted databases (as described before) for the *addSpecies* **dada2** function. 
The script was developed for and with the purpose of parsing the 16S rRNA gene **GTDB** and **RefSeq-RDP** databases for the *addSpecies* **dada2** function, 
but it should function well with any other database formatted for dada2  (may be after some adjustments!). 

The database formatted for the *addSpecies* dada2 function differs from the general database formatted for the *assignTaxonomy* dada2 function, in the sense that only contais the **seq_id** followed by the Species name, i.e., generic name/genus and specific name/species, like this:

    >Seq_id Genus species
    ATGCATCGGCTAGACTGATGACGTGCTGCTAGTTAGCTGCCTAGTCGTAGCTAGCTAGCTGTAGCATGCATCGTGCTCGATCGTAGCT

Therefore, the only difference of the database formatted for the *addSpecies* function compared with the database formatted for the *assignTaxonomy* function, relies on the number of taxonomic levels, the former it has only one - Species level - and the latter could contain seven levels, i.e., Kingdom-Species. Thus, from the database formatted for the *assignTaxonomy* function it is possible to parse it to use it for the *addSpecies* function. 

 <br>

 ### *parsing_DB_dada2_spp_assign.py* python script

 The *parsing_DB_dada2_spp_assign.py* python script that I made to parse a database formatted for the *assignTaxonomy* function into a compatible database formatted for the *addSpecies* function relies on the following:
  
  + [Python](https://www.python.org) 3.6.5
  + [Biopython](https://biopython.org) 1.72

#### How *parsing_DB_dada2_spp_assign.py* works?

The *parsing_DB_dada2_spp_assign.py* python script can be inspected inside the **script** folder, in this GitHub repository. 

Download it and type:

    parsing_DB_dada2_spp_assign.py -h

This will give as output how the script functions: 

    usage: parsing_DB_dada2_spp_assign.py [-h] [-i INPUT] [-o OUTPUT] [-f]

    Parsing assignTaxonomy() dada2 formatted databases for addSpecies() assignment

    optional arguments:
      -h, --help  show this help message and exit
      -i INPUT    Input database file path to parse
      -o OUTPUT   Output database file path parsed
      -f          Trim/Filter database

The script takes as input **-i** the **Input database file path to parse**, i.e., a database formatted for the *assignTaxonomy* function (in fasta format), and it gives as output **-o** the **Output database file path parsed**, i.e., a database formatted for the *addSpecies* function. The **-f** argument is optional, and if **True** the **Output database file path parsed** is also trimmed/filter: all the fasta sequences with length < half-average sequence length as well as any ambigous fasta sequence, i.e., N > 0, are discarded. 





 #### A little explanation about the code   

First, the **DB** is imported and, then, it iterates over the **DB** to retrieve two lists: one of fasta headers (**fasta_headers**) and other of fasta sequences (**fasta_seqs**).

    DB = list(SeqIO.parse(args.input, "fasta"))

    fasta_headers = []
    fasta_seqs = []

    for seq in DB:
	    fasta_head = seq.id
	    fasta_seq = seq.seq
	    fasta_seqs.append(fasta_seq)
	    fasta_headers.append(fasta_head)

Then, the main piece of code of this script:

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


This piece of code iterates over the list of fasta headers (**fasta_headers**) to retrieve the species name and sequence ID. It does it by splitting each fasta header (that corresponds to the full taxonomic path) by **;** character. If at Species level there is only the sequence ID and there is not a species names, the species name given to such as sequence is **unclassified species**. If at Species level there is not a species name neither a sequence ID, the species name given to such as sequence is **unclassified species** and the sequence ID is **UNKNOWN_** plus the nº of that sequence in the database. In addition, characters such as **(**, **)**, **_**, **'** or **"** are discarded from the species name. Then the species name of each sequence is appended to the **species_name** list and the sequence ID to the **fasta_id**.


The code below joins the sequence ID (**fasta_id**), species name (**species_name**) and fasta sequence (**fasta_seqs**). This database is parsed for the *addSpecies* function. However, it contains unclassified species. 

    species_fasta = []

    for spp_fasta in range(len(species_name)): 
	    species_fasta.append(SeqRecord(fasta_seqs[spp_fasta], id=fasta_id[spp_fasta], description=species_name[spp_fasta]))


Since the aim of using *addSpecies* function is to get the species name, **unclassified species** names that were given to those sequences without any description at Species level, are not welcome. Therefore, the following piece of code separates those sequences from the remaining. 


    unclassified_spp = []
    classified_spp_cld = []

    for sp in species_fasta:
	    if sp.description == 'unclassified species':
		    unclassified_spp.append(sp)
	    else: 
		    classified_spp_cld.append(sp)

Still there is sequences with unclassified notations such as environmental clades, clones or whatsoever. Usually these names are written in upper case and alphanumeric (**UBA10521**) or in lower case (**uncultured nanoarchaeote**). Therefore the following piece of code tries to retrieve classified species names based on the fact that species names start by a upper case letter followed by a lower case letter. 

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


If the option **-f** is passed, all the fasta sequences with length < half-average sequence length as well as any ambigous fasta sequence, i.e., N > 0, are discarded from the remaining sequences.


    if args.f is True:
	
	seq_len = []
	
	for len_seq in classified_spp:
		seq_len.append(len(len_seq.seq))

	trim_seq_len = sum(seq_len)/len(seq_len) * 1/2

	classified_spp_trimmed = []

	for len_seq in classified_spp:
		if (len(len_seq.seq) > trim_seq_len) and (len_seq.seq.count('N') == 0): 
			classified_spp_trimmed.append(len_seq)


#### Example: GTDB and RefSeq-RDP DBs

The [**GTDB** and **RefSeq-RDP**](https://zenodo.org/record/2541239#.XM2G8SOZPOQ) databases were downloaded and decompressed. Then, each database was parsed as follows: 

    # GTDB
    python3 parsing_DB_dada2_spp_assign.py -i GTDB_bac-arc_ssu_r86.fa -o GTDB_dada2 -f

    # RefSeq-RDP
    python3 parsing_DB_dada2_spp_assign.py -i RefSeq-RDPv2_16S_species.fa -o RefSeq-RDP_dada2 -f

For each database were generated 6 files (compressed inside the folder **DBs**). 

<br>

Starting by the GTDB:


| Output Databases         | Description           | No. of Sequences  |
| ------------- |:-------------:| -----:|
| GTDB_dada2_ALL.fa      | all | 21559 |
|   GTDB_dada2_CLASS_SPP_CLD.fa   |   all species and clades/clones    |  15192  |
| GTDB_dada2_UNCLASS_SPP.fa |   unclassified species   |   6367  |
| GTDB_dada2_CLASS_SPP.fa |   only classified species    |  13163   |
| GTDB_dada2_UNCLASS_CLD.fa |    clades/clones   |   2029  |
| GTDB_dada2_TRIM_CLASS_SPP.fa |   only classified spp. trimmed/filtered    |   12442  |

<br> 

Now for the RefSeq-RDP:


| Output Databases         | Description           | No. of Sequences  |
| ------------- |:-------------:| -----:|
| RefSeq-RDP_dada2_ALL.fa      | all | 15336 |
|   RefSeq-RDP_dada2_CLASS_SPP_CLD.fa   |   all species and clades/clones    |  12828  |
| RefSeq-RDP_dada2_UNCLASS_SPP.fa |   unclassified species   |   2508  |
| RefSeq-RDP_dada2_CLASS_SPP.fa |   only classified species    |  12803   |
| RefSeq-RDP_dada2_UNCLASS_CLD.fa |    clades/clones   |  25   |
| RefSeq-RDP_dada2_TRIM_CLASS_SPP.fa |   only classified spp. trimmed/filtered    |   10899  |

<br>

Personally, I recommend the use of **GTDB_dada2_TRIM_CLASS_SPP.fa** (compressed inside **DBs/GTDB/** directory) and of **RefSeq-RDP_dada2_TRIM_CLASS_SPP.fa** (compressed inside **DBs/RefSeq-RDP/** directory) to use as input of *addSpecies* dada2 function. 





  