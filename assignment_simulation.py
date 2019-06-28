#!/usr/bin/env python3

import random
import sys
import string


# init some constants
CONTIG_NUM = 100
MIN_CONTIG_LEN = 500
MAX_CONTIG_LEN = 2000
READ_LEN = 100

NUC = ['A','C','G','T']
QUAL = ['A','B','C','D','E','F','G','H','I']
SYM = list(string.hexdigits)



class Contig:
        
    def __init__(self, num):
        # create a header and a random sequence for a contig     
        self.header = f'>contig{num}'
        self.seq = ''.join(random.choices(NUC, k=random.randrange(MIN_CONTIG_LEN, MAX_CONTIG_LEN, 1)))
        
    def __repr__(self):
        return(f'{self.header}\n{self.seq}\n')
        

class Genome:
    
    def __init__(self):
        # genome is a list of contigs
        self.contigs = list()
        
        # init the contigs with an incrementing number
        curr_num = 0
        for i in range(1, CONTIG_NUM + 1):
            self.contigs.append(Contig(num=curr_num))
            curr_num += 1
        

class Fastq:
    
    def __init__(self, genome, MN):
        # fastq file is a set of reads
        self.reads = set()
        self.genome = genome
        
        # the number of genome-derived reads is taken from the matrikelnummer
        read_num = int(str(MN)[-4:])
        
        # generate reads from the genome
        for i in range(1, read_num + 1):
            cont = self.genome.contigs[random.randrange(0, 99, 1)] # pick random contig
            read_start = random.randrange(0, 350, 1) # pick random start within that contig
            read_end = read_start + 100
            read = cont.seq[read_start:read_end]
            
            entry = gen_fastq_entry(seq=read) # generate the whole fastq entry
            self.reads.add(entry)
            
        # additionally add a random number of non-genome-derived reads
        for i in range(1, random.randrange(10, 500, 1)):
            read = ''.join(random.choices(NUC, k=READ_LEN))
            entry = gen_fastq_entry(seq=read)
            self.reads.add(entry)
        

def write_file(filename, obj):
    # writes contigs or fastq entries to a file
    with open(filename, 'w') as outfile:
        [outfile.write(str(i)) for i in obj]

def gen_fastq_entry(seq):
    # takes a sequence and creates the other components of a fastq entry
    h = ''.join(random.choices(SYM, k=30)) # header consists of random characters
    head = '@' + h
    comm = '+' + h
    q = ''.join(random.choices(QUAL, k=READ_LEN)) # quality is also randomized - hq only
    return(f'{head}\n{seq}\n{comm}\n{q}\n')

    
# add the real matrikelnummern and some fake ones, to make the students glob again
real_MNs = ['01345185', '01267766', '01260019','1446986','51804732','1345177','1657223','11838220','1306668','1206041','11837680','11732015','11812630','11724735','1309902','1300649','1307172','1225181','1276775','1440875']
fake_MNs = [str(random.randrange(1000000, 1999999, 1)) for i in range(1, 50)]


for i in real_MNs + fake_MNs:
    # generate and write a genome for every matrikelnummer
    student_genome = Genome()
    write_file(filename=f'{i}.fasta', obj=student_genome.contigs)
    
    # generate sequencing data for every matrikelnummer
    student_reads = Fastq(genome=student_genome, MN=i)
    write_file(filename=f'{i}.fastq', obj=student_reads.reads)
    
    

