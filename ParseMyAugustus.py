#! /usr/bin/env python3

# augustus file parsing


class Augfile:

    def __init__(self, name):
        self.name = name
        self.contig = name.split('/')[1].split('.')[0]
        self.genes = [] # list of genes in this augustus file

    def read_genes(self):
        # reads all genes into the genes list
        genenum = 0
        raw_gene = ''
        with open(self.name, 'r') as infile:
            for line in infile:
                if line.startswith('# start gene '):
                    # start grabbing lines at every line of '#start gene'
                    genenum = line.rstrip('\n').split('# start gene ')[1][1:]

                    nline = infile.readline()

                    while nline.startswith('# end gene') is False:
                        raw_gene += nline
                        nline = infile.readline()
                    # until the line '# end gene' appears
                    # create a Gene object with raw lines - parse what we need later
                    gene = Gene(raw=raw_gene, contig=self.contig, num=genenum)

                    # append to gene list of augustus file and reset the gene content
                    self.genes.append(gene)
                    raw_gene = ''


class Gene:

    def __init__(self, raw, contig, num):
        # init with contig and gene number
        # raw contains everything in augfile between '# start gene' and '# end gene'
        self.contig = contig
        self.num = num
        self.raw = raw

    def getSeq(self):
        # grab the predicted protein seq and clean up
        raw_split = self.raw.split('protein sequence = [')[1]
        raw_split2 = raw_split.split(']\n# Evidence')[0]
        prot_seq = raw_split2.replace('\n# ', '')
        self.seq = prot_seq

    def grabFeatures(self):
        # parse 'raw' into lines that contain features
        # all lines that are not commented out and not empty
        rawList = self.raw.split('\n')
        feat = [i for i in rawList if i.startswith('#') is False and len(i) is not 0]
        self.features = [Feature(line=f, contig=self.contig, geneNum=self.num) for f in feat]

    def writeSeq(self):
        # create a fasta header and write the protein sequence 
        head = f'>{self.contig}_{self.num}'
        print(head)
        print(self.seq)

    def writeFeatures(self):
        # modify the last field of each gtf entry in the gene
        # write gtf section for the whole gene, incl some padding
        print(f'### {self.contig}_{self.num}')
        for f in self.features:
            f.modField()
            print(f.line)


class Feature:

    def __init__(self, line, contig, geneNum):
        # parse one line in the augustus gtf file
        l = line.split('\t')

        self.geneID = f'{contig}_{geneNum}'
        self.contig, self.source, self.feat, self.start, self.end = l[0], l[1], l[2], l[3], l[4]
        self.score, self.strand, self.frame, self.comment = l[5], l[6], l[7], l[8]

    def modField(self):
        # change the last field of the gtf entry depending on feature type
        if self.feat == 'gene':
            self.commentNew = self.geneID
        elif self.feat == 'transcript':
            tN = self.comment.split('.')[-1]
            self.commentNew = f'{self.geneID}.{tN}'
        elif self.feat in ['start_codon', 'stop_codon', 'CDS', 'intron']:
            tID = self.comment.split(';')[0].split(' ')[-1].split(".")[1].rstrip('"')
            tID_new = f'{self.geneID}.{tID}'
            self.commentNew = f'transcript_id "{tID_new}"; gene_id "{self.geneID}";'
            
        else:
            print(f'warning: unknown feature in gene - {self.geneID} {self.feat}')

        # construct new gtf entry with modified last field
        self.line = '\t'.join([self.contig, self.source, self.feat, self.start,
                              self.end, self.score, self.strand, self.frame, self.commentNew])
        
        