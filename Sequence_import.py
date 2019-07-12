import re

def read_fasta(file):
    '''Reads multiple sequences on a fasta file. Returns a list of the sequences'''
    seqs=[]
    seq=""
    lines=open(file).readlines()
    limit=len(lines)-1
    for i in range(limit):
        if lines[i][0] == '>':
            j=i+1
            while j<=limit and lines[j][0]!='>':
                seq+=lines[j].strip()
                j+=1
            seqs.append(seq)
            seq=""
    return seqs

def read_fasta_titles(file):
    '''Reads the titles of multiple sequences on a fasta file. Returns a list of the titles. '''
    seqs=[]
    seq=""
    lines=open(file).readlines()
    limit=len(lines)-1
    for i in range(limit):
        if lines[i][0] == '>':
            seq=lines[i].strip()
            seqs.append(seq)
            seq=""
    return seqs



def read_fasta_specie(file):
    '''Reads the species of the sequence indated between [ ] '''
    seqs=[]
    seq=""
    lines=open(file).readlines()
    limit=len(lines)-1
    for i in range(limit):
        if lines[i][0] == '>':
            seq=lines[i].strip()
            seqs.append(seq)
            seq=""
    species=[]
    for title in seqs:
        sp=re.findall(r"\[(.*?)\]", title)
        species+=sp
    return species
    
def export_fasta(file, seqs, seqsID=None):
    file=open("%s.txt" % (file), "w" )
    if seqsID != None:
        for i in range(len(seqsID)):
            _=file.write("> %s \n" % (seqsID[i]))
            _=file.write(seqs[i] + "\n\n")
    else:
        for i in range(len(seqs)):
            _=file.write("> Sequence %s \n" % (i+1))
            _=file.write(seqs[i] + "\n\n")
    file.close()
