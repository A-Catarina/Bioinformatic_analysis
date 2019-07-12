from Sequence_import import *
from MySeq import MySeq
from SubstMatrix import SubstMatrix
from PairwiseAlignment import PairwiseAlignment
from MultipleAlignment import MultipleAlignment
from UPGMA import UPGMA
from MyBlast import MyBlast
from MyGraph import MyGraph
import re


print('''\n                                                   Algoritmos para Bioinformática
                                                            2018/2019                      
_____________________________________________________________________________________________________________________________________________''')

#Sequence Import
##Source Sequence

ref_seq=read_fasta("source.fasta") #list of reference sequences
ref_sp=read_fasta_specie("source.fasta") #list of reference sequence's species

##Library
lib_seq=read_fasta("seqdump.txt") # list of database sequences
lib_sp=read_fasta_specie("seqdump.txt") # List of database sequence's species

#Blast

print("\n1. Selection of the 10 most similar sequences to the reference sequence using blast-like algorithm.")
print("\n This may take a while.")

def top_n_blast(blast_object,ref_seq, sp, n):
    ''' Finds the top_n most similar sequences from the database to the reference sequence. It was implemented so that n can be defiened by the user.'''
    similar_seqs=[]
    similar_seqs_sp=[]
    while len(similar_seqs) != n:
        bestAlign=blast_object.bestAlignment(ref_seq[0])
        if sp[bestAlign[4]] not in similar_seqs_sp:
            similar_seqs.append(blast_object.db[bestAlign[4]])
            similar_seqs_sp.append(sp[bestAlign[4]])
            blast_object.removeSequenceDB(bestAlign[4]) 
            del sp[bestAlign[4]]
        else:
            blast_object.removeSequenceDB(bestAlign[4]) 
            del sp[bestAlign[4]]
    return similar_seqs, similar_seqs_sp

def n_similar_seqs(ref_seq, ref_sp, lib_seq, lib_sp,n):
    '''Returns a list with the n most similar sequences to a reference sequences, that included'''
    #elimination from library sequences of the same specie as the reference sequence and adding the last to the library
    new_lib=[]
    new_lib_sp=[]
    for i in range(len(lib_seq)):
        for j in range(len(ref_seq)):
            if lib_sp[i]!=ref_sp[j]:
                new_lib.append(lib_seq[i])
                new_lib_sp.append(lib_sp[i]) 
    new_lib+=ref_seq
    new_lib_sp+=ref_sp          
    blast=MyBlast()
    #constructing the db for the MyBlast object
    for i in new_lib :
        blast.addSequenceDB(i) 
    most_similar_seq, most_similar_seq_sp=top_n_blast(blast,ref_seq, new_lib_sp,n)
    return most_similar_seq, most_similar_seq_sp

# a total of 11 sequences are selected, one corresponds to the reference sequence and the other refer to the 10 most similar
most_similar_seq, most_similar_seq_sp=n_similar_seqs(ref_seq, ref_sp, lib_seq, lib_sp,11) 

print("\n Sequences referreing the the following species were selected:\n")
for i in most_similar_seq_sp:
    print(" ",i)
input("\n Press Enter to continue. \n____________________________________________________________________________________________________________________________________________")


#Multiple Alignment

print("\n2. Multiple alignment of the reference sequence with the 10 sequences found.")
print("\n This may take a while.\n")

def m_align(seqs, seqsID):
    '''Performs and returns the multiple alignment of the selected sequences'''
    for i in range(len(seqs)):
        seqs[i]=MySeq(seqs[i], 'PROTEIN')
    for i in range(len(seqsID)): #Substitute spaces for "_" on the species names
        seqsID[i]=re.sub(" ", "_", seqsID[i])
    sm = SubstMatrix()
    sm.read_submat_file("blosum62.mat", "\t")
    aseq = PairwiseAlignment(sm, -8)
    ma = MultipleAlignment(seqs, aseq)
    al=ma.align_consensus()
    return al

m_al=m_align(most_similar_seq, most_similar_seq_sp)

def print_m_align(m_align, max_len, seqsID):
    fixed_width=max(len(i) for i in seqsID)
    x=0
    while x<max(len(i) for i in m_align):
        for i in range(len(seqsID)):
            res=""
            for j in range(x, x+max_len):
                if j > (max(len(i) for i in m_align)-1):
                    break
                else:
                    res+=m_align[i][j]
            print(("%"+str(fixed_width)+"s"+" | "+ res) % (seqsID[i]))
        print('\n')
        x+=max_len

print_m_align(m_al, 80, most_similar_seq_sp)

export_fasta("m_alignment", m_al, most_similar_seq_sp)

input("\n Press Enter to continue. \n____________________________________________________________________________________________________________________________________________")

    
##Phylogenetic Tree

print("\n3. Construction of the phylogenetic tree. \n\n Results  will be exported in newick format to the tree_viz.txt file.\n")
print("\n This may take a while. \n")
    
def phylo_tree(seqs,seqsID=None):
    '''Constructs, print and return the phylogenetic tree considering the sequences preveously selected'''
    sm = SubstMatrix()    
    sm.read_submat_file("blosum62.mat", "\t")
    alseq = PairwiseAlignment(sm, -8)    
    up  = UPGMA(seqs, alseq)    
    tree = up.run()    
    if seqsID==None:
        tree.print_tree()
    else:
        tree.print_tree(seqsID)
    return tree


def tree2Newick(tree, seqsID=None):
    '''Performs the convertion of the phylogenetic tree the the newick format'''
    path=")"
    if seqsID==None:
        if tree.value >= 0:
            path="%s" % (tree.value) + path
        elif tree.left.value >=0 and tree.right.value>=0:
            path="(%s , %s " % (tree.left.value, tree.right.value)+path
        elif tree.left.value>=0:
            path=",%s " % (tree.left.value) + path
            path="("+tree2Newick(tree.right) + path
        elif tree.right.value>=0:
            path=",%s" % (tree.right.value) + path
            path="("+tree2Newick(tree.left) +  path
        else:
            if (tree.right != None):
                path=tree2Newick(tree.right) + path
            if (tree.left != None):
                path="("+ tree2Newick(tree.left) + ","+ path
    else:
        if tree.value >= 0:
                path="%s" % (seqsID[tree.value]) + path
        elif tree.left.value >=0 and tree.right.value>=0:
            path="(%s , %s " % (seqsID[tree.left.value], seqsID[tree.right.value])+path
        elif tree.left.value>=0:
            path=",%s " % (seqsID[tree.left.value]) + path
            path="("+tree2Newick(tree.right, seqsID) + path
        elif tree.right.value>=0:
            path=",%s" % (seqsID[tree.right.value]) + path
            path="("+tree2Newick(tree.left,seqsID) +  path
        else:
            if (tree.right != None):
                path=tree2Newick(tree.right,seqsID) + path
            if (tree.left != None):
                path="("+ tree2Newick(tree.left,seqsID) + ","+ path
    return path 
        

def export_tree(tree):
    '''Tree export to the file tree_viz.txt'''
    file=open("tree_viz.txt","w")
    _=file.write(tree)
    file.close
    
tree=phylo_tree(most_similar_seq, most_similar_seq_sp)

newick_t=tree2Newick(tree,most_similar_seq_sp)
newick_t+=";"

export_tree(newick_t)

input("\n Press Enter to continue. \n____________________________________________________________________________________________________________________________________________")


##Graph
print("\n4. Construction of a graph.\n\n Results will be exported in dot format to the graph_viz.txt file.")
print("\n This may take a while. \n")

def create_graph (seqs, dist_cut):
    '''Construction of the graph considering the sequences preveously selected'''
    sm = SubstMatrix()    
    sm.read_submat_file("blosum62.mat", "\t")
    alseq = PairwiseAlignment(sm, -8)    
    up  = UPGMA(seqs, alseq)    
    graph=MyGraph()
    for i in range(up.matdist.num_rows()):
        for j in range(0,i):
            if up.matdist[i][j] < dist_cut:
                graph.add_edge(i,j)
                graph.add_edge(j,i)
    return graph

def export_graph(graph,seqsID):
    '''Export graph in dot format to the file graph_viz.txt'''
    used=[]
    file=open("graph_viz.txt","w")
    _=file.write("graph G { \n")
    for i in graph.graph.keys():
        used.append(i)
        for j in graph.graph[i]:
            if j not in used:
             _=file.write("\t" +str(seqsID[i])+" -- " + str(seqsID[j])+ ";\n")
    _=file.write("}")
    file.close()


def print_graph(self, seqsID=None):
    ''' Prints the content of the graph as adjacency list '''
    if seqsID==None:
        for v in self.graph.keys():
            print (v, " -> ", self.graph[v])
    else:
        res=self.graph
        for k in res.keys():
            for v in res[k]:
                v=seqsID[v]
            k=seqsID[k]
        for v in res.keys():
            print (res[v], " -> ", res[v])

graph=create_graph(most_similar_seq, 20)

graph.print_graph()
    
export_graph(graph, most_similar_seq_sp)

