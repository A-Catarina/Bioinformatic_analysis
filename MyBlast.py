# -*- coding: utf-8 -*-

class MyBlast:
    '''
    Classe para matrizes de pontos
    '''

    def __init__(self, filename = None, w = 3):
        '''
        Construtor
        '''
        if filename is not None:
            self.readDatabase(filename)
        else:
            self.db = []
        self.w = w
        self.map = None

    def readDatabase(self, filename):
        """From file with sequences line by line read the sequences to a list"""
        seqs=[]
        seq=""
        lines=open(filename).readlines()
        limit=len(lines)-1
        for i in range(limit):
            while i<=limit:
                seq+=lines[i].strip()
            seqs.append(seq)
            i+=1    
            seq=""
        return seqs

    
    def addSequenceDB(self, seq):
        """Add an extra sequence to DB"""
        self.db.append(seq)

    def removeSequenceDB(self, index):
        del self.db[index]
        
    def buildMap (self, query):
        dic={}
        for i in range(len(query)-self.w+1):
            if query[i:i+self.w] in dic.keys():
                dic[query[i:i+self.w]].append(i)
            else:
                dic[query[i:i+self.w]]=[i]
        return dic
            
    def getHits (self, seq, query):
        m=self.buildMap(query)
        res=[]
        for i in range(len(seq)-self.w+1):
            subseq=seq[i:i+self.w]
            if subseq in m:
                l=m[subseq]
                for ind in l:
                    res.append((ind,i))
        return res

    
    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0       
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]: 
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk
    
        k = 0
        matbw = 0   
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]: 
                matbw+=1
                bestk = k+1
            k+=1       
        size += bestk
        
        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        
    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best
 
 
    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]  
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res

    def get_hits(self,seq, m,w):
        res=[]
        for i in range(len(seq)-w+1):
            subseq=seq[i:i+w]
            for k in m.keys():
                msm=0
                while msm<2:
                    for j in range(w):
                        if subseq[j]!=k[j]:
                            msm+=1
                if msm<2:
                    l=m[k]
                    for ind in l:
                        res.append((ind,i))
        return res

def test2():
    mb = MyBlast()
    mb.addSequenceDB("cgtgcactgtacgtgactgatctgtact")
    mb.addSequenceDB("gcgtgaggxtgcgtagxttacgt")
    mb.addSequenceDB("agtcgatgtcgatgctgaccta")
    query2 = "cgacgacgacgacgaatgatg"
    m=mb.buildMap(query2)
    a=mb.get_hits(mb.db[0],m, mb.w)
    print(a)
    r = mb.bestAlignment(query2)
    print(r)

test2()
print('a')
