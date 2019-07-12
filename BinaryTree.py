class BinaryTree:

    def __init__(self, val, dist = 0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

        
    def get_cluster(self):
        leaves=[]
        if self.value>=0:
            leaves.append(self.value)
        else:
            if self.left!=None :
                leaves.extend(self.left.get_cluster())
            if self.right!=None:
                leaves.extend(self.right.get_cluster())
        return leaves

    def print_tree(self,seqsID=None):
        if seqsID==None:
            self.print_tree_rec(0, "Root")
        else:
            self.print_tree_rec(0, "Root",seqsID)

    
    def print_tree_rec (self, level, side,seqsID=None):
        tabs = ""
        for i in range(level): tabs += "\t"
        if self.value >= 0:
            if seqsID==None:
                print(tabs, side, " - value:", self.value)
            else: 
                print(tabs, side, " - value:", seqsID[self.value])
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if (self.left != None): 
                self.left.print_tree_rec(level+1, "Left", seqsID)
            if (self.right != None): 
                self.right.print_tree_rec(level+1, "Right", seqsID)
               
    
    def size(self):
        ''' size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves'''
        numleaves = 0
        numnodes = 0
        if self.value >= 0:
            numleaves = 1
        else: 
            if (self.left != None):
                resl = self.left.size()
            else: resl = (0,0)
            if (self.right != None):  
                resr = self.right.size() 
            else: resr = (0,0)
            numnodes += (resl[0] + resr[0] + 1)
            numleaves += (resl[1] + resr[1])
        return numnodes, numleaves

    def exists_leaf(self, leafnum):
        if self.value>=0:
            if self.value==leafnum:
                return True
        else:    
            if self.left!=None:
                resl= self.left.exists_leaf(leafnum)
                if resl == True : return True
            if self.right!=None:
                resr=self.right.exists_leaf(leafnum)
                if resr==True : return True
        return False

    def common_ancestor(self, leaf1, leaf2):
        ''' Return simplest tree that contains leaf1, leaf2'''
        if self.value >= 0: return None
        if self.left.exists_leaf(leaf1):
            if self.left.exists_leaf(leaf2):
                return self.left.common_ancestor(leaf1,leaf2)
            if self.right.exists_leaf(leaf2):
                return self
            return None
        if self.right.exists_leaf(leaf1):
            if self.right.exists_leaf(leaf2):
                return self.right.common_ancestor(leaf1,leaf2)
            if self.left.exists_leaf(leaf2):
                return self
        return None

    def distance_leaves(self, leafnum1, leafnum2):
       ''' distance between leafnum1 and leafnum2 using the common ancestor function.
       '''
       ca=self.common_ancestor(leafnum1, leafnum2)
       return 2* ca.distance

def test():              
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    e = BinaryTree(-1, 2.0, b, c)
    f = BinaryTree(-1, 1.5, d, a)
    g = BinaryTree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())

    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1,4).print_tree()
    print(g.distance_leaves(1,4))
    print(g.distance_leaves(1,2))

if __name__ == '__main__':
    test()
