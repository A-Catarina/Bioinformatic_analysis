## Graph represented as adjacency list using a dictionary
## keys are vertices
## values of the dictionary represent the list of adjacent vertices of the key node

class MyGraph:
    def __init__(self, g = {}):
        ''' Constructor - takes dictionary to fill the graph as input; default is empty dictionary '''
        self.graph = g    
    def print_graph(self):
        ''' Prints the content of the graph as adjacency list '''
        for v in self.graph.keys():
            print (v, " -> ", self.graph[v])
        
    ## get basic info
    def get_nodes(self):
        ''' Returns list of nodes in the graph '''
        return list(self.graph.keys())  
    def get_edges(self): 
        ''' Returns edges in the graph as a list of tuples (origin, destination) '''
        edges=[]
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v,d))
        return edges    
    def size(self):
        ''' Returns size of the graph : number of nodes, number of edges '''
        return len(self.get_nodes(), len(self.get_edges()))
      
    ## add nodes and edges    
    def add_vertex(self, v):
        ''' Add a vertex to the graph; tests if vertex exists not adding if it does '''
        if v not in self.graph:
            self.graph[v]=[]        
    def add_edge(self, o, d):
        ''' Add edge to the graph; if vertices do not exist, they are added to the graph ''' 
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    ## successors, predecessors, adjacent nodes
        
    def get_successors(self, v):
        return list(self.graph[v])     # needed to avoid list being overwritten of result of the function is used           
             
    def get_predecessors(self, v):
        res = []
        for k in self.graph.keys(): 
            if v in self.graph[k]: 
                res.append(k)
        return res
    
    def get_adjacents(self, v):
        suc = self.get_successors(v)
        pred = self.get_predecessors(v)
        res = pred
        for p in suc: 
            if p not in res: res.append(p)
        return res
        
    ## degrees    
    
    def out_degree(self, v):
        return len(self.get_successors(v))
    
    def in_degree(self, v):
        return len(self.get_predecessors(v))
        
    def degree(self, v):
        len(self.get_adjacents(v))
        
    def all_degrees(self, deg_type = "inout"):
        ''' Computes the degree (of a given type) for all nodes.
        deg_type can be "in", "out", or "inout" '''
        degs = {}
        for v in self.graph.keys():
            if deg_type == "out" or deg_type == "inout":
                degs[v] = len(self.graph[v])
            else: degs[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs[d] + 1
        return degs
    
    def highest_degrees(self, all_deg= None, deg_type = "inout", top= 10):
        if all_deg is None: 
            all_deg = self.all_degrees(deg_type)
        ord_deg = sorted(list(all_deg.items()), key=lambda x : x[1], reverse = True)
        return list(map(lambda x:x[0], ord_deg[:top]))
        
    
    ## topological metrics over degrees

    def mean_degree(self, deg_type = "inout"):
        ''' average degree of all nodes: sum of all degrees divided by number of nodes'''
        degs=self.all_degrees(deg_type)
        return sum(degs.values())/float(len(degs))


        
    def prob_degree(self, deg_type = "inout"):
        # count the number of occurrences of each degree in the network and derive its frequencies
        deg=self.all_degrees(deg_type)
        prob={}
        for k in deg.keys():
            if deg[k] in prob.keys():
                prob[deg[k]]+=1
            else:
                prob[deg[k]]=1
        for i in prob:
            prob[i]/=float(len(deg))
        return prob

    
    
    ## BFS and DFS searches    
    
    def reachable_bfs(self, v):
        l = [v]   # list of nodes to be handled
        res = []  # list of nodes to return the result
        while len(l) > 0:
            node = l.pop(0)  # implements a queue: LILO
            if node != v: res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res
        
    def reachable_dfs(self, v):
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0) # implements a stack: 
            if node != v: res.append(node)
            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(s, elem)
                    s += 1
        return res    
    
    def distance(self, s, d):
        if s == d: return 0
        l = [(s,0)]
        visited = [s]
        while len(l) > 0:
            node, dist = l.pop(0)
            for elem in self.graph[node]:
                if elem == d: return dist + 1
                elif elem not in visited: 
                    l.append((elem,dist+1))
                    visited.append(elem)
        return None
        
    def shortest_path(self, s, d):
        if s == d: return 0
        l = [(s,[])]
        visited = [s]
        while len(l) > 0:
            node, preds = l.pop(0)
            for elem in self.graph[node]:
                if elem == d: return preds+[node,elem]
                elif elem not in visited: 
                    l.append((elem,preds+[node]))
                    visited.append(elem)
        return None

                
    ## clustering
        
    def clustering_coef(self, v):
        adjs = self.get_adjacents(v)
        if len(adjs) <=1: return 0.0
        # calculate the number of links of the adjacent nodes 
        ligs = 0
        # compare pairwisely if nodes in this list are connected between them
        for i in adjs:
            for j in adjs:
                if i != j:
                    # check if i and j are connected to each other; if yes increment counter of links
                    if j in self-graph[i] or i in self.graph[j]:
                        links+=1
        return float(ligs)/(len(adjs)*(len(adjs)-1))
        
    def all_clustering_coefs(self):
        # go through all the nodes and calculate its cc
        # put those in a dictionary and return
        ccs={}
        for k in self.graph.keys():
            ccs[k]=self.clustering_coef(k)
        return ccs

        
    def mean_clustering_coef(self):
        # get all the clustering coefficients
        # and return the mean of all ccs
        degs=self.all_clustering_coefs()
        return sum(degs.values())/float(len(degs))

    
if __name__ == "__main__":
    gr = MyGraph()
    gr.add_vertex(1)
    gr.add_vertex(2)
    gr.add_vertex(3)
    gr.add_vertex(4)
    gr.add_edge(1,2)
    gr.add_edge(2,3)
    gr.add_edge(3,2)
    gr.add_edge(3,4)
    gr.add_edge(4,2)
    gr.print_graph()
    #print(gr.size())
    
    # print (gr.get_successors(2))
    # print (gr.get_predecessors(2))
    # print (gr.get_adjacents(2))
    # 
    # print (gr.in_degree(2))
    # print (gr.out_degree(2))
    # print (gr.degree(2))
    # 
    #print(gr.all_degrees("inout"))
    #print(gr.all_degrees("in"))
    #print(gr.all_degrees("out"))
    # 
    # gr2 = MyGraph({1:[2,3,4], 2:[5,6],3:[6,8],4:[8],5:[7],6:[],7:[],8:[]})
    # print(gr2.reachable_bfs(1))
    # print(gr2.reachable_dfs(1))
    # 
    # print(gr2.distance(1,7))
    # print(gr2.shortest_path(1,7))
    # print(gr2.distance(1,8))
    # print(gr2.shortest_path(1,8))
    # print(gr2.distance(6,1))
    # print(gr2.shortest_path(6,1))
    # 
    # print(gr2.reachable_with_dist(1))

    # print(gr.mean_degree())
    # print(gr.prob_degree())
    # print(gr.mean_distances())
    #print (gr.clustering_coef(1))
    #print (gr.clustering_coef(2))
