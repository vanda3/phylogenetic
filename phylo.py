from aux import blosum62
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster
import math
import sys
import pygraphviz as pgv

dist=[]
accum=0
labels=[]

class prot:
    def __init__(self, name, sequence):
        self.name=name
        self.sequence=sequence

class Node:
    def __init__(self, idx, left, right):
        self.nr=-1
        self.left_nr=-1
        self.right_nr=-1
        self.idx=idx
        self.n_leaves=0
        self.left=left
        self.right=right
        self.label=""
        self.sequence=[]
        self.isLeaf=False
        self.distance=0
        self.height=0
    def leaf(self, label, sequence):
        self.label=label
        self.sequence=sequence
        self.isLeaf=True
        self.n_leaves=1
    def node(self, distance):
        if(self.right!=None):
            self.n_leaves+=self.right.n_leaves
        if(self.left!=None):
            self.n_leaves+=self.left.n_leaves
        self.distance=truncate(distance,2)
        if(self.isLeaf==False):
            self.sequence.append(self.left.sequence)
            self.sequence.append(self.right.sequence)

def fasta(name,debug):
    fam={}
    name=str(name)+".txt"
    file = open(name,mode='r+')
    sequence=""
    header=""
    i=0
    for line in file:
        line=line.strip()
        if (line.startswith(">")):
            if(sequence!=""):
                info=prot(header,sequence)
                if(debug):
                    print("prot[",i,"-",header,"]: ",sequence)
                fam[i]=info
                i+=1
            header=line.strip(">")
            sequence=""
        else:
            sequence+=line
    if(debug):
        print("prot[",i,"-",header,"]: ",sequence)
    info=prot(header,sequence)
    fam[i]=info
    return fam

def debugger(node,level):
    if node == None:
        return
    print("level: ", level, "  idx: ", node.idx, " distance: ", node.distance)
    print("Sequence: ", node.sequence)
    debugger(node.left, level+1)
    debugger(node.right, level+1)

def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n


def printTree(g, node):
    if(node.left!=None):
        if(node.left.isLeaf):
            g.add_node(node.left.idx,label=node.left.label)
            g.add_edge(node.idx,node.left.idx,label=str(node.left.distance),weight=(node.left.distance/5+1), penwidth=2)
        else:
            g.add_node(node.left.idx,label="", height="0.1", width="0.1")
            g.add_edge(node.idx,node.left.idx,label=str(node.left.distance),weight=(node.left.distance/5+1), penwidth=2)
        printTree(g, node.left)
    if(node.right!=None):
        if(node.right.isLeaf):
            g.add_node(node.right.idx,label=node.right.label)
            g.add_edge(node.idx,node.right.idx,label=str(node.right.distance),weight=(node.left.distance/5+1), penwidth=2)
        else:
            g.add_node(node.right.idx,label="", height="0.1", width="0.1")
            g.add_edge(node.idx,node.right.idx,label=str(node.right.distance),weight=(node.left.distance/5+1), penwidth=2)
        printTree(g, node.right)


def distances(A, B):
    count = 0
    for i in range(0, len(A)):
        a = A[i]
        b = B[i]
        if(a != '-' and b != '-'):
            if((a, b) in blosum62 and blosum62[(a, b)] < 0):
                count += 1
            if((b,a) in blosum62 and blosum62[(b, a)] < 0):
                count += 1
        if(a=='-' or b=='-'):
            count+=1
    return count

def Qmatrix(distMatrix,leaves):
    Q={}
    n=len(leaves)
    for a in leaves:
        for b in leaves:
            if (a!=b and (a,b) in distMatrix):
                count1=0
                count2=0
                for c in leaves:
                    if ((a,c) in distMatrix):
                        count1+=distMatrix[(a,c)]
                    if ((c,b) in distMatrix):
                        count2+=distMatrix[(c,b)]
                Q[(a,b)]=(n-2)*distMatrix[(a,b)]-count1-count2
    return Q
 
def joining(distMatrix, sequences, debug):
    leaves={}
    nodes=[]
    if(debug):
        print()
    for i in range(0,len(sequences)):
        node=Node(str(i), None, None)
        node.leaf(sequences[i].name,sequences[i].sequence,)
        leaves[str(i)]=node
        nodes.append(str(i))
    m=1
    while len(leaves) > 2:
        n=len(leaves)
        Q=Qmatrix(distMatrix,leaves)
        #find the two clusters the have the minimun distance
        if(debug):
            print("/////////// CYCLE: ",m)
        minimum = sys.maxsize # high value
        min_id1=-1
        min_id2=-1
        for a in leaves:
            for b in leaves:
                if (a,b) in Q:
                    if(a!=b and minimum > Q[(a,b)]):
                        minimum = Q[(a,b)]
                        min_id1 = a
                        min_id2 = b
        if(debug):
            print("MIN DISTANCE:", minimum,"   ",min_id1, "+",min_id2)
        new_id=str(min_id1)+","+str(min_id2)
        nodes.append(new_id)
        #update distance values
        for key in leaves:
            if(key!=min_id1 and key!=min_id2):
                # ((a,b),c)=((a,c)|C_a|+(b,c)|C_b|)/(|C_a| + |C_b|)
                distMatrix[(new_id, key)] = (1/2)*(distMatrix[(min_id1, key)] + distMatrix[(min_id2, key)] - distMatrix[(min_id1, min_id2)])
                distMatrix[(key, new_id)] = distMatrix[(new_id, key)]
                if(debug):
                    print("DIST (",new_id,",",key,") = ", distMatrix[(key,new_id)])
        sum1=0
        sum2=0
        for key in leaves:
            if(key!=min_id1):
                sum1+=distMatrix[(key,min_id1)]
            if(key!=min_id2):
                sum2+=distMatrix[(key,min_id2)]
        distMatrix[(min_id1,new_id)]=(1/2)*distMatrix[(min_id1,min_id2)]+(1/(2*(n-2)))*(sum1-sum2)
        distMatrix[(new_id,min_id1)]=distMatrix[(min_id1,new_id)]
        if(debug):
            print("DIST (",new_id,",",min_id1,") = ", distMatrix[(new_id,min_id1)])
        distMatrix[(min_id2,new_id)]=(1/2)*distMatrix[(min_id1,min_id2)]+(1/(2*(n-2)))*(sum2-sum1)    
        distMatrix[(new_id,min_id2)]=distMatrix[(min_id2,new_id)]
        if(debug):
            print("DIST (",new_id,",",min_id2,") = ", distMatrix[(new_id,min_id2)])
        leaves[min_id1].node(distMatrix[(min_id1,new_id)])
        leaves[min_id2].node(distMatrix[(min_id2,new_id)])
        node=Node(new_id,leaves[min_id1],leaves[min_id2])
        leaves[new_id]=node
        
        del leaves[min_id1]
        del leaves[min_id2]

        m+=1
        if(debug):
            print()
    #find last two clusters
    keys=list(leaves.keys())
    new_id=str(keys[0])+str(keys[1])
    node1=leaves[keys[0]]
    node2=leaves[keys[1]]
    node=Node(new_id,node1,node2)
    node.node(distMatrix[(keys[0],keys[1])]+1)
    leaves[new_id]=node
    #debugger(node,0)
    return node

def upgma(distMatrix,sequences,debug):
    global dist
    labels=[]
    accum=0
    nr=len(sequences)
    leaves={}
    if(debug):
        print()
    for i in range(0,len(sequences)):
        node=Node(str(i), None, None)
        node.leaf(sequences[i].name,sequences[i].sequence)
        node.nr=i
        labels.append(node.label)
        leaves[str(i)]=node

    m=1
    while len(leaves) > 2:
        #find the two clusters the have the minimun distance
        if(debug):
            print("/////////// CYCLE: ",m)
        minimum = sys.maxsize # high value
        min_id1=-1
        min_id2=-1
        for a in leaves:
            for b in leaves:
                if (a,b) in distMatrix:
                    if(a!=b and minimum > distMatrix[(a,b)]):
                        minimum = distMatrix[(a,b)]
                        min_id1 = a
                        min_id2 = b
        if(debug):
            print("MIN DISTANCE:", minimum,"   ",min_id1, "+",min_id2)
        new_id=str(min_id1)+","+str(min_id2)
        leaves[min_id1].node(minimum)
        leaves[min_id2].node(minimum)
        #update distance values
        for key in leaves:
            if(key!=min_id1 and key!=min_id2):
                # ((a,b),c)=((a,c)|C_a|+(b,c)|C_b|)/(|C_a| + |C_b|)
                distMatrix[(new_id, key)] = (distMatrix[(min_id1, key)]* len(leaves[min_id1].sequence)   +distMatrix[(min_id2, key)] * len(leaves[min_id2].sequence))/(len(leaves[min_id1].sequence) + len(leaves[min_id2].sequence))
                distMatrix[(key, new_id)] = distMatrix[(new_id, key)]
                if(debug):
                    print("DIST (",new_id,",",key,") = ", distMatrix[(key,new_id)])
        
        node=Node(new_id,leaves[min_id1],leaves[min_id2])
        node.nr=nr
        node.node(minimum)
        accum+=minimum
        dist.append([leaves[min_id1].nr,leaves[min_id2].nr,float(accum),node.n_leaves])
        nr+=1
        leaves[new_id]=node
        del leaves[min_id1]
        del leaves[min_id2]

        m+=1
        if(debug):
            print()

    #find last two clusters
    keys=list(leaves.keys())
    new_id=keys[0]+","+keys[1]
    minimum = (distMatrix[(keys[0], keys[1])]* len(leaves[keys[0]].sequence)   +distMatrix[(keys[0], keys[1])] * len(leaves[keys[1]].sequence))/(len(leaves[keys[0]].sequence) + len(leaves[keys[1]].sequence))
    leaves[keys[0]].node(minimum)
    leaves[keys[1]].node(minimum)
    node=Node(new_id,leaves[keys[0]],leaves[keys[1]])
    node.nr=nr
    node.node(minimum)
    accum+=minimum
    dist.append([leaves[keys[0]].nr,leaves[keys[1]].nr,float(accum),node.n_leaves])
    leaves[new_id]=node

    #debugger(node,0)

    #root
    return node,labels

if __name__ == '__main__':
    isProtein=True
    distMatrix = {}
    debug=False
    print("Which family will you evaluate?")
    name=str(input())
    fam=fasta(name,debug)

    if(debug):
        print("/////////// CYCLE: ",0)

    for i in range(0,len(fam)):
        for j in range(i+1,len(fam)):
            if(i!=j):
                A=fam[i].sequence
                B=fam[j].sequence
                distMatrix[(str(i),str(j))] = distances(A, B)
                # simetrico
                distMatrix[(str(j),str(i))] = distMatrix[(str(i),str(j))]
                if(debug):
                    print("DIST (",i,",",j,") = ", distMatrix[(str(i),str(j))])
    
    node1,labels=upgma(distMatrix,fam,debug)
    node2=joining(distMatrix,fam,debug)

    G2=pgv.AGraph()
    G2.add_node(node2.idx,label="", color="red", height="0.2", width="0.2")
    printTree(G2, node2)
    G2.outputorder="edgesfirst"  
    G2.layout(prog="neato")
    G2.draw('nj.png')
    ytdist = np.array(dist)
    scipy.cluster.hierarchy.dendrogram(ytdist, leaf_font_size=4, orientation="left", labels=labels)
    plt.savefig('upgma.png')