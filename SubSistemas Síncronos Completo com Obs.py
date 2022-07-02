from deslab import *
import networkx as nx
import matplotlib as plt
import numpy as np
from tabulate import tabulate
syms ('x01 x11 x21 x31 x41 x02 x12 x22 x32 x03 x13 x23 x33 1 2 3')
Seq1 = ['100','000','111','100','000','001','100','000','111']
Seq4 = ['100','110','100','000','001']
Seq2 = ['1000','0110','1001','1111','1000','0110','1011','1111','1001']
Seq3 = ['10000','10101','10001','10011','10000','10101','10001','00011','10011','10100','00001','10011','10101']
Seqfinal = [Seq1,Seq4]

Seq5 = ['100','110','001','000','100','000']
Seq6 = ['100','111','000']
Seqsimples = [Seq5,Seq6]
table =[('sigma',' ')]
def seqbreak (Seq,IOs):     
    Nsys = len(IOs)
    Seqn = []
    for i in range(Nsys):
        h = 0
        Se =[]
        for k in range(len(Seq)):
            S = list(Seq[k])
            for j in range(len(S)):
                if j+1 in IOs[i]:
                    S[j] = S[j]
                else:
                    S[j] = '-'
            S = [''.join(str(e) for e in S)]
            if len(Se) == 0:
                Se = Se + S
                Sx = tuple(['x'+str(i)+str(h)] + [''.join(str(e) for e in S)])
                Seqn = Seqn + [Sx]
            if S[0] != Se[h]:
                Se = Se + S
                Sx = tuple(['x'+str(i)+str(h+1)] + [''.join(str(e) for e in S)])
                h = h+1
                Seqn = Seqn + [Sx]
    L1 = []
    L2 = []
    L3 = []
    for i in range(len(Seqn)):
        if Seqn[i][0][1] == '0':
            L1 = L1 + [Seqn[i]]
        if Seqn[i][0][1] == '1':
            L2 = L2 + [Seqn[i]]
        if Seqn[i][0][1] == '2':
            L3 = L3 + [Seqn[i]]
##    print('L1='+str(L1),'L2='+str(L2),'L3='+str(L3),'L4='+str(L4),'L5='+str(L5))
    return(L1,L2,L3)

def ksequence (Seq,IOs,k):
    Tfinal = []
    Seqfinal = []
    for i in range(len(Seq)):
        L1,L2,L3 = seqbreak(Seq[i],IOs)
        L = [L1,L2,L3]
        Seqk = []
        T = []
        for i in range(len(L)):
            for p in range(len(L[i])):      
                if p == 0:              #Estado Inicial
                    x0 = L[i][0][1]
                    X0 = L[i][0][1]
                    for l in range(k-1):
                        x0 = x0 + '|' + X0
                    Init = x0
                    S = tuple(['x'+str(i)+'0'] + [x0])
                    Seqk = Seqk + [S]
                xn = x0.split('|')
                if p > 0:               #Estados 1 até inf
                    for j in range(k-1):
                        xn[j]=xn[j+1]
                    xn[k-1] = L[i][p][1]
                    xnp = xn[0]
                    for u in range(len(xn)-1):
                        xp = xn[u+1]
                        xnp = xnp +  '|' + xp
                    T = T + [(x0,'sigma',xnp)]
                    S = tuple(['x'+str(i)+str(p)] + [xnp])
                    Seqk = Seqk + [S]
                    x0 = xnp
            if xnp != Init:
                T = T + [(xnp,'sigma',Init)]
        Tfinal = Tfinal + T
        Seqfinal = Seqfinal + Seqk
    Out1,Out2,Out3 = [],[],[]
    for i in range(len(Seqfinal)):
        if Seqfinal[i][0][1] == '0':
            Out1 = Out1 + [Seqfinal[i][1]]
        if Seqfinal[i][0][1] == '1':
            Out2 = Out2 + [Seqfinal[i][1]]
        if Seqfinal[i][0][1] == '2':
            Out3 = Out3 + [Seqfinal[i][1]]
    Out1 = list(dict.fromkeys(Out1))
    Out2 = list(dict.fromkeys(Out2))
    Out3 = list(dict.fromkeys(Out3))
    T1,T2,T3 = [],[],[]
    for i in range(len(Tfinal)):
        if Tfinal[i][0] in Out1:
            T1 = T1 + [Tfinal[i]]
        if Tfinal[i][0] in Out2:
            T2 = T2 + [Tfinal[i]]
        if Tfinal[i][0] in Out3:
            T3 = T3 + [Tfinal[i]]
    T1 = list(dict.fromkeys(T1))
    T2 = list(dict.fromkeys(T2))
    T3 = list(dict.fromkeys(T3))
    if len(Out1) > 0:
        x0 = Out1[0]
        G1= fsa(Out1,['sigma'],T1,x0,[],table, name = 'Subsistema 1')
        G1.setgraphic(style ='observer',direction = 'UD')
        #draw (G1)
    if len(Out2) > 0:
        x0 = Out2[0]
        G2= fsa(Out2,['sigma'],T2,x0,[],table, name = 'Subsistema 2')
        G2.setgraphic(style ='observer',direction = 'UD')
        #draw (G2)
        T1f,T2f,E1,E2= [],[],[],[]
        for i in range(len(T1)):
            T = (T1[i][0],T1[i][2].split('|')[k-1],T1[i][2])
            T1f = T1f + [T]
            x0 = Out1[0]
            E1 = E1 + [T1[i][2].split('|')[k-1]]
        Go1 = fsa(Out1,E1,T1f,x0,[],table)
        Go1.setgraphic(style ='observer',direction = 'UD')
        Gobs1 = observer(Go1,E1)
        #draw(Gobs1)
        for i in range(len(T2)):
            T = (T2[i][0],T2[i][2].split('|')[k-1],T2[i][2])
            T2f = T2f + [T]
            x0 = Out2[0]
            E2 = E2 + [T2[i][2].split('|')[k-1]]
        Go2 = fsa(Out2,E2,T2f,x0,[],table)
        Go2.setgraphic(style ='observer',direction = 'UD')
        Gobs2 = observer(Go2,E2)
        #draw(Gobs2)
    if len(Out3) > 0:
        x0 = Out3[0]
        G3= fsa(Out3,['sigma'],T3,x0,[],table, name = 'Subsistema 3')
        G3.setgraphic(style ='observer',direction = 'UD')
        draw (G3)
    return (Out1,Out2,Out3,T1,T2,T3)

def TransformAut(Seq,IOs,k):
    Out1,Out2,Out3,T1,T2,T3 = ksequence (Seq,IOs,k)
    L1,LL1,L2,LL2,L3,LL3=[],[],[],[],[],[]
    if len(Out1)>0:
        for i in range(len(Out1)):
            O = Out1[i].split('|')
            LL1 = LL1 + [tuple(['x'+str(i)] + [Out1[i]])]
            L1 = L1 + [tuple(['x'+str(i)] + [O[k-1]])]
        for j in range(len(LL1)):
            m = 0
            while m < len(T1):
                if T1[m][0] == LL1[j][1]:
                    Tl = list(T1[m])
                    Tl[0] = LL1[j][0]+'|'+LL1[j][1].split('|')[k-1]
                    T1.pop(m)
                    T1 = T1 + [tuple(Tl)]
                    m=0
                if T1[m][2] == LL1[j][1]:
                    Tl = list(T1[m])
                    Tl[2] = LL1[j][0]+'|'+LL1[j][1].split('|')[k-1]
                    T1.pop(m)
                    T1 = T1 + [tuple(Tl)]
                    m=0
                if T1[m][0] != LL1[j][1] and T1[m][2] != LL1[j][1]:
                    m = m+1
    if len(Out2)>0:
        for i in range(len(Out2)):
            O = Out2[i].split('|')
            LL2 = LL2 + [tuple(['y'+str(i)] + [Out2[i]])]
            L2 = L2 + [tuple(['y'+str(i)] + [O[k-1]])]
        for j in range(len(LL2)):
            m = 0
            while m < len(T2):
                if T2[m][0] == LL2[j][1]:
                    Tl = list(T2[m])
                    Tl[0] = LL2[j][0]+'|'+LL2[j][1].split('|')[k-1]
                    T2.pop(m)
                    T2 = T2 + [tuple(Tl)]
                    m=0
                if T2[m][2] == LL2[j][1]:
                    Tl = list(T2[m])
                    Tl[2] = LL2[j][0]+'|'+LL2[j][1].split('|')[k-1]
                    T2.pop(m)
                    T2 = T2 + [tuple(Tl)]
                    m=0
                if T2[m][0] != LL2[j][1] and T2[m][2] != LL2[j][1]:
                    m = m+1
    if len(Out3)>0:
        for i in range(len(Out3)):
            O = Out3[i].split('|')
            LL3 = LL3 + [tuple(['z'+str(i)] + [Out3[i]])]
            L3 = L3 + [tuple(['z'+str(i)] + [O[k-1]])]
        for j in range(len(LL3)):
            m = 0
            while m < len(T3):
                if T3[m][0] == LL3[j][1]:
                    Tl = list(T3[m])
                    Tl[0] = LL3[j][0]+'|'+LL3[j][1].split('|')[k-1]
                    T3.pop(m)
                    T3 = T3 + [tuple(Tl)]
                    m=0
                if T3[m][2] == LL3[j][1]:
                    Tl = list(T3[m])
                    Tl[2] = LL3[j][0]+'|'+LL3[j][1].split('|')[k-1]
                    T3.pop(m)
                    T3 = T3 + [tuple(Tl)]
                    m=0
                if T3[m][0] != LL3[j][1] and T3[m][2] != LL3[j][1]:
                    m = m+1
    return(LL1,LL2,LL3,T1,T2,T3)
            

def AutomatoSintese (Seq,IOs,k):
    LL1,LL2,LL3,T1,T2,T3 = TransformAut(Seq,IOs,k)
    comum = VerifierCommom(IOs)
    T1comum = []
    T2comum = []
    T1part = []
    T2part = []
    Tn = []
    X = []
    E = []
    X0 = []
    Tfinal = []
    if len(T2)>0:
        for i in range(len(T1)):   
            if [T1[i][0].split('|')[1][c] for c in comum] != [T1[i][2].split('|')[1][c] for c in comum]:
                T1comum = T1comum + [T1[i]]
            else:
                T1part = T1part + [T1[i]]        
        for i in range(len(T2)):
            if [T2[i][0].split('|')[1][c] for c in comum] != [T2[i][2].split('|')[1][c] for c in comum]:
                T2comum = T2comum + [T2[i]]
            else:
                T2part = T2part + [T2[i]]
##        print(T1comum)
##        print(T2comum)
        for i in range(len(T1comum)):                                       #Transições comuns
            for j in range(len(T2comum)):
                Pre1 = [T1comum[i][0].split('|')[1][c] for c in comum]
                #print('Pre1='+ str(Pre1))
                Pre2 = [T2comum[j][0].split('|')[1][c] for c in comum]
                #print('Pre2='+ str(Pre2))
                Pos1 = [T1comum[i][2].split('|')[1][c] for c in comum]
                #print('Pos1='+str(Pos1))
                Pos2 = [T2comum[j][2].split('|')[1][c] for c in comum]
                #print('Pos2='+str(Pos2))
                if  Pre1 == Pre2 and Pos1 == Pos2:
                    Tp = JoinFunc(T1comum[i],T2comum[j])                    #Transições comuns sempre andam ambos
                    Tn = Tn + [Tp]
##                    if k > 1:
##                        if 'x0' not in Tp[2] and 'y0' not in Tp[2]:
##                            Tn = Tn + [Tp]
##                        if 'x0,y0' in Tp[2]:
##                            Tn = Tn + [Tp]
##                    else:
##                       Tn = Tn + [Tp] 
        for i in range(len(T1part)):                                        #Transições particular Anda T1 e mantém T2
            for j in range(len(T2)):
                Pre1 = [T1part[i][0].split('|')[1][c] for c in comum]
                Pre2 = [T2[j][0].split('|')[1][c] for c in comum]
                Pos1 = [T1part[i][2].split('|')[1][c] for c in comum]
                Pos2 = [T2[j][2].split('|')[1][c] for c in comum]
                if  Pre1 == Pre2:
                    Tp1 = JoinFunc(T1part[i],(T2[j][0],'sigma',T2[j][0]))
                    Tn = Tn + [Tp1]
                    #print('Tp1='+str(Tp1))
##                    if k > 1:
##                        if 'x0' not in Tp1[2]:
##                            Tn = Tn + [Tp1]
##                    else:
##                        Tn = Tn + [Tp1]
        for i in range(len(T1)):                                            #Transições particular Anda T2 e mantém T1
            for j in range(len(T2part)):
                Pre1 = [T1[i][0].split('|')[1][c] for c in comum]
                Pre2 = [T2part[j][0].split('|')[1][c] for c in comum]
                Pos1 = [T1[i][2].split('|')[1][c] for c in comum]
                Pos2 = [T2part[j][2].split('|')[1][c] for c in comum]
                if  Pre2 == Pre1:
                    Tp2 = JoinFunc((T1[i][0],'sigma',T1[i][0]),T2part[j])
                    Tn = Tn + [Tp2]
                    #print('Tp2='+str(Tp2))
##                    if k > 1:
##                        if 'y0' not in Tp2[2]:
##                            Tn = Tn + [Tp2]
##                    else:
##                        Tn = Tn + [Tp2]
        for i in range(len(T1part)):                                        #Transições particular Anda T1 e T2
            for j in range(len(T2part)):
                Pre1 = [T1part[i][0].split('|')[1][c] for c in comum]
                Pre2 = [T2part[j][0].split('|')[1][c] for c in comum]
                Pos1 = [T1part[i][2].split('|')[1][c] for c in comum]
                Pos2 = [T2part[j][2].split('|')[1][c] for c in comum]
                if  Pre1 == Pre2 and Pos1 == Pos2:
                    Tp3 = JoinFunc(T1part[i],T2part[j])
                    Tn = Tn + [Tp3]
                    #print('Tp3='+str(Tp3))
##                    if k > 1:
##                        if 'x0' not in Tp3[2] and 'y0' not in Tp3[2]:
##                            Tn = Tn + [Tp3]
##                        if 'x0,y0' in Tp3[2]:
##                            Tn = Tn + [Tp3]
##                    else:
##                        Tn = Tn + [Tp3]
        for i in range(len(Tn)):
            if Tn[i][0].split('|')[1] != Tn[i][2].split('|')[1]:
                Tfinal = Tfinal+[Tn[i]]
        Tfinal = list(dict.fromkeys(Tfinal))
        for i in range(len(Tfinal)):
            X = X + [Tfinal[i][0]]+ [Tfinal[i][2]]
            E = E + [Tfinal[i][1]]
            if 'x0,y0' in Tfinal[i][0]:
                X0 = X0 + [Tfinal[i][0]]
            if 'x0,y0' in Tfinal[i][2]:
                X0 = X0 + [Tfinal[i][2]]
        X = list(dict.fromkeys(X))
        E = list(dict.fromkeys(E))
        X0 = list(dict.fromkeys(X0))
    else:
        #print(T1)
        Tfinal = []
        for i in range(len(T1)):
            T = (T1[i][0],T1[i][2].split('|')[1],T1[i][2])
            Tfinal = Tfinal + [T]
        #print(Tfinal)
        for i in range(len(Tfinal)):
            X = X + [Tfinal[i][0]]+ [Tfinal[i][2]]
            E = E + [Tfinal[i][1]]
            if 'x0' in Tfinal[i][0]:
                X0 = X0 + [Tfinal[i][0]]
            if 'x0' in Tfinal[i][2]:
                X0 = X0 + [Tfinal[i][2]]
        X = list(dict.fromkeys(X))
        E = list(dict.fromkeys(E))
        X0 = list(dict.fromkeys(X0))
##    print(X)
##    print(X0)
##    print(Tfinal)
    Gn = fsa(X,E,Tfinal,X0,[],table, name = 'Sistema Total k='+str(k))
    Gn.setgraphic(style ='rectangle',direction = 'UD')
    Gn = ac(Gn)
    Tfinal = transitions(Gn)
    #draw (Gn)
    Gobs = observer(Gn,E)
    #draw (Gobs)
    return(Gn,Tfinal)

def VerifierCommom(IOs):
    commom = []
    for i in range(len(IOs)):
        for k in range(len(IOs)):
            for j in range(len(IOs[i])):
                if IOs[i][j] in IOs[k] and k!=i:
                    commom = commom + [IOs[i][j]-1]
    commom = list(dict.fromkeys(commom))
    return(commom)

               
def JoinFunc(T1,T2):
    Pre1 = T1[0].split('|')[1]
    Pre2 = T2[0].split('|')[1]
    Pos1 = T1[2].split('|')[1]
    Pos2 = T2[2].split('|')[1]
    tamanho = len(Pre1)
    Pre = ['-']*tamanho
    Pos = ['-']*tamanho
    for i in range(tamanho):
        if Pre1[i] == Pre2[i]:
            Pre[i] = Pre1[i]
        if Pre1[i] == '-' and Pre2[i] != '-':
            Pre[i] = Pre2[i]
        if Pre1[i] != '-' and Pre2[i] == '-':
            Pre[i] = Pre1[i]
        if Pre1[i] != '-' and Pre2[i] != '-' and Pre1[i] != Pre2[i]:
            Pre[i] = 'c'
    for i in range(tamanho):
        if Pos1[i] == Pos2[i]:
            Pos[i] = Pos1[i]
        if Pos1[i] == '-' and Pos2[i] != '-':
            Pos[i] = Pos2[i]
        if Pos1[i] != '-' and Pos2[i] == '-':
            Pos[i] = Pos1[i]
        if Pos1[i] != '-' and Pos2[i] != '-' and Pos1[i] != Pos2[i]:
            Pos[i] = 'c'
    Pre = ''.join(Pre)
    Pos = ''.join(Pos)
    Tn = (T1[0].split('|')[0]+','+T2[0].split('|')[0]+'|'+str(Pre),str(Pos),T1[2].split('|')[0]+','+T2[2].split('|')[0]+'|'+str(Pos))
    #print(Tn)
    return (Tn)

def IOtot (IOs):
    IOtot = []
    for i in range(len(IOs)):
        for j in range(len(IOs[i])):
            if IOs[i][j] not in IOtot:
                IOtot = IOtot + [IOs[i][j]]
    return(IOtot)
    
def Lexc (Seq,IOs,k,n):
    Gcomp,Tcomp = AutomatoSintese (Seq,IOs,k)
    IOt = IOtot (IOs)
    print(IOt)
    Gtot,Ttot = AutomatoSintese (Seq,[IOt],k)
    gcomp = Gcomp.Graph
    gtot = Gtot.Graph
    Tot1,Tot2,Totfinal = [],[],[]
    Nodescomp = list(gcomp.nodes)
    Nodestot = list(gtot.nodes)
    for i in range(len(Nodescomp)):
        if 'x0,y0' in Nodescomp[i]:
            linha1 = i
    for i in range(len(Nodestot)):
        if 'x0' in Nodestot[i]:
            linha2 = i
    for c in range(1,n+1):
        Acomp = nx.to_numpy_matrix(gcomp)
        Acomp = Acomp**c
        Seqs1 = int(Acomp[linha1].sum())
        Tot1 = Tot1 + ['n = ' + str(c) + ':' + str(Seqs1)]
        Atot = nx.to_numpy_matrix(gtot)
        Atot = Atot**c
        Seqs2 = int(Atot[linha2].sum())
        Tot2 = Tot2 + ['n =' + str(c) + ':' + str(Seqs2)]
        #print(Seqs1-Seqs2)
        Totfinal = Totfinal + ['n = ' + str(c) + ' : ' + str(Seqs1 - Seqs2)]
    #print('Modelo Monolítico: Número de Estados =' + str(len(Gtot.X))+'\n                   Número de Transições=' + str(len(Ttot)))
    #print('Modelo Composto:   Número de Estados =' + str(len(Gcomp.X))+'\n                   Número de Transições=' + str(len(Tcomp)))
    #print('Modelo Monolítico \n' +str (tabulate(Tot2)))
    #print('Modelo Composto \n' + str(tabulate(Tot1)))
    print('Linguagem em Excesso \n' + str(tabulate(Totfinal)))
    #return (Totfinal)
