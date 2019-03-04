import pickle
import numpy as np
#1st key:CommuneId
#2nd key:AgeBin
#0:Homme ; 1:Femme
#0:Maison;1:Appartement;2:Autre
    # 11 : Before 1949
    # 12 : 1949-1974
    # 13 : 1975-1981
    # 14 : 1982-1989
    # 15 : 1990-1998
    # 21,22,23,24,25,26,27,28,29,30,31,32: 1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010
    # 99: Under construction (i assume 2011)    


ageBins=[0, 15, 20, 25, 40, 55, 60,65, 80]
ageBinsReduced=[0,4,25,65]
ageBuild=[11,12,13,14,15,21,22,23,24,25,26,27,28,29,30,31,32,99]

f=open('BTX_TD_AGEB_2008.dat','r')
data=pickle.load(f)
f.close()

data2d={}
for com in data.keys():
   data2d[com]=np.zeros((9,18))
   for i,ageBin in enumerate(sorted(data[com].keys())):
       for j,buildAge in enumerate(range(18)):
           data2d[com][i,j]=int(np.array(data[com][ageBin][0][0][j])+np.array(data[com][ageBin][0][1][j])+np.array(data[com][ageBin][0][2][j])\
                           +np.array(data[com][ageBin][1][0][j])+np.array(data[com][ageBin][1][1][j])+np.array(data[com][ageBin][1][2][j]))


data2dr={}
for com in data2d.keys():
    data2dr[com]=np.zeros((4,18))
    data2dr[com][0,:]=np.array([round(i) for i in 0.25*data2d[com][1,:]])
    data2dr[com][1,:]=np.array([round(i) for i in 0.75*data2d[com][1,:]])+data2d[com][2,:]
    data2dr[com][2,:]=data2d[com][3,:]+data2d[com][4,:]+data2d[com][5,:]+data2d[com][6,:]
    data2dr[com][3,:]=data2d[com][7,:]+data2d[com][8,:]

