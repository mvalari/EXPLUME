import pickle
import numpy as np

ageBins=[]#=[0, 15, 20, 25, 40, 55, 65, 80]
cnsDates=[]
comIds=[]

#Read once to get comIds, ageBins, and dates
f=open('BTX_TD_AGEB_2008.ins','r')
for i,ln in enumerate(f):
 if i>0:
   comId,sex,age,type,ageB,cnt=ln.split()
   if comId not in comIds:comIds.append(comId)
   if age not in ageBins:ageBins.append(age)
   if ageB not in cnsDates:cnsDates.append(ageB)


data2d={}
for com in comIds:
    data2d[com]={}
    data2d[com]=np.zeros((len(ageBins),len(cnsDates)))

ageBins=sorted(ageBins)
cnsDates=sorted(cnsDates)

f.seek(0)
for i,ln in enumerate(f):
 if i>0:
   comId,sex,age,type,ageB,cnt=ln.split()
   
   data2d[comId][ageBins.index(age),cnsDates.index(ageB)]=data2d[comId][ageBins.index(age),cnsDates.index(ageB)]+int(cnt)
f.close()
data2dr={}
for com in data2d.keys():
    data2dr[com]=np.zeros((4,18))
    data2dr[com][0,:]=np.array([round(i) for i in 0.25*data2d[com][1,:]])
    data2dr[com][1,:]=np.array([round(i) for i in 0.75*data2d[com][1,:]])
    data2dr[com][2,:]=data2d[com][2,:]+data2d[com][3,:]+data2d[com][4,:]+data2d[com][5,:]
    data2dr[com][3,:]=data2d[com][6,:]+data2d[com][7,:]


fout=open('BUILDING_AGE_MVAL.ins','w')

for com in data2dr.keys():
   for date in range(len(cnsDates)):
       fout.write(com +' ;  '+str(int(data2dr[com][0,date])) +' ;   '+str(int(data2dr[com][1,date])) +' ;   '+str(int(data2dr[com][2,date])) +' ;   '+str(int(data2dr[com][3,date])) +' ;   '+str(cnsDates[date])+'\n')

fout.close()

