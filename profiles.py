import numpy as np
from bisect import bisect_left
from random import uniform
import sys,cPickle

def strati_sampling(data_POP1,data_POP5,data_ACT2,virt_idf,virt_com,samp_dir):
                  
    sample=dict(); pp=dict(); data=dict(); k=-1; pop_idf=np.sum(data_POP1.values()); mod=0; obs=0
        
    #data_POP1: 1st dimension the gender, 2nd the 101 age groups
    #data_POP5: 1st dimension the gender, 2nd the 11 age groups, 3rd the 6 activity groups (active working, unemployed, retired, pupils/students, at_home, other inactive)
    #data_ACT2: 1st dimension the gender, 2nd the 11 age groups, 3rd salarie/non-salarie, 4th part-time/full-time
        
    #gen_id=0: Male
    #gen_id=1: Female
        
    #age_id=0: <=3
    #age_id=1: 4-24
    #age_id=2: 25-64
    #age_id=3: >=65
    
    #act_id=0: Infants, children before pre-primary school (0-3)
    #act_id=1: pupils, students without salary (4-18 for pupils and >18 for students)    
    #act_id=2: active working
    #act_id=3: active not working (unemployed)
    #act_id=4: inactive (retired, at_home, other)

    #con_id=0: no contract
    #con_id=1: Full-time
    #con_id=2: Part-time
        
    for com in [i for i in data_POP5.keys() if i<>95633]:
                        
        data[com]=np.zeros((2,4,5,3)); sample[com]=0.; tmp=np.sum(data_ACT2[com],axis=2) # i sum the salarie/non salarie that i dont need
                
        #sample[com]=int(round(float(virt_idf)*np.sum(data_POP1[com])/pop_idf,0))
        sample[com]=virt_com[com]
        
        for g in [0,1]: # i create a unified dataset from the data_POP1,data_POP5 and data_ACT2 files
        
            data_age=[sum(data_POP5[com][g,0:2]),sum(data_POP5[com][g,2:10]),data_POP5[com][g,10]]; data_con=[sum(tmp[g,0:2]),sum(tmp[g,2:10]),tmp[g,10]] 
            
            for a in range(3): # 3 age groups (4-24, 25-64, >=65)
                data[com][g,a+1,1,0]=data_age[a][3]; data[com][g,a+1,3,0]=data_age[a][1]; data[com][g,a+1,4,0]=data_age[a][2]+data_age[a][4]+data_age[a][5]
                data[com][g,a+1,2,1]=data_con[a][0]; data[com][g,a+1,2,2]=data_con[a][1] 

            data[com][g,0,0,0]=sum(data_POP1[com][g,0:4]); data[com][g,1,1,0]+=sum(data_POP1[com][g,4:15]) # i add the 0-3 and 3-14 age groups from the data_POP1 file
            
        data[com]=np.rint(data[com])                       
        for g in [i for i in [0,1] if np.sum(sample[com])>0]: # 1st strata (gender)               
            strata_lvl1=np.rint(sample[com]*np.sum(data[com][g])/np.sum(data[com]))
                        
            for age in [i for i in range(4) if np.sum(data[com][g])>0.]: # 2nd strata (age)                
                strata_lvl2=np.rint(strata_lvl1*np.sum(data[com][g,age])/np.sum(data[com][g]))            
                  
                for act in [i for i in range(5) if np.sum(data[com][g,age])>0.]: # 3rd strata (activity)                
                    strata_lvl3=np.rint(strata_lvl2*np.sum(data[com][g,age,act])/np.sum(data[com][g,age]))
                    
                    for con in [i for i in range(3) if np.sum(data[com][g,age,act,i])>0.]: # 4th strata (contract type)                                                
                        strata_lvl4=np.rint(strata_lvl3*np.sum(data[com][g,age,act,con])/np.sum(data[com][g,age,act]))
                        
                        for p in range(int(strata_lvl4)): k+=1; pp[k]=[com,g,age,act,con] # create the profile
                                                                                             
                        #if g==0 and age==2 and act==2 and con==1: w=str(com)+' '+str(round(data[com][g,age,act,con+1]/np.sum(data[com]),5))+' '+str(round(float(strata_lvl4)/sample[com],5)); print >> f1,w
                    
    st=open('data/data_strati.dat','wb'); cPickle.dump(data,st,-1); st.close(); return pp
    
def workplace(data_NAV1,data_DET_,pp): # Selects the workplace and schooling location (BTX_TD_NAV1_2009.xls, BTX_FM_DET_2009.xls)
          
    num_wp=len(data_NAV1[data_NAV1.keys()[0]][0][0]); cor_age={0:1,1:1,2:2,3:2,4:2,5:2,6:2,7:2,8:2,9:2,10:3}; data_DET_[78606][0]=[100,0,100]
   
    # id=1: Work in the commune they live
    # id=2: Work in another commune of the same department
    # id=3: Work in another department
      
    for k,v in pp.items(): # Here we loop all individuals

        wp_tot=0; id_selection=0; home=v[0]; gender=v[1]; age=v[2]; activity=v[3]; group_age=cor_age[age]; # The age in the profile is an age group. In the activity database age groups are divided per 5yr from (15 to 65)

        # for pupils/students i use the BTX_FM_DET_2009.xls file. in the TOTAL sheet it gives the number of people excersising schooling activities at the commune of residence. 
        # i assume that if not the person goes to id=2 and that there are not id=3 in this case.
        if activity==1: wp_p=[data_DET_[home][0][0]/data_DET_[home][0][2]]; id_selection=selection(wp_p,1); # pupils/students
                                                                
        elif activity==2: # only people that work
              
           if np.sum(data_NAV1[home][gender,group_age])>0.: wp_p=[data_NAV1[home][gender,group_age,w]/np.sum(data_NAV1[home][gender,group_age]) for w in xrange(num_wp)]; id_selection=selection(wp_p,1) 
           else: id_selection=1 # rare cases where contr_tot=0 then we use id=1
                
           if id_selection>=4: id_selection=3 # those who work abroad or in another region i use them as woking in another department
           
        pp[k].append(id_selection) # Add the sixth pp column, the workplace id
                
    return pp

def going_to(data_NAV1,data_DET_,data_DTR_,pp): # Selects the workplace commune (for the workers) (BTX_FM_DET_2009.xls) or the school place (for the pupils/students) (BTX_FM_DTR_2009.xls)

    list_deps=[75,77,78,91,92,93,94,95]
       
    # Valid only for pupils, students, people that work. The movements database covers only people that are going to other than their home commune.                
    for k,v in pp.items(): # Here we loop all individuals

        home=v[0]; activity=v[3]; workplace=v[5]; id_final=0; home_dep=int(home/1000.)        
        
        if activity in [1,2] and workplace in [2,3]:
           
           if activity==1: curr_dic=data_DET_ # pupils/students
           if activity==2: curr_dic=data_DTR_ # people that work           
           
           # there are no data for that home id. INSEE has data only data when a specific commune to commune flux is larger than 100
           if home not in [i for i in curr_dic.keys()]: id_final=1 # we write in the xls file a value of 1 so it will be treated in the diaries                                                                                        
           else: # there are data on where people commute so perform a stochastic selection

              tmp=[i for i in curr_dic[home] if int(i/1000.) in list_deps] # we just exclude the communes that are located outside IdF
                            
              if workplace==2: tmp=[i for i in tmp if int(i/1000.)==home_dep] # i select the communes in the same department
              else : tmp=[i for i in tmp if int(i/1000.)<>home_dep] # i select the communes in the other departments              
                 
              if tmp==[]: pp[k].append(1); continue  # cases where all the valid fluxes are less than <100 people
                                                                      
              if activity==1: tot=curr_dic[home][0][1] # the total population of people excersising schooling activities
              else: tot=np.sum(data_NAV1[home],axis=(0,1))[workplace-1] # this the total population of the commune corresponding to the working population of workplace 2 or 3
                                                                                                                 
              gt_p=[curr_dic[home][dest][0]/tot for dest in tmp]; gt_p.append(1.-sum(gt_p)); id_selection=selection(gt_p,0) # the remaining % represents the outward flux for <=100 destinations
                                                                           
              if id_selection==len(tmp): id_final=1
              else: id_final=tmp[id_selection]
                                                                
        elif workplace==1: id_final=home
        
        pp[k].append(id_final) # Add the seventh pp column, the workplace/school commune id
            
    return pp

def transport(data_NAV2,pp): # Selects the transport means to workplace (only for people that work) (BTX_TD_NAV2A_2009.xls)
     
    num_tr=len(data_NAV2[data_NAV2.keys()[0]][0])

    # id=1: they don't move (it means they work at home or stay at home)
    # id=2: On foot
    # id=3: On 2-wheels (bikes and bicycles)
    # id=4: Car
    # id=5: Public transport (the selection of which PT type is done inside the diaries)
    
    for k,v in pp.items(): # Here we loop all individuals

        id_selection=0; home=v[0]; gender=v[1]; activity=v[3]; workplace=v[5]; home_dep=int(int(home)/1000.)
        
        if activity==2: # only people that have a working activity
                   
           if np.sum(data_NAV2[home][gender,workplace-1])>0.:
              tr_p=np.array([data_NAV2[home][gender,workplace-1,t]/np.sum(data_NAV2[home][gender,workplace-1]) for t in xrange(num_tr)]) # % of each workplace group in the total pop of that commune, gender and age group
                          
	          # If the person is working far away then it cannot be in id=1 or go to work on foot, 2-wheels so i zero these percentages and transfer them to the car/PT modes.
              # Also if the selection is PT but there are no PT stops in the communes involved this will be handled inside the diaries (so the profile for those people are not correct here)
              if workplace>=2: tmp=sum(tr_p[0:3])/sum(tr_p[3::]); tr_p=[0.,0.,0.,tr_p[3]*(1.+tmp),tr_p[4]*(1.+tmp)]                                 
              id_selection=selection(tr_p,1)

              if id_selection==1: id_selection=2 # these are few cases and they complicate the diaries. i assume that these people walk on-foot
             
           else: id_selection=4 # to adress some very very few cases where total populations zero for some gender/workplace combinations
                                                               
        pp[k].append(id_selection) # Add the eighth pp column, the transportation type
        
    return pp

def houseage(data_AGEB,pp): # Selects the house age (BUILDING_AGE.ins)
     	                                        
    for k,v in pp.items(): home=v[0]; age=v[2]; dt=data_AGEB[home][age]; a_p=[a/np.sum(dt) for a in dt]; id_selection=selection(a_p,1); pp[k].append(id_selection) # adds the ninth pp column, the house age group
                              
    return pp

def officeage(pp): # Selects the office age. the same will be used for schools

    # statistics per department on the share of buildings in the 3 groups: <=1974, 1975-2004, >=2005 
    a_p={75:[0.76,0.22,0.02],77:[0.56,0.41,0.03],78:[0.55,0.42,0.03],91:[0.32,0.58,0.10],92:[0.55,0.41,0.04],93:[0.39,0.59,0.02],94:[0.58,0.39,0.03],95:[0.45,0.50,0.05]}

    for k,v in pp.items(): office_dep=int(int(v[0])/1000.); id_selection=selection(a_p[office_dep],1); pp[k].append(id_selection)
        		
    return pp    

#################################################################
##### I dont use the information below for the time being #######
#################################################################

def house(data_PRINC1,pp): # Selects the house type and number of rooms in the house (BTX_TD_PRINC1_2009.xls)
     
    num_h=len(data_PRINC1[data_PRINC1.keys()[0]][0]); num_r=len(data_PRINC1[data_PRINC1.keys()[0]][0][0]); list_ages_in_house_type=[0,20,25,40,55,65,80]

    # id=1: house
    # id=2: apartment
    # id=3: other (shelters, hotels etc)

    for k,v in pp.items(): # Here we loop all individuals, v[1]: home, v[2]: age

        h_p=[]; r_p=[]; h_tot=0; r_tot=0.; id_selection=0; home=v[1]; age=v[2]

        # The age in the profile (pp) is the exact number e.g. 46 years old. But in the house type database age groups are divided in 7 groups ([0,20,25,40,55,65,80])
        for index_a,a in enumerate(list_ages_in_house_type[0:6]): # The 80 number is not iterated
            if age>=a and age<list_ages_in_house_type[index_a+1]: group_age=index_a # in data_PRINC1 the "age" dimension of the array has integers e.g. element 0 is age <20
            if age>=80: group_age=6 # element 6 is age>=80
            
        # 1st lets try to select a house type for this age group
        for h in xrange(num_h): h_tot=h_tot+sum(data_PRINC1[home][group_age][h][:]) # Sums all house types and all room number groups within this age group
        
        # In some group ages (mostly for <20) the data on the PRINC1 file do not match those of the population. In PRINC1 almost all communes have zero number of houses where people <20age live.
        # As a result a house id cannot be selected for them. When this happens we target the PRINC1 statistics of the following age groups
        while h_tot==0.:              
              if group_age<6: group_age=group_age+1
              if group_age==6: group_age=group_age-1 # for the >80 we select the previous group
              for h in xrange(num_h): h_tot=h_tot+sum(data_PRINC1[home][group_age][h][:]) # Re-do the calculation for the updated age xrange
              if h_tot<>0.: break   
        
        try: 
            for h in xrange(num_h): h_p.append(sum(data_PRINC1[home][group_age][h][:])/h_tot) # % of each house type group in the total pop of that commune,age group
            id_selection=selection(h_p,1); id_house_type=id_selection
        except: ZeroDivisionError; id_house_type=-9999 # if h_tot=0 then we set this individual as invalid

        pp[k].append(id_selection) # Add the ninth pp column, the house type
                   
        if pp[k][-1] not in [1,2,3]: print h_p,pp[k][-1]; sys.exit(0)          
                   
        # 2nd lets try to select a room number group for the above selected house type for this age group           
        r_tot=sum(data_PRINC1[home][group_age][id_house_type-1][:]) # Sums all #rooms groups within this commune,age group and house type
        try:
            for r in xrange(num_r): r_p.append(data_PRINC1[home][group_age][id_house_type-1][r]/r_tot) # % of each house type group in the total pop of that commune,age group
            id_selection=selection(r_p,1)
        except: ZeroDivisionError; id_selection=-9999 # if r_tot=0 then we set this individual as invalid
            
        pp[k].append(id_selection)
                          
    return pp

def surface(data_PRINC26,pp): # Selects the surface of the house

    # The fit function used above
    def poly_fit(x, c0, c1, c2, c3): return c3*x**3 + c2*x**2 + c1*x + c0    

    from numpy import polynomial as P

    for k,v in pp.items(): pp[k].append(0)
    return pp # i deactivate this for the time beeing (i dont use it and it is slow)
    
    num_s=len(data_PRINC26[data_PRINC26.keys()[0]]) # group types according to surface
        
    # id=1: <40 m2
    # id=2: 40-100 m2
    # id=3: >100 m2

    for k,v in pp.items(): # Here we loop all individuals, v[1]: home, v[8]: house type, v[9]: # rooms

        home=v[1];house=v[8];rooms=v[9]; s_p=[];count=[];s_tot=0;id_selection=-9999;tol=0.03 # tolerance to the matching of the fitting                
        s_tot=sum(data_PRINC26[home][house-1][:]) # Sums all surface types within this house type group
        
        try:
           for s in xrange(num_s): s_p.append(data_PRINC26[home][house-1][s]/s_tot); count.append(data_PRINC26[home][house-1][s]) # % of each surface type group in the total pop of that commune, house type group                
           id_selection=selection(s_p,1)
              
           # Some of the statistics are strange, e.g. 1 room houses to be >100m2. These are limited cases and the model is stochastic but lets fix that issue
           if rooms in [5,6] and id_selection==1: id_selection=3 # 5 or 6 room houses that were selected as <40m2
           if rooms in [1,2] and id_selection==3: id_selection=1 # 1 or 2 room houses that were selected as >100m2

           # We fit a curve to the 3 numbers (house surface groups): X is the surface, Y is the # of houses
           x = [10.,11.,40.,100.,151.] # I assume that the number for <40m2 goes to x=11 (smallest houses in France), 40-100m2: x=40 and >100m2: x=101
           y = [0.,count[0],count[1],count[2],0.] # the zeros means that for <11m2 and >150m2 also there are no houses with thes surfaces
           c, stats = P.polynomial.polyfit(x,y,3,full=True) # its a 3rd degree polyfit

           # scan all surface values between the xrange group selected by id_selection and find the exact surface.
           # Here we have to respect also the # of rooms. The random generator might give a house of 20m2 but this cannot be a 4 room house.
           avg_room_size=(x[id_selection+1]-x[id_selection])/rooms; start_surface=int(x[id_selection])+(rooms-1)*avg_room_size; end_surface=int(x[id_selection])+rooms*avg_room_size
           if rooms in [1,2] and id_selection==2: end_surface=30*rooms+20 # Control: i assume that e.g. 1 room house >40m2 cannot reach all the way up to 100m2
              
           r1=poly_fit(start_surface,c[0],c[1],c[2],c[3]); r2=poly_fit(end_surface,c[0],c[1],c[2],c[3]); tmp=[r1,r2]
           if min(tmp)<0: tmp=[0.,max(tmp)] # Due to the fitting sometimes r1 or r2 might be negative.

           # Now select a Y random float number (# of houses) between r1 and r2 and try to find in which X (house surface) it belongs
           rn=uniform(min(tmp),max(tmp))
                                          
           for s in np.arange(start_surface,end_surface,0.1):                                    
               if poly_fit(s,c[0],c[1],c[2],c[3])>rn-tol*rn and poly_fit(s,c[0],c[1],c[2],c[3])<rn+tol*rn: actual_surf=round(s,1); break
               else: actual_surf=50 # to adress some times a case where you havean error f
                                                                   
        except: ZeroDivisionError; id_selection=-9999 # if s_tot=0 then we set this individual as invalid               
                                    
        pp[k].append(actual_surf)
                          
    return pp

def heating(data_PRINC3,pp): # Selects the house heating type

    for k,v in pp.items(): pp[k].append(0)
    return pp

    num_h=len(data_PRINC3[data_PRINC3.keys()[0]]) # group types according to heating type

    # id=1: district heating
    # id=2: natural gas
    # id=3: oil, mazut
    # id=4: electricity
    # id=5: gas from a tank
    # id=6: other

    for k,v in pp.items(): # Here we loop all individuals, v[1]: home, v[8]: house type+ #rooms

        h_p=[]; h_tot=0;id_selection=-9999; home=v[1];house=v[8]; h_tot=sum(data_PRINC3[home][house-1][:]) # Sums all heating type groups this house type group
        
        try:
           for h in xrange(num_h): h_p.append(data_PRINC3[home][house-1][h]/h_tot) # % of each surface type group in the total pop of that commune, house type group
           id_selection=selection(h_p,1)
        except: ZeroDivisionError; id_selection=-9999 # if h_tot=0 then we set this individual as invalid
                                    
        pp[k].append(id_selection)
                          
    return pp
	
def create_profile_file(pp,update,samp_dir):
    
    import os.path 
    
    try: os.makedirs('samples')
    except: pass 
    try: os.makedirs(samp_dir)
    except: pass 

    age=dict(); hdr=''; f1=open(samp_dir+'/profiles.ins', "w"); list_d=['id',"home","gender","age","activity","contract","workplace","going_to","transport","houseage","officeage"]
    
    if update==1: list_d.append("moves")
    
    for d in list_d: hdr+=d+' ' # The header
    print >> f1,hdr
    
    for key,v in pp.items():
        val=''; val=str(int(key)+1)+' '
        for i in xrange(len(list_d)-1): val+=str(v[i])+' '
                  
        age[key]=[]; age[key].append(v[0]); age[key].append(v[2]);  print >> f1,val # in the age.dat i put home and age
                                
    f=open(samp_dir+'/age.dat','wb'); cPickle.dump(age,f,-1); f.close(); f=open(samp_dir+'/profiles.dat','wb'); cPickle.dump(pp,f,-1); f.close(); f1.close()
    
def selection(perc_list,write_id):

    sum_p=0.; id_selection=0; sorted_list=[0.]           
    
    for index_p,p in enumerate(perc_list): sum_p+=p; sorted_list.append(sum_p) #in index_p=0, sum_p is zero and it adds this as the 1st element of sorted_list
        
    if len(sorted_list)>2 and 0.999 < sum(sorted_list) < 1.001: sorted_list[len(sorted_list)-1]=1. # Some times the sum gives a slightly lower than 1 number (e.g. 0.9999999998)
    elif sorted_list[len(sorted_list)-1]<>1.: sorted_list.append(1.) # if the array has only 0. and one % value then the above if will mistakenly remove the % and place an 1.
        
    tau = uniform(0,1)                    
    for index_p,p in enumerate(sorted_list[0:len(sorted_list)-1]): # The last item is not iterated its p=1.0           
        if tau>=p and tau<sorted_list[index_p+1]: id_selection=index_p+write_id; break
        
    return id_selection   
