from read_data import read_communes_coords,read_squales
from profiles import create_profile_file
from expo import get_progress_bar,get_i_j
from visual import plot_stats_diaries,get_projection
from shapely.geometry import Point, LineString, Polygon
from collections import defaultdict
from bisect import bisect_left
from random import shuffle,choice,uniform
from decimal import Decimal
from sets import Set
from os.path import join
import math,os,cPickle,sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def select_diary(stops,speeds2grid,data_diaries,pp,weekend,samp_dir,grid_coords,dom,corresp,grid_polyg,dx,dy,write_asc_dia,static):

    per_act_pop=[0.,0.,0.]; workers,school=0,0; mot_pool=[[0,0,0,0,0,0] for _ in xrange(len(pp))]; res=int(math.sqrt((grid_polyg[1][0]).area)*100); proj=get_projection(dx,dy,grid_coords,res)
    acts=[int(pp[k][3]) for k in pp.keys()]; IO_metro=read_squales(); stk_str=['','_STATIC']; communes=dict(); communes=read_communes_coords() # read the IdF.txt file with the communes centroid coords                                                        
    
    if weekend==0: print ''; wkd='wday'; print 'preparing calculation...'             
    else: data_diaries=change_weekend_statistics(data_diaries); wkd='wend'
    
    # At this point we need the total population for actives, pupil/students, inactive. This is needed for the probability of movement in the check_for_movement routine. I add the 0-3 group to inactives
    st='data/data_strati.dat'; data_strati=cPickle.load(open(st,'rb'))            
    for c in data_strati.keys(): per_act_pop[0]+=np.sum(data_strati[c][:,:,2:4,:]); per_act_pop[1]+=np.sum(data_strati[c][:,:,1,:]); per_act_pop[2]+=np.sum(data_strati[c][:,:,4,:])+np.sum(data_strati[c][:,:,0,:])
    
    for k in xrange(len(pp)): # work and school is also valid in the weekend
        activity=int(pp[k][3]) # in mot_pool motives are: 'Work','Prof','School','Market','Fun','Personal'        
        if activity<>1: mot_pool[k][2]=1 # invalid school for not pupil/student
        if activity<>2: mot_pool[k][0]=1; mot_pool[k][1]=1 # invalid work/professional for all working         
        if activity==1: school+=1
        elif activity==2: workers+=1
        
    if weekend==0: print ''; print 'your sample has: '+str(workers)+' workers and '+str(school)+' pupils/students'; print ''
         
    flux_deps,flux_mode_pt=get_fluxes(data_diaries) # construct the fluxes dictionaries
              
    # here i start creating the diaries for each hour of the day
    ind=diary(data_diaries,pp,communes,per_act_pop,mot_pool,stops,weekend,flux_deps,flux_mode_pt,acts,speeds2grid,corresp,samp_dir,dom,grid_polyg,dx,dy,IO_metro,proj,static)
    write_diaries(ind,pp,weekend,samp_dir,write_asc_dia); f=open(samp_dir+'/diaries_'+wkd+stk_str[static]+'.dat','wb'); cPickle.dump(ind,f,-1); f.close()    
	    
def create_polygons(grid_coords,dx,dy,dom):
    
    plt.ioff(); fig=plt.figure(); ax = fig.add_subplot(111)

    LL=get_cell_coord(1,1,dx,grid_coords); LR=get_cell_coord(dx+1,1,dx,grid_coords); ax.set_xlim(LL[0]-1.,LR[0]+1.) # get the coords of corner cells of the domain
    UR=get_cell_coord(dx+1,dy+1,dx,grid_coords); ax.set_ylim(LL[1]-1.,UR[1]+1.)
    
    grid_polyg=dict()
    for y in xrange(dy):
        for x in xrange(dx): grid_polyg[y*dx+x+1]=[]; points=find_corners(x+1,y+1,grid_coords,dx,dy); grid_polyg[y*dx+x+1].append(Polygon(points)); poly=patches.Polygon(points,fc='none',lw=0.2); ax.add_patch(poly)
        
    f=open('data/depts_idf.prn','r'); deps=dict(); points=[]        
    for l in f: # get the perimeter points of each department
        if '@' in l:  	   
            if len(points)>0: deps[dp].append(points); points=[]
            dp=l.replace(" ","")[1:3]
            if dp<>'10': deps[dp]=[]     
        else: lon,lat=l.strip().split(' ')[0],l.strip().split(' ')[-1]; tmp=[float(lon),float(lat)]; points.append(tmp)
    
    for d in deps.keys(): poly=patches.Polygon(deps[d][0],fc='none',lw=0.2); ax.add_patch(poly)
        
    fig.savefig('domains/'+dom+'/'+dom+'.pdf', dpi=300); plt.close(fig); st=open('domains/'+dom+'/grid_polygons.dat','wb'); cPickle.dump(grid_polyg,st,-1); st.close()   

    return grid_polyg

def get_cell_coord(cell_x,cell_y,dx,grid_coords):
    
    cell_id=(cell_y-1)*dx+cell_x; coord=grid_coords[cell_id] # the grid_coords dictionary of id_cells starts from 1 NOT 0.
    
    return coord

def find_corners(x,y,grid_coords,dx,dy):
 
    points=[]
    id_cell=(y-1)*(dx+1)+x; x1=grid_coords[id_cell][0]; y1=grid_coords[id_cell][1]; p=[x1,y1]; points.append(p) # LL corner 
    id_cell=(y-1)*(dx+1)+(x+1); x1=grid_coords[id_cell][0]; y1=grid_coords[id_cell][1]; p=[x1,y1]; points.append(p) # LR corner 
    id_cell=y*(dx+1)+(x+1); x1=grid_coords[id_cell][0]; y1=grid_coords[id_cell][1]; p=[x1,y1]; points.append(p) # UR corner    
    id_cell=y*(dx+1)+x; x1=grid_coords[id_cell][0]; y1=grid_coords[id_cell][1]; p=[x1,y1]; points.append(p) # UL corner   
                      
    return points 

def comm2cells(dx,dom,grid_polyg,nb,pop2com,pop,coms,indx,c2c):

    previous=0; done=0; cor=dict(); ranges=int(nb/4)+10

    for c in coms: # c is the commune and pop_cell the index of the pop grid cell found inside that commune
                
        p=[]; pop_p=[]; cor[c]=np.zeros((ranges,3)) # id_cell (from 1 NOT 0), lon, lat (centroid)
        for pop_cell in pop2com[c]: pop_cell_x=pop_cell[1]; pop_cell_y=pop_cell[0]; p.append(pop['pop'][pop_cell_x][pop_cell_y]) # the pop 1km cells move from left to right and from north to south           
        pop_p=[cl/sum(p) for cl in p] # the pop_p contains the % of each pop cell to the total population of all the cells under a commune   
                
        for k_r in xrange(ranges):
                           
            # now i have the id of the population grid in which the person is located (curr_lon,lat refers to the centroids of the population grid cells)              
            try: id_selection=selection(pop_p,0); wpc=pop2com[c][id_selection]; curr_lon,curr_lat=pop['lon'][wpc[1]][wpc[0]],pop['lat'][wpc[1]][wpc[0]] 
            except: pass          
            
            point=Point(curr_lon,curr_lat) # create the point of the pop 1km grid cell (based on the centroid)
            for cell in c2c[c]: # the last col and row of the grid is omitted           
                i,j=get_i_j(cell,dx) # chimere grid cell and location inside that grid (the location of the centroid of the pop grid cell)                
                if (grid_polyg[(j-1)*dx+i][0]).contains(point): cor[c][k_r,0],cor[c][k_r,1],cor[c][k_r,2]=cell,round(curr_lon,3),round(curr_lat,3); break
            
            # covers some rare cases where there is not correspondence
            if cor[c][k_r,0]==0.: cell=choice(c2c[c]); i,j=get_i_j(cell,dx); gr=grid_polyg[(j-1)*dx+i][0]; cn=gr.centroid; cor[c][k_r,0],cor[c][k_r,1],cor[c][k_r,2]=cell,round(cn.x,3),round(cn.y,3)
                        
        if indx==0: done+=1; previous=get_progress_bar(done,len(coms),previous)
                
    f=open('domains/'+dom+'/comm2cells_'+str(indx)+'.dat','wb'); cPickle.dump(cor,f,-1); f.close()    
    
def get_speeds(roads,grid_polyg,dom):
    
    speeds2grid=dict() # contains the grid cells as keys and a speed
        
    for r in roads: # all road segments          
        point=Point(roads[r][0],roads[r][1]) # create the point (based on the Xa,Ya of the network file)
        for c in [i for i in grid_polyg.keys() if (grid_polyg[i][0]).contains(point)]:             
            if c not in speeds2grid.keys(): speeds2grid[c]=[[],[]]
            speeds2grid[c][0].append(roads[r][3]); speeds2grid[c][1].append(roads[r][4]) # appends the fluxes and speeds of all the roads found in this grid cell
     
    st=open('domains/'+dom+'/speeds2grid.dat','wb'); cPickle.dump(speeds2grid,st,-1); st.close()   
              
    return speeds2grid
    
def change_weekend_statistics(data_diaries):

    changes_nb=[0.7,0.6,0.72]; changes_motive=[0.113,0.2,0.14,1.1,1.73,0.44]; coeff=[0.0047,0.0047,0.007,0.007,0.0094,0.012,0.012,0.035,0.047,0.054,0.07,0.07,0.065,0.087,0.091,0.088,0.098,0.07,0.065,0.044,0.023,0.014,0.012,0.0102]    
    changes_dis=[[1. for i in xrange(8)],[1. for i in xrange(8)],[1.2 for i in xrange(8)],[1.3 for i in xrange(8)],[1. for i in xrange(8)]]
        
    # changes in total movements for the weekend, this is based on the total number of movements per category, it is not hourly distributed
    # EGT, December 2006: Les deplacements de fin de semaine, page 5 (2001 data)    
    for new_act in xrange(2):  # new_act is actives, in-actives, students
        tmp=0.; coeff1=[]
        for h in xrange(24): tmp+=data_diaries['Nb'][new_act][h]; coeff1.append(coeff[h]); coeff1.append(coeff[h])
        for h in xrange(24): data_diaries['Nb'][new_act][h]=round(tmp*coeff[h]*changes_nb[new_act],1) # i redistribute within the day based on the plot of page 8 in EGT janvier 2013

    # changes in movements per hour and motive for the weekend, EGT, December 2006: Les deplacements de fin de semaine (2001 data)
    for m in xrange(6):
        tmp=0.; tmp=sum(data_diaries['obs_mot'][m][h] for h in xrange(24))            
        for h in xrange(24): data_diaries['obs_mot'][m][h]=round(tmp*coeff1[h]*changes_motive[m],1)
            
    # the reductions comes from the EGT, December 2006, page 11 (2001 data) and they correspond to changes in the average time not the distance.
    # the average time actually changes significantly only for foot and bicycles not for the others. I assume that the increase in time is due to increase of the distance which is not entirely true
    for t in xrange(5): # motives are: Work,Proff,School,Market,Fun,Personal
        for d in xrange(8): data_diaries['distance'][t][d]=data_diaries['distance'][t][d]*changes_dis[t][d] # modes of transport: public, car, foot, bicycle, bike 

    return data_diaries
        
def get_fluxes(data_diaries):
        
    flux_mode_pt=[[],[],[],[]]; flux_deps={'75':[[],[],[],[]],'92':[[],[],[],[]],'93':[[],[],[],[]],'94':[[],[],[],[]],'77':[[],[],[],[]],'78':[[],[],[],[]],'91':[[],[],[],[]],'95':[[],[],[],[]]}
    deps=['75','77','78','91','92','93','94','95']
        
    for index_st,st in enumerate(data_diaries['public']): 
        for index_end in xrange(4): flux_mode_pt[index_st].append([float(x) for x in st[index_end].split(',')]) # the data for bus,metro,RER are separated by a comma                       
                                
    for end in data_diaries['flux']:
        for index_st in xrange(8):
            tmp=[float(x) for x in end[index_st].split(',')]
            for t in xrange(4): flux_deps[deps[index_st]][t].append(tmp[t]) # public, car, foot, 2-wheels
                
    return flux_deps,flux_mode_pt
       
def diary(data_diaries,pp,communes,per_act_pop,mot_pool,stops,weekend,flux_deps,flux_mode_pt,acts,speeds2grid,corresp,samp_dir,dom,grid_polyg,dx,dy,IO_metro,proj,static):

    trans_str=['Home','Car','Foot','Bicycle','Bike','Bus','Metro','RER','tram']; motive_str=['Work','Prof','School','Market','Fun','Personal','at_Home','to_Home']    
    deps=[75,77,78,91,92,93,94,95]; nb=len(pp); cs=[0 for i in xrange(nb)]; wk_str={0:'wday',1:'wend'}; rem_path=[[] for k in xrange(nb)]; 
    temp=[60,75,105,135,165,195,225,255,285,315,345,375]; tmp1=[0.25,0.13,0.05,0.04,0.05,0.04,0.04,0.03,0.04,0.1,0.1,0.13] # time that a pre-primary school child passes at child care (in minutes) and fraction of the pop at each time selection    
    mot={'at_Home':[1,60],'Work':[6,61],'Prof':[6,62],'Personal':[4,63],'Market':[0,64],'Fun':[0,65],'School':[7,66],'Car':[21,21],'Bus':[22,22],'RER':[23,23],'Metro':[24,24],'tram':[25,25],'Foot':[31,31],'Bicycle':[32,32],'Bike':[33,33]}
    age_groups=[0,15,25,35,45,55,65,75,150]; wr_trv_traj=[[[] for h in xrange(24)] for k in xrange(nb)]; time_curr_dest=[[[' '] for _ in xrange(24)] for _ in xrange(nb)]
    ind=[[[0 for _ in xrange(7)] for _ in xrange(24)] for _ in xrange(nb)]; make_video=[]; valid_cells=[[] for _ in xrange(nb)]; stk_str=['','_STATIC']
    mv_per_mot=[sum(data_diaries['Nb_mot'][m][h] for m in xrange(6)) for h in range(24)]
    
    def get_indents(hr):              
        if hr<=9: o=1; o1=2; o2=10; o3=15
        else: o=0; o1=0; o2=12; o3=17    
        return o,o1,o2,o3
         
    def get_mot_per(hr,tot): mot_p_start=[data_diaries['Nb_mot'][m][hr]/tot for m in xrange(6)]; return mot_p_start               
    
    for h in xrange(24):
                
        indends=[[0,0,0,0] for i in range(3)]
        for t in [-2,-1,0]:
            try: o,o1,o2,o3=get_indents(h+t); indends[abs(t)]=[o,o1,o2,o3]
            except: pass 
                
        dist_travel=0.; allow_check=1; pep_mv=[0.,0.,0.,0.]; 
        if weekend==0: print 'hour '+str(h)
        
        for k in xrange(nb): # Here we loop all individuals
                
            c=0.; home=pp[k][0]; k_r=int((k/4)//1); age=pp[k][2]; activity=int(pp[k][3]); mot_pool[k][3]=0
            # we create a 24 hour array for each person. then in each hour we will put 5 elements that change every hour. These elements are: 
            # 1) remaining time for travelling
            # 2) remaining time for activity
            # 3) transport type used (1: car, 2: foot, 3: bicycle, 4: bike, 5: bus, 6: metro, 7: RER, 8: tram)
            # 4) motive (Work, Professional business, School, Market, Recreation, Personal business, at home, to home)
            # 5) destination commune
            # 6) How many moves were performed so far
            # 7) the distance travelled for this movement
            # These will tell the current hour if there is some remaining time for travelling, activity from the previous hour. The 7 elements are written with the fill_pr_hr subroutine
            
            # before determining what to do we read any remaining period from the previous hour.
            if h>0 and cs[k] in [1,3,4]:
               time_travel,tmACT_spend=ind[k][h-1][0],ind[k][h-1][1]; trans=trans_str.index(ind[k][h-1][2]); motive=motive_str.index(ind[k][h-1][3]); dest_com=ind[k][h-1][4]; dist_travel=ind[k][h-1][6]
               if ind[k][h-1][6]<>0.: dist_travel=0. # so the distance will be written (in the txt files) only in the beginning of the travelling
               case=write_ind(ind,k,h,trans,time_travel,tmACT_spend,motive,home,dest_com,1,dist_travel); cs[k]=case # we re-determine the case
                   
               if time_travel<>0.: wr_trv_traj[k][h].append(rem_path[k]); rem_path[k]=[]
               else: wr_trv_traj[int(k)][h].append([])             
                             
            # in the first hour or in h>0 and case=0,5,10 we want to know how the individual will proceed through the day. 10 is equivalent to 0 here 
            # In cases 1,3,4 this is not valid as the person already has a task for the whole duration of the current hour.
            if h==0 or (h>0 and cs[k] in [0,5,10]):

               # regular case where we must see what the person will do next            
               prev_tr_tm=0.; num_moves=0; passed=0; home_dep=str(int(int(home)/1000.)); traj=[]; cc=int(home); mtv=ind[k][h][3]
    
               if h>0: cc=int(ind[k][h-1][4]); num_moves=ind[k][h-1][5]; prev_tr_tm=ind[k][h-1][0]+ind[k][h-1][1] # get how many movements were performed so far by that individual
                              
               # check for movement. This is valid for all individuals at h=0 or if the person is at home. In the latter case the person can continue to be at home or move
               on_the_move,pep_mv,allow_check=check_for_movement(data_diaries,activity,h,per_act_pop,home_dep,ind,pep_mv,allow_check,k,weekend,nb,acts,static)
                              
               # There are cases in which a person stays at home. The 2nd case is valid REGARDLESS of the movement check:
               # 1) if the person in the previous hour was at home (cs[10]) then it will stay at home
               # 2) the person was at home and we are in hour 23, so we don't initiate a new movement
               if (cs[k]==0 and (on_the_move==2 or h==23)) or cs[k]==10: dest_com=home; motive=6; trans,time_travel,tmACT_spend,dist_travel=0,0,0,0 # so it will go to case 0  
               else: # It will enter here if the person moves and cases are 0 or 5.
                                    
                  # motive=7 is that the person is going back home. This happens when:
                  # 1) Between the hours 22:00 to 00:00 when a person went out to do something then it comes back home, it cannot exercise another activity.               
                  # 2) the person was out for an activity and the movement check came negative                                         
                  if on_the_move==2 or (h in [22,23] and cs[k]==5): motive=7                                         
                  else: # the person is moving but with what kind of motive?
                                                                           
                     # work and school is allowed only once per day. if some activity in not valid for a group or already exercised them i transfer the % to the other activities
                     red1=[]; red2=[]; sum_red1=0.; sum_red2=0.; mot_p=get_mot_per(h,mv_per_mot[h])
                                                                                                         
                     if h>=17 and activity in [0,1]: mot_pool[k][2]=1 # no school activities after 17:00 for ages<=14
                     if h in [21,22,23,0,1,2,3,4,5,6,7]: mot_pool[k][3]=1 # the market is closed between 21:00 and 7:00

                     for i in [j for j in xrange(6) if mot_pool[k][j]==1]: red1.append(i); sum_red1+=mot_p[i]; mot_p[i]=0. # activities already exercised or invalid for that person             
                     for i in [j for j in xrange(6) if j not in red1]: red2.append(i); sum_red2+=mot_p[i] # activities NOT exercised yet          
                
                     if sum_red2<>0.:
                       for i in [j for j in xrange(len(mot_p)) if j in red2]: mot_p[i]+=(mot_p[i]/sum_red2)*sum_red1 # reconstruct the mot_p array
                       motive=selection(mot_p,0)
                       if motive in [0,2]: mot_pool[k][motive]=1 # Work (motive=0), Proff, School, Market, Fun, Personal (motive=5)                             
                     else: motive=7 # if sum_red2=0 it means that this person cannot exercise any activity at this hour; e.g some activities are done and the others are not allowed
                                          
                  # get destination, distance, transport type
                  trsp_type,dst_prf,workplace=int(pp[k][7]),int(pp[k][6]),int(pp[k][5]); dest_com,dist_travel,trans=get_dest_dist_trans(cc,data_diaries,trsp_type,dst_prf,workplace,communes,home,deps,motive,home_dep,h,mtv,stops,flux_deps,flux_mode_pt) 
                                                                        
                  # based on the above i want to know how much time is needed to perform a movement and the time spend on that activity
                  crs1,crs2=corresp[cc][k_r],corresp[dest_com][k_r]; time_travel,passed,traj=travel_time(deps,data_diaries,trans,dist_travel,dest_com,prev_tr_tm,cc,speeds2grid,crs1,crs2,grid_polyg,dx,dy)   
                  
                  if motive==7: tmACT_spend=0 # motive=7 means the person goes back home to finish the day                                        
                  elif activity==0 and motive==2: id_selection=selection(tmp1,0); tmACT_spend=temp[id_selection] # for 0-3 time spend in child care (BUDGET-ESPACE-TEMPS-ACTIVITES DES ENFANTS: LIEUX DE LOISIRS ET DE GARDE)
                  else: 
                     indx=deps.index(int(home_dep)); tmACT_spend=int(round(data_diaries['time_act'][motive][indx],0)) # i take the time from the excel input file
                     if motive==0 and int(pp[k][4])==2: tmACT_spend=int(np.rint(tmACT_spend/2.)) # if the person has a part-time contract then the person works half hours.                  
          
                  tm_tot=prev_tr_tm+time_travel+tmACT_spend  # to ensure that an activity will not stop at e.g. 20:00, it helps with the writing
                  if Decimal(float((tm_tot))/60.)._isinteger() or Decimal(float((tm_tot))/59.)._isinteger(): tmACT_spend=max(0,tmACT_spend-5)
                                       
               # write the dictionary and keep in which case the person is assigned 
               if str(ind[k][h][len(ind[k][h])-1])[6:11]<>'23.59':
                  case=write_ind(ind,k,h,trans,time_travel,tmACT_spend,motive,home,dest_com,0,dist_travel); cs[k]=case                  
                  if traj<>[]: wr_trv_traj[k][h].append(traj[0:passed]); rem_path[k]=traj[passed::] # cells from 1 NOT 0    
               
            # get the time spend in the activity and current location for all activities of all hours of that individual: 
            # element 0 is the time spend for that activity in this hour
            # element 1 is the id_cell of the current location (starting from 1)
            # element 2 is the id_cell of the destination location (for those who travel)
            # element 3 is the I/O case
            # element 4 is the actual I/O
            # element 5 the start commune
            # element 6 the end commune			
            # element 7 is the motive                          
            for a,act in enumerate(ind[k][h][7::]): # all the different activities performed in this hour            
            
                if a>0: time_curr_dest[k][h].append(' ') # it was initialized only for one action                
                inde=indends[0]; 
                
                try: spend=str(int(act[9-inde[1]:11-inde[1]])-int(act[3-inde[0]:5-inde[0]])+1); c+=int(spend)                
                except: print act
                
                for m in [i for i in mot.keys() if i in act]: IO_case=mot[m][0]; mtv=mot[m][1]      
                if mtv in [64,65]: IO_case=get_indoor_outdoor_act(act,age) # this is for market or recreation: here it will return case 4 (indoors) or 5 (outdoors)
                                
                if IO_case<=7: curr_loc=str(act)[o2:o3]; dest_loc=curr_loc; dest_cell=int(corresp[int(dest_loc)][k_r,0]) # the person is static (indoors or outdoors)
                else: # the person is moving; find the beginning and destination of this person    

                    # the destination location is stored in the 4th element of the ind array of this hour. cells start from 1    
                    dest_loc=int(ind[k][h][4]); dest_cell=int(corresp[dest_loc][k_r,0]) 
                    
                    if a==0 and (h==0 or (h==1 and '(' in str(ind[k][0][7]))): curr_loc=home                                                          
                    elif a==0 and h>0:     
                        inde=indends[1]
                                               
                        if 'T(' not in str(ind[k][h-1][-1]): curr_loc=str(ind[k][h-1][-1])[inde[2]:inde[3]] # the person began the hour with a movement and during the previous hour was not moving                                                                                  
                        elif len(ind[k][h-1])>8: curr_loc=str(ind[k][h-1][-2])[inde[2]:inde[3]] # the current commune is in the previous hour, 2 acts before the end                                                                                                                      
                        elif len(ind[k][h-1])==8: 
                             inde=indends[2]; curr_loc=str(ind[k][h-2][7])[inde[2]:inde[3]] # the current commune is 2 hours before (last act)
                             if '(' in curr_loc: curr_loc=str(ind[k][h-2][-1])[inde[2]:inde[3]] # some bugged case
                                                         
                    else: curr_loc=str(ind[k][h][7+a-1])[inde[2]:inde[3]] # the current commune is in the previous act of this hour                                                                             
                    
                    if wk_str[weekend]=='wday' and int(curr_loc)<>int(dest_loc): # to make the video
                       p1=proj(corresp[int(curr_loc)][k_r,1],corresp[int(curr_loc)][k_r,2]); p2=proj(corresp[int(dest_loc)][k_r,1],corresp[int(dest_loc)][k_r,2]); make_video.append([[p1[0],p1[1]],[p2[0],p2[1]]])
                                
                I_O=0.; curr_cell=int(corresp[int(curr_loc)][k_r,0])                 
                if IO_case not in [1,4,6,7]: I_O=get_IO_ratio(IO_case,IO_metro,curr_loc,dest_loc) # IOs only for transport. Indoors is taken inside exposure from the SIREN database	                                               
                if c<=60.: time_curr_dest[k][h][a]=','.join([spend,str(curr_cell),str(dest_cell),str(IO_case),str(I_O),str(curr_loc),str(dest_loc),str(mtv)])
                else: del ind[k][h][-1]; del time_curr_dest[k][h][-1] # this is to fix some bug on the diaries where in one hour you get more than 60m
                                
            # i want to keep in an array (it is placed after the 24th hour) the valid cells for each individual. this will save time in the exp module
            valid_cells[k].extend((int(curr_cell),int(dest_cell),traj))
            if h==23: final=[i for i in map(int,list(Set(np.hstack(valid_cells[k])))) if i>0]; time_curr_dest[k].append(final) # i keep the unique cells
        
    f=open(samp_dir+'/trajectories_'+wk_str[weekend]+'_'+dom+stk_str[static]+'.dat','wb'); cPickle.dump(wr_trv_traj,f,-1); f.close(); 
    f=open(samp_dir+'/input_expo_'+wk_str[weekend]+stk_str[static]+'.dat','wb'); cPickle.dump(time_curr_dest,f,-1); f.close()
    if wk_str[weekend]=='wday': f=open(samp_dir+'/patches.dat','wb'); cPickle.dump(make_video,f,-1); f.close()
    
    return ind

def get_indoor_outdoor_act(act,age): # age here is the group defined in the profiles

    per=[48,13,39,54,16,9,11,8,14,7,16]; a=1; case=4; act_per_age=[0,90,77,53] # act_per_age are those that exercise a sport

    # specifically for market we want to know if it is outdoors. I use IdF statistics (Les cahiers de l'Enquete Globale de Transport-Les deplacements pour achats - Julliet 2006) 
    # e.g. that 14% of movements are in marches (open market mostly)
    if 'Market' in act: 
         out_p=[0.14]; out=selection(out_p,0)
         if out==0: case=5
             
    elif 'Fun' in act: # I base the selection on statistics of IdF (Note Rapide).
        
        # the act_per_age is the % to do sports per age group. Otherwise they exercise something indoors (e.g cinema)
        act_p=[act_per_age[age]/100.]; s=selection(act_p,0)
            
        if s==0: # if they do sports then we seek if it is outdoors
        
           # we scan to see which sport activity the individual exercises. The % in per array correspond to:
           # marche, gym, velo, natation, course a pied, foot, musculation, relaxation, sports collectif, danse, tennis, raquettes, combat, petanque, roller, golf, equitation, autres, ski, nautices, peche
           while a==1: 
                 for index_p,p in enumerate(per): act_p=[p/100.]; a=selection(act_p,0); which_act=index_p      
           if which_act not in [1,3,6]: case=5 # outdoor activities

    return case # returns case 5 (outdoors) or 4 (indoors)

def get_IO_ratio(case,IO_metro,curr_loc,dest_loc): # this is for the transport IO ratio
        
    IO=[]; curr_dep=int(int(curr_loc)/1000.); dest_dep=int(int(dest_loc)/1000.); tmp=''   

    if case==23: # RER                
        IO=[0.,0.,0.,0.] # O3-->waiting, travelling / PM2.5-->waiting, travelling
        if curr_dep==75: IO[0]=0.; IO[2]=IO_metro[1] # people waiting inside
        else: IO[0]=1.; IO[2]=1. # otherwise the people wait outside (IO=1)
        if curr_dep==75 and dest_dep==75: IO[1]=0.; tau=uniform(3.2,5.4); IO[3]=round(tau,1) # this is the travelling IO
        else: IO[1]=0.; tau=uniform(2.9,3.2); IO[3]=round(tau,1) 
        return '/'.join(map(str, IO)) 
    
    #case21: car; case22: bus; case23: RER; case24: metro; case25: tram; case5,31: foot; case32: bicycle; case33: motorbike  
    def cs5(): IO.append(1.); tau=uniform(1.,2.); IO.append(round(tau,1)); return IO # foot
    def cs6(): IO.append(1.); tau=uniform(2.3,3.4); IO.append(round(tau,1)); return IO # 2-wheels
    def cs2(): IO.append(0.); tau1=uniform(0.9,2.); tau2=uniform(0.9,2.1); tau3=uniform(0.9,3.3); IO.extend((round(tau1,1),round(tau2,1),round(tau3,1))); return IO # car (ozone, pm rural, pm periferique, pm paris)    
    
    # for PT i have 4 numbers: O3-->waiting, travelling / PM2.5-->waiting, travelling
    def cs3(): IO.extend((1.,0.)); tau=uniform(1.7,3.7); IO.extend((1.,round(tau,1))); return IO # bus    
    def cs4(): IO.extend((0.,0.)); tau=uniform(5.5,8.5); IO.extend((IO_metro[1],round(tau,1))); return IO # metro (average from L1 and L14)
    def cs7(): IO.extend((1.,0.)); tau=uniform(1.,1.9); IO.extend((1.,round(tau,1))); return IO # tram
                       
    options = {5:cs5,21:cs2,22:cs3,24:cs4,25:cs7,31:cs5,32:cs6,33:cs6}; options[case]()
        
    return '/'.join(map(str, IO)) # combines the IOs in a string separated by '/' e.g. 0./0.8/3.2/4.3 (half for ozone, half for PM2.5)
    
def check_for_movement(data_diaries,activity,h,per_act_pop,home_dep,ind,pep_mv,allow_check,k,weekend,nb_individuals,activities,static):
    
    trans_str=['Car','Foot','Bicycle','Bike','Bus','Metro','RER','tram']; pep_mv=[0,0,0]; scale=sum(per_act_pop)/nb_individuals; on_the_move=2
    new_act_dic={2:0,3:0,1:1,4:2,0:2}; new_act=new_act_dic[activity] # activity=2,3: actives, activity=1: school/student, activity=4: inactive, activity=0: 0-3 (i use the variation of inactives)
    
    # we seek whereas this active, inactive, pupil/student will be on the move or not. We have each hour the number of journeys for the 3 categories. 
    # when we calculate the possibility of movement in each hour we use the statistics on the number of people on the move in IdF. From these people we have to exclude
    # those that started their move in the previous hour because we want the new movements. To find those that started their movement previously i use the output of the 
    # model in the previous hour and scale it to the total population
    if allow_check==1: # it will enter here only for one individual in each hour
    
       # here i want to calculate the people that moved in the previous hour and they continue to move in the current hour
       for k in xrange(nb_individuals):  
           for a in xrange(len(ind[k][h-1])-7): # number of activities performed in the previous hour
               if any(x in ind[k][h-1][a+7] for x in trans_str): new_act_tmp=new_act_dic[activities[k]]; pep_mv[new_act_tmp]+=1                                   
       allow_check=0 # we don't want this to happen for every individual   
    
    if static<>1: p=max(data_diaries['Nb'][new_act][h]-pep_mv[new_act]*scale,0); perc=min(p/per_act_pop[new_act],0.99); on_the_move=selection([perc],1)        

    if activity==0 and h in [23,0,1,2,3,4,5,6]: on_the_move=2 # children 0-3 don't move between 23:00 and 7:00
    
    return on_the_move,pep_mv,allow_check

def get_dest_dist_trans(curr_com,data_diaries,transport_type,dest_in_prof,workplace,communes,home,deps,motive,home_dep,h,mtv,stops,flux_deps,flux_mode_pt):
        
    # CASE1: i know the destination commune e.g., the person is working or goes to school and we already know from INSEE where do they go or if the person returns home
    # work follows the RED allowed path, school the BLUE path
    if (motive in [0,2] and dest_in_prof>1) or motive==7: # the dest_in_prof is the destination commune in the profiles.ins
            
       dest_com=int(dest_in_prof)
       if motive==7: dest_com=int(home)
              
       trans=get_transport(data_diaries,home_dep,transport_type,motive,mtv,curr_com,dest_com)
       if trans==0: trans=get_PT(home,curr_com,dest_com,stops,communes,flux_mode_pt,h)
       
       dist_travel=get_distance(trans,communes,dest_com,curr_com)
       
    # CASE2: we have some information about the destination (only for the work motive and that it is an intra or inter-departement movement). This is the YELLOW path
    elif motive==0 and workplace>=2: # in this case on-foot/bicycle is not valid. this exception is already included in the profiles

       chg={0:0,1:1,2:2,3:3,4:3}; trans=get_transport(data_diaries,home_dep,transport_type,motive,mtv,curr_com,0); tr=chg[trans]; dest=1
            
       if workplace==3: flx_p=[flux_deps[home_dep][tr][i]/sum(flux_deps[home_dep][tr][2::]) for i in range(2,10)]; dest=selection(flx_p,2) # inter-department (workplace=3)
   
       dest_com,dist_travel=get_dest_com(communes,data_diaries,home,home_dep,dest,deps,trans,curr_com)
       if trans==0: trans=get_PT(home,curr_com,dest_com,stops,communes,flux_mode_pt,h) # now we know the destination commune so we can translate public transport to its modes
       
    # CASE3: absolutely no information about the destination. This is the GREEN path
    # selects a destination commune, Intra-communal (dest=0), Intra-department, Paris, 92, 93, 94, 77, 78, 91, 95 (dest=9)
    # even if a person is currently in some other than home department i still get the selection based on the home commune because the statistics
    # are in respect to the home commune regardless from where the person is located
    else:

       chg={0:0,1:1,2:2,3:3,4:3}; trans=get_transport(data_diaries,home_dep,transport_type,motive,mtv,curr_com,0); tr=chg[trans]
       
       if trans in [2,3]: dest=0 # on-foot it is intra-communal
       else: flx_p=[flux_deps[home_dep][tr][i]/sum(flux_deps[home_dep][tr]) for i in xrange(10)]; dest=selection(flx_p,0)                          
                  
       dest_com,dist_travel=get_dest_com(communes,data_diaries,home,home_dep,dest,deps,trans,curr_com)
       if trans==0: trans=get_PT(home,curr_com,dest_com,stops,communes,flux_mode_pt,h) # now we know the destination commune so we can translate public transport to its modes
       
    return dest_com,dist_travel,trans

def get_transport(data_diaries,home_dep,transport_type,motive,mtv,curr_com,dest_com):
   
    deps=[75,92,93,94,77,78,91,95]; transform={2:2,3:3,4:1,5:0}; motive_str=['Work','Prof','School','Market','Fun','Personal'] # Work (motive=0), Proff, School, Market, Fun, Personal (motive=5) 

    #if a person returns at home from e.g. work then i have to add him in the movements of the "work" motive
    if motive==7:
       try: motive=motive_str.index(mtv)
       except: motive=0 # covers some rare cases with an unexplained error
        
    if motive==0: # for work we also know the transport type
       
       # tranport_type<>0 is in the case that the person goes to work and we know the means of transport from the profiles.ins file. But in the profiles.ins the trans id is different from the one here.                
       trans=transform[int(transport_type)] # transport_type is the element in the profiles.ins
       
       # here we select if it is bike or bicycle. the selection is based on the fact that 41% of the movements are bicycles and 59% bikes (statistics for IdF)
       if trans==3: tmp=[0.41,0.59]; trans=selection(tmp,3)
              
    # all other motives but work   
    # Enquete globale transport: Motorisation et usage de la voiture en Ile-de-France (Janvier 2013). These are the distributions by motive. Only for the shopping (market) i use the report titled
    # Les cahiers de l'Enquete Globale de Transport-Les deplacements pour achats (Juliet 2006)
    # Work,Proff,School,Market,Fun,Personal (inside each element: 0 public, 1 car, 2 foot, 3 bicycle, 4 motorbike)        
    else: trans_p=[float(x)/100. for x in data_diaries['transport'][motive][deps.index(int(home_dep))].split(',')]; tmp=sum(trans_p); trans_p.insert(0,round(1.-tmp,2)); trans=selection(trans_p,0)            
       
    if dest_com>0 and curr_com<>dest_com and trans in [2,3]: # to cover cases were foot,bicycle is chosen but the person is going to a different commune. i choose car for those cases
       trans_p=[float(x)/100. for x in data_diaries['transport'][motive][deps.index(int(home_dep))].split(',')]; tmp=sum(trans_p); trans_p.insert(0,round(1.-tmp,2))      
       tmp=sum(trans_p[2::])/sum(trans_p[0:2]); trans_p=[trans_p[0]*(1.+tmp),trans_p[1]*(1.+tmp),0.,0.,0.]; trans=selection(trans_p,0)
        
    return trans

def get_PT(home,curr_com,dest_com,stops,communes,flux_mode_pt,h): # public transport but which one? changes trans=0 to: 5 bus, 6 metro, 7 RER, 8 tram
        
    # the first case is if the person will go by PT based on INSEE but there is not stop it the commune. I assume a bus because anyway i know 
    # that the map of the stops is incomplete for the outer communes of IdF.
    if str(curr_com) not in stops.keys() or str(dest_com) not in stops.keys(): trans=5
    else:    
       # even if a person is currently in some other than home department i still get the selection based on the home commune because the statistics
       # are in respect to the home commune regardless from where the person is located in the current hour
       home_id=communes[int(home)][0]-1; dest_id=communes[int(dest_com)][0]-1 # get the id of the element of the home commune (0 paris, 1 coeur, 2 central, 3 autres)

       # i reformulate for tram, now it is bus,metro,RER,tram, this is based on the fact that 92.7% of the movements are with metro and 7.3% tram (statistics for IdF)
       data=[i for i in flux_mode_pt[dest_id][home_id]]; data.append(data[2]); data[3]=0.073*data[1]; data[1]=0.927*data[1]                              
       pub_p=[data[t]/sum(data) for t in xrange(4)]; pub_p=check_for_valid_pt_modes(curr_com,dest_com,stops,pub_p,h); trans=selection(pub_p,5)    

    return trans
       
def check_for_valid_pt_modes(curr_com,dest_com,stops,p_p,h):

    tmp2=[]; tmp1=[]; tmp=[]; modes=['Bus','metro','rer','tram']; pub_p=[0.,0.,0.,0.]
        
    for tr in stops[str(curr_com)]:       
        if tr not in stops[str(dest_com)]: tmp.append(modes.index(tr)) # if e.g. metro is not valid in tmp it will include the index 2

    if int(curr_com)==int(dest_com): tmp.append(2) # we cannot use RER to go to the same commune it is too short distance       
    if h in xrange(1,5): tmp.extend((1,2,3)) # between 1:00 and 5:00 the metro, RER, tram are not operational       
    
    tmp = list(set(tmp)) # keep the unique items    
    for index_p,p in enumerate(p_p): 
        if index_p not in tmp: tmp1.append(p) # in tmp1 we keep the % of the valid modes
        else: tmp2.append(p) # in tmp2 we keep the % of the invalid modes

    for index_p,p in enumerate(p_p): 
        if index_p not in tmp: pub_p[index_p]=p+(p/sum(tmp1))*sum(tmp2)
    
    if sum(pub_p)==0.: pub_p=[1.,0.,0.,0.] # if nothing is valid i choose bus because i know that the bus network is anyway incomplete (rare cases mostly in the GC)
        
    return pub_p    
    
def get_dest_com(communes,data_diaries,home,home_dep,dest,deps,trans,curr_com): 
        
    curr_dep=int(int(curr_com)/1000); tolerance=[3,15,15,15,6,6,6,15] # how many tries it will make until it increases the distance for the search loop

    if dest==0: dest_com=curr_com; dist_travel=get_distance(trans,communes,dest_com,curr_com) # will stay in the same commune                
    
    # if the person has to go outside the current commune but in the same department by bicycle i assume travelling to one of the adjacent communes (stored in communes[curr_com][4])
    elif dest==1 and trans==3: dest_com=choice(communes[curr_com][4]); dist_travel=get_distance(trans,communes,dest_com,curr_com)                  
    else:
    
        # to select a destination commune we will be based on the statistics on the per movement distance of the residents in each department 
        tries=0.; tmp1=[]; indx=deps.index(int(home_dep)); avg_dst=data_diaries['distance'][trans][indx] # the distance travelled per movement. I input it here as the maximum distance that the person can travel.    

        if dest==1: dest_dep=curr_dep # same department, dst_tol is used below see comments           
        else: dest_dep=deps[dest-2] # same other department (will give one of the: 92, 93, 94, 77, 78, 91, 95)
        
        for key in [i for i in communes.keys() if int(int(i)/1000.)==int(dest_dep) and i<>int(curr_com)]: tmp1.append(int(key)) # creates tmp with all the communes of the destination department
        shuffle(tmp1)

        for dest_com in tmp1: # search for a community that is inside those limits, i assume equal possibility for a person to go to each one of the communes inside the department
            dist_travel=get_distance(trans,communes,dest_com,curr_com)	
            if 0. < dist_travel <= avg_dst: break
            else:
               tries+=1
               if tries==tolerance[indx]: avg_dst*=1.1; tries=0 # if the person goes far away then the default average distance will never get a hit, so i gradually increase it.
               
    return dest_com,dist_travel

def get_distance(trans,communes,dest_com,curr_com):

    rng={0:[1.9,2.],1:[1.9,2.],2:[0.2,0.4],3:[0.8,2.],4:[1.9,2.],5:[1.9,2.],6:[1.9,2.],7:[1.9,2.],8:[1.9,2.]} # 1 car, 2 foot, 3 bicycle, 4 motorbike, 5 bus, 6 metro, 7 RER, 8 tram
    
    lon_home=float(communes[int(curr_com)][1]); lat_home=float(communes[int(curr_com)][2]); surf=float(communes[int(curr_com)][3]) # get the lat,lon,surface of the home commune
    lon_dest=float(communes[int(dest_com)][1]); lat_dest=float(communes[int(dest_com)][2]) # get the lat,lon of the destination commune
    
    # here we have 2 cases:
    # 1) the destination commune is other than the current location. In that case i get the distance based on the centroid coordinates of the start/end commune
    if int(curr_com)<>int(dest_com):

       degrees_to_radians = math.pi/180.; phi1 = (90. - lat_home)*degrees_to_radians; phi2 = (90. - lat_dest)*degrees_to_radians # phi = 90 - latitude
       theta1 = lon_home*degrees_to_radians; theta2 = lon_dest*degrees_to_radians # theta = longitude
       cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + math.cos(phi1)*math.cos(phi2)); dist_travel = math.acos(cos)*6373. # this will export the distance in Km
       
    # 2) the travelling is within the commune. In this case i get the surface of the commune, calculate the radius (spherical assumption)
    else: r=math.sqrt(surf/math.pi); dist_travel=uniform(rng[trans][0]*r,rng[trans][1]*r)       
              
    return dist_travel # its in Km
    
def travel_time(deps,data_diaries,trans,dist_travel,dest_com,previous_act_tr_time,curr_com,speeds2grid,crs1,crs2,grid_polyg,dx,dy):
                
    cc,lon1,lat1,ds,lon2,lat2=int(crs1[0]),crs1[1],crs1[2],int(crs2[0]),crs2[1],crs2[2]; c=0; cells_passed=0; time_travel=0.; curr_dep=int(int(curr_com)/1000.); indx_curr=deps.index(int(curr_dep))
    what_to_search=[cc+1,cc-1,cc-dx,cc+1-dx,cc-1-dx,cc+dx,cc+1+dx,cc-1+dx]; travel_traj=[cc]; j=int(math.sqrt(dx**2+dy**2))/2; trajectory=LineString([[lon1,lat1],[lon2,lat2]])
    try: speed=data_diaries['speed'][trans-1][indx_curr]
    except: speed=data_diaries['speed'][0][0] # some bugged rare cases
    
    if cc==ds: # 1 car, 2 foot, 3 bicycle, 4 motorbike, 5 bus, 6 metro, 7 RER, 8 tram
       if trans not in [1,4,5]: # other than car,bus,motorbike i get the speed from the diaries_input.xls file
          dest_dep=int(int(dest_com)/1000.); indx_dest=deps.index(int(dest_dep)); speed_dest=data_diaries['speed'][trans-1][indx_dest]; speed_av=(speed+speed_dest)/2.
          time_travel+=int(round((dist_travel/speed_av)*60.,2)) # i make the crude assumption that the moving speed id the average of the two speeds
       else: # bus,car,motorbike
         if cc in speeds2grid.keys(): s_p=[float(s)/sum(speeds2grid[cc][0]) for s in speeds2grid[cc][0]]; id_selection=selection(s_p,0); speed=speeds2grid[cc][1][id_selection]                                
         time_travel+=int(round((dist_travel/speed)*60.,2))
    
    else: # different start/end cell
    
       # after the starting cell has been added i can continue to fill the path until i reach the ending cell. The c<j is to avoid the endless loop 
       # in rare cases where the des_cell is not added in the travel_traj array
       while ds not in travel_traj and c<j: # travel traj cells start from 1 NOT 0
             c+=1             
             for cell in what_to_search: # what_to_search is the id of the grid cells shorted (starting from 1)
                 # when i find an intersecting polygon then i reset the array of cells to search for the next loop. the new array is the neighbouring cells. 
                 # the break statement breaks the for loop and goes back to the while loop		              
                 if cell not in travel_traj and cell in grid_polyg.keys() and trajectory.intersects(grid_polyg[cell][0]):
                    travel_traj.append(cell); what_to_search=[cell+1,cell-1,cell-dx,cell+1-dx,cell-1-dx,cell+dx,cell+1+dx,cell-1+dx]; break
       if travel_traj[-1]<>ds: travel_traj.append(ds)
  
       # here we could do a simple division to calculate the time since the speed is the same in every cell but we want to know how many cells are passed in this hour
       for c in travel_traj:            
           if trans in [1,4,5]: # car,bus,motobike
              # there are some cells that belong to the path but there is not a road segment in there so there is not speed to select. i use the method of average department speed
              if c in speeds2grid.keys(): s_p=[float(s)/sum(speeds2grid[c][0]) for s in speeds2grid[c][0]]; id_selection=selection(s_p,0); speed=speeds2grid[c][1][id_selection] 

           dist_per_cell=dist_travel/float(len(travel_traj)); time_travel+=round((dist_per_cell/speed)*60.,2)
           if time_travel<=(60.-previous_act_tr_time+1):               
              cells_passed+=1
              if cells_passed==len(travel_traj): break # with this i calculate how many cells will be passed in the travelling in the current hour
             
       time_travel=int(time_travel)           
       if time_travel+previous_act_tr_time>=60.: cells_passed=max(min(cells_passed,len(travel_traj)-1),1) # to be sure that at least one cell will be at the current hour ot left for the next hour of travelling       

    if trans in [5,6,7,8]: time_travel+=4 # for the PT modes i assume a 4 min waiting in the dock or bus stop
    time_travel=max(time_travel,5); time_travel=min(time_travel,90) # minimum of 5 mins and maximum of 90 mins
    if trans==2: time_travel=min(time_travel,30) # on foot maximum of 30 mins walking (arbitrary selection)
       
    tm_tot=previous_act_tr_time+time_travel # previous_act_tr_time is the time already spend by an activity or travelling  
    if Decimal(float((tm_tot))/60.)._isinteger() or Decimal(float((tm_tot))/59.)._isinteger(): time_travel=time_travel-2            
       
    return time_travel,cells_passed,travel_traj
         
def write_ind(ind,k,h,trans,time_travel,tmACT_spend,motive,home,dest_com,rem,dist_travel):

    def fill_pr_hr(tr,ac,ts,mt,ds,mv,dst,k,h,ind): ind[k][h][0]=int(tr); ind[k][h][1]=ac; ind[k][h][2]=ts; ind[k][h][3]=mt; ind[k][h][4]=ds; ind[k][h][5]=mv; ind[k][h][6]=round(dst,1)
    def convert_time_to_string(tm):
        if tm==0: final='00'
        elif tm<10: final='0'+str(tm)
        else: final=str(tm)    
        return final
    
    s='.59';s1='.00 '; trans_str=['Home','Car','Foot','Bicycle','Bike','Bus','Metro','RER','tram']; motive_str=['Work','Prof','School','Market','Fun','Personal','at_Home','to_Home']
    dest_com=str(int(dest_com)); tm_ACT_prev=0; mot=str(motive)
    
    # if rem=1 we are in the case where this subroutine just writes the remaining times of the on-going task
    if rem==1: # in this case tmACT_spend,time_travel are the residuals from the previous task
                      
       if int(time_travel)==0:
          if int(tmACT_spend)>60: c=4
          elif int(tmACT_spend)<=60: c=5         
       else: # if time_travel<>0 it means that there is residual travelling time from the previous hour. After the travelling the activity will continue                        
          if int(time_travel)>60: c=1
          elif int(time_travel)<60:
             c=3
             if motive==7: c=10 # the person travelled and returned home
             
       tmTRV_str=convert_time_to_string(int(time_travel-1)); tmTRV_str1=convert_time_to_string(int(time_travel)); tmACT_str=convert_time_to_string(int(tmACT_spend-1)) # make the strings for the files
       num_moves=ind[k][h-1][5] # if the current hour is full with something then we just repait the moves from the previous hour

    else: # in this case tmACT_spend, time_travel are the new times calculated for the new task. This is where it goes in h=0 which can only start with cases 0,1,3
              
       add_minute=-1
       if h>0: tm_ACT_prev=ind[k][h-1][1]
       if tm_ACT_prev==0.: add_minute=0 # this time is what is already written in the ind for this hour and it is a number <60
                 
       if int(time_travel)==0 and int(tmACT_spend)==0: c=0 # the only case in rem=0 when there is no travelling time
       else: # if time_travel<>0 it means that there is travelling to be done. After the travelling some activity will continue
          if int(tm_ACT_prev)+int(time_travel)>60: c=1
          elif int(tm_ACT_prev)+int(time_travel)<60:
             c=3 # the person travelled and continues an activity in the next hour
             if motive==7: c=10 # the person travelled and returned home
                        
       if h==0:
          if c==0: num_moves=0
          else: num_moves=1
       else:             
          if c==0: num_moves=ind[k][h-1][5]
          else: num_moves=ind[k][h-1][5]+1

       if h>0 and num_moves>1 and tm_ACT_prev<>0: i=convert_time_to_string(int(tm_ACT_prev)); s1='.'+i+' ' # if the prev activity ended e.g. at 3:51 then the travel will continue from 3:52
          
       tmTRV_str=convert_time_to_string(int(tm_ACT_prev-1+time_travel+add_minute)); tmTRV_str1=convert_time_to_string(int(tm_ACT_prev-1+time_travel+add_minute+1))

    def case0(): # in this case the individual is at home, this is valid for all hours
        hr=str(h)+s1+str(h)+s; wrt=hr+' '+str(int(home))+' '+'at_Home'; ind[k][h].append(wrt)
        fill_pr_hr(0,0,trans_str[trans],motive_str[motive],str(int(home)),num_moves,dist_travel,k,h,ind)

    def case1(): # in this case the travelling will extend the next hour, this is valid for all hours 
        hr=str(h)+s1+str(h)+s; wrt=hr+' '+'T('+mot+')'+' '+trans_str[trans]; ind[k][h].append(wrt) # write the travelling (it will take the whole hour)
        tmTRV_rem=int(tm_ACT_prev+time_travel)-60; tmACT_rem=int(tmACT_spend); fill_pr_hr(tmTRV_rem,tmACT_rem,trans_str[trans],motive_str[motive],dest_com,num_moves,dist_travel,k,h,ind)
                                             
    def case3(): # in this case the travelling ended in the current hour but the activity will continue in the next hour, this is valid for all hours    
        hr=str(h)+s1+str(h)+'.'+tmTRV_str; wrt=hr+' '+'T('+mot+')'+' '+trans_str[trans]; fill_pr_hr(0,0,trans_str[trans],motive_str[motive],dest_com,num_moves,dist_travel,k,h,ind); ind[k][h].append(wrt) # write the travelling        
        hr=str(h)+'.'+tmTRV_str1+' '+str(h)+s; wrt=hr+' '+dest_com+' '+motive_str[motive]; tmACT_rem=int(tmACT_spend)-(60-int(tm_ACT_prev+time_travel)); ind[k][h][1]=tmACT_rem; ind[k][h].append(wrt)

    def case4(): # in this case we have only activity which extends to the next hour, this is valid only for h>0        
        hr=str(h)+s1+str(h)+s; wrt=hr+' '+dest_com+' '+motive_str[motive]; ind[k][h].append(wrt)
        tmACT_rem=int(tmACT_spend)-60; fill_pr_hr(0.,tmACT_rem,trans_str[trans],motive_str[motive],dest_com,num_moves,dist_travel,k,h,ind)

    def case5(): # in this case we have only activity which ends on the current hour, this is valid only for h>0
        hr=str(h)+s1+str(h)+'.'+tmACT_str; wrt=hr+' '+dest_com+' '+motive_str[motive]; ind[k][h].append(wrt)
        fill_pr_hr(0,0,trans_str[trans],motive_str[motive],dest_com,num_moves,dist_travel,k,h,ind)

    def case10(): # in this case the travelling ended in the current hour and the person has returned home.
        hr=str(h)+s1+str(h)+'.'+tmTRV_str; wrt=hr+' '+'T('+mot+')'+' '+trans_str[trans]; fill_pr_hr(0,0,trans_str[trans],motive_str[6],dest_com,num_moves,dist_travel,k,h,ind); ind[k][h].append(wrt) # write the travelling
        hr=str(h)+'.'+tmTRV_str1+' '+str(h)+s; wrt=hr+' '+dest_com+' '+motive_str[6]; tmACT_rem=0.; ind[k][h][1]=tmACT_rem; ind[k][h].append(wrt)                        
      
    options = {0: case0, 1: case1, 3: case3, 4: case4, 5: case5, 10: case10}; options[c]()
      
    return c
    
def write_diaries(ind,pp,weekend,samp_dir,write_asc_dia):
 
    if weekend==0: print 'writing output...'
    trans_str=['Car','Foot','Bicycle','Bike','Bus','Metro','RER','tram']; wk_string=['weekday','weekend']

    if write_asc_dia==1:
    
       if not os.path.isdir(samp_dir+'/diaries/'): os.makedirs(samp_dir+'/diaries')
       if not os.path.isdir(samp_dir+'/diaries/'+wk_string[weekend]): os.mkdir(samp_dir+'/diaries/'+wk_string[weekend])

       # write the diary files for each individual
       print ''; print 'ascii files will be written for each individual...'	    
       for k in xrange(len(pp)):
           f = open(samp_dir+'/diaries/'+wk_string[weekend]+'/Individual_no'+str(k+1)+'.txt', "w")        
           num_moves=ind[k][23][5] # how many moves has the person performed in the day
           head='age='+str(int(pp[k][2]))+' home='+str(int(pp[k][0]))+' activity='+str(int(pp[k][3]))+' moves='+str(int(num_moves))+'/'+str(int(pp[k][-1])); print >> f,head # print some header with basic profile info        
           for h in xrange(24): 
               acts=len(ind[k][h])-7 # how many difference "activities" does the person have  
               for a in xrange(acts): 
                   if ind[k][h][6]<>0. and any(t in ind[k][h][a+7] for t in trans_str): print >> f,ind[k][h][a+7]+' ('+str(round(ind[k][h][6],2))+')'
                   else: print >> f,ind[k][h][a+7]
           f.close()
   
def get_statistics(ind,data_diaries,plt_stats,weekend,samp_dir,pp): 
        	
    nb=len(pp); ii=0.; jj=0.; p_above_5=0; cnt=[[0 for _ in range(nb)] for _ in range(6)]; tt=0; deps=[75,92,93,94,91,95,77,78]; school=0; work=0
    parisians=0; p=0; dmp=0; tot_time=np.zeros((6,8)); tot_dst=dict(); tot_mov=dict(); moves=np.zeros((24,6,8))
    motive_str=['Work','Prof','School','Market','Fun','Personal']; wk_string=['weekday','weekend']; mot_str=['wrk','pro','sco','mar','fun','per']
    obs_mov=['18.7','14.1','15.2','15.7','16.0','20.4']; tot_mov_ab3=0; p_above_3=0; trans_str=['Car','Foot','Bicycle','Bike','Bus','Metro','RER','tram']; time_act=np.zeros((7))  
    avg_dst=['6.1','0.4','2.0','6.1','-','-','-','-']; avg_tim=['23','12','19','22','-','-','-','-']      
    obs_trv={75:[['16','15','5','4'],['13','51','3','3'],['1','57','2','1'],['4','79','2','1'],['10','46','4','2'],['8','51','3','2']], \
             92:[['37','11','4','4'],['33','41','2','3'],['11','56','0','1'],['26','63','2','1'],['31','44','2','1'],['32','42','3','1']], \
             93:[['38','10','2','2'],['43','32','1','1'],['8','67','1','0'],['35','48','1','1'],['37','40','1','1'],['34','32','2','3']], \
             94:[['41','11','1','4'],['40','36','0','3'],['10','60','1','1'],['36','54','1','1'],['36','42','1','2'],['37','35','1','1']], \
             77:[['68','6','1','2'],['67','22','0','1'],['27','44','0','0'],['68','28','1','0'],['55','38','1','1'],['52','30','1','1']], \
             78:[['58','6','2','3'],['55','28','0','3'],['22','48','4','0'],['58','35','0','1'],['54','34','3','0'],['60','26','0','1']], \
             91:[['61','7','1','3'],['66','25','0','1'],['25','43','2','0'],['62','30','1','1'],['55','34','2','1'],['64','19','1','1']], \
             95:[['54','8','2','2'],['55','28','1','2'],['17','56','2','0'],['54','39','2','0'],['47','42','2','1'],['49','32','0','0']]}
        
    for d in deps: tot_dst[d]=np.zeros((6,8)); tot_mov[d]=np.zeros((6,8))
    
    print ''; print 'calculating statistics ('+wk_string[weekend]+')'; print ''; all_mot=motive_str+['at_Home']

    for k in [i for i in xrange(nb) if pp[i][3]==1]: school+=1
    for k in [i for i in xrange(nb) if pp[i][3]==2]: work+=1
           
    for h in xrange(24):
        
        o=0; o1=0; o2=-1
        if h<=9: o=1; o1=2; o2=1
        
        for k in xrange(len(ind)):

            age=pp[k][2]; home=int(pp[k][0]); home_dep=int(int(home)/1000.) # here the age dict has the home sector in the 1st array element
                                                                                                      
            for a in xrange(len(ind[k][h])-7):
                
                diff=int(str(ind[k][h][a+7])[9-o1:11-o1])-int(str(ind[k][h][a+7])[3-o:5-o]) # time spend in activity

                if any(x in ind[k][h][a+7] for x in all_mot): found=[x for x in all_mot if x in ind[k][h][a+7]]; m=all_mot.index(found[0]); time_act[m]+=diff
                                   
                if any(x in ind[k][h][a+7] for x in trans_str):
                                                
                   found=[x for x in trans_str if x in ind[k][h][a+7]]; index_t=trans_str.index(found[0]); mtv=int(str(ind[k][h][a+7])[13-o2:14-o2]) # the motive is inside the T(...) sub-string
                                                                                 
                   if mtv==7: #if a person returns at home from e.g. work then i have to add him in the work motive movements
                   
                      if a==0: tmp=str(ind[k][h-1]) # search for the motive in the previous hour
                      else: tmp=str(ind[k][h]) # search for the motive in the previous activity                      
                      try: found=[x for x in motive_str if x in tmp]; mtv=motive_str.index(found[0]) 
                      except: tmp=str(ind[k][h-2]); found=[x for x in motive_str if x in tmp]; mtv=motive_str.index(found[0]) # search for the motive two hour before
                                      
                   cnt[mtv][k]+=1; moves[h,mtv,index_t]+=1; tot_time[mtv,index_t]+=diff; tot_dst[home_dep][mtv,index_t]+=ind[k][h][6]; tot_mov[home_dep][mtv,index_t]+=1
                   if age>=1: tot_mov_ab3+=1
                               
    to_div=[work,work,school,nb,nb,nb,nb]    
    for m in range(7): print 'average time spend for '+all_mot[m]+': '+str(round(time_act[m]/(float(to_div[m])*60.),1))+'h'    
    print ''    
        
    for m in range(6): print 'people in the sample moving at least once for '+motive_str[m]+': '+str(len([i for i in cnt[m] if i>0])); tt+=sum(cnt[m])
    print ''    
    
    for m in range(6): print 'number of movements for '+motive_str[m]+': '+str(sum(cnt[m]))+' ('+str(round(100.*sum(cnt[m])/tt,1))+'% - '+obs_mov[m]+'%)'
    print 'total movements: '+str(tt)
    
    print ''; dst=np.zeros((len(deps))); mov=np.zeros((len(deps)))
    for in_d,d in enumerate(deps): mov[in_d]=np.sum(tot_mov[d][:,:]); dst[in_d]=np.sum(tot_dst[d][:,:]); print 'distance travelled per movement in '+str(d)+'='+str(round(dst[in_d]/mov[in_d],1))+' Km'

    ii=sum(dst[1:4]); jj=sum(mov[1:4]); print 'average distance travelled per movement in Petit Couronne='+str(round(ii/jj,1))+' Km (3.6 km)'
    ii=sum(dst[4::]); jj=sum(mov[4::]); print 'average distance travelled per movement in Grand Couronne='+str(round(ii/jj,1))+' Km (5.7 km)'
    print 'average distance travelled per movement in IdF='+str(round(sum(dst[:])/sum(mov),1))+' Km (4.5 km)'
    
    print ''
    for k in [i for i in pp.keys() if pp[i][2]>=1]: p_above_3+=1 # count the people aged>=4
    print 'average number of movements per person aged above 3 in IdF='+str(round(tot_mov_ab3/float(p_above_3),2))+' (3.5)'; print ''
                
    print '                              '+' Car '+'        Foot '+'       Bicycle '+'  Bike ' 
    for d in obs_trv.keys():
        for m in xrange(6):
            mm=''; base='% of movements in '+str(d)+' by '+mot_str[m]+'='; tt=sum(tot_mov[d][m,:])
            for t in range(4): mm+=str(round(100.*tot_mov[d][m,t]/tt,1))+'% ('+obs_trv[d][m][t]+'%)'+' '
            print base+mm
        print ''

    tim_IdF=np.sum(tot_time, axis=0); dst_IdF=np.zeros((8)); mov_IdF=tot_mov[75]+tot_mov[92]+tot_mov[93]+tot_mov[94]+tot_mov[91]+tot_mov[95]+tot_mov[77]+tot_mov[78]
    for t in xrange(8): dst_IdF[t]+=sum(tot_dst[d][:,t])        
        
        
    for t in range(8): tt=sum(mov_IdF[:,t]); print 'average time/distance spend per movement in '+trans_str[t]+'='+str(round(tim_IdF[t]/tt,1))+' min ('+avg_tim[t]+' min) / '+str(round(dst_IdF[t]/tt,1))+' Km ('+avg_dst[t]+' km)'
    print 'average daily time spend for movements (all modes)'+'='+str(round(sum(tim_IdF)/p_above_3,1))+' minutes'
    
    # count the persons that had not moved at all during the day (i count only the age>=5 because the statistics are for those only)        
    for k in [i for i in xrange(len(ind)) if ind[i][23][5]==0 and pp[i][2]>=1]: p+=1
    for k in [i for i in xrange(len(ind)) if ind[i][23][5]==0 and pp[i][2]>=1 and int(pp[i][0]/1000.)==75]: dmp+=1
    for k in [i for i in pp.keys() if pp[i][2]>=1 and int(pp[i][0]/1000.)==75]: parisians+=1
      
    print ''; print 'people above 3 not moving for the whole duration of the day='+str(round(float(p)/float(nb),3)*100.)+'% (6.9%)'
    print 'Parisians above 3 not moving for the whole duration of the day='+str(round(float(dmp)/float(parisians),3)*100.)+'% (5.1%)'
	
    if weekend==0 and plt_stats==1: plot_stats_diaries(moves,nb,samp_dir,data_diaries)

def selection(perc_list,write_id):

    sum_p=0.; id_selection=0; sorted_list=[0.]           
    
    for index_p,p in enumerate(perc_list): sum_p+=p; sorted_list.append(sum_p) #in index_p=0, sum_p is zero and it adds this as the 1st element of sorted_list
        
    if len(sorted_list)>2 and 0.999 < sum(sorted_list) < 1.001: sorted_list[len(sorted_list)-1]=1. # Some times the sum gives a slightly lower than 1 number (e.g. 0.9999999998)
    elif sorted_list[len(sorted_list)-1]<>1.: sorted_list.append(1.) # if the array has only 0. and one % value then the above if will mistakenly remove the % and place an 1.
        
    tau = uniform(0,1)                    
    for index_p,p in enumerate(sorted_list[0:len(sorted_list)-1]): # The last item is not iterated its p=1.0           
        if tau>=p and tau<sorted_list[index_p+1]: id_selection=index_p+write_id; break
        
    return id_selection