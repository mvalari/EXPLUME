import os,cPickle,sys,math,matplotlib
import pylab as plt
import numpy as np
from expo import get_progress_bar
from sets import Set
from dateutil.easter import *
from scipy.signal import butter,lfilter
from scipy.optimize import curve_fit
from string import lowercase
    
def images_for_video(lines,proj,samp_dir,procc):

    import matplotlib.patches as patches
    from multiprocessing import Process as prll
	
    try: os.makedirs(samp_dir+'/trajectories/')
    except: pass

    def split_for_multi(total):    
        procc=10; chunks=[]; chunk=int((len(total)/float(procc))//1)                    
        for index_d,d in enumerate(total[0::chunk]): chunks.append([]); chunks[index_d]=total[index_d*chunk:(index_d+1)*chunk]                
        if len(chunks)>procc:
           diff=len(chunks)-procc; rem=np.hstack(chunks[-diff::])  
           for i in range(diff): del chunks[-1]        
           chunks.append(rem.tolist())        
        return chunks    
    
    def video(pat,step,chunk):
        fig=plt.figure(figsize=(10,10))
        for indx in chunk:	
            ax=fig.add_subplot(111); proj.readshapefile('data/communes_polyline','idf'); proj.readshapefile('data/Departments_region','dp',linewidth=1.5,color='k',antialiased=1)        
            for p in pat[0:indx*step+1]: ax.add_patch(p)    
            plt.savefig(samp_dir+'/trajectories/trj'+str(indx*step)+'.jpg',dpi=100)        
    
    print len(lines)
    
    sys.exit(0)
    
    threads=[]; step=30; print ''; print 'saving trajectory images...'; pat=[patches.Polygon(i,lw=0.8,color='r') for i in lines]; chunks=split_for_multi(range(int(len(pat)/step)))
        
    for chunk in chunks: p=prll(target=video,args=[pat,step,chunk]); threads.append(p)
    for thread in threads: thread.start()
    for thread in threads: thread.join()
            
def plot_SIREN(data,conc_period):
              
    s_str={0:[1,2,12],1:[3,4,5],2:[6,7,8],3:[9,10,11],4:range(1,13)}; max_tp=[5,3,3]; IO_f='IO'; seas_str={1:'_DJF',2:'_MAM',3:'_JJA',4:'_SON',0:'_annual'}; pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']
    if conc_period>2008: IO_f='IO_'+str(conc_period)
    
    IO=np.zeros((2,5,3,5,21)) # pols,seasons,BT,ageB,bins
    for indx,s in enumerate([4,0,1,2,3]): # i transform the distribution from monthly to seasonal (and the annual)
        for m in s_str[s]: 
            for d in [75,77,78,91,92,93,94,95]: IO[:,indx,:,:,:]+=data[d][:,m-1,:,:,:]
            
    for indx_pol,pol in enumerate(['O3','PM25']):    
        for indx_bl,bl in enumerate(['Houses','Offices','Schools']):
            fig=plt.figure(figsize=(180,30*max_tp[indx_bl])); plt.subplots_adjust(wspace=0.4,hspace=0.4)
        
            for s in range(5):                
                for tp in range(max_tp[indx_bl]):
        
                    ax=fig.add_subplot(max_tp[indx_bl],5,5*tp+s+1); ax.set_xlabel('I/O ratio ',fontsize=120,labelpad=50); c=np.linspace(0,1,21); dt=IO[indx_pol,s,indx_bl,tp]
                    freq=dt/sum(dt); ax.axis([0,1.,0,1.1*max(freq)]); xtickPos=c[0::2]; plt.setp(ax,xticks=xtickPos); ax.set_xticklabels(c[0::2],size=70,rotation=-45,ha='left')
                    plt.tick_params(axis='y', which='major', labelsize=70, pad=25); plt.tick_params(axis='x', which='major', color='white', width=4, length=5, pad=20)                
                    plt.suptitle(pol_str[indx_pol]+' I/O ratios frequency distributions', fontsize=250); ax.bar(left=c,height=freq,width=0.05,color='Chartreuse', alpha=0.5,linewidth=2)
                    ax.set_ylabel('normalized frequency',fontsize=120,labelpad=50)
                    
                    bins=np.linspace(0.,0.05*len(dt),len(dt)); tmp1=dt.astype(int); tmp=np.repeat(bins,tmp1); md=np.median(tmp); mn=np.mean(tmp)                    
                    ax.annotate('mean = '+str(round(mn,2)), xy=(0.02,0.97), xycoords='axes fraction', fontsize=110, ha='left',va='top')
                    
                    if tp==0: plt.title(seas_str[s][1::], fontsize=200, y=1.12)
                                                            
            fig.savefig('data/IO_SIREN/'+IO_f+'_'+pol+'_'+bl+'.pdf', dpi=300); plt.close(fig)
        
def plot_profiles(pp,samp_dir,building_stock):

    c=[[0,0,0,0,0],[0,0,0]]; cl=['g','r','b','yellow','orange']; lg=[['cat1','cat2','cat3','cat4','cat5'],['cat1','cat2','cat3']]; cl1=[[cl[1],cl[4],cl[3],cl[2],cl[0]],[cl[1],cl[3],cl[2]]]    
    
    # plot resident house
    for k in pp.keys(): ha=pp[k][8]-1; c[0][ha]+=1; oa=pp[k][9]-1; c[1][oa]+=1
	            
    ttl='% of population living'; plt.ioff(); fig=plt.figure(figsize=(45,20)); fig.suptitle(ttl, fontsize=150); matplotlib.rcParams['font.size']=60; plt.subplots_adjust(wspace=0.4, hspace=0.2)
    
    for indx_bl,bl in enumerate(['houses','offices']):
    
        data=np.array(c[indx_bl])/float(len(pp)); ax=fig.add_subplot(1,2,indx_bl+1); wedge=ax.pie(data, labels=lg[indx_bl], colors=cl1[indx_bl], autopct='%1.1f%%', shadow=True)				
        for w in wedge[0]: w.set_edgecolor('white')                       
                    					
    plt.ion(); fig.savefig(samp_dir+'/building_age_'+str(building_stock)+'.pdf', dpi=100); plt.close(fig)
                                    
def plot_stats_diaries(moves,nb_individuals,samp_dir,data_diaries): # moves contains number of movements per hour (24), motive (6) and transport type (8)
    
    pop=[2234,1313,1408,1208,1562,1516,1319,1169]; obs_tot=[60,30,15,15,15,90,285,1000,1700,800,650,700,710,720,650,750,1150,1500,1350,1000,550,225,150,130]; moves*=sum(pop[:])*1000./nb_individuals
    motive_str=['commuting','professional','school','shopping','recreation','personal']; trans_str=['car','on-foot','public transport']  
    
    try: os.makedirs('samples'); os.makedirs(samp_dir)
    except: pass 
                                                                 
    plt.ioff(); fig=plt.figure(figsize=(50,30))       
    ax = fig.add_subplot(111); plt.tick_params(labelsize=60); dtm=np.sum(moves, axis=(1,2)); dto=np.asarray(obs_tot)*1000.; ax.set_xlabel('hour',fontsize=60); ax.set_ylabel('people moving',fontsize=60)
    ax.grid(True,linestyle='-',color='0.75'); t=np.array([i+1 for i in xrange(24)]); ax.set_xlim(1,24); ax.plot(t,dtm,'r-o',lw=4.,label='model'); ax.plot(t,dto,'b-o',lw=4.,label='obs')            
    ax.set_ylim(0.,max(max(dtm),max(dto))*1.15); plt.legend(loc='best',prop={'size':71},labelspacing=0.5,shadow='true'); plt.ion()        
    fig.savefig(samp_dir+'/diurnal.pdf', dpi=300); plt.close(fig)
    
    # plot per motive        
    plt.ioff(); fig=plt.figure(figsize=(30,20)); fig.suptitle('based on a sample of '+str(nb_individuals)+' individuals', fontsize=50)   
    for m in [0,1,2,3,4,5]:
        panel=231+m; ax = fig.add_subplot(panel); ax.set_xlabel('hour',fontsize=20); ax.set_ylabel('people in movement (in thousands)',fontsize=20); plt.title(motive_str[m], fontsize=30)
        dtm=np.sum(moves,axis=2)[:,m]; dto=np.asarray(data_diaries['obs_mot'][m])*1000.
        ax.grid(True,linestyle='-',color='0.75'); t = np.array([i+1 for i in xrange(24)]); ax.plot(t,dtm,'r-o',lw=1.5,label='model'); ax.plot(t,dto,'b-o',lw=1.5,label='obs')
        ax.set_xlim(1,24); plt.tick_params(labelsize=20); plt.legend(loc='best',prop={'size':23},labelspacing=0.5,shadow='true'); plt.ion()
    fig.savefig(samp_dir+'/diurnal_per_motive.pdf', dpi=300); plt.close(fig)   
        
    # plot per transport type                    
    plt.ioff(); fig=plt.figure(figsize=(30,10)); plt.suptitle('based on a sample of '+str(nb_individuals)+' individuals', fontsize=50)
    for t in [0,1,2]:
        panel=131+t; ax = fig.add_subplot(panel); ax.set_xlabel('hour',fontsize=20); ax.set_ylabel('people in movement (in thousands)',fontsize=20); plt.title(trans_str[t], fontsize=30)
        dtm=np.sum(moves,axis=1)[:,m]; dto=np.asarray(data_diaries['obs_trans'][t])*1000.
        ax.grid(True,linestyle='-',color='0.75'); t = np.array([i+1 for i in xrange(24)]); ax.plot(t,dtm,'r-o',lw=1.5,label='model'); ax.plot(t,dto,'b-o',lw=1.5,label='obs')
        ax.set_xlim(1,24); plt.tick_params(labelsize=20); plt.legend(loc='best',prop={'size':23},labelspacing=0.5,shadow='true'); plt.ion()
    fig.savefig(samp_dir+'/diurnal_per_transport.pdf', dpi=300); plt.close(fig)
	
def plot_stats_expo(samp_dir,out_dir,path_to_plots,dom,period_stats,procc,path_conc,building_stock,chim_label,static,run_label):

    import fnmatch
    from expo import get_grid_dims,split_for_multi
    from os import path,listdir,makedirs
    from datetime import datetime,timedelta
        	            
    st=open('domains/'+dom+'/grid_polygons.dat','rb'); grid_polyg=cPickle.load(st); st.close(); surf=(grid_polyg[1][0]).area; res=int(math.sqrt(surf)*100); pp_file=samp_dir+'/profiles.dat'
    st='domains/'+dom+'/COORDS_'+dom+'.dat'; grid_coords=cPickle.load(open(st,'rb')); dx,dy=get_grid_dims(grid_coords); bin_for_REF=0.1; print ''; print 'plotting...'; bld_npg=''
    if int(building_stock)>2008: pp_file=samp_dir+'/profiles_'+str(building_stock)+'.dat'
    
    st=open(pp_file,'rb'); pp=cPickle.load(st); st.close(); nb=len(pp)
    
    if static==1: run_label+='_STATIC'; out_dir+=run_label 
    if int(building_stock)>2008: bld_npg='_BLD'+str(building_stock)
    plot_dir=path_to_plots+'/ind_'+str(nb)+'_'+run_label; plot_REF=path_to_plots+'/ind_'+str(nb)+'_REF'
    for D in ['path_to_plots','plot_dir','plot_REF']: exec("if not path.isdir(%s): makedirs(%s)" %(D,D))
        
    # open files
    f=open(out_dir+'/exp_flux.dat','rb'); epf=cPickle.load(f); f.close(); f=open(out_dir+'/exp_per_ind.dat','rb'); epi=cPickle.load(f); f.close()
    f=open(out_dir+'/exp_avg.dat','rb'); epd=cPickle.load(f); f.close(); f=open(out_dir+'/exp_per_bin.dat','rb'); epb=cPickle.load(f); f.close()
    f=open(out_dir+'/exp_per_bin_per_day.dat','rb'); bpd=cPickle.load(f); f.close(); f=open(out_dir+'/exp_per_cell.dat','rb'); epc=cPickle.load(f); f.close()
    f=open(samp_dir+'/exp_per_bin_REF.dat','rb'); epbR=cPickle.load(f); f.close(); f=open(samp_dir+'/exp_per_bin_per_day_REF.dat','rb'); bpdR=cPickle.load(f); f.close()
    f=open(samp_dir+'/exp_avg_REF.dat','rb'); epdR=cPickle.load(f); f.close(); f=open(samp_dir+'/exp_per_ind_REF.dat','rb'); epiR=cPickle.load(f); f.close()
    
    # create the final list with the dates that the user has selected with period_stats. The valid date is taken from the file name not the conc files
    exp_flist=sorted([f for f in listdir(out_dir) if fnmatch.fnmatch(f, 'exp_'+dom+'*.dat')],key=str.lower); flist_all=[]
    start_date=datetime.strptime(exp_flist[0][5+len(dom):15+len(dom)],'%Y-%m-%d'); end_date=datetime.strptime(exp_flist[-1][5+len(dom):15+len(dom)],'%Y-%m-%d'); delta=end_date-start_date    
    for index_d,d in enumerate(xrange(delta.days+1)): dt=str(start_date+timedelta(days=index_d))[0:10]; flist_all.append(dt)    
    if (period_stats[0] or period_stats[1]) not in flist_all: print 'the dates provided in the expo.par file are not in the simulation output in '+out_dir; sys.exit(1)
    
    print 'there are '+str(len(flist_all))+' exposure day(s) in the '+out_dir+' directory'; flist=flist_all[flist_all.index(period_stats[0]):flist_all.index(period_stats[1])+1]  
    print 'the user has selected a sub-period of '+str(len(flist))+' day(s) ('+flist[0]+' to '+flist[-1]+')'; dates_procc=split_for_multi(procc,flist); proj=get_projection(dx,dy,grid_coords,res)     
    cf=path_conc+'/'+chim_label+'/conc_avg_'+flist[0]+'_'+flist[-1]+'.dat'; cday=path_conc+'/'+chim_label+'/conc_avg_per_day_'+flist[0]+'_'+flist[-1]+'.dat'; f=open(cf,'rb'); conc_avg=cPickle.load(f)
    
    # plot exposure per area of occurance
    #print ''; print 'plotting fluxes...'; plot_flux(epf,plot_dir)
        
    # export exposure per various groups
    #print 'export average stats...'; avg_stats(plot_dir,out_dir,nb,flist,bld_npg,samp_dir); export_time(epi,bld_npg,plot_dir,samp_dir,flist,nb); export_exp(epi,bld_npg,plot_dir,samp_dir,flist,nb)
    
    # box plots
    #print 'plotting boxplots...'; plot_box(epi,plot_dir,samp_dir,bld_npg); #plot_box(epiR,plot_REF,samp_dir,bld_npg)
        
    # time series
    #print 'plotting time series...'; plot_timeseries(plot_dir,flist,epd)

    # correlation
    #print 'plotting correlation...'; correlation(plot_dir,bpd,epd,epdR,epi,epiR,cday,flist,run_label)
          
    # low-pass time-series filter
    if len(flist)<1825: print 'your time-series consists of less than a total of 5 years, low pass filter will not be applied...'
    else: print 'applying low-pass filter...'; plot_filter(plot_dir,flist,bpd,epd,epi,0); plot_filter(plot_REF,flist,bpdR,epdR,epiR,1); plot_filter_ep(plot_dir,flist,epb,epbR,bpd,bpdR,cday,run_label)
                            
    # distributio/pie/PDF plots average exposure in the whole period
    print 'plotting period and seasonal average exposure distribution plots...'
    con_md_s=[[[] for i in range(5)],[[] for i in range(5)]]
    for s in [i for i in range(5) if np.hstack(epb[i]['tot'][0])<>[]]: con_md_s[0][s],con_md_s[1][s]=get_modes(epb[s]['tot'],nb) # get the modes per season and period                       
    mode_cntr=plot_distr(epb,plot_dir,nb,con_md_s,0,bin_for_REF)    
    for s in [i for i in range(5) if np.hstack(epb[i]['tot'][0])<>[]]: data=epb[s]; plot_distr_prof(s,plot_dir,nb,data,mode_cntr,0,bin_for_REF,epi); plot_cumul(s,plot_dir,nb,data); #plot_pie(s,plot_dir,'n/d',data,0,con_md_s,[0])
        
    # distribution plots for the REF case
    #con_md_s=[[[] for i in range(5)],[[] for i in range(5)]]
    #for s in [i for i in range(5) if epbR[i]['tot'][0]<>[]]: con_md_s[0][s],con_md_s[1][s]=get_modes(epbR[s]['tot'],nb) # get the modes per season and period    
    #mode_cntr=plot_distr(epbR,plot_REF,nb,con_md_s,1,bin_for_REF)     
    #for s in [i for i in range(5) if epbR[i]['tot'][0]<>[]]: data=epbR[s]; plot_distr_prof(s,plot_REF,nb,data,mode_cntr,1,bin_for_REF,epiR)
    
    # get modes for each day    
    try: f=open(out_dir+'/modes_per_day_'+flist[0]+'_'+flist[-1]+'.dat','rb'); con_md_d=cPickle.load(f); f.close()
    except:
       print 'calculating modes per day...'; con_md_d=[[[] for i in range(len(flist))] for _ in range(2)]; done=0; previous=0
       for index_day in range(len(flist)): con_md_d[0][index_day],con_md_d[1][index_day]=get_modes(bpd[0][index_day],nb); done+=1.; previous=get_progress_bar(done,len(flist),previous)
       f=open(out_dir+'/modes_per_day_'+flist[0]+'_'+flist[-1]+'.dat','wb'); cPickle.dump(con_md_d,f,-1); f.close()
        
    for index_pol,pol in enumerate('O3 PM25'.split()):
        print ''; c=np.array([0. for i in range(100)]); print pol; print '-------'
        for md in con_md_d[index_pol]: c[len(md)]+=1                                        
        for indx,per in enumerate(c):
           if per/sum(c)>0.02: print str(round(100.*per/sum(c),0))+'% of days have '+str(indx+1)+' mode(s)'
                      
    # maps of average period exposure    
    print ''; print 'plotting period average exposure maps...'    
    data,data1,data2=np.zeros((dx,dy,4)),np.zeros((dx,dy,4)),np.zeros((dx,dy,4)); X,Y=proj(epc[:,:,0],epc[:,:,1]); ephc=exp_to_home(epi,nb,dom,samp_dir,building_stock,dx,dy,epc[:,:,0:2])
    for s in [i for i in range(5) if np.sum(epc[:,:,2+i])<>0.]:
        for index_j,j in enumerate([0,1,2+s,7+s]): data[:,:,index_j]=epc[:,:,j]; data1[:,:,index_j]=conc_avg[:,:,j]; data2[:,:,index_j]=ephc[:,:,j]
        plot_map(s,'con',data1,cf,plot_dir,proj,X,Y); plot_map(s,'n/d',data,cf,plot_dir,proj,X,Y); plot_map(s,'home',data2,cf,plot_dir,proj,X,Y)
        
def exp_to_home(epi,nb,dom,samp_dir,building_stock,dx,dy,coords): # put the total exposure of each individual to their home cell
        
    from expo import get_i_j
    
    dt_p=[0,0]; dt_t=[0,0]; bld=samp_dir+'/profiles'; tmp_p,tmp_t,ephc=np.zeros((dx,dy,10)),np.zeros((dx,dy,10)),np.zeros((dx,dy,12)); np.seterr(invalid='ignore'); epi[np.isnan(epi)]=0.
        
    for p in [0,1]: dt_p[p]=epi[:,0,p,-1,:]+epi[:,1,p,-1,:]; dt_t[p]=epi[:,0,-1,-1,:]+epi[:,1,-1,-1,:] # sum weekday/weekend for exposure and time
               
    if int(building_stock)>2008: bld=samp_dir+'/profiles_'+str(building_stock)    
    st=open(bld+'.dat','rb'); pp=cPickle.load(st); st.close(); st=open('domains/'+dom+'/comm2cells.dat','rb'); corresp=cPickle.load(st); st.close() 
    st=open('domains/'+dom+'/idf_cells.dat','rb'); idf_cells=cPickle.load(st); st.close() 

    for k in xrange(nb):
        home=pp[k][0]; k_r=int((k/4)//1); home_cell=int(corresp[int(home)][k_r,0]); ii,jj=get_i_j(home_cell,dx)
        if home_cell not in idf_cells: continue
        tmp_p[ii-1,jj-1,0:5]+=dt_p[0][:,k]; tmp_p[ii-1,jj-1,5::]+=dt_p[1][:,k]; tmp_t[ii-1,jj-1,0:5]+=dt_t[0][:,k]; tmp_t[ii-1,jj-1,5::]+=dt_t[1][:,k]            
    
    ephc[:,:,2::]=tmp_p/tmp_t; ephc[np.isnan(ephc)]=0.; ephc[:,:,0:2]=coords; corresp.clear(); pp.clear(); return ephc

def plot_flux(data,plot_dir):
    
    # data here is exp_flux with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes (6), home area (3), exposure area (3)
    
    lbl=['work / school','total exposure']; ttl=['annual','DJF','MAM','JJA','SON']; pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; lg=['Paris','suburban','rural']; W=0.5   
    cl = ['y','r','gold']; dt=np.zeros((5,3,6,3,3)); dt=np.sum(data,axis=1); plt.ioff(); fig=plt.figure(figsize=(120,220)); plt.subplots_adjust(wspace=0.3,hspace=0.4)
    fl=[['P','S','R','P','S','R'],['S','P','S','S','P','S'],['R','R','P','R','R','P']]; ttl1=['Paris','suburbs','rural','Paris','suburbs','rural']; np.seterr(divide='ignore', invalid='ignore')
    
    dt1=np.zeros((5,3,7,3,3)); dt1[:,:,[0,1,2,3,4,5],:,:]=dt; dt1[:,:,-1,:,:]=np.sum(dt,axis=2); x=np.linspace(2,3*len(lbl),len(lbl)); c=np.hstack([[i-1.,i,i+1.] for i in x])
    new=np.zeros((5,2,3,3*len(lbl))); tot=np.sum(dt1,axis=-1); ss=['np.zeros((6))','new[s,indx_pol,0]','new[s,indx_pol,0]+new[s,indx_pol,1]']; new1=np.zeros((5,2,3,3*len(lbl)))
    dt2=np.zeros((5,2,7,3,3)); dt2[:,0,:,:,:]=dt1[:,0,:,:,:]/dt1[:,-1,:,:,:]*2.; dt2[:,1,:,:,:]=dt1[:,1,:,:,:]/dt1[:,-1,:,:,:]*2.    
    for i in range(3): dt1[:,:,:,i,0]/=tot[:,:,:,i]/100.; dt1[:,:,:,i,1]/=tot[:,:,:,i]/100.; dt1[:,:,:,i,2]/=tot[:,:,:,i]/100.
        
    for s in range(5):
        for ind in range(2):
            new[s,ind,0,0]=dt1[s,ind,1,0,0]; new[s,ind,0,1]=dt1[s,ind,1,1,1]; new[s,ind,0,2]=dt1[s,ind,1,2,2]; new[s,ind,0,3]=dt1[s,ind,6,0,0]; new[s,ind,0,4]=dt1[s,ind,6,1,1]; new[s,ind,0,5]=dt1[s,ind,6,2,2]
            new[s,ind,1,0]=dt1[s,ind,1,0,1]; new[s,ind,1,1]=dt1[s,ind,1,1,0]; new[s,ind,1,2]=dt1[s,ind,1,2,1]; new[s,ind,1,3]=dt1[s,ind,6,0,1]; new[s,ind,1,4]=dt1[s,ind,6,1,0]; new[s,ind,1,5]=dt1[s,ind,6,2,1]
            new[s,ind,2,0]=dt1[s,ind,1,0,2]; new[s,ind,2,1]=dt1[s,ind,1,1,2]; new[s,ind,2,2]=dt1[s,ind,1,2,0]; new[s,ind,2,3]=dt1[s,ind,6,0,2]; new[s,ind,2,4]=dt1[s,ind,6,1,2]; new[s,ind,2,5]=dt1[s,ind,6,2,0]            

            new1[s,ind,0,0]=dt2[s,ind,1,0,0]; new1[s,ind,0,1]=dt2[s,ind,1,1,1]; new1[s,ind,0,2]=dt2[s,ind,1,2,2]; new1[s,ind,0,3]=dt2[s,ind,6,0,0]; new1[s,ind,0,4]=dt2[s,ind,6,1,1]; new1[s,ind,0,5]=dt2[s,ind,6,2,2]
            new1[s,ind,1,0]=dt2[s,ind,1,0,1]; new1[s,ind,1,1]=dt2[s,ind,1,1,0]; new1[s,ind,1,2]=dt2[s,ind,1,2,1]; new1[s,ind,1,3]=dt2[s,ind,6,0,1]; new1[s,ind,1,4]=dt2[s,ind,6,1,0]; new1[s,ind,1,5]=dt2[s,ind,6,2,1]
            new1[s,ind,2,0]=dt2[s,ind,1,0,2]; new1[s,ind,2,1]=dt2[s,ind,1,1,2]; new1[s,ind,2,2]=dt2[s,ind,1,2,0]; new1[s,ind,2,3]=dt2[s,ind,6,0,2]; new1[s,ind,2,4]=dt2[s,ind,6,1,2]; new1[s,ind,2,5]=dt2[s,ind,6,2,0]            
            
    for indx_pol,pol in enumerate(['O3','PM25']):                                                      
        for s in range(5): # annual and the 4 seasons
            ax=fig.add_subplot(5,2,2*s+indx_pol+1); plt.title(ttl[s]+' '+pol_str[indx_pol],fontsize=170,y=1.1); ax.axis([0.5,3*len(lbl)+2,0,100]); ax.set_axisbelow(True)
            xtickNames=plt.setp(ax,xticklabels=ttl1); plt.setp(xtickNames, fontsize=110, y=-0.01, ha='center', rotation=-20); plt.setp(ax,xticks=c+W/2.)
            plt.tick_params(axis='y',labelsize=90,pad=35); ax.set_ylabel('% of exposure',fontsize=140,labelpad=40)
            
            for index_l,l in enumerate(lbl): ax.text(x[index_l]+W/2., 103, l, fontsize=90, ha='center', va='center', fontweight='bold')
            
            for i in range(3):
                exec("ax.bar(c, new[s,indx_pol,%s], W, color=cl[%s], bottom=%s, edgecolor='none'); bars = ax.patches" %(i,i,ss[i])); kk=-1; exec("V=%s" %(ss[i]))                                                                
                for bar,value in zip(bars,new[s,indx_pol,i]): 
                    kk+=1; H=bar.get_height(); X=bar.get_x()
                    if value>=5:                    
                       ax.text(X+W/2., V[kk]+value/2., str(int(round(value,0)))+'%', fontsize=70, ha='center',va='center')                                      
                       ax.text(X-0.1, V[kk]+value/2., fl[i][kk], fontsize=70, ha='center', va='center')
                       ax.text(X+W/2., V[kk]+value/2.-4, str(round(new1[s,indx_pol,i,kk],1)), fontsize=70, ha='center',va='center')
                          
    fig.savefig(plot_dir+'/exposure_flux.pdf', dpi=300); plt.close(fig)   
        
def correlation(plot_dir,bpd,epd,epdR,epi,epiR,cday,dates,run_label):

    from numpy import percentile as pctl
    from scipy.stats import pearsonr as PC
    from scipy.stats import spearmanr as SC

    # bpd here is exp_per_bin_per_day with dimensions: index group (7: 'tot','P','PC','GC','1974','2005','2012'), index day, pollutants (2), bins    
    # epd here is exp_avg with dimensions: days, pollutants (2), modes+total (7), conc (daily average, 25th, 75th)
    # epdR here is exp_avg_REF with dimensions: days, pollutants (2)  
    # epi here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
    # epiR here is exp_per_ind_REF with dimensions: season (5), pollutants,time (3), nb individuals
    # cday is: pollutants (2), department+IdF (9), days
                    
    season={1:0,2:0,12:0,3:1,4:1,5:1,6:2,7:2,8:2,9:3,10:3,11:3}; dts=np.array([season[int(day[5:7])] for day in dates]); indx=np.where(dts>=0)[0]    
    f=open(cday,'rb'); amb=cPickle.load(f); f.close(); np.seterr(divide='ignore', invalid='ignore'); lbl=['house','work/school','other indoor','outdoor','car','public transport','total']
    unit=[' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']; f=open(plot_dir+'/correlation.txt',"w"); epi=np.sum(epi,axis=1); grps=['P','PC','GC','1974','2005','2012']
    
    print >> f,'co-founders (O3 with PM2.5)'; print >> f,''; print >> f,'time-series:'    
    o3=np.array(amb[0][-1])[indx]; pm=np.array(amb[1][-1])[indx]; pr=round(PC(o3,pm)[0],2); sp=round(SC(o3,pm)[0],2); print >> f,'AMB',pr,sp    
    o3=epdR[:,0][indx]; pm=epdR[:,1][indx]; pr=round(PC(o3,pm)[0],2); sp=round(SC(o3,pm)[0],2); print >> f,'REF',pr,sp    
    o3=epd[:,0,-1,0][indx]; pm=epd[:,1,-1,0][indx]; pr=round(PC(o3,pm)[0],2); sp=round(SC(o3,pm)[0],2); print >> f,run_label,pr,sp; print >> f,''; print >> f,'intra-population (99th percentiles):'    
        
    t=epiR[0,-1]; o3=epiR[0,0]/t; pm=epiR[0,1]/t; pr=round(PC(o3,pm)[0],2); sp=round(SC(o3,pm)[0],2)    
    per=pctl(o3,99); pm=pm[np.where(o3>=per)[0]]; o3=o3[o3>=per]; pr1=round(PC(o3,pm)[0],2); sp1=round(SC(o3,pm)[0],2); print >> f,'REF',pr,sp,'('+str(pr1),str(sp1)+')'
    t=epi[0,-1,-1]; o3=epi[0,0,-1]/t; pm=epi[0,1,-1]/t; pr=round(PC(o3,pm)[0],2); sp=round(SC(o3,pm)[0],2)                        
    per=pctl(o3,99); pm=pm[np.where(o3>=per)[0]]; o3=o3[o3>=per]; pr1=round(PC(o3,pm)[0],2); sp1=round(SC(o3,pm)[0],2); print >> f,run_label,pr,sp,'('+str(pr1),str(sp1)+')'
        
    for indx_pol,pol in enumerate(['O3','PM25']):

        print >> f,''; print >> f, pol; print >> f,''; print >> f,'time-series of groups with ambient concentrations:'; amb_paris=amb[indx_pol][0]
        for g in range(6):
            expo=np.zeros((len(dates)))        
            for d in xrange(len(dates)): tmp=bpd[g+1][d][indx_pol]; exp=np.dot(np.linspace(0.1,0.1*len(tmp),len(tmp)),tmp)/sum(tmp); expo[d]=exp        
            print >> f,grps[g],round(PC(expo,amb_paris)[0],2),round(SC(expo,amb_paris)[0],2)  
            
        if epd.shape[0]==len(dates): 
           print >> f,''; print >> f,'time-series of MEs with ambient concentrations:'; conc=amb[indx_pol][-1]
           for l in range(7): expo=epd[:,indx_pol,l,0]; print >> f,lbl[l],round(PC(expo,conc)[0],2),round(SC(expo,conc)[0],2)                              
        else: print 'the dates inside the exp_avg.dat file are different from the ones selected in the expo.par file. skipping...'
        
        E_REF=epdR[:,indx_pol][indx]; E_sce=epd[:,indx_pol,-1,0][indx]; pr=round(PC(E_REF,E_sce)[0],2); sp=round(SC(E_REF,E_sce)[0],2); print >> f ,''; print >> f, 'time-series with REF:',pr,sp
        
        print >> f ,''; print >> f, 'intra-population of MEs with total:'; expo=[epi[0,indx_pol,l]/epi[0,-1,l] for l in range(7)]            
        for l in range(6): ex,tot=expo[l],expo[-1]; valid=~np.isnan(ex); ex=ex[valid]; tot=tot[valid]; pr=round(PC(ex,tot)[0],2); sp=round(SC(ex,tot)[0],2); print >> f,lbl[l]+':',pr,sp
                
        # correlation of total exposure between the current scenario and the REF
        E_sce=expo[-1]; t_REF=epiR[0,-1]; E_REF=epiR[0,indx_pol]/t_REF; valid=~np.isnan(E_sce); E_sce,E_REF=E_sce[valid],E_REF[valid]; pr=round(PC(E_REF,E_sce)[0],2); sp=round(SC(E_REF,E_sce)[0],2)               
        print >> f ,''; print >> f, 'intra-population of total with REF:',pr,sp
                
        #plt.ioff(); fig=plt.figure(figsize=(20,20)); plt.subplots_adjust(bottom=0.2, wspace=0.3,hspace=0.3); plt.tick_params(axis='both',which='major',labelsize=50,pad=35)
        #ax = fig.add_subplot(1,1,1); ax.grid(True,linestyle='-',color='0.75'); ax.plot(E_sce,E_REF,'bo',markersize=1.7, markeredgecolor='none')       
        #ax.annotate('p='+str(sp), xy=(0,0),xycoords='axes fraction',xytext=(0.98,0.98),textcoords='axes fraction',size=80,ha='right',va='top')        
        #ax.set_xlabel(run_label+' '+unit[indx_pol],fontsize=60,labelpad=20); ax.set_ylabel('REF '+unit[indx_pol],fontsize=60,labelpad=20)            
        #ax.set_xlim(min(E_sce),max(E_sce)); ax.set_ylim(min(E_REF),max(E_REF)); fig.savefig(plot_dir+'/exposure_correl_'+pol+'.pdf', dpi=300); plt.close(fig)
        
def export_time(data,bld_npg,plot_dir,samp_dir,dates,nb):

    from datetime import datetime
    from datetime import timedelta as TD

    # data here is exp_per_ind with dimensions: season (annual+4seasons), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
 
    basic={'age':['<4','4-24','25-64','>64'],'res':['P','PC','GC'],'ac':['infants','education','working','unemployed','retired']}; npg=cPickle.load(open(samp_dir+'/npg_wday'+bld_npg+'.dat','rb')); years=[]
    f=open(plot_dir+'/time.txt',"w"); wrt=str(len(dates))+' days ('+dates[0]+' to '+dates[-1]+')'; print >> f,wrt; print >> f,'time spend in MEs'; wk=['weekdays','weekends','total']
    ttl=['group','house','work / school','other_indoor','outdoor','car','public','total']; title = '{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}'.format(*ttl); c=0; holidays=[]
    
    for day in [i for i in dates if datetime.strptime(i,'%Y-%m-%d').strftime("%A") in ['samedi','dimanche','Saturday','Sunday']]: c+=1 # count weekend days   
    for day in [i for i in dates if i[5::] in ['01-01','05-01','05-08','07-14','08-15','11-01','11-11','12-25']]: c+=1 # add bank holidays as weekends
    years=[int(day[0:4]) for day in dates]; years=Set(years) # get the unique years within the dates
    for y in years: east=easter(y)+TD(days=1); d=TD(days=1); holidays.append(str(east+d)); d=TD(days=38); holidays.append(str(east+d)); d=TD(days=49); holidays.append(str(east+d))        
    for day in [i for i in dates if i in holidays]: c+=1 # add Easter Monday, Ascension and Pentecote as weekends    
    tot_days=[len(dates)-c,c,len(dates)] # weekday days, weekend days, all days
        
    for w in range(3):
        try: new=data[:,w,:,:,:]
        except: new=data[:,0,:,:,:]+data[:,1,:,:,:]
        print >> f,''; print >> f,'-------------------------'+wk[w]+'-------------------------'; print >> f,''; print >> f,title
        
        for group in basic.keys():        
            print >> f,''; print >> f,group
            for sub_grp in basic[group]: 
                indx=[i for i,x in enumerate(npg[group]) if x==sub_grp]; wrt=[sub_grp]                          
                for m in range(7): new1=new[0,-1,m]; tmp=str(round(sum(new1[indx])/(60.*tot_days[w]*len(indx)),2)); wrt.append(tmp)        
                line = '{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}'.format(*wrt); print >> f,line
                 
def export_exp(data,bld_npg,plot_dir,samp_dir,dates,nb):

    # data here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
        
    grps={0:[['infants','education','working','unemployed','retired'],-1,'total'],1:[['cat1','cat2','cat3','cat4','cat5'],0,'home'],2:[['cat1','cat2','cat3'],1,'work'],3:[['cat1','cat2','cat3'],2,'other_ind']}
    groups=['ac','hage_cat','oage_cat','oage_cat']; s_str=['ANN','DJF','MAM','JJA','SON']; npg=cPickle.load(open(samp_dir+'/npg_wday'+bld_npg+'.dat','rb')); new=data[:,0,:,:,:]+data[:,1,:,:,:]
    f=open(plot_dir+'/exposure.txt',"w"); wrt=str(len(dates))+' days ('+dates[0]+' to '+dates[-1]+')'; print >> f,wrt; wrt='population of '+str(nb)+' individuals'; print >> f,wrt; print >> f,''
                
    for indx_pol,pol in enumerate(['O3','PM25']): 
        print >> f,''; print >> f,pol
        
        for indx_grp,grp in enumerate(groups):                
            print >> f,''; ttl=''        
            for sub_grp in grps[indx_grp][0]: ttl+=sub_grp+'/'
            ttl=ttl[0:-1]; indx=[]; m=grps[indx_grp][1]; wrt=grp+' '+grps[indx_grp][2]+' ('+ttl+')'; print >> f,wrt; print >> f,''
            for sub_grp in grps[indx_grp][0]: indx.append([i for i,x in enumerate(npg[grp]) if x==sub_grp])                        
                        
            for s in range(5):
                dt=[]; new1=new[s,indx_pol,m]; base=s_str[s]+': '; wrt=''                                
                for i in range(len(grps[indx_grp][0])): tmp=sum(new1[indx[i]]); dt=int(tmp/(60.*1e6)); wrt+=str(dt)+' '                                        
                print >> f,base+wrt
                       
def avg_stats(plot_dir,out_dir,nb,dates,bld_npg,samp_dir): # write some general statistics

    # epd here is exp_avg with dimensions: days, pollutants (2), modes+total (7), conc (daily average, 25th, 75th)
    # epi here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
          
    f=open(out_dir+'/exp_per_dep.dat','rb'); epdep=cPickle.load(f); f.close(); f=open(out_dir+'/exp_per_ind.dat','rb'); epi=cPickle.load(f); f.close(); season=['annual','DJF','MAM','JJA','SON']
    f=open(out_dir+'/exp_avg.dat','rb'); epd=cPickle.load(f); f.close(); lbl=['house','work','other_ind','outdoor','car','public']; ttl=''
    npg=cPickle.load(open(samp_dir+'/npg_wday'+bld_npg+'.dat','rb')); new=epi[:,0,:,-1]+epi[:,1,:,-1] 
    for l in lbl: ttl+=l+'/'
    
    f=open(plot_dir+'/stats.txt',"w"); wrt=str(len(dates))+' days ('+dates[0]+' to '+dates[-1]+')'; print >> f,wrt; wrt='population of '+str(nb)+' individuals'; print >> f,wrt; ttl=ttl[0:-1]
    for indx,pol in enumerate('O3 PM25'.split()): 
        
        a=epd[:,indx,0,0]; print >> f,''; print >> f,pol; wrt='date with the max total exposure conc.: '+dates[np.where(a==max(a))[0][0]]+' ('+str(round(max(a),1))+')'; print >> f,wrt; print >> f,''        
        for s in range(5): wrt=season[s]+' average exposure: '+ str(round(np.sum(epi[s,:,indx,-1]/np.sum(epi[s,:,-1,-1])),1)); print >> f,wrt
        print >> f,''
        
        for indx_w,w in enumerate(['weekday','weekend']): 
            print >> f,w
            for s in range(5): 
                base= season[s]+' % of exposure in '+ttl+': '; wrt=''
                for i in range(6): wrt+=str(round(100.*np.sum(epi[s,indx_w,indx,i])/np.sum(epi[s,indx_w,indx,-1]),1))+'% '
                print >> f,base+wrt
            print >> f,''            
                
        for d in [75,'PC','GC']: 
            wrt='over '+str(d)+': '
            for s in range(5): wrt+=str(epdep[d][indx][s])+'  '
            print >> f,wrt

        print >> f,''; tmp=[epdep[92][indx][0],epdep[93][indx][0],epdep[94][indx][0]]; print >> f,'std for PC deps=',round(np.std(tmp),2)
        tmp=[epdep[77][indx][0],epdep[78][indx][0],epdep[91][indx][0],epdep[95][indx][0]]; print >> f,'std for GC deps=',round(np.std(tmp),2)
        tmp=[epdep[d][indx][0] for d in [91,92,93,94,95,77,78]]; print >> f,'std for PC,GC deps=',round(np.std(tmp),2); print >> f,''
                          
        for grp in ['P','PC','GC']:         
            valid=[i for i,x in enumerate(npg['res']) if x==grp]; wrt='for the inhabitants of '+grp+'  '
            for s in range(5): C=np.sum(new[s,indx][valid]); time=np.sum(new[s,-1][valid]); wrt+=str(round(C/time,1))+'  '                                                       
            print >> f,wrt
        
def get_modes(data,nb):

    from scipy import signal
    from scipy.signal import find_peaks_cwt as pks
        
    con_modes=[]; step=0.1
    for index_pol in [0,1]:
                
        dt=np.asarray(data[index_pol]); min_con,max_con=get_min_max_plot(dt,nb,0.001,step); tmp=np.asarray([i for i in dt[int(round(min_con/step,0)):int(round(max_con/step,0))+1]/nb])
        fit=[[],[]]; valid_peaks=[]; valid_lows=[]; find_peaks = pks(tmp, np.arange(1,5)); find_lows=signal.argrelextrema(tmp, np.less_equal,order=5); per=[tmp[i] for i in find_peaks]; todel=[]
                                                                                            
        for index_p,p in enumerate(per): # if there are 2 peaks one next to the other we keep only one
            if index_p>0 and p==per[index_p-1]: todel.append(index_p-1)
        peaks=[round(min_con+(i*step),1) for index_i,i in enumerate(find_peaks) if index_i not in todel]
                                     
        per=[tmp[i] for i in find_lows[0]]; todel=[]                                                                                               
        for index_p,p in enumerate(per): # if there are 2 dips one next to the other we keep only one
            if index_p>0 and p==per[index_p-1]: todel.append(index_p-1)               
        lows=[round(min_con+(i*step),1) for index_i,i in enumerate(find_lows[0]) if index_i not in todel] # contains the boundaries as well
            
        for index_l,l in enumerate(lows[0:-1]): # check that each peak found is between two dips
            pp=[]; pp1=[]
            for p in [i for i in peaks if l<i<lows[index_l+1]]: pp.append(p); pp1.append(dt[int(round(p/step,0))])                                
            if pp<>[]: valid_peaks.append(pp[pp1.index(max(pp1))]) # if between 2 dips we found 2 peaks we keep the highest one
                
        for index_p,p in enumerate(valid_peaks[0:-1]): # check that each low found is between two peaks
            ll=[]; ll1=[]
            for l in [i for i in lows if p<i<valid_peaks[index_p+1]]: ll.append(l); ll1.append(dt[int(round(l/step,0))])
            if ll<>[]: valid_lows.append(ll[ll1.index(min(ll1))]) # if between 2 peaks we found 2 dips we keep the lowest one
                
        todel=[]; valid_lows=[min_con]+valid_lows+[max_con] # add the boundaries which are always lows
                              
        for mode in range(len(valid_peaks)): # Gaussian fit of the modes
            start,end=int(round(valid_lows[mode]/step,0)),int(round(valid_lows[mode+1]/step,0)); x=[round(step*i+0.1,1) for i in range(start,end+1)]; y=(dt[start:end+1]).tolist()
            if mode>0: del x[0]; del y[0]
                                                
            x=np.asarray(x); y=np.asarray(y); mn=np.dot(x,y)/sum(y); std=np.sqrt(np.dot(y,(x-mn)**2)/sum(y)); 
            try: coeff,matrix=curve_fit(gau+ss,x,y,p0=[1.,mn,std]); fit[0].append(coeff[1]); fit[1].append(abs(coeff[2])) # 0 holds the means and 1 holds the std of the modes
            except: fit[0].append(mn); fit[1].append(abs(std)) # the fit fails (rare cases) so i use the calculated mean and std
                                    
        for mode in range(len(valid_peaks)-1): # test for bi-modality (http://www.jstor.org/stable/3087302)
            r=(fit[1][mode]/fit[1][mode+1])**2; S=np.sqrt(-2.+3.*r+3.*r**2.-2.*r**3.+2.*(1.-r+r**2.)**(1.5))/(r+np.sqrt(r))                        
            if abs(fit[0][mode]-fit[0][mode+1])<=S*(fit[1][mode]+fit[1][mode+1]) and len(valid_lows[1:-1])>1: todel.append(mode+1) # i don't merge if it is to leave 0 modes        
                       
        valid_lows=[i for index_i,i in enumerate(valid_lows) if index_i not in todel]; con_modes.append(valid_lows[1:-1])
            
    return con_modes[0],con_modes[1]

def plot_distr(data,plot_dir,nb,con_modes,case,bin_for_REF):

    # data here is exp_per_bin (or exp_per_bin_REF) with dimensions: season (5), groups string (18), pollutants (2)
            
    def gauss(x, *p): A, mu, sigma = p; return A*np.exp((-(x-mu)**2.)/(2.*sigma**2.))
    
    plt.ioff(); fig=plt.figure(figsize=(120,200)); plt.subplots_adjust(wspace=0.3,hspace=0.5); unit=[' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']; ttl=['annual','DJF','MAM','JJA','SON']; out_str=''
    pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; mode_cntr=[[[] for _ in range(5)] for _ in range(2)]; step=0.1
    if case==1: step=bin_for_REF
        
    for indx_pol,pol in enumerate(['O3','PM25']):    
        for s in [i for i in range(5) if np.hstack(data[i]['tot'][indx_pol])<>[]]: # annual and the 4 seasons
   
            stp=int(round(step/0.1,0)); tmp=data[s]['tot'][indx_pol]; new_data=np.array([np.sum(tmp[i:i+stp]) for i in range(0,len(tmp)+stp,stp)])
            min_con,max_con=get_min_max_plot(new_data,nb,0.001,step); dt=[(float(i)/float(nb))*100. for i in new_data[int(round(min_con/step,0)):int(round(max_con/step,0))+1]]    
            ax=fig.add_subplot(5,2,2*s+1+indx_pol); ax.set_ylabel('% population',fontsize=150,labelpad=100); ax.set_xlabel('concentration exposure '+unit[indx_pol],fontsize=150,labelpad=40)   
            c = [round(min_con+float(i)*step,1) for i in xrange(0,int(round((max_con-min_con)/step,0)+1))]; plt.title(ttl[s]+' '+pol_str[indx_pol], fontsize=200, y=1.05)            
            ax.axis([min(c), max(c)+step,0,1.1*max(dt)]); xtickPos=c[0::int(round(1./step,0))]; plt.setp(ax,xticks=xtickPos); ax.set_xticklabels(c[0::int(round(1./step,0))],size=70,rotation=-45,ha='left')
            plt.tick_params(axis='y', which='major', labelsize=100, pad=25); plt.tick_params(axis='x', which='major', color='white', width=4, length=5, pad=20)  
            ax.bar(left=c,height=dt,width=step,color='Chartreuse', alpha=0.5,linewidth=2)
                                    
            if len(con_modes[indx_pol][s]) in [1,2] and case==0: #fit Gaussian curves (not for the REF or IO1 cases)                      
               bounds=[min_con]+con_modes[indx_pol][s]+[max_con]
               if len(con_modes[indx_pol][s])==2: cl=['g','b','r']; mode_str={0:'low',1:'medium',2:'high'}
               elif len(con_modes[indx_pol][s])==3: cl=['g','b','r','purple']; mode_str={0:'low',1:'medium',2:'medium-high',3:'high'}
               else: cl=['g','r']; mode_str={0:'low',1:'high'}               
                              
               for m in range(len(bounds)-1):                                                            
                   x=np.arange(bounds[m],bounds[m+1]+step,step); start,end=int(round(bounds[m]/step,0)),int(round(bounds[m+1]/step,0)); y1=new_data[start:end+1]                                                         
                   if str(x[-1])<>str(bounds[m+1]): x=[i for i in x[0:-1]]                   
                   y=[(float(i)/float(nb))*100. for i in y1]; mn=np.dot(x,y)/sum(y); std=np.sqrt(np.dot(y,(x-mn)**2)/sum(y)); coeff,matrix=curve_fit(gauss,x,y,p0=[1.,mn,std])                   
                   x=np.arange(min_con-step,max_con+2*step,step); y_fit=gauss(x,coeff[0],coeff[1],coeff[2]); ax.plot(x,y_fit,color=cl[m],lw=10.,label=mode_str[m]); mode_cntr[indx_pol][s].append(coeff[1])               
                   ax.plot([coeff[1],coeff[1]],[0.,0.99*coeff[0]],color=cl[m],ls='-',lw=20.) # draw the mean of the mode
               
               plt.legend(loc='best',prop={'size':120},labelspacing=0.5,shadow=True,fancybox=True)
            
    fig.savefig(plot_dir+'/exposure_dist'+out_str+'.pdf', dpi=300); plt.close(fig); return mode_cntr
        
def plot_box(data,plot_dir,samp_dir,bld_npg):
    
    # data here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
    # epiR here is exp_per_ind_REF with dimensions: season (5), pollutants,time (3), nb individuals
    
    from matplotlib.patches import Polygon,Patch
    
    lbl=['house','work / school','other_indoor','outdoor','car','public','total']; ttl=['annual','DJF','MAM','JJA','SON']; unit=[' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']; yax_str=[]     
    npg=cPickle.load(open(samp_dir+'/npg_wday'+bld_npg+'.dat','rb')); basic_groups=['tot','ac','hage_cat','sens']; pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']
    groups={'tot':['tot'],'ac':['infants','education','working','unemployed','retired'],'hage_cat':['cat1','cat2','cat3','cat4','cat5'],'sens':['sens1','sens2']}    
    grp_lbls={'ac':['infants','education','working','unemployed','retired'],'hage_cat':['cat1','cat2','cat3','cat4','cat5'],'sens':['prof1','prof2','population']}
    np.seterr(divide='ignore', invalid='ignore'); new=[0,0]; boxColors = ['royalblue','c','m','r','darkkhaki','g','chocolate']
      
    for p in [0,1]:  
       try: new[p]=(data[:,0,p,:,:]+data[:,1,p,:,:])/(data[:,0,-1,:,:]+data[:,1,-1,:,:]) # sum weekday/weekend and get the exposure conc for each individual
       except: new[p]=data[:,p,:]/data[:,-1,:] # get the exposure conc for each individual in the REF case
       
    for indx_pol,pol in enumerate(['O3','PM25']):     
        plt.ioff(); fig=plt.figure(figsize=(220,40*len(groups))); plt.subplots_adjust(wspace=0.4,hspace=0.4)
    
        for indx_grp,grp in enumerate(basic_groups):
                        
            indx=[range(new[0].shape[-1])]; pos=range(1,new[0].shape[1]+1); wd=0.7; sz=70; ttl_y=1.04; p=range(new[0].shape[1]); Z=[-2.5,-1.7,-0.9,0.,0.9,1.7,2.5]
            l_max=7; lg=[]; st=l_max/2.; p0,p1,p2=0.5,2.,1.; mx=[60,120]
            if grp=='hage_cat': p=[0]; l_max=1; Z=[0.]; st=1.; p0,p1,p2=1.,1.,0. # we will use indoors for the building age group                                   
            if grp=='sens' or plot_dir[-3::]=='REF': p=[6]; l_max=1; Z=[0.]; st=1.; p0,p1,p2=1.,1.,0. # we will use the total 
            if plot_dir[-3::]=='REF': pos=[1]           
            lg=[Patch(color=boxColors[l], label=lbl[l]) for l in p]
                
            if grp<>'tot': # get the ids of individuals in each group in question             
               curr_np=npg[grp]; indx=[]; pos=[]; wd=0.3; sz=35; ttl_y=1.09
               for sub_grp in groups[grp]: indx.append([i for i,x in enumerate(curr_np) if x==sub_grp])
               if grp=='sens': indx.append(range(new[0].shape[-1])) # to add the exposure of the whole population
               if grp in ['hage_cat','sens']: sz=70; mx=[30,30]; lg=[Patch(color=boxColors[0], label=lbl[l]) for l in p]
               for i in range(0,len(indx)): pos+=[p0+st*(p1*i+p2)+j for j in Z] # pos controls the positions of boxes in the x-axis                                       
            
            for s in range(5): # annual and the 4 seasons
                dt=[]; ax=fig.add_subplot(len(groups),5,5*indx_grp+s+1); plt.title(ttl[s]+' '+pol_str[indx_pol],fontsize=170,y=ttl_y); medians = range(l_max*len(indx))                               
                for indx_sub_grp in range(len(indx)):                                 
                    for l in [i for i in p if plot_dir[-3::]<>'REF']: tmp=new[indx_pol][s,l]; tmp=tmp[indx[indx_sub_grp]]; tmp=tmp[~np.isnan(tmp)]; dt.append(tmp)
                    if plot_dir[-3::]=='REF': tmp=new[indx_pol][s]; tmp=tmp[indx[indx_sub_grp]]; tmp=tmp[~np.isnan(tmp)]; dt.append(tmp); mx=[50,30]
                                                                                    
                pl=ax.boxplot(dt,sym='+',positions=pos, widths=0.5-0.01*l_max, showfliers=False); plt.setp(pl['boxes'],color='black'); plt.setp(pl['whiskers'],color='black')
                ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5); ax.set_xlim(0.5,l_max*len(indx)+0.5); ax.set_axisbelow(True); ax.set_ylim(0, mx[indx_pol])
                xtickNames=plt.setp(ax,xticklabels=lbl); plt.setp(xtickNames, rotation=-35, fontsize=70, y=-0.01, ha='left'); plt.tick_params(axis='y',which='major',labelsize=100,pad=35)           
                ax.set_ylabel('concentration'+unit[indx_pol],fontsize=120); ax.yaxis.labelpad=90; top=mx[indx_pol]

                if grp<>'tot': 
                   xtickPos=[st+i for i in range(0,l_max*len(indx),l_max)]; plt.setp(ax,xticks=xtickPos); xtickNames=plt.setp(ax,xticklabels=grp_lbls[grp])                                      
                   plt.setp(xtickNames, rotation=-20, fontsize=100, y=-0.03, ha='center')                   
                   leg=plt.legend(handles=lg,prop={'size':35},labelspacing=0.15,fancybox=True, bbox_to_anchor=(0.,1.,1.,.1), loc=3, ncol=l_max, mode="expand", borderaxespad=0.)
                   frm=leg.get_frame(); frm.set_linewidth(5); frm.set_edgecolor('grey'); frm.set_facecolor('SeaGreen'); frm.set_alpha(alpha=0.3)
                
                for wsk in pl['whiskers']: wsk.set(linestyle='--',lw=7)            
                for cap in pl['caps']: cap.set(lw=8)

                for l in range(l_max*len(indx)):
                                                                          
                    y=l; box=pl['boxes'][l]; boxX = []; boxY = []; box.set(lw=6); X=dt[l]; X = X[~np.isinf(X)]
                    if grp<>'tot': y=l%l_max
                                    
                    for j in range(5): boxX.append(box.get_xdata()[j]); boxY.append(box.get_ydata()[j]); boxCoords = zip(boxX,boxY)      
                    boxPolygon = Polygon(boxCoords, facecolor=boxColors[y]); ax.add_patch(boxPolygon); med=pl['medians'][l]; medianX=[]; medianY=[]               
                    for j in range(2): medianX.append(med.get_xdata()[j]); medianY.append(med.get_ydata()[j]); plt.plot(medianX, medianY, 'k'); 

                    if X.tolist()==[] or np.sum(X)==0.: upperLabels = 'n/a'
                    else: plt.plot([np.average(med.get_xdata())], [np.average(X)], color='yellow', marker='*', markeredgecolor='k', markersize=60); upperLabels = str(round(medianY[0],1))
                    
                    ax.text(pos[l], top-(top*0.04), upperLabels, horizontalalignment='center', size=sz, weight='bold', color=boxColors[y], rotation=45 )
                        
        fig.savefig(plot_dir+'/exposure_box_'+pol+'.pdf', dpi=300); plt.close(fig) 
        
def plot_timeseries(plot_dir,dates,data):
  
    # data here is exp_avg with dimensions: days, pollutants (2), modes+total (7), conc (daily average, 25th, 75th)

    if data.shape[0]<>len(dates): print 'the dates inside the exp_avg.dat file are different from the ones selected in the expo.par file. skipping...'; return plot_dir
    
    plt.ioff(); fig=plt.figure(figsize=(60,40)); plt.subplots_adjust(left=None,bottom=None,right=None,top=None, wspace=0.3, hspace=0.3); unit=['','',' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']
    for index_pol,pol in enumerate(['O3','PM25']):      
        avg_exp=data[:,index_pol,-1,0]; min_err=data[:,index_pol,-1,1]; max_err=data[:,index_pol,-1,2] # percentiles           
        ax=fig.add_subplot(2,1,index_pol+1); ax.set_ylabel('concentration exposure'+unit[index_pol+2],fontsize=50)                          
        ax.errorbar(xrange(len(dates)),avg_exp,yerr=[min_err,max_err],linestyle="None",marker="None",color="green",mew=2); ax.plot(xrange(len(dates)),avg_exp,'r-',lw=2,label=pol)                         
        plt.xlim(-1,len(dates)); xtickPos=[i for i in range(0,len(dates)+1,int(len(dates)/10))]; plt.setp(ax,xticks=xtickPos);                
        ax.set_xticklabels([i for i in dates[0::(int(len(dates)/10))]],ha='left',size=30,rotation=-35); ax.yaxis.labelpad=40                
        plt.tick_params(axis='y', which='major', labelsize=50, pad=20); plt.legend(loc='best',prop={'size':80},labelspacing=0.5,shadow='true'); plt.ion()
           
    fig.savefig(plot_dir+'/exposure_tm.pdf', dpi=300); plt.close(fig) 
    
def plot_filter(plot_dir,dates,data,data1,data2,case):
 
    # data here is exp_per_bin_per_day with dimensions: index group (7: 'tot','P','PC','GC','1974','2005','2012'), index day, pollutants (2), bins
    # data here is exp_per_bin_per_day_REF with dimensions: groups (2: population, Parisians), index day, pollutants (2), bins
    # data1 here is the exp_avg with dimensions: days, pollutants (2), modes+total (7), conc (daily average, 25th, 75th)
    # data1 here is the exp_avg_REF with dimensions: days, pollutants (2)
    # data2 here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
    # data2 here is exp_per_ind_REF with dimensions: season (5), pollutants,time (3), nb individuals
        
    pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; unit=['','',' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']; order = 1; fl=365
    mode_str={0:'low exp.',1:'high exp.'}; line_str={0:'g-',1:'r-'}; line_fit_str={0:'g',1:'r'}; yax_str=['population fraction','concentration']; md_str=['low','high']
                
    for index_pol,pol in enumerate(['O3','PM25']):                              
                                              
        plt.ioff(); fig=plt.figure(figsize=(40,40)); plt.subplots_adjust(wspace=0.2,hspace=0.4); signal_freq = len(dates); low_freq = signal_freq/float(fl)
        
        if case==0: tot=data1[:,index_pol,-1,0]; t=np.sum(data2[0,:,-1,6,:]); CUT=np.sum(data2[0,:,index_pol,6,:])/t # scenario case
        else: tot=data1[:,index_pol]; t=np.sum(data2[0,-1,:]); CUT=np.sum(data2[0,index_pol,:])/t # REF case
                
        ydata=get_data_per_mode(data[0],CUT,index_pol,dates) # CUT limit represents the average annual exposure                
        
        for tp in [0,1]: # against population fraction and against concentration exposure.
            ax=fig.add_subplot(2,1,1+tp); ax.set_ylabel(yax_str[tp]+unit[index_pol+2*tp],fontsize=70); avg=[]; m_str=''; s_str=''; get_mx=[]
            plt.xlim(fl,len(dates)); xtickPos=[i for i in range(fl,len(dates)+1,int(len(dates)/10.))]; plt.setp(ax,xticks=xtickPos); ax.yaxis.labelpad=40                                
            ax.set_xticklabels(dates[fl::(int(len(dates)/10))],ha='left',size=30,rotation=-35); plt.tick_params(axis='y',which='major',labelsize=50,pad=20)
            ax.annotate(lowercase[tp:tp+1]+')', xy=(0,0), xycoords='axes fraction', fontsize=130, xytext=(-250, 1150), textcoords='offset points',ha='right',va='top')
            if tp==0: plt.title('annual '+pol_str[index_pol], fontsize=130, y=1.15)        
            for d in xtickPos[1:-1]: ax.plot([d,d],[0.,1000.],color='black',ls='--',lw=5.)
            if index_pol==1: # WHO limit line
               ax.plot([0.,100000.],[10.,10.],color='black',ls=':',lw=5.); ax.annotate('WHO limit value',xy=(0,0),xycoords='axes fraction',xytext=(3650,10.5),textcoords='data',size=30,ha='right',va='bottom')
               if tp==1: get_mx=[8.]
            
            for mode in xrange(len(ydata)): # elements 0,2,4... are the pop fraction and 1,3,5... the concentration (per mode)
                dt=np.array(ydata[mode][tp::2]); dt[np.isnan(dt)]=0.; ma=filter(dt, low_freq, signal_freq, order); avg.append(round(np.mean(ma[fl::]),2-tp))
                ax.plot(range(0,len(dates)),ma,line_fit_str[mode],lw=5,label=mode_str[mode]); m_str+=mode_str[mode]+'='+str(avg[mode])+'  '; get_mx.extend((max(ma[fl::]),min(ma[fl::])))
                if tp==1: s_str+='std_'+md_str[mode]+'='+str(round(np.std(dt),1))+' '            
            
            if tp==1: 
               ma1 = filter(np.asarray(tot), low_freq, signal_freq, order); ax.plot(range(0,len(dates)),ma1,'gold',linestyle='--',lw=5,label='total exp.')
               m_str+='total exp.='+str(round(np.mean(ma1[fl::]),1))+'  '; get_mx.extend((max(ma1[fl::]),min(ma1[fl::])))
               ax.annotate(s_str,xy=(0,0),xycoords='axes fraction',xytext=(0.99,0.01),textcoords='axes fraction',size=60,ha='right',va='bottom')          
                        
            plt.ylim(min(get_mx)*0.9,max(get_mx)*1.1); ax.annotate(m_str,xy=(0,0),xycoords='axes fraction',xytext=(15,970),textcoords='offset points',size=60)            
            plt.legend(loc='upper left',prop={'size':50},labelspacing=0.2,shadow='true',fancybox=True)
                      
        fig.savefig(plot_dir+'/exposure_ma_'+pol+'_annual.pdf', dpi=300); plt.close(fig)
        
def plot_filter_ep(plot_dir,dates,epb,epbR,bpd,bpdR,cday,run_label): # these plots are representative for Paris

    # exp_per_bin (or exp_per_bin_REF) with dimensions: season (5), groups string (18), pollutants (2) 
    # exp_per_bin_per_day with dimensions: index group (7: 'tot','P','PC','GC','1974','2005','2012'), index day, pollutants (2), bins
    # exp_per_bin_per_day_REF with dimensions: groups (2: population, Parisians), index day, pollutants (2), bins
    # cday is: pollutants (2), departments+Idf (9), days

    pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; unit=[' (ppb)',r' ($\mu$'+'g/m'+r'$^3$)']; order = 1; fl=365; signal_freq = len(dates); low_freq = signal_freq/float(fl)
    plt.ioff(); fig=plt.figure(figsize=(90,40)); plt.subplots_adjust(wspace=0.3, hspace=0.4); f=open(cday,'rb'); CA=cPickle.load(f); f.close(); label=['BASE','REF']
    if len(dates)<1825: print 'your time-series consists of less than a total of 5 years, low pass filter will not be applied...'; return
            
    for index_pol in [0,1]:
        
        for plot in [0,1]: # ambient with REF modes / ambient with scenario modes
                
            if plot==0: data=bpd[1]; tmp=epb[0]['P'][index_pol] # scenario
            else: data=bpdR[1]; tmp=epbR[0]['P'][index_pol] # REF case
          
            CUT=np.dot(np.linspace(0.1,0.1*len(tmp),len(tmp)),tmp)/sum(tmp); ydata=get_data_per_mode(data,CUT,index_pol,dates)
            dt1=np.array(ydata[0][1::2]); dt1[np.isnan(dt1)]=0.; dt2=np.array(ydata[1][1::2]); dt2[np.isnan(dt2)]=0.
            
            ax=fig.add_subplot(2,2,index_pol+2*plot+1); ax.set_ylabel('concentration '+unit[index_pol],fontsize=70); plt.xlim(fl,len(dates)); xtickPos=[i for i in range(fl,len(dates)+1,int(len(dates)/10.))]
            plt.setp(ax,xticks=xtickPos); ax.yaxis.labelpad=40; ax.set_xticklabels(dates[fl::(int(len(dates)/10))],ha='left',size=30,rotation=-35); plt.tick_params(axis='y',which='major',labelsize=50,pad=20)                    
            ax.annotate(lowercase[plot:plot+1]+')',xy=(0,0),xycoords='axes fraction',xytext=(-0.08,1.07),textcoords='axes fraction',size=130)             
            amb=CA[index_pol][0]; plt.title(label[plot]+' vs. '+'ambient ('+pol_str[index_pol]+')', fontsize=120, y=1.05)
            for d in xtickPos[1:-1]: ax.plot([d,d],[0.,1000.],color='black',ls='--',lw=5.)
                
            ma = filter(amb, low_freq, signal_freq, order); ax.plot(range(0,len(dates)),ma,'blue',linestyle='--',lw=4,label='ambient')
            ma1 = filter(dt1, low_freq, signal_freq, order); ax.plot(range(0,len(dates)),ma1,'green',linestyle='--',lw=4,label='low')
            ma2 = filter(dt2, low_freq, signal_freq, order); ax.plot(range(0,len(dates)),ma2,'red',linestyle='--',lw=4,label='high')                                                                                              
            get_mx=[max(ma[fl::]),min(ma[fl::]),max(ma1[fl::]),min(ma1[fl::]),max(ma2[fl::]),min(ma2[fl::])]
                                                
            plt.ylim(min(get_mx)*0.9,max(get_mx)*1.1); plt.legend(loc='upper left',prop={'size':50},labelspacing=0.2,shadow='true',fancybox=True)
                      
    fig.savefig(plot_dir+'/exposure_ma.pdf', dpi=300); plt.close(fig)
       
def plot_cumul(s,plot_dir,nb,data):

    # data here is exp_per_bin with dimensions: season (5), groups string (18), pollutants (2)

    from numpy.polynomial import polynomial as P
    
    pl=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; unit=['ppb',r'$\mu$'+'g/m'+r'$^3$']; seas_str={0:'_annual',1:'_DJF',2:'_MAM',3:'_JJA',4:'_SON'}; nr=1; pop_groups=['tot']; nc=2
    plt.ioff(); fig=plt.figure(figsize=(45*nc,40*nr)); plt.subplots_adjust(bottom=0.2, wspace=0.3, hspace=0.3); y=range(0,101); y1=[10.,25.,50.,75.,90.]; step=0.1
                   
    for index_pol,pol in enumerate('O3 PM25'.split()):
        for gr in pop_groups:
                                                                                
            dt=data[gr][index_pol]; min_c,max_c=get_min_max_plot(dt,nb,0.001,step); n_dt=dt[int(round(min_c/step,0)):int(round((max_c/step)+1,0))]; x=[]; x1=[]; tmp=[]
                                               
            for indx_b,b in enumerate(n_dt): # each 0.1 bin
                for k in range(int(b)): tmp.append(round(min_c+step/2.+indx_b*step,1)) # each individual
                        
            for pr in y: x.append(np.percentile(tmp,pr)) # percentiles per bins of 2
                                                           
            ax = fig.add_subplot(nr,nc,1+index_pol); ax.grid(True,linestyle='-',color='0.75'); plt.title(seas_str[s][1::]+' '+pl[index_pol],fontsize=180,y=1.04); ax.axis([min(x)+step,max(x)-step,0,100])
            ax.set_ylabel('percentile',fontsize=140,labelpad=40); ax.set_xlabel('concentration exposure ('+unit[index_pol]+')',fontsize=120,labelpad=50)            
            plt.tick_params(axis='y', which='major', labelsize=80, pad=30); plt.tick_params(axis='x', which='major', color='white', width=4, length=5, labelsize=80, pad=20)                           
            c = [round(min_c+step*float(i),1) for i in xrange(len(n_dt)+1)][0::int(round(1./step,0))]            
            xtickPos=c; plt.setp(ax,xticks=xtickPos); ax.set_xticklabels(c,size=60,rotation=-45,ha='left'); smooth=P.polyfit(x,y,deg=4); x_new = np.arange(min(x),max(x), 0.2); y_new=[]
                                         
            for x in x_new: y_new.append(smooth[-1]*x**4.+smooth[-2]*x**3.+smooth[-3]*x**2.+smooth[-4]*x+smooth[-5])                              
            ax.plot(x_new,y_new,linewidth=8,linestyle='--',color='black'); ii=[str(i) for i in x_new]                        

            for j in y1: ind=min(range(len(y_new)), key=lambda i: abs(y_new[i]-j)); x1.append(round(x_new[ind],1)) # get the the 10,25,50,75,90,95 percentiles                                                             
            for i in range(len(y1)): ax.scatter(x1[i],y_new[ii.index(str(x1[i]))],marker='o', c='b', s=500*(i+1), edgecolors='none', label=y1[i]) # plot the 10,25,50,75,90,95 percentiles
                        
    plt.ion(); fig.savefig(plot_dir+'/exposure_CD'+seas_str[s]+'.pdf', dpi=100); plt.close(fig)
            
def plot_distr_prof(s,plot_dir,nb,data,mode_cntr,case,bin_for_REF,data1):

    # data here is exp_per_bin (or exp_per_bin_REF) with dimensions: season (5), groups string (18), pollutants (2)
    # data1 here is exp_per_ind with dimensions: season (5), weekday/weekend (2), pollutants,time (3), modes+total (7), nb individuals
    # data1 here is exp_per_ind_REF with dimensions: season (5), pollutants,time (3), nb individuals

    from scipy.interpolate import spline
    
    unit=['ppb',r'$\mu$'+'g/m'+r'$^3$']; cl=['g','r','b','yellow','Chartreuse','LimeGreen','SpringGreen','SeaGreen','gold','orange','purple']; nc=2; step=0.1
    pop_groups=[['P','PC','GC'],['1974','2005','2012'],['working','unemployed+retired'],['<4','4-24','25-64','>64']]
    lg=[['urban','suburban','rural'],['pre-1974','1975-2005','post-2005'],['actives','unemployed+retired'],['<4','4-24','25-64','>64']]; nr=len(pop_groups)
    cl=[[cl[0],cl[2],cl[1]],[cl[1],cl[2],cl[0]],[cl[2],cl[1]],[cl[0],cl[2],cl[9],cl[-1],cl[1]]]; seas_str={0:'_annual',1:'_DJF',2:'_MAM',3:'_JJA',4:'_SON'}
    ttl1=['geographical region of residence','construction year of residence','professional activity','population age']; pl=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']; plt.ioff(); fig=plt.figure(figsize=(40*nc,22*nr))
    lbl_out='exposure_nb'; ttl='exposure distributions for population sub-groups'; plt.subplots_adjust(wspace=0.3, hspace=0.5); fig.suptitle(ttl, fontsize=150)
            
    for index_pol,pol in enumerate('O3 PM25'.split()): 
        stp=int(round(step/0.1,0)); tmp=data['tot'][index_pol]; new_data=np.array([np.sum(tmp[i:i+stp]) for i in range(0,len(tmp)+stp,stp)]); min_con,max_con=get_min_max_plot(new_data,nb,0.001,step)

        if case==1: mode_cntr=[[[] for _ in range(5)] for _ in range(2)]; step=bin_for_REF; avg=np.sum(data1[s,index_pol])/np.sum(data1[s,-1]); min_con,max_con=get_min_max_plot(new_data,nb,0.,step)
        else: avg=np.sum(data1[s,:,index_pol,-1])/np.sum(data1[s,:,-1,-1]) # the average exposure
        
        for indx_sg,sub_group in enumerate(pop_groups): # the basic groups
        		        
            l=[]; l_str=''; new_max=0.; indx=2*indx_sg+index_pol; p=[[] for _ in sub_group]             
            ax=fig.add_subplot(nr,nc,2*indx_sg+1+index_pol); ax.set_ylabel('population fraction',fontsize=85,labelpad=50); ax.set_xlabel('exposure concentration ('+unit[index_pol]+')',fontsize=70,labelpad=30)                               
            c = [round(min_con+step*float(i),1) for i in xrange(0,int(round((max_con-min_con)/step,0)+1))]; plt.title(ttl1[indx_sg], fontsize=100, y=1.05)
            xtickPos=c[0::int(round(1./step,0))]; plt.setp(ax,xticks=xtickPos); ax.set_xticklabels(c[0::int(round(1./step,0))],size=45,rotation=-45,ha='left'); plt.tick_params(axis='y',which='both', labelsize=45, pad=20)                               
            ax.grid(True, linestyle='-', which='both', color='lightgrey',alpha=0.5); plt.tick_params(which='both', color='black', width=2, length=12); c_new=np.linspace(min(c),max(c),300)
            if mode_cntr[index_pol][s]==[]: ax.plot([avg,avg],[0.,30.],color='black',ls='--',lw=10.) # if there are no modes draw the average of the distribution
            for m in mode_cntr[index_pol][s]: ax.plot([m,m],[0.,30.],color='black',ls='--',lw=10.) # draw the mean of the mode
            ax.annotate(lowercase[indx:indx+1]+')', xy=(0,0), xycoords='axes fraction', fontsize=130, xytext=(-300, 1050), textcoords='offset points',ha='right',va='top')
            if indx_sg==0: ax.annotate(pl[index_pol], xy=(0,0), xycoords='axes fraction', fontsize=130, xytext=(1000, 1200), textcoords='offset points',ha='center',va='center')
                        
            for index_g,g in enumerate(sub_group):
                        
                if '+' in g: pos=g.index('+'); tmp=data[g[0:pos]][index_pol]+data[g[pos+1::]][index_pol] # if i need to aggregate two groups
                else: tmp=data[g][index_pol]                       
            
                l_str+='p['+str(index_g)+'][0],'; new_data=np.array([np.sum(tmp[i:i+stp]) for i in range(0,len(tmp)+stp,stp)])            
                dt=[(float(i)/nb)*100. for i in new_data[int(round(min_con/step,0)):int(round(max_con/step,0))+1]]                                                                                                                                       
                smooth=spline(c,dt,c_new); p[index_g]=ax.plot(c_new,smooth,linewidth=7,linestyle='-',color=cl[indx_sg][index_g],alpha=0.8)
                if max(dt)>new_max: new_max=max(dt); ax.axis([min(c),max(c),0,1.1*new_max])                    
            
            l_str=l_str[0:len(l_str)-1]; exec("l=[%s]" %l_str); leg=plt.legend(l,lg[indx_sg],loc='best',prop={'size':65},labelspacing=0.2)
            frm=leg.get_frame(); frm.set_linewidth(5); frm.set_edgecolor('grey'); frm.set_facecolor('SeaGreen'); frm.set_alpha(alpha=0.3)
            
    plt.ion(); fig.savefig(plot_dir+'/'+lbl_out+seas_str[s]+'.pdf', dpi=100); plt.close(fig)
                          
def plot_pie(s,plot_dir,dates,data,indx_chunk,con_modes,dates_procc):

    # data here is exp_per_bin with dimensions: season (5), groups string (18), pollutants (2)
  
    cl=['g','r','b','yellow','Chartreuse','LimeGreen','SpringGreen','SeaGreen','gold','orange','purple']; pos={3:[0.78,0.5,0.21],2:[0.7,0.26]}; seas_str={0:'',1:'_DJF',2:'_MAM',3:'_JJA',4:'_SON'}
    pop_groups=[['P','PC','GC'],['P-->P','PC-->PC','GC-->GC','GC-->P','PC-->P'],['<4','4-24','25-64','>64'],['infants','education','working','unemployed','retired'],['1974','2005','2012']]                
    lbl_out='exposure_pie'; ttl='% of population exposed in the exposure modes per group type'; nc=len(pop_groups); l='s'; msg=seas_str[s][1::]+' period' 
        
    for index_pol,pol in enumerate('O3 PM25'.split()): 
        if con_modes[index_pol][s]==[]: continue
        n_modes=len(con_modes[index_pol][s])+1; md=['low','medium','high']; plt.ioff(); fig=plt.figure(figsize=(170,15*(n_modes+1))); fig.suptitle(ttl, fontsize=150); 
        matplotlib.rcParams['font.size']=70; plt.subplots_adjust(wspace=0.8, hspace=0.2)
        if n_modes==2: md={0:'low',1:'high'}            
            
        for index_sub_group,sub_group in enumerate(pop_groups): # the basic groups
            ydata=get_data_per_mode(s,2,sub_group,data,con_modes,index_pol,['n/d'])
                                
            for mode in xrange(n_modes): 
                
                dt=ydata[mode][0::2]; ax=fig.add_subplot(len(ydata),nc,index_sub_group+nc*mode+1); fig.text(0.02,pos[n_modes][mode],md[mode],fontsize=120); todel=[]                                        
                lg=[['downtown','suburbs','rural'],['Paris->Paris','suburban->suburban','rural->rural','rural->Paris','Suburban->Paris'],['<4','4-24','25-64','>64'],['actives','inactives'],['<1974','1975-2005','>2005']]
                cl1=[[cl[0],cl[2],cl[1]],[cl[0],cl[2],cl[1],cl[9],cl[10]],[cl[0],cl[2],cl[3],cl[4],cl[-1]],[cl[0],cl[1]],[cl[1],cl[2],cl[0]]]; lg1=lg; cl2=cl1
                    
                for indx,d in enumerate(dt):
                    if d<=0.05 or math.isnan(d): todel.append(lg1[index_sub_group][indx])
                lg1[index_sub_group],cl2[index_sub_group]=[i for i in lg1[index_sub_group] if i not in todel],[i for i in cl2[index_sub_group] if i not in todel]
					
                if lg1[index_sub_group]<>[] and sum(dt)>0.:					
                   dt=[d for d in dt if d>0.05]; explode=[0 for i in range(len(dt))]; explode[dt.index(max(dt))]=0.1                                             
                   wedge=ax.pie(dt, explode=explode, labels=lg1[index_sub_group], colors=cl2[index_sub_group], autopct='%1.1f%%', shadow=True)					
                   for w in wedge[0]: w.set_edgecolor('white')                       
                else: ax.text(0.45, 0.45, 'n/a', fontsize=100)
					
        plt.ion(); fig.savefig(plot_dir+'/'+lbl_out+'_'+pol+seas_str[s]+'.pdf', dpi=100); plt.close(fig)
		
def plot_map(s,type,data,cf,plot_dir,proj,X,Y):

    seas_str={0:'_annual',1:'_DJF',2:'_MAM',3:'_JJA',4:'_SON'}; ttl_str={0:'average period',1:'winter',2:'spring',3:'summer',4:'autumn'}
    unit=['ppb',r'$\mu$'+'g/m'+r'$^3$']; plot_ttl=['','EXPO-CHIMERE']; done=0; previous=0; fl_max=1; pol_str=['O'+r'$_3$','PM'+r'$_2$'+r'$_.$'+r'$_5$']
    
    for p in [i for i in range(2) if type<>'con']: # remove outlier cells     
        new=data[:,:,2+p]; new1=new[new>0.] # select non-zero cells                                                             
        try: 
           outliers=Set(new1[abs(new1-np.mean(new1))>5*np.std(new1)])
           for o in outliers: new[new==o]=np.mean(new1) # find cells with values>4*std and replace the value with domain avg
           if outliers.tolist()<>[]: print 'WARNING: found some outlier cells in '+ttl_str[s]+' exposure'
        except: pass            
            
    if type in ['n/d','home']: # generate the delta with the concentrations of chimere
       st=open(cf); data_avg=cPickle.load(st); fl_max=2; fl_name='exposure'; ttl1=ttl_str[s]+' exposure'; dataf=data; new_data=np.zeros((dataf.shape[0],dataf.shape[1],6))                      
       new_data[:,:,0:4]=dataf[:,:,0:4]; new_data[:,:,4]=dataf[:,:,2]-data_avg[:,:,2+s]; new_data[:,:,5]=dataf[:,:,3]-data_avg[:,:,7+s]; dataf=new_data; dataf[dataf==0.]=-999.    
       if type=='home': fl_name='exposure_home'
        
       for i in xrange(dataf.shape[0]):
           for j in xrange(dataf.shape[1]):
               if dataf[i,j,2]==0.: dataf[i,j,4]=0.
               if dataf[i,j,3]==0.: dataf[i,j,5]=0.                                        
                      
    elif type=='con': fl_name='concentration'; ttl1=ttl_str[s]+' concentration'; dataf=data; dataf[dataf==0.]=-999.
                     
    for index_pol,pol in enumerate('O3,PM25'.split(',')):
        
        if (fl_max==1 and index_pol==0) or fl_max==2: plt.ioff(); fig=plt.figure(figsize=(10,20)); plt.subplots_adjust(wspace=0.1,hspace=0.1); fig.suptitle(ttl1, fontsize=25)
              
        for fl in xrange(fl_max): # one is the exposure and one the difference exp-chimere
            
            cmap = plt.get_cmap('jet',20); negs=1; tmp=dataf[:,:,index_pol+2+2*fl]; min_val,max_val=round(min(tmp[tmp<>-999.]),1),round(max(tmp[tmp<>-999.]),1)
                        
            if type in ['n/d','home']: indx=fl; ttl=plot_ttl[fl]+' '+pol_str[index_pol]
            else: indx=index_pol; ttl=pol
                                                              
            step,increase_step=0.5,0.5; bins=0; panel=210+indx+1; ax=fig.add_subplot(panel); ttl=seas_str[s][1::]+' '+ttl
            if max_val-min_val>=15.:                                 
               if math.trunc(min_val)%2==1.: min_val=math.trunc(min_val)+1.                                
               step,increase_step=2.,1.
                   
            elif max_val-min_val>=5.: step,increase_step=1.,1.
            if max_val-min_val<=2.: step,increase_step=0.2,0.2
                
            if math.trunc(min_val)==math.floor(max_val): bounds=[min_val+step]+[round(min_val+step+i*step,2) for i in xrange(1,100) if min_val+step+i*step<max_val]; bins=len(bounds)
                                
            while bins not in range(3,13):
                                           
                  bounds=[]										   
                  for i in xrange(1,100):                          
                                            
                      if min_val>0. and math.trunc(min_val)+i*step>math.trunc(max_val): bins=len(bounds); break
                      elif min_val<0. and math.trunc(min_val)+i*step>math.floor(max_val): bins=len(bounds); break                       
                          
                      if min_val<0. and i==1: bounds=[math.trunc(min_val)]
                      bounds.append(float(math.trunc(min_val)+i*step))
                          
                  if bins not in range(3,13) or (bins in range(3,13) and min_val<0. and max_val>0.): step+=increase_step
                  if min_val<0. and max_val>0. and 0. not in bounds: bins=0                                      
                                            
            if type in ['n/d','home'] and fl==1: # only for the differences with CHIMERE
               if min_val<0. and max_val<=0.: cmap = plt.get_cmap('cool',13); negs=8
               elif min_val>=0. and max_val>0.: cmap = plt.get_cmap('hot_r',13); negs=8
               elif min_val<0. and max_val>0.: negs=len([i for i in bounds if i<0])+1
                   
            bounds=[bounds[0]-2*step]+[bounds[0]-step]+bounds+[bounds[-1]+step]; cmaplist=[cmap(i) for i in [0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19]] # the last negative is cmap(11)                                               
            if cmap.N<>20: cmaplist = [cmap(i) for i in xrange(cmap.N)]
                                
            new_cmaplist=cmaplist[(10-negs):(10-negs+len(bounds))]; new_cmaplist[0]=(1.,1.,1.,1.)                         
            cmap = cmap.from_list('Custom cmap',new_cmaplist,len(new_cmaplist)); norm = matplotlib.colors.BoundaryNorm(bounds,cmap.N)                                
            L=proj.pcolormesh(X,Y,dataf[:,:,index_pol+2+2*fl],cmap=cmap); L.set_clim([bounds[0],round(max_val+step,1)]); L=proj.drawmapscale(1.8,48.17,0.,0.,50,barstyle='fancy',fontsize=13)
            proj.readshapefile('data/communes_polyline','idf'); proj.readshapefile('data/Departments_region','dp',linewidth=2.,color='k',antialiased=1); plt.title(ttl, fontsize=35, y=1.02)               
            cbar = plt.colorbar(fraction=0.1,pad=0.04,shrink=1.,cmap=cmap,norm=norm,orientation='horizontal',spacing='uniform',ticks=bounds[1:-1],boundaries=bounds,extend='both')                
            cbar.ax.tick_params(labelsize=14); cbar.set_label(unit[index_pol],size=24)
    
        if type in ['n/d','home']: plt.ion(); plt.savefig(plot_dir+'/map_'+fl_name+'_'+pol+seas_str[s]+'.pdf',dpi=100)
   
    if type=='con': plt.savefig(plot_dir+'/map_'+fl_name+seas_str[s]+'.pdf',dpi=100)        
    elif type not in ['n/d','home']: plt.ion(); plt.savefig(plot_dir+'/map_'+fl_name+'.pdf',dpi=100); plt.close(fig) # for the daily
        
def get_data_per_mode(data,CUT,index_pol,dates):
      
    # for the moving average plot data is exp_per_bin_per_day with dimensions: index group (7: 'tot','P','PC','GC','1974','2005','2012'), index day, pollutants (2)
    # for the moving average plot in REF data is exp_per_bin_per_day_REF with dimensions: groups (2: population, Parisians), index day, pollutants (2), bins

    stp=0.1; final=[[],[]]; av=[0.,round(CUT,1),1000.]; np.seterr(invalid='ignore')
                                       
    for index_d,day in enumerate(dates):
                        
        dt=data[index_d][index_pol].astype(int); bins=np.linspace(0.1,0.1*len(dt),len(dt)); tmp=np.repeat(bins,dt); median=np.median(tmp)
        nb=sum(dt); new=np.zeros((10000)); new[0:len(dt)]=dt; md=[np.mean(tmp[tmp<median]),np.mean(tmp[tmp>=median])]
                                
        # append % of population and concentration of exposure                                            
        for mode in range(2): start,end=int(round(av[mode]/stp,0)),int(round(av[mode+1]/stp,0)); tot=sum(new[start:end]); final[mode].extend((tot/nb,md[mode])) 
                    
    return final
		
def get_min_max_plot(dat,nb_individuals,tol,step):
        
    for c in [i for i in dat if i<>0 and i>tol*nb_individuals]: min_con=((dat.tolist().index(c))*step); break    
    for index_c,c in enumerate(dat[::-1]): # this is the base of the max con. e.g max_con=18 represents the range from 18 to 18.2
        if c<>0 and c>tol*nb_individuals: max_con=(len(dat)-index_c-1)*step; break

    return min_con,max_con
    
def get_projection(dx,dy,GC,res):
    
    from diaries import get_cell_coord
    from shapely.geometry import Polygon
    from mpl_toolkits.basemap import Basemap,cm
    
    # get the coords of corner cells of the domain
    LL,LR,UR,UL=get_cell_coord(1,1,dx,GC),get_cell_coord(dx,1,dx,GC),get_cell_coord(dx,dy,dx,GC),get_cell_coord(1,dy,dx,GC); pL=Polygon([LL,LR,UR,UL]); lon0,lat0=pL.centroid.x,pL.centroid.y
    proj=Basemap(projection='lcc',lat_0=lat0,lon_0=lon0,width=(dx-1)*res*1000,height=(dy-1)*res*1000); return proj
    
def filter(dt, low_freq, signal_freq, order): b, a = butter_lowpass(low_freq, signal_freq, order=order); y = lfilter(b, a, dt); return y    
def butter_lowpass(low_freq, signal_freq, order): nyq = 0.5 * signal_freq; normal_cutoff = low_freq/nyq; b, a = butter(order, normal_cutoff, btype='low', analog=False); return b, a