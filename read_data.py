import xlrd as xls
import numpy as np
import netCDF4,cPickle,os,fnmatch,sys,multiprocessing
from os.path import isfile,join,isdir
from os import listdir
from collections import defaultdict

def cells2com(dom): # this file contains the communes and with which cells each commune interacts (cells start from 1)

    c2c=dict(); f=open('data/cells2comm_'+dom+'.prn')    
    for line in f:
        dt=[]; dt.append([str(n) for n in line.strip().split(' ')]); dt1=[]; dt1.append([k for k in [i for i in dt[0] if i<>'']])
        
        # create a dictionary with keys are the communes and the arrays include the cells inside that commune.
        if int(dt1[0][0]) not in c2c.keys(): c2c[int(dt1[0][0])]=[]        
        c2c[int(dt1[0][0])].append(int(dt1[0][1]))
        
    f=open('domains/'+dom+'/cells2comm_'+dom+'.dat','wb'); cPickle.dump(c2c,f,-1); f.close()    
        
    return c2c
    
def read_build_stats():

    build=dict()
    for d in [75,77,78,91,92,93,94,95]: build[d]=[[] for _ in range(2)] # build type (home,office),age groups (<74,74-05,>05)
    
    try: xls_file=xls.open_workbook('data/build_stats.xls',on_demand=True); sheet_nm=xls_file.sheet_names(); num_sheets=len(sheet_nm)
    except: print 'the file build_stats.xls was not found in the data directory.'; sys.exit(0)
    
    for index_sheet in xrange(num_sheets): # sheets are homes,offices                        
        data=xls_file.sheet_by_name(sheet_nm[index_sheet]); num_cells,num_rows=data.ncols-1,data.nrows-1; curr_cell = 0        
        while curr_cell < num_cells: # the 3 age groups for buildings
              curr_cell+=1; curr_row=0; cell_value=0
              while curr_row < num_rows: curr_row+=1; d=int(data.cell_value(curr_row,0)); cell_value=data.cell_value(curr_row,curr_cell); build[d][index_sheet].append(int(cell_value))
        
    # houses of the 1st category are split in 2: 75% becomes case1 and 25% case2, case2 then becomes case3 and case3 becomes case4. The 0 is for the >2012 buildings
    for d in [75,77,78,91,92,93,94,95]: cs1=int(0.75*build[d][0][0]); cs2=int(0.25*build[d][0][0]); cs3=build[d][0][1]; cs4=build[d][0][2]; build[d][0]=[cs1,cs2,cs3,cs4,0.]
        
    return build
    
def read_SIREN(conc_period):

    import glob
    from datetime import datetime,timedelta
    from sets import Set
    from numpy import searchsorted as srch
    from collections import Counter as cnt
    
    nbl=[131040,44140]; IO_indoor=dict(); max_tp={'Maison':5,'Ecoles':3,'Bureaux':3}; st_date=[1,10]; years=range(1991,2001) 
    IO_file='IO_indoor.dat'; print 'compiling SIREN I/O database...please wait...'; cs=np.linspace(0,1,21); lcs=len(cs)+1
    if conc_period>2008: years=range(2046,2056); IO_file='IO_indoor_'+str(conc_period)+'.dat'
                  
    for d in [75,77,78,91,92,93,94,95]:    
        IO_indoor[d]=np.zeros((2,12,3,5,21))    

        for indx_bl,bl in enumerate('Maison Bureaux Ecoles'. split(' ')):
            for tp in range(max_tp[bl]): # type based on age
                tmp=[[[] for _ in range(2)] for _ in range(12)]; bl_n=bl               
                
                for yr in years:
                    if indx_bl>0 and tp==0: bl_n='Bureaux-Ecoles'
                    dirs=glob.glob('data/IO_SIREN/'+str(d)+'/'+str(yr)+'/'+bl_n+' '+str(tp+1)+'*')                                     
                                                                                                           
                    for indx,fl in enumerate(dirs): # jan to sep and oct to dec                       
                        print 'processing '+fl; f=open(glob.glob(fl+'/'+'*')[0]); start_date=datetime(yr, st_date[int(fl[-3])-1], 1, 1, 0); t=-3            
                        
                        for ln,line in enumerate(f):
                            dt=[]; dt.append([str(n) for n in line.strip().split('\t')]); dt1=[]; dt1.append([k.replace(',', '.') for k in [i for i in dt[0] if i<>'']])                        
                            if nbl[indx]>ln>0: t+=3; curr_date=start_date+timedelta(minutes=t); mn=curr_date.month; tmp[mn-1][0].append(float(dt1[0][-7])); tmp[mn-1][1].append(float(dt1[0][-4]))
                        f.close()                       
                
                for m in range(12): # i do not take into account values above 1                 
                    for p in range(2): indx=srch(cs,tmp[m][p]); C=cnt(indx); X=np.array([C[j] for j in range(1,lcs)]); IO_indoor[d][p,m,indx_bl,tp]=X
                    	                                                                                                                                                                                   
    f=open('data/IO_SIREN/'+IO_file,'wb'); cPickle.dump(IO_indoor,f,-1); f.close()
	    
    return IO_indoor	

def network_stops():

    stops=dict(); f=open('data/pt_net.prn')    
    for line in f:
        dt=[]; dt.append([str(n) for n in line.strip().split(' ')]); dt1=[]; dt1.append([k for k in [i for i in dt[0] if i<>'']])

        # create a dictionary with keys are the communes and the arrays include the type of stop included. 
        # the prn file also has the coord of each stop, i do not utilize that information
        if dt1[0][4] not in stops.keys(): stops[dt1[0][4]]=[]
        if str(dt1[0][3]) not in stops[dt1[0][4]]: stops[dt1[0][4]].append(dt1[0][3])
        
    return stops

def read_road_network():
  
    roads=dict(); f=open('data/road_network.prn') # read the complete road network    
    for index_line,line in enumerate([i for i in f if i[0:2]<>'Xa']):         
        dt=[]
	dt.append([str(n) for n in line.strip().split(' ')])
	dt1=[]
	dt1.append([k for k in [i for i in dt[0] if i<>'']])
        roads[index_line]=[]
	roads[index_line].extend((float(dt1[0][0]),float(dt1[0][1]),int(dt1[0][4]),int(dt1[0][5]),int(dt1[0][6])))
    return roads

def read_BP():
     
    BP=dict(); f=open('data/BP.prn') # the coords of the periferique    
    for index_line,line in enumerate(f):         
        dt=[]; dt.append([str(n) for n in line.strip().split(' ')]); dt1=[]; dt1.append([k for k in [i for i in dt[0] if i<>'']])
        BP[index_line]=[]; BP[index_line].extend((float(dt1[0][0]),float(dt1[0][1]))) 
 
    return BP

def read_squales(): # notice: there are no measurements for ozone, only NO2 and PM2.5
    
    squales=defaultdict(dict); months_dic=dict; sq=dict(); IO_metro=dict(); obs_tmp=defaultdict(dict); obs=dict(); pol_list=['NO2','PM10']; mg2ppb=[1.92,1.]
    for m in '01 02 03 04 05 06 07 08 09 10 11 12'.split(): sq[m]=[[],[]]; obs[m]=[[],[]]; IO_metro[m]=[[],[]]

    months_dic={'janvier':'01','fevrier':'02','mars':'03','avril':'04','mai':'05','juin':'06','juillet':'07','aout':'08','septembre':'09','octobre':'10','novembre':'11','decembre':'12'}
  
    f=open('data/squales.prn')    
    for line in f:
        dt=[]; dt.append([str(n) for n in line.strip().split(' ')])
        
        if dt[0][0].isdigit() and len(dt[0])>1:
           day=int(dt[0][0])
           
           if day in [i for i in xrange(1,32)]: # isolate only the lines with data
              if day<=9: day_str='0'+str(day)
              else: day_str=str(day)
 	      mn_str=dt[0][1]; mn=months_dic[mn_str]; date=str(day_str)+'/'+str(mn)+'/2013'; h=str(dt[0][3])[0:2]
        
	      if h[1]==':': h='0'+h[0]
	      if date not in squales.keys() or h not in squales[date].keys(): squales[date][h]=[[] for _ in xrange(2)] # for the 2 pollutants
	      for index_p,p in enumerate(pol_list):
           
                if dt[0][7+index_p].isdigit(): squales[date][h][index_p].append(int(dt[0][7+index_p]))                   
                elif dt[0][7+index_p][0]=='<': squales[date][h][index_p].append(int(dt[0][7+index_p][1])) 

    for m in '01 02 03 04 05 06 07 08 09 10 11 12'.split(): # average per month
        tmp=[[],[]]
        for date in [i for i in squales.keys() if i[3:5]==str(m)]: 
            for h in squales[date].keys(): tmp[0].append(squales[date][h][0]); tmp[1].append(squales[date][h][1])     
        for index_pol in xrange(2): temp=tmp[index_pol]; sq[m][index_pol]=round(np.mean(np.hstack(temp))/mg2ppb[index_pol],2) 
                                                         
    f.close; squales.clear()
    
    f=open('data/obs.prn')  
    for line in f:
        dt=[]; dt.append([str(n) for n in line.strip().split(' ')])
                
        if dt[0][0][0].isdigit():
           
           dt1=[]
           for e in [i for i in dt[0] if i<>'']: dt1.append(e) # remove the '' elements
           
           h=dt1[1]
           if int(dt1[1])<=9: h='0'+dt1[1]
        
           date_tmp=[]; date_tmp.append([str(n) for n in dt1[0].strip().split('/')]); mn=date_tmp[0][1]; year=date_tmp[0][2]
           
           if int(date_tmp[0][1])<=9: mn='0'+date_tmp[0][1]
           
           date=str(day_str)+'/'+str(mn)+'/'+str(year)
           
           if date not in obs_tmp.keys() or h not in obs_tmp[date].keys(): obs_tmp[date][h]=[[] for _ in xrange(2)]
           
           for index_pol in [i for i in xrange(2) if dt1[i+2]<>'n/d']: obs_tmp[date][h][index_pol].append(float(dt1[index_pol+2]))  
           
    for m in '01 02 03 04 05 06 07 08 09 10 11 12'.split(): # average per month
        tmp=[[],[]]
        for date in [i for i in obs_tmp.keys() if i[3:5]==str(m)]: 
            for h in obs_tmp[date].keys(): tmp[0].append(obs_tmp[date][h][0]); tmp[1].append(obs_tmp[date][h][1])            
        for index_pol in xrange(2): temp=tmp[index_pol]; obs[m][index_pol]=round(np.mean(np.hstack(temp))/mg2ppb[index_pol],2) 
           
    f.close(); obs_tmp.clear(); tmp=[]; tmp1=[]  

    # for ozone i use an inverse relationship for the IO: e.g IO for NO2=2 means that ozone is half inside e.g. IO for ozone equal to 0.5
    for m in '01 02 03 04 05 06 07 08 09 10 11 12'.split(): IO_metro[m][0]=1./(sq[m][0]/obs[m][0]); IO_metro[m][1]=sq[m][1]/obs[m][1]  # for the 2 pollutants 
    for m in '01 02 03 04 05 06 07 08 09 10 11 12'.split(): tmp.append(IO_metro[m][0]); tmp1.append(IO_metro[m][1])
    IO_metro=[round(np.mean(tmp),1),round(np.mean(tmp1),1)] # i average the monthly IO into a single annual value
        
    return IO_metro
    
def read_tunnels(dom):

    tunnels=dict(); line=0; f = open('data/tunnels.txt', "r")	       
  
    for l in f.readlines():
        line+=1; ID,lon1,lat1,lon2,lat2 = l.split() # coords of the start and end of the tunnel (the ID is some internal coding) 
        if line>1 and int(ID) not in tunnels.keys(): tunnels[int(ID)]=[]; tunnels[int(ID)].append(lon1); tunnels[int(ID)].append(lat1); tunnels[int(ID)].append(lon2); tunnels[int(ID)].append(lat2)
  
    f.close()

    return tunnels
    
def read_conc(path,dom,chim_label):
 
    outdir='domains/'+dom+'/conc/'+chim_label; tmp=[i for i in xrange(23,1000,24)]; conc=dict() # O3 in ppb and PM25 in microgrammes/m3
  
    try: os.makedirs('domains')
    except: pass 
    try: os.makedirs('domains/'+dom)
    except: pass 
    try: os.makedirs('domains/'+dom+'/conc')
    except: pass 
    try: os.makedirs('domains/'+dom+'/conc/'+chim_label)
    except: pass
 
    # gets all the files in the concentration output dir and the dimensions of the domain   
    if not os.path.isdir(path): print ''; print 'the concentration path '+path+' does not exist.'; sys.exit(1)
        
    dlist = [path]+[join(path,d) for d in listdir(path) if isdir(join(path,d))] # all the directories in the path provided        
    for dir in dlist:        
        print 'we are in '+dir; flist=[f for f in listdir(dir) if fnmatch.fnmatch(f, 'out*.nc')]; num_files = len(flist); print 'Found '+str(num_files)+' netcdf files'               
        if num_files==0: continue
        nc=netCDF4.Dataset(join(dir,flist[0]),'r'); lon=nc.variables['lon'][:,:]; lat=nc.variables['lat'][:,:]; dims=lon.shape; dx,dy=dims[1],dims[0]; nc.close() # get domain dimensions     

        grid_coords=get_domain_coords(dx,dy,lon,lat,dom) # get the domain coordinates of all cells (LL corner of each cell)
        
        for f in flist: # for all files in the directory
              
            conc['O3']=[[] for _ in xrange(24)]; conc['PM25']=[[] for _ in xrange(24)]; nc=netCDF4.Dataset(join(dir,f),'r'); date=nc.variables['Times'][:]; print 'reading '+f
            for h in xrange(len(date)-1): # i skip the last hour of each file cause its the 00:00 of the next file

                yr=str(date[h])[2:3]+str(date[h])[6:7]+str(date[h])[10:11]+str(date[h])[14:15]; mn=str(date[h])[22:23]+str(date[h])[26:27]; d=str(date[h])[34:35]+str(date[h])[38:39]
                hr=int(str(date[h])[46:47]+str(date[h])[50:51]); dt=yr+'-'+mn+'-'+d
                        
                for pol in 'O3 PM25'.split():                    
                    var=nc.variables[pol][h,0] # i get the 1st layer concentrations only         

                    for y in xrange(dy):
                        for x in xrange(dx): val=var[y,x]; exec("conc['%s'][hr].append(round(val,4))" %pol)                      
                        
                if h in tmp: fo=open(outdir+'/conc_'+dt+'.dat','wb'); cPickle.dump(conc,fo,-1); fo.close(); conc['O3']=[[] for _ in xrange(24)]; conc['PM25']=[[] for _ in xrange(24)]
                      
            nc.close(); conc.clear()
        
    return grid_coords

def get_domain_coords(dx,dy,lon,lat,dom):

    from shapely.geometry import Point, Polygon
    from mpl_toolkits.basemap import Basemap 
    from diaries import create_polygons

    LL=[lon[0,0],lat[0,0]]; LR=[lon[0,dx-1],lat[0,dx-1]]; UL=[lon[dy-1,0],lat[dy-1,0]]; UR=[lon[dy-1,dx-1],lat[dy-1,dx-1]]; points=[]; points.extend((LL,LR,UR,UL))    
    pL=Polygon(points) # create a polygon of domain boundaries
    lon0=pL.centroid.x;lat0=pL.centroid.y # get its centroid which is actually the centroid of the domain
    proj=Basemap(projection='lcc',lat_0=lat0,lon_0=lon0,width=dx*4000.,height=dy*4000.) # create the projection
    X_center,Y_center=proj(lon,lat); X_corner=X_center-2000.; Y_corner=Y_center-2000.; X,Y=proj(X_corner,Y_corner,inverse=True) # get corners and transform back to degrees

    # create a dictionary with the lat,lons of the grid
    id_cell=0; grid_coords=dict()
    for y in xrange(dy):
        for x in xrange(dx): id_cell+=1; grid_coords[id_cell]=[]; grid_coords[id_cell].append(X[y,x]), grid_coords[id_cell].append(Y[y,x])   
    f=open('domains/'+dom+'/COORDS_'+dom+'.dat','wb'); cPickle.dump(grid_coords,f,-1); f.close() # these are the LL corners of the cells, cells start from 1 NOT 0  

    # now i create a fake row, column in the grid in order to help create the grid_polygon dictionary
    tmp1=np.empty((dy))
    for y in xrange(dy): add=X_corner[y,-1]+4000.; tmp1[y]=add 
    c = np.hstack((X_corner, np.atleast_2d(tmp1).T)); X_corner=c

    tmp1=np.empty((dx)); tmp1=X_corner[-1]+(X_corner[-1]-X_corner[-2]); c=np.vstack((X_corner,tmp1)); X_corner=c

    tmp1=np.empty((dy))
    for y in xrange(dy): add=Y_corner[y,-1]+(Y_corner[y,-1]-Y_corner[y,-2]); tmp1[y]=add
    c = np.hstack((Y_corner, np.atleast_2d(tmp1).T)); Y_corner=c

    tmp1=np.empty((dx+1))
    for x in xrange(dx+1): add=Y_corner[-1,x]+4000.; tmp1[x]=add
    c = np.vstack((Y_corner,tmp1)); Y_corner=c
 
    grid_coords.clear(); id_cell=0; X,Y=proj(X_corner,Y_corner,inverse=True) # get corners and transform back to degrees
    for y in xrange(dy+1):
        for x in xrange(dx+1): id_cell+=1; grid_coords[id_cell]=[]; grid_coords[id_cell].append(X[y,x]), grid_coords[id_cell].append(Y[y,x])   

    if not os.path.isfile('domains/'+dom+'/grid_polygons.dat'): print ''; print 'generating vector form of the grid...'; grid_polyg=create_polygons(grid_coords,dx,dy,dom) 

    grid_coords.clear(); st=open('domains/'+str(dom)+'/COORDS_'+str(dom)+'.dat','rb'); grid_coords=cPickle.load(st); st.close() 

    return grid_coords # keys start from 1

def read_communes_coords():

    communes=dict(); neig=dict()

    f = open('data/IdF.txt', "r"); f1 = open('data/neighbours.txt', "r")

    # this dictionary contains the idf communes as keys and the neighbouring commune in the array
    for l in f1.readlines(): 
        start,neighbour= l.split(" ")
        if int(start) not in neig.keys(): neig[int(start)]=[]        
        neig[int(start)].append(int(neighbour))

    for l in f.readlines():
        CODE,ID,lon,lat,SURF = l.split() 
        if int(CODE) not in communes.keys(): 
           communes[int(CODE)]=[]; communes[int(CODE)].extend((int(ID),lon,lat,SURF,neig[int(CODE)]))

    f.close();f1.close()

    return communes

def read_diaries(label):

    data_diaries=dict()

    #open the diaries datafile
    try: xls_file=xls.open_workbook('data/diaries_input_'+label+'.xls',on_demand=True); sheet_nm=xls_file.sheet_names(); num_sheets=len(sheet_nm)
    except: print 'the file diaries_input_'+label+'.xls was not found in the data directory.'; sys.exit(1)
    
    for index_sheet in [i for i in xrange(num_sheets) if sheet_nm[i] not in ['Diaries']]:
            
        max_dim=2 # The dic contains a 2D array, 1st dimension is the travels by public transport and the second by car (for per activity tabs)
        if sheet_nm[index_sheet] in 'moves'.split(): max_dim=1 # The dic contains a 1d array
        elif sheet_nm[index_sheet] in 'Nb obs_trans'.split(): max_dim=3 # The dic contains a 3d array, the number of people in movement for actives, inactive, students        
        elif sheet_nm[index_sheet] in 'public'.split(): max_dim=4 # The dic contains a 4d array
        elif sheet_nm[index_sheet] in 'distance'.split(): max_dim=5 # The dic contains a 5d array, distance travelled per department
        elif sheet_nm[index_sheet] in 'time_act transport obs_mot Nb_mot'.split(): max_dim=6 # the time spend in each activity type        
        elif sheet_nm[index_sheet] in 'speed'.split(): max_dim=8 # the time spend per movement for the intercommunale movements, also the per day moves by car, public, others
        elif sheet_nm[index_sheet] in 'flux'.split(): max_dim=10 # the fluxes of movements amongst deps
        
        if sheet_nm[index_sheet] not in data_diaries.keys(): data_diaries[sheet_nm[index_sheet]]=[[] for _ in xrange(max_dim)] # create the dictionary with the raw data

        # Get the max number of cols,rows (starts from 1)
        exec('data%s=xls_file.sheet_by_name(sheet_nm[index_sheet])' %index_sheet)
        exec('num_cells=data%s.ncols-1' %index_sheet); exec('num_rows=data%s.nrows-1' %index_sheet) 

        curr_cell = 0
        while curr_cell < num_cells: 
              curr_cell += 1; curr_row = 0
              while curr_row < num_rows:
                    curr_row += 1      
                    exec('cell_value = data%s.cell_value(curr_row,curr_cell)' %index_sheet) # rows,columns
                    data_diaries[sheet_nm[index_sheet]][curr_cell-1].append(cell_value)

    f=open('data/diaries_input_'+label+'.dat','wb'); cPickle.dump(data_diaries,f,-1); f.close()
    
    return data_diaries

def read_demographics(fl,list_deps):
        
    ind=fl[7:11] # e.g. POP1, POP5 
    if ind=='PRIN': ind=fl[7:13] # It will give ind e.g. PRINC1
    if ind=='PRINC2': ind=fl[7:14] # It will give ind e.g. PRINC24
    if ind[0:1]<>'D': exec('data_%s=dict()' %ind) # 1d key dict (all but the movements statistics)
    if ind[0:1]=='D': exec('data_%s=defaultdict(dict)' %ind ) # 2d key dict (for the movements statistics)
    
    #open the xls INSEE datafile
    xls_file=xls.open_workbook('data/'+fl,on_demand=True); num_sheets=2
    if ind[0:4]<>'PRIN' and ind[0:4]<>'DTR_' and ind[0:4]<>'DET_': data1=xls_file.sheet_by_name('Communes'); data2=xls_file.sheet_by_name('Arrondissements municipaux') # Contains the metropolitan areas
    if ind[0:4]=='PRIN': data1=xls_file.sheet_by_name('COM'); data2=xls_file.sheet_by_name('ARM') # The buildings statistics
    if ind[0:1]=='D': data1=xls_file.sheet_by_name('TOTAL'); data2=xls_file.sheet_by_name('FLUX>=100') # The movements statistics
    
    for index_sheet in xrange(1,num_sheets+1):  
    
        keep_gender_male=[];keep_house_type=[];keep_move_type=[]
        
        # Get the max number of cols,rows (starts from 1)        
        exec('num_cells=data%s.ncols-1' %index_sheet); exec('num_rows=data%s.nrows-1' %index_sheet) 

        curr_row = -1
        while curr_row < num_rows:
              curr_row+=1	  
     
              exec('cell_value = data%s.cell_value(curr_row,0)' %index_sheet) # rows,columns
                 
              if cell_value == 'CODGEO':                                                    
                 start_row = curr_row+1 #1 line after the 'CODGEO header'
                          
                 curr_cell = -1
                 while curr_cell < num_cells: 
                       curr_cell += 1
                       exec('cell_value = data%s.cell_value(curr_row,curr_cell)' %index_sheet) # rows,columns
                       if 'SEXE1' in cell_value: keep_gender_male.append(curr_cell) # For the person statistics
                       if 'TYPLR' in cell_value: keep_house_type.append(curr_cell) # For the house statistics
                       if cell_value[0:4] in ['DCLT','DCET']: keep_move_type.append(curr_cell) # For the movements statistics (FLUX>=100 sheet)
                       if cell_value[0:4]=='C09_': keep_move_type.append(curr_cell-2) # For the movements statistics (TOTAL sheet)
                 break
                
        for curr_row in xrange(start_row,num_rows+1):

            got_fin_id=0; exec('cell_value = data%s.cell_value(curr_row,0)' %index_sheet) # rows,columns
            if ind[0:1]=='D' and index_sheet==2: exec('cell_value_fin = data%s.cell_value(curr_row,2)' %index_sheet)
            if curr_row>=start_row and (cell_value[0:2] not in ['2A','2B']) and cell_value<>'75056' and (int(cell_value)/1000 in list_deps): # Keep the department ids that you are interested in     
               got_com_id=int(cell_value) # Stores the keys of the dictionary, here keys are the community id                
               if ind[0:1]=='D' and index_sheet==2 and cell_value_fin[0:2] not in ['AL','BE','LU','MO','SU','ZZ']: got_fin_id=int(cell_value_fin) # Final destination commune id           
               
               if ind=='POP1': # Basic population statistics
                  if got_com_id not in data_POP1.keys(): data_POP1[got_com_id]=np.zeros((2,101)) # The dic contains a 2D array, 1st dimension is the gender and then 101 age groups
                  exec('for indx,curr_cell in enumerate(keep_gender_male): data_POP1[got_com_id][0,indx]=data%s.cell_value(curr_row, curr_cell)' %index_sheet)                                            
                  exec('for indx,curr_cell in enumerate([i+1 for i in keep_gender_male]): data_POP1[got_com_id][1,indx]=data%s.cell_value(curr_row, curr_cell)' %index_sheet)                                            

               elif ind[0:1]=='D': # Movements statistics (DET or DTR)
                  exec('data_%s[got_com_id][got_fin_id]=[]' %ind) # The dic contains a 1D array
                  exec('for curr_cell in keep_move_type: data_%s[got_com_id][got_fin_id].append(data%s.cell_value(curr_row, curr_cell+2))' %(ind,index_sheet))          
                  
               elif ind=='PRINC1': # Statistics on the number of rooms in the houselhold

                  # House type groups in the PRINC1 file 
                  # 1 : maisons
                  # 2 : appartements
                  # 3 : autres
                  # Age groups in the PRINC1 file are <20,20-25,25-40,40-55,55-65,65-80,>80
                                                      
                  # The dic data_PRINC1 contains a 3D array, 1st dimension the age (7 groups), 2nd dimension the house type (3 groups), 3rd dimension the pop per tot room number of their house
                  if got_com_id not in data_PRINC1.keys(): data_PRINC1[got_com_id]=[[[] for _ in xrange(3)] for _ in xrange(14)]
                  for index_cell,curr_cell in enumerate(keep_house_type[::42]): #every 42 columns the house type group changes                  
                      for a in xrange(7): # The age groups
                          for n in xrange(6): exec('data_PRINC1[got_com_id][a][index_cell].append(data%s.cell_value(curr_row, curr_cell+n*7+a))' %index_sheet)

               elif ind=='PRINC26': # Statistics on the surface of the house

                  # House type groups in the PRINC2 file 
                  # 1 : moins de 40 m2
                  # 2 : de 40 m2 a moins de 100 m2
                  # 3 : 100 m2 ou plus
                                                      
                  # The dic data_PRINC26 contains a 2D array, 1st dimension the house type (3 groups), 2nd dimension the pop per surface type
                  if got_com_id not in data_PRINC26.keys(): data_PRINC26[got_com_id]=[[] for _ in xrange(3)]
                  for index_cell,curr_cell in enumerate(keep_house_type[::15]): #every 15 columns the house type group changes                
                      for s in xrange(3): # 3 surface groups
                          sm=0.
                          for j in xrange(5): # This is another level in the xls file which we dont use (STOCD)
                              exec('sm=sm+data%s.cell_value(curr_row, curr_cell+j*3+s)' %index_sheet)                             
                          data_PRINC26[got_com_id][index_cell].append(sm)

               elif ind=='PRINC3': # Statistics on the heating type

                  # House type groups in the PRINC30M file 
                  # 1 : chauffage urbain
                  # 2 : gaz de ville ou de reseau
                  # 3 : fioul (mazout)
                  # 4 : electricite
                  # 5 : gaz en bouteilles ou en citerne
                  # 6 : autre
                                                      
                  # The dic data_PRINC3 contains a 2D array, 1st dimension the house type (3 groups), 2nd dimension the pop per heating type
                  if got_com_id not in data_PRINC3.keys(): data_PRINC3[got_com_id]=[[] for _ in xrange(6)]
                  for index_cell,curr_cell in enumerate(keep_house_type[::30]): #every 30 columns the house type group changes                
                      for h in xrange(6): # 6 heating type groups
                          sm=0.
                          for j in xrange(5): # This is another level in the xls file which we dont use (STOCD)
                              exec('sm=sm+data%s.cell_value(curr_row, curr_cell+j*6+h)' %index_sheet)                             
                          data_PRINC3[got_com_id][index_cell].append(sm)
                      
               elif ind=='ACT2': # Statistics on the contract type of the invividual
                                    
                  # Create a 4d array dictionary: 1st dimension the gender, 2nd dimension age, 3rd salarie/non-salarie, 4rd dimension the pop of the 2 types of work (full-time, part-time)
                  if got_com_id not in data_ACT2.keys(): data_ACT2[got_com_id]=np.zeros((2,11,2,2))
                  for index_cell,curr_cell in enumerate(keep_gender_male[0:22:2]): #every 2 columns the age group changes and we have the 2 columns (salarie, non salaries)                                    
                      for s in xrange(2): # salarie
                          for tp in xrange(2): exec('data_ACT2[got_com_id][0,index_cell,s,tp]=data%s.cell_value(curr_row, curr_cell+tp*44+s)' %index_sheet)

                  for index_cell,curr_cell in enumerate([i+11*2 for i in keep_gender_male][0:22:2]): #every 2 columns the age group changes and we have the 2 columns (salarie, non salaries)                      
                      for s in xrange(2):
                          for tp in xrange(2): exec('data_ACT2[got_com_id][1,index_cell,s,tp]=data%s.cell_value(curr_row, curr_cell+tp*44+s)' %index_sheet)
                
               elif ind=='NAV1': # Statistics on the active population and workplace

                  # Activity groups in the NAV1 file   
                  # 1 : commune de residence
                  # 2 : autre commune du departement de residence
                  # 3 : autre departement de la region de residence
                  # 4 : autre region en France metropolitaine
                  # 5 : autre (Dom, Com, etranger)
                                    
                  # The dic data_NAV1 contains a 3D array, 1st dimension the gender, 2nd dimension the 11 age groups, 3rd dimension the pop of each of the 5 categories above
                  if got_com_id not in data_NAV1.keys(): data_NAV1[got_com_id]=np.zeros((2,11,5))
                  for index_cell,curr_cell in enumerate(keep_gender_male[::5]): #every 5 columns the age group changes and we have the 6 different activities of the population                      
                      for a in xrange(5): exec('data_NAV1[got_com_id][0,index_cell,a]=data%s.cell_value(curr_row, curr_cell+a)' %index_sheet)                                        
                  for index_cell,curr_cell in enumerate([i+11*5 for i in keep_gender_male][::5]): 
                      for a in xrange(5): exec('data_NAV1[got_com_id][1,index_cell,a]=data%s.cell_value(curr_row, curr_cell+a)' %index_sheet)

               elif ind=='NAV2': # Statistics on the active population and the mode of transport

                  # Activity groups in the NAV2 file                   
                  # 1 : pas de transport
                  # 2 : marche a pied
                  # 3 : deux roues
                  # 4 : voiture, camion, fourgonnette
                  # 5 : transports en commun
                                    
                  # The dic data_NAV2 contains a 3D array, 1st dimension the gender, end the workplace group, 3rd dimension the pop of each of the 5 categories above
                  if got_com_id not in data_NAV2.keys(): data_NAV2[got_com_id]=np.zeros((2,5,5))
                  for index_cell,curr_cell in enumerate(keep_gender_male[::5]):
                      for a in xrange(5): exec('data_NAV2[got_com_id][0,index_cell,a]=data%s.cell_value(curr_row, keep_gender_male[0]+index_cell+5*a)' %index_sheet)
                  for index_cell,curr_cell in enumerate([i+5*5 for i in keep_gender_male][::5]): 
                      for a in xrange(5): exec('data_NAV2[got_com_id][1,index_cell,a]=data%s.cell_value(curr_row, keep_gender_male[0]+index_cell+5*a)' %index_sheet)
                      
               elif ind=='POP5': # Statistics on the inactive population

                  # Activity groups in the POP5 file                  
                  # 11 : actifs ayant un emploi
                  # 12 : chomeurs
                  # 21 : retraites ou preretraites
                  # 22 : eleves, etudiants, stagiaires non remuneres
                  # 24 : femmes ou hommes au foyer
                  # 26 : autres inactifs
                   
                  # The dic data_POP5 contains a 3D array, 1st dimension the gender, 2nd dimension the 11 age groups (15-65 years), 3rd dimension the pop of each of the 6 activity groups
                  if got_com_id not in data_POP5.keys(): data_POP5[got_com_id]=np.zeros((2,11,6))      
                  for index_cell,curr_cell in enumerate(keep_gender_male[::6]): #every 6 columns the age group changes and we have the 6 different activities of the population       
                      for a in xrange(6): exec('data_POP5[got_com_id][0,index_cell,a]=data%s.cell_value(curr_row, curr_cell+a)' %index_sheet)
                                            
                  for index_cell,curr_cell in enumerate([i+11*6 for i in keep_gender_male][::6]): # females                       
                      for a in xrange(6): exec('data_POP5[got_com_id][1,index_cell,a]=data%s.cell_value(curr_row, curr_cell+a)' %index_sheet)
                      
    f=open('data/'+ fl[0:len(fl)-4]+ '.dat','wb'); exec('cPickle.dump(data_%s,f,-1)' %ind); f.close()
    
def read_age_build():

    # 1st column holds the home commune, the 2nd, 3rd, 4th, 5th columns the number of buildings per age group (<3, 4 to 24, 25 to 64, >65). the 6th is the building age group:
    
    # 11 : Before 1949
    # 12 : 1949-1974
    # 13 : 1975-1981
    # 14 : 1982-1989
    # 15 : 1990-1998
    # 21,22,23,24,25,26,27,28,29,30,31,32: 1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010
    # 99: Under construction (i assume 2011)    
    
    data_age_build=dict(); ch_grp_age={11:1,12:1,13:2,14:2,15:2,21:2,22:2,23:2,24:2,25:2,26:2,27:2,28:3,29:3,30:3,31:3,32:3,99:3}; f=open('data/BUILDING_AGE.ins')
    
    for line in f:
        dt=[int(n) for n in line.strip().split(';')]; com=dt[0]
                
        if com not in data_age_build.keys(): data_age_build[com]=np.zeros((4,4)) # i use 4 age building groups (<1974 no rehab, <1974 rehab, 1975-2005, >2005)
        
        for a in range(4): data_age_build[com][a,ch_grp_age[dt[-1]]]+=dt[a+1] # here i put all the <1974 buildings in dimension 2, position 1
    
    for com in data_age_build.keys(): 
    
        data_age_build[com][:,0]=0.25*data_age_build[com][:,1]; data_age_build[com][:,1]=0.75*data_age_build[com][:,1] # I split the <1974 to rehabilitated (75%) and not (25%)
        
        for a in [i for i in range(4) if np.sum(data_age_build[com][i]==0.)]: data_age_build[com][a,0]=1.
                        
    f.close(); return data_age_build
                  
def read_definition_file():
   
   f=open('expo.par')
   for line in f:
   
       if '#' in line or not line.strip(): continue
       else: eq=line.index("="); param_str=line[0:eq]; param_nb=line[eq+1::]; exec(line)
              
       if param_str in 'exp_start_date,exp_end_date,stats_start_date,stats_end_date'.split(",") and len(param_nb)<12: print 'ERROR: the start/end dates are incorrect. Please check your configuration file.'; sys.exit(1)
       if param_str=='procc' and int(param_nb)>multiprocessing.cpu_count(): print 'ERROR: your system has '+str(multiprocessing.cpu_count())+' processors and you have selected '+param_nb; sys.exit(1)    
       if param_str in ['building_stock','conc_period'] and int(param_nb) not in [2008,2050]: print 'the only dates supported for the future scenarios are 2008 and 2050'; sys.exit(0)
         
   return label,chim_label,run_label,dom,chimere_exp,static,procc,traj_images,building_stock,conc_period,allow_steps_prf,allow_steps_dia,allow_steps_con,allow_steps_exp,write_asc_dia, \
   path_to_conc,path_for_output,path_to_plots,calc_stats_dia,calc_stats_exp,plot_stats_exp,exp_start_date,exp_end_date,stats_start_date,stats_end_date
