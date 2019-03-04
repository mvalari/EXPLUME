f=open('expo.par')
for line in f:

       if '#' in line or not line.strip(): continue
       else: eq=line.index("="); param_str=line[0:eq]; param_nb=line[eq+1::]; exec(line)

       if param_str in 'exp_start_date,exp_end_date,stats_start_date,stats_end_date'.split(",") and len(param_nb)<12: print 'ERROR: the start/end dates are incorrect. Please check your configuration file.'; sys.exit(1)
       if param_str=='procc' and int(param_nb)>multiprocessing.cpu_count(): print 'ERROR: your system has '+str(multiprocessing.cpu_count())+' processors and you have selected '+param_nb; sys.exit(1)
       if param_str in ['building_stock','conc_period'] and int(param_nb) not in [2008,2050]: print 'the only dates supported for the future scenarios are 2008 and 2050'; sys.exit(0)

   return label,chim_label,run_label,dom,chimere_exp,static,procc,traj_images,building_stock,conc_period,allow_steps_prf,allow_steps_dia,allow_steps_con,allow_steps_exp,write_asc_dia, \
   path_to_conc,path_for_output,path_to_plots,calc_stats_dia,calc_stats_exp,plot_stats_exp,exp_start_date,exp_end_date,stats_start_date,stats_end_date

