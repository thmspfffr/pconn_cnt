%% PLOT SENSOR-LEVEL DFA RESULTS
% pconn_cnt_sens_dfa_plot

clear

v_cnt         = 20;
v_res         = 19;

SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 23 24];

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';

% load ~/pconn_cnt/matlab/sensorselection_new.mat
load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
addpath /home/tpfeffer/pconn/matlab

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

% subsel = 1:length(SUBJLIST);


%%


addpath ~/pcbi/
ord       = pconn_randomization;

% subsel = 1:length(SUBJLIST);
clear dfa_cnt_all dfa_res_all tmp cvar_cnt_all pow_cnt_all cvar_res_all pow_res_all var_cnt_all par

SUBJLIST =  [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

for ifoi = 1:4
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      tmp_cnt = pcbi_cnt(isubj);

      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));

      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));

      dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2); 
      var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2); 
      cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      pow_cnt_all(:,isubj,m,ifoi)  = nanmean(par.pow,2); clear par


      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));

      dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);    
      pow_res_all(:,isubj,m,ifoi)  = nanmean(par.pow,2);
      
    end
  end
end

% cnt = cnt(SUBJLIST,:);
  dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
  var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
  cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));
  pow_cnt_all  = double(pow_cnt_all(:,SUBJLIST,:,:));

  dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
  var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
  cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
  pow_res_all  = double(pow_res_all(:,SUBJLIST,:,:));

%% INDIV SUBJECTS
% SUBJ = [1:5; 6:10; 11:15; 16:19];

str = 'dfa';

for ifoi = 1 : 5

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

for isubj = 1 : 1
h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')

subj = SUBJLIST(isubj);

subplot(1,3,1)

par = nanmean(squeeze(all_par_res(:,isubj,:)),3);

pars = [];
figure; set(gcf,'color','white');
pars.scale=[0.75 0.9];
pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.resolution = 300;
showfield_colormap(par(:,1),sa.locs_2D,pars);
% showfield_colormap(a,sa.locs_2D,pars);

colormap(hot)

subplot(1,3,2)

par = nanmean(squeeze(all_par_cnt(:,isubj,:)),3);
par = nanmean(nanmean(squeeze(all_par_cnt(:,:,:)),3),2);
% figure; set(gcf,'color','white');
% pars.scale=[0.5 0.75];
pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.resolution = 300;
showfield_colormap(par(:,1),sa.locs_2D,pars);

subplot(1,3,3)

par = nanmean(squeeze(all_par_cnt(:,isubj,:)-all_par_res(:,isubj,:)),3);

% figure; set(gcf,'color','white');
% pars.scale=[0.5 0.75];
pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.resolution = 300;
showfield_colormap(par(:,1),sa.locs_2D,pars);

print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_task-rest_s%d_f%d_v%d.jpg',str,subj,ifoi,v_cnt))%%

close
end
end
error('!!')



%% PHARMA COMPARISON & TASK VS REST
% -------------------------------------------

ifoi = 3;
str = 'dfa';

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

contrasts = [2 1; 3 1; 2 3];
dfa_d2 = nanmean(nanmean(all_par_cnt,3),2)-nanmean(nanmean(all_par_res,3),2);

for icontr = 1 : 1

  % PHARMA CONTRAST DURING CNT
  dfa_d1 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2)),2);

  % CONTRAST CNT VS REST

  figure; set(gcf,'color','white');

  r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
  pars.scale=[-r r];

  pars.cbar = 0;
  pars.markersize = 0;
  pars.linewidth = 9;
  pars.resolution = 300;

  showfield_colormap(dfa_d1,sa.locs_2D,pars);
  colormap(parula)
  print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))

end

showfield_colormap(dfa_d2,sa.locs_2D,pars);
colormap(parula)

print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_task-rest_f%d_v%d.jpg',str,ifoi,v_cnt))

%% COMPARE MEASURES 
% -------------------------------------

ifoi = 1;

for isubj = 1 : 18

  p_dfa  = nanmean(dfa_res_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_res_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(:,isubj,:,ifoi),3);
  
  r_res(1,isubj) = corr(p_dfa,p_var);
  r_res(2,isubj) = corr(p_dfa,p_cvar);
  r_res(3,isubj) = corr(p_cvar,p_var);
  
  p_dfa  = nanmean(dfa_cnt_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(:,isubj,:,ifoi),3);
  
  r_cnt(1,isubj) = corr(p_dfa,p_var);
  r_cnt(2,isubj) = corr(p_dfa,p_cvar);
  r_cnt(3,isubj) = corr(p_cvar,p_var);
  
  

end

for i = 1 : 3
  
  [~,p_res(i)] = ttest(r_res(:,i));
 	[~,p_cnt(i)] = ttest(r_cnt(:,i));
  
end

%%
ifoi = 1;

clear r_res r_cnt

for isens = 1 : 268
  
  p_dfa  = nanmean(dfa_res_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_res_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(isens,:,:,ifoi),3);
  
  r_res(1,isens) = corr(p_dfa',p_var');
  r_res(2,isens) = corr(p_dfa',p_cvar');
  r_res(3,isens) = corr(p_cvar',p_var');
  
  p_dfa  = nanmean(dfa_cnt_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(isens,:,:,ifoi),3);
  
  r_cnt(1,isens) = corr(p_dfa',p_var');
  r_cnt(2,isens) = corr(p_dfa',p_cvar');
  r_cnt(3,isens) = corr(p_cvar',p_var');
  
end

figure; set(gcf,'color','white');
for i = 1 : 3
  
  subplot(3,2,i*2-1)

  pars.scale=[-1 1];
  pars.cbar = 0;
  pars.markersize = 0;
  pars.linewidth = 9;
  pars.resolution = 300;

  showfield_colormap(r_res(i,:),sa.locs_2D,pars);
  colormap(parula)

  subplot(3,2,i*2)
  
  showfield_colormap(r_cnt(i,:),sa.locs_2D,pars);
  colormap(parula)
  
end

%%