%% PLOT SENSOR-LEVEL DFA RESULTS
% pconn_cnt_sens_dfa_plot

clear
v_res = 21;
v = 19;
SUBJLIST =  [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';
load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
addpath /home/tpfeffer/pconn/matlab

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m


addpath ~/pcbi/
%%
ord       = pconn_randomization;

for ifoi = 1: 5;
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      tmp_cnt = pcbi_cnt(isubj);

      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));

      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));

      dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2); clear tmp
      var_cnt_all(:,isubj,m,ifoi)  = double(nanmean(par.var,2)); clear tmp
      cvar_cnt_all(:,isubj,m,ifoi) = double(nanmean(par.cvar,2)); clear tmp

      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));

      dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_res_all(:,isubj,m,ifoi)  = double(nanmean(par.var,2));
      cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);    

    end
  end
end

cnt          = cnt(SUBJLIST,:);
dfa_cnt_all  = dfa_cnt_all(:,SUBJLIST,:,:);
var_cnt_all  = var_cnt_all(:,SUBJLIST,:,:);
cvar_cnt_all = cvar_cnt_all(:,SUBJLIST,:,:);

dfa_res_all  = dfa_res_all(:,SUBJLIST,:,:);
var_res_all  = var_res_all(:,SUBJLIST,:,:);
cvar_res_all = cvar_res_all(:,SUBJLIST,:,:);


%% INDIV SUBJECTS

str = 'dfa';


for ifoi = 1 : 5;

if strcmp(str,'var')
  d1 = double(var_cnt_all(:,:,2,ifoi));
  d2 = double(var_res_all(:,:,1,ifoi));
elseif strcmp(str,'cvar')
  d1 = double(cvar_cnt_all(:,:,2,ifoi));
  d2 = double(cvar_res_all(:,:,1,ifoi));
elseif strcmp(str,'dfa')
  d1 = double(dfa_cnt_all(:,:,2,ifoi));
  d2 = double(dfa_res_all(:,:,1,ifoi));
end

 

for isubj = 1:18
  
subj = SUBJLIST(isubj);

h=figure; set(h,'color','white'); 

pars.cbar = 0;
pars.markersize = 0;
pars.linewidth = 4;
pars.resolution = 300;

% ---------------------------------
% PLOT REST
% ---------------------------------

subplot(1,3,1)

par = d1(:,isubj);

showfield_colormap(par,sa.locs_2D,pars);
colormap(hot)

% ---------------------------------
% PLOT TASK
% ---------------------------------

subplot(1,3,2)

par = d2(:,isubj);

showfield_colormap(par,sa.locs_2D,pars);
colormap(hot)

% ---------------------------------
% PLOT DIFFERENCE
% ---------------------------------

subplot(1,3,3)

par = d1(:,isubj)-d2(:,isubj);

showfield_colormap(par,sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_task-rest_s%d_f%d_v%d.jpg',str,subj,ifoi,v))%%

end
end
error('!!')

% ---------------------------------
%% PLOT PHARMA CONTRAST
% ---------------------------------
contrasts = [2 1; 3 1; 2 3];

ifoi = 3;
str  = 'dfa'; 

clear d1 d2

for icontr = 1 : 3

  if strcmp(str,'var')
    d1 = nanmean(var_cnt_all(:,:,contrasts(icontr,1),ifoi)-var_cnt_all(:,:,contrasts(icontr,2),ifoi),2);
  elseif strcmp(str,'cvar')
    d1 = nanmean(cvar_cnt_all(:,:,contrasts(icontr,1),ifoi)-cvar_cnt_all(:,:,contrasts(icontr,2),ifoi),2);
  elseif strcmp(str,'dfa')
    d1 = nanmean(dfa_cnt_all(:,:,contrasts(icontr,1),ifoi)-dfa_cnt_all(:,:,contrasts(icontr,2),ifoi),2);
  end


  load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v));
  senssel = find(stats.mask);

  pars            = [];
  pars.markersel  = senssel;
  pars.markersize = 20;
  pars.markersize = 0;
  pars.linewidth  = 9;
  pars.resolution = 300;
  pars.cbar       = 0;

  % ---------------------------------
  % PHARMA CONTRAST DURING CNT
  % ---------------------------------

  r = max([abs(min(d1)) abs(max(d1))]);
  pars.scale= [-r r];

  figure; set(gcf,'color','white');
  showfield_colormap(double(d1),sa.locs_2D,pars);
  colormap(cmap)
  print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_c%d_s%d_f%d_v%d.jpg',str,icontr,isubj,ifoi,v))

end

%%
% ---------------------------------
% CONTRAST TASK VS REST
% ---------------------------------
str = 'cvar';

for  ifoi = 1 : 5

if strcmp(str,'var')
  d2 = double(nanmean(nanmean(var_cnt_all(:,:,:,ifoi),3),2)-nanmean(nanmean(var_res_all(:,:,:,ifoi),3),2));
elseif strcmp(str,'cvar')
  d2 = double(nanmean(nanmean(cvar_cnt_all(:,:,:,ifoi),3),2)-nanmean(nanmean(cvar_res_all(:,:,:,ifoi),3),2));
elseif strcmp(str,'dfa')
  d2 = double(nanmean(nanmean(dfa_cnt_all(:,:,:,ifoi),3),2)-nanmean(nanmean(dfa_res_all(:,:,:,ifoi),3),2));
end

load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_stat_clusterstat_task-rest_f%d_v%d.mat',str,ifoi,v));
senssel = find(stats.mask);

pars            = [];
pars.markersel  = senssel;
pars.markersize = 20;
pars.linewidth  = 9;
pars.resolution = 300;
pars.cbar       = 0;

r = max([abs(min(d2)) abs(max(d2))]);
pars.scale = [-r r];

showfield_colormap(d2,sa.locs_2D,pars);

print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_%s_topo_task-rest_stats_f%d_v%d.jpg',str,ifoi,v))

% ---------------------------------

close all

end








%% OTHER STUFF



% 
% 
% % HIGH VS LOW COUNT
% idx_hi=nanmean(cnt,2)>prctile(nanmean(cnt,2),50);
% idx_lo=nanmean(cnt,2)<prctile(nanmean(cnt,2),50);
% 
% d = nanmean(var_cnt_all(:,idx_hi,:),3)-nanmean(var_cnt_all(:,idx_lo,:),3);
% d = nanmean(d,2);
% 
% figure; set(gcf,'color','white');
% 
% pars.markersize = 0;
% pars.linewidth = 9;
% pars.resolution = 300;
% % pars.scale = [-0.05 0.05];
% %   pars.senssel = senssel;
% 
% showfield_colormap(d,sa.locs_2D,pars);
% 
% print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_topo_hicnt_vs_lowvnt_f%d_v%d.jpg',ifoi,v))
% 
% %% CORRELATE EACH SENSORY WITH COUNTS
% clear r
% 
% all_cnt = nanmean(cnt,2);
% all_dfa = nanmean(var_cnt_all,3);
% 
% for ichan = 1 : size(var_cnt_all,1)
%   
%   tmp = corrcoef(all_dfa(ichan,:),all_cnt');
%   r(ichan) = tmp(1,2);
%   
% end
% 
% figure; set(gcf,'color','white');
% 
% pars.markersize = 0;
% pars.linewidth = 9;
% pars.resolution = 300;
% pars.scale = [-0.75 0.75];
% %   pars.senssel = senssel;
% 
% showfield_colormap(r,sa.locs_2D,pars);
% 
% print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_topo_corr_f%d_v%d.jpg',ifoi,v))
% 
% 
% %%
% %
% for isubj = 1:18
%   figure; set(gcf,'color','white');
%   pars.cbar = 0;
% 
%   % PHARMA CONTRAST DURING CNT
%   d = nanmean(var_cnt_all(:,isubj,2)-var_cnt_all(:,isubj,1),2);
% 
%   range = max([abs(min(d)) abs(max(d))]);
% 
%   pars.scale= [-range range];
%   pars.markersize = 0;
%   pars.linewidth = 9;
%   pars.resolution = 300;
% 
%   showfield_colormap(double(d),sa.locs_2D,pars);
% %   colormap(winter)
%   print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_dfa_topo_a-p_s%d_f%d_v%d.jpg',isubj,ifoi,v))
% end