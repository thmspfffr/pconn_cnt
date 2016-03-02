%% COMPUTE DFA SENSOR LEVEL TOPO STATISTICS
% pconn_cnt_sens_dfa_clusterstat

% Implement statistical test as described in Nichols & Holmes, 2001

% (1) Single threshold permutation test
% (2) Cluster-based permutation test

% tpfeffer | thms.pfffr@gmail.com | 05-05-15

clear

% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
v_cnt         = 19;
v_res = 21;
% v_out     = 2;
SUBJLIST =  [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
str = 'cvar';
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v         = 2;
% v_out     = 3;
% SUBJLIST  = [2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
% --------------------------------------------------------

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';

load ~/pconn/matlab/sensorselection.mat
% load /home/tpfeffer/pconn_cnt/proc/src/pconn_sa_s9_m1_b1_v1.mat

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
addpath ~/Documents/MATLAB/plot2svg_20120915/plot2svg_20120915/
addpath ~/pcbi/
%
addpath /home/tpfeffer/pconn/matlab
ord       = pconn_randomization;

subsel = 1:length(SUBJLIST);
clear dfa_all idx

%% CLUSTER-BASED PERMUTATION TEST
% FOR WITHIN SUBJECTS DESIGNS
ord       = pconn_randomization;

% subsel = 1:length(SUBJLIST);
clear dfa_cnt_all dfa_res_all tmp
load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
%
clear dfa_cnt_all dfa_res_all tmp


for ifoi = 1:5
  for isubj = SUBJLIST
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      tmp_cnt = pcbi_cnt(isubj);
      
      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
      
      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
      
      dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2); clear par
      
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
      
      dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      
    end
  end
end

% cnt = cnt(SUBJLIST,:);
dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));

dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
%% COMPARE PHARACOLOGICAL CONDITIONS

if strcmp(str,'dfa')
  par_cnt = dfa_cnt_all;
  par_res = dfa_res_all;
elseif strcmp(str,'var')
  par_cnt = var_cnt_all;
  par_res = var_res_all;
elseif strcmp(str,'cvar')
  par_cnt = cvar_cnt_all;
  par_res = cvar_res_all;
end

load '/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s4_m1_b1_v6.mat'
clear data_hi
% -------------------------------------------------------
% GET RID OF MISSING SENSORS
% -------------------------------------------------------
cfg             = [];
cfg.method      = 'template';
cfg.layout      = 'CTF275';
n       = ft_prepare_neighbours(cfg);


FOI = [1 2 3 4 5 6 7 8 9 10 11];
% figure; set(gcf,'color','white');
contrast = [2 1; 3 1; 2 3];

for ifoi = 1 : 5
  
  for icontr = 1 : 3
    
    dat     = [par_cnt(:,:,contrast(icontr,1),ifoi) par_cnt(:,:,contrast(icontr,2),ifoi)];
    
    data_low.dimord                  = 'subj_chan_freq_time';
    data_low.time                   = [1 2];
    data_low.freq                   = [1 2];
    data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
    
    cfg                  = [];
    cfg.channel          = 'all';
    cfg.latency          = [1 1];
    cfg.frequency        = [1 1];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.computeprob      = 'yes';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.025;
    cfg.clusterstatistic = 'maxsum';
    cfg.clustertail      = 0; %1 = right
    cfg.tail             = 0; %1 = right
    cfg.alpha            = 0.025;
    cfg.minnbchan        = 2;
    cfg.numrandomization = 10000;
    cfg.avgovertime      = 'yes';
    cfg.avgoverfreq      = 'yes';
    
    % specifies with which sensors other sensors can form clusters
    cfg_neighb.method           = 'template';
    cfg_neighb.template         = 'CTF275_neighb.mat';
    cfg_neighb.feedback         = 'no';
    cfg.neighbours = n;
    
    n_subj = length(SUBJLIST);
    design = zeros(2,2*n_subj);
    design(1,:) = repmat(1:n_subj,1,2);
    design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
    
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    
    [stats] = ft_freqstatistics(cfg, data_low);
    
    save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt),'stats');
    
  end
end

%% COMPARE TASK VS REST

for ifoi = 1 : 5
  
  dat     = [nanmean(par_cnt(:,:,:,ifoi),3) nanmean(par_res(:,:,:,ifoi),3)];
  
  data_low.dimord                  = 'subj_chan_freq_time';
  data_low.time                   = [1 2];
  data_low.freq                   = [1 2];
  data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
  
  cfg                  = [];
  cfg.channel          = 'all';
  cfg.latency          = [1 1];
  cfg.frequency        = [1 1];
  cfg.method           = 'montecarlo';
  cfg.statistic        = 'depsamplesT';
  cfg.computeprob      = 'yes';
  cfg.correctm         = 'cluster';
  cfg.clusteralpha     = 0.025;
  cfg.clusterstatistic = 'maxsum';
  cfg.clustertail      = 0; %1 = right
  cfg.tail             = 0; %1 = right
  cfg.alpha            = 0.025;
  cfg.minnbchan        = 2;
  cfg.numrandomization = 10000;
  cfg.avgovertime      = 'yes';
  cfg.avgoverfreq      = 'yes';
  
  % specifies with which sensors other sensors can form clusters
  cfg_neighb.method           = 'template';
  cfg_neighb.template         = 'CTF275_neighb.mat';
  cfg_neighb.feedback         = 'no';
  cfg.neighbours = n;
  
  n_subj = length(SUBJLIST);
  design = zeros(2,2*n_subj);
  design(1,:) = repmat(1:n_subj,1,2);
  design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
  
  cfg.design   = design;
  cfg.uvar     = 1;
  cfg.ivar     = 2;
  
  [stats] = ft_freqstatistics(cfg, data_low);
  
  save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_clusterstat_cnt-rest_f%d_v%d.mat',str,ifoi,v_cnt),'stats');
  
end

%% COMPARE ATOMOX(TASK-REST) VS. PLACEBO(TASK-REST)
%
for ifoi = 1 : 5
dat     = [par_cnt(:,:,2,ifoi)-par_cnt(:,:,1,ifoi) par_res(:,:,2,ifoi)-par_res(:,:,1,ifoi)];

data_low.dimord                  = 'subj_chan_freq_time';
data_low.time                   = [1 2];
data_low.freq                   = [1 2];
data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);

cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [1 1];
cfg.frequency        = [1 1];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; %1 = right
cfg.tail             = 0; %1 = right
cfg.alpha            = 0.025;
cfg.minnbchan        = 2;
cfg.numrandomization = 10000;
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';

% specifies with which sensors other sensors can form clusters
cfg_neighb.method           = 'template';
cfg_neighb.template         = 'CTF275_neighb.mat';
cfg_neighb.feedback         = 'no';
cfg.neighbours = n;

n_subj = length(SUBJLIST);
design = zeros(2,2*n_subj);
design(1,:) = repmat(1:n_subj,1,2);
design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stats] = ft_freqstatistics(cfg, data_low);

  save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_clusterstat_diffdiff_f%d_v%d.mat',str,ifoi,v_cnt),'stats');


end

error('!')


%% BEHAVIOR


ifoi = 5;
contr = 1;

clear im tmp dfa dfa_cnt_all dfa_res_all

for isubj = SUBJLIST
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    
    tmp_cnt = pcbi_cnt(isubj);
    
    cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
    
    for iblock = 1 : 2
      
      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,im,ifoi,v));
      
      tmp(:,iblock) = dfa.MarkerValues; clear dfa
      
    end
    
    dfa_cnt_all(:,isubj,im) = nanmean(tmp,2); clear tmp
    
    load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    dfa_res_all(:,isubj,im) = nanmean(par.dfa,2);
    
  end
end

load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_dfa_stat_clusterstat_contr%d_f%d_v%d.mat',icontr,ifoi,v));

cnt = cnt(SUBJLIST,:);
dfa_cnt_all = dfa_cnt_all(:,SUBJLIST,:);
dfa_res_all = dfa_res_all(:,SUBJLIST,:);

%%


clc

% WHOLE BRAIN CORRELATION

a     = squeeze(nanmean(dfa_cnt_all,1));
slp   = pconn_regress(a(:)',cnt(:));
[r,p] = corrcoef(a(:),cnt(:));

figure; set(gcf,'color','white'); hold on
line([min(a(:)) max(a(:))],[slp(2)*min(a(:))+slp(1) slp(2)*max(a(:))+slp(1)],'linewidth',7.5)
scatter(a(:),cnt(:),200,'facebolor','r','markeredgecolor','w');
set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025],'fontsize',24);
title(sprintf('Whole-brain: r = %.2f | p = %.2f',r(1,2),p(1,2)))
xlabel('Difference DFA'); ylabel('Difference switches');

axis([min(a(:))-0.1*min(a(:)) max(a(:))+0.1*max(a(:)) min(cnt(:))-2*min(cnt(:)) max(cnt(:))+0.2*max(cnt(:))]);

print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_behav_wb_f%d_v%d.eps',ifoi,v))

if sum(stats.mask)==0
  error('No significant clusters found!');
end

% CHANGE OF SWITCHES AFTER ATOMOX

d_dfa = squeeze(nanmean(dfa_cnt_all(stats.mask,:,2),1)'-nanmean(dfa_cnt_all(stats.mask,:,1),1)');
d_cnt = cnt(:,2)-cnt(:,1);
slp   = pconn_regress(d_dfa',d_cnt);

figure; set(gcf,'color','white'); hold on
line([min(d_dfa) max(d_dfa)],[slp(2)*min(d_dfa)+slp(1) slp(2)*max(d_dfa)+slp(1)],'linewidth',7.5)
scatter(d_dfa,d_cnt,200,'facebolor','r','markeredgecolor','w')
set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025]);
axis([-0.05 0.15 -35 55]);

xlabel('Difference DFA');ylabel('Difference switchs');

[r,p]=corrcoef(d_dfa,d_cnt);
title(sprintf('Cluster: r = %.2f | p = %.2f',r(1,2),p(1,2)))

print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_behav_pharmachange_f%d_c%d_v%d.eps',ifoi,icontr,v))













