%% COMPUTE DFA SENSOR LEVEL TOPO STATISTICS
% pconn_cnt_sens_var_clusterstat

% Implement statistical test as described in Nichols & Holmes, 2001

% (1) Single threshold permutation test
% (2) Cluster-based permutation test

% tpfeffer | thms.pfffr@gmail.com | 05-05-15

clear


% --------------------------------------------------------
% --------------------------------------------------------

str           = 'cvar';
v_cnt         = 19;
v_res         = 21;
SUBJLIST      = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

% --------------------------------------------------------
% --------------------------------------------------------

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/pconn/matlab
addpath ~/pcbi/
outdir    = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir   = '/home/tpfeffer/pconn_cnt/proc/plots/';


run  ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
load ~/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
load ~/pconn/matlab/sensorselection.mat

ord       = pconn_randomization;

%% PERFORM CLUSTER BASED PERMUITATON TEST

% DEFINE VARIABLE OF INTEREST
fprintf('Processing %s ...\n',str)

for ifoi = 1:4
  
  fprintf('Processing freq %d ...\n',ifoi)
  
  clear im tmp dfa var_cnt_all var_res_all
  
  for isubj = SUBJLIST
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      tmp_cnt = pcbi_cnt(isubj);
      
      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
      
      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
      
      var_cnt_all(:,isubj,m,ifoi)  = double(nanmean(par.var,2)); 
      cvar_cnt_all(:,isubj,m,ifoi) = double(nanmean(par.cvar,2)); 
      dfa_cnt_all(:,isubj,m,ifoi)  = double(nanmean(par.dfa,2)); clear par
      
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
      
      var_res_all(:,isubj,m,ifoi)  = double(nanmean(par.var,2));
      cvar_res_all(:,isubj,m,ifoi) = double(nanmean(par.cvar,2));
      dfa_res_all(:,isubj,m,ifoi)  = double(nanmean(par.dfa,2)); clear par
      
    end
  end
  
  cnt = cnt(SUBJLIST,:);
  
  var_cnt_all = var_cnt_all(:,SUBJLIST,:,:);
  var_res_all = var_res_all(:,SUBJLIST,:,:);
  
  cvar_cnt_all = cvar_cnt_all(:,SUBJLIST,:,:);
  cvar_res_all = cvar_res_all(:,SUBJLIST,:,:);
  
  dfa_cnt_all = dfa_cnt_all(:,SUBJLIST,:,:);
  dfa_res_all = dfa_res_all(:,SUBJLIST,:,:);
  
  
  %% COMPARE PHARACOLOGICAL CONDITIONS
  
	if strcmp(str,'var')
    d1 = var_cnt_all(:,:,2,ifoi);
    d2 = var_cnt_all(:,:,1,ifoi);
  elseif strcmp(str,'cvar')
    d1 = cvar_cnt_all(:,:,2,ifoi);
    d2 = cvar_cnt_all(:,:,1,ifoi);
  elseif strcmp(str,'dfa')
    d1 = dfa_cnt_all(:,:,2,ifoi);
    d2 = dfa_cnt_all(:,:,1,ifoi);
  end

  load '/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s4_m1_b1_v4.mat'
  clear data_hi
  % -------------------------------------------------------
  % GET RID OF MISSING SENSORS
  % -------------------------------------------------------
  cfg             = [];
  cfg.method      = 'template';
  cfg.layout      = 'CTF275';
  n       = ft_prepare_neighbours(cfg);
  
  FOI       = [1 2 3 4 5];
  contrast  = [2 1; 3 1; 2 3];
  
  for icontr = 1 : 3
    
    dat     = [d1 d2];
    
    data_low.dimord                 = 'subj_chan_freq_time';
    data_low.time                   = [1 2];
    data_low.freq                   = [1 2];
    data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
    
    cfg                   = [];
    cfg.channel           = 'all';
    cfg.latency           = [1 1];
    cfg.frequency         = [1 1];
    cfg.method            = 'montecarlo';
    cfg.statistic         = 'depsamplesT';
    cfg.computeprob       = 'yes';
    cfg.correctm          = 'cluster';
    cfg.clusteralpha      = 0.025;
    cfg.clusterstatistic  = 'maxsum';
    cfg.clustertail       = 0; %1 = right
    cfg.tail              = 0; %1 = right
    cfg.alpha             = 0.025;
    cfg.minnbchan         = 2;
    cfg.numrandomization  = 10000;
    cfg.avgovertime       = 'yes';
    cfg.avgoverfreq       = 'yes';
    
    % specifies with which sensors other sensors can form clusters
    cfg_neighb.method    	= 'template';
    cfg_neighb.template  	= 'CTF275_neighb.mat';
    cfg_neighb.feedback 	= 'no';
    cfg.neighbours        = n;
    
    n_subj = length(SUBJLIST);
    design = zeros(2,2*n_subj);
    design(1,:) = repmat(1:n_subj,1,2);
    design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
    
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    
    [stats] = ft_freqstatistics(cfg, data_low);
    
    save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,FOI(ifoi),v_cnt),'stats');
    
  end
  
  %% COMPARE TASK VS REST
  
  if strcmp(str,'var')
    d1 = squeeze(nanmean(var_cnt_all(:,:,:,ifoi),3));
    d2 = squeeze(nanmean(var_res_all(:,:,:,ifoi),3));
  elseif strcmp(str,'cvar')
    d1 = squeeze(nanmean(cvar_cnt_all(:,:,:,ifoi),3));
    d2 = squeeze(nanmean(cvar_res_all(:,:,:,ifoi),3));
  elseif strcmp(str,'dfa')
    d1 = squeeze(nanmean(dfa_cnt_all(:,:,:,ifoi),3));
    d2 = squeeze(nanmean(dfa_res_all(:,:,:,ifoi),3));
  end
  
  dat     = [d1 d2];
  
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
    
  save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_%s_stat_clusterstat_task-rest_f%d_v%d.mat',str,FOI(ifoi),v_cnt),'stats');
  
end
  

%% OTHER CODE BEHAVIOR MOSTLY







%   %% COMPARE ATOMOX(TASK-REST) VS. PLACEBO(TASK-REST)
%   
%   dat     = [var_cnt_all(:,:,2)-var_res_all(:,:,2) var_cnt_all(:,:,1)-var_res_all(:,:,1)];
%   
%   data_low.dimord                  = 'subj_chan_freq_time';
%   data_low.time                   = [1 2];
%   data_low.freq                   = [1 2];
%   data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
%   
%   cfg                  = [];
%   cfg.channel          = 'all';
%   cfg.latency          = [1 1];
%   cfg.frequency        = [1 1];
%   cfg.method           = 'montecarlo';
%   cfg.statistic        = 'depsamplesT';
%   cfg.computeprob      = 'yes';
%   cfg.correctm         = 'cluster';
%   cfg.clusteralpha     = 0.025;
%   cfg.clusterstatistic = 'maxsum';
%   cfg.clustertail      = 0; %1 = right
%   cfg.tail             = 0; %1 = right
%   cfg.alpha            = 0.025;
%   cfg.minnbchan        = 2;
%   cfg.numrandomization = 10000;
%   cfg.avgovertime      = 'yes';
%   cfg.avgoverfreq      = 'yes';
%   
%   % specifies with which sensors other sensors can form clusters
%   cfg_neighb.method           = 'template';
%   cfg_neighb.template         = 'CTF275_neighb.mat';
%   cfg_neighb.feedback         = 'no';
%   cfg.neighbours = n;
%   
%   n_subj = length(SUBJLIST);
%   design = zeros(2,2*n_subj);
%   design(1,:) = repmat(1:n_subj,1,2);
%   design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
%   
%   cfg.design   = design;
%   cfg.uvar     = 1;
%   cfg.ivar     = 2;
%   
%   [stats] = ft_freqstatistics(cfg, data_low);
  
%   dfa2plot = nanmean(var_cnt_all(:,:,2)-var_res_all(:,:,2),2)-nanmean(var_cnt_all(:,:,1)-var_res_all(:,:,1),2);%.*stats.mask;
%   %
%   pars            = [];
%   pars.cbar       = 0;
%   % pars.scale      = [-0.05 0.05];
%   pars.markersize = 0;
%   pars.linewidth  = 9;
%   pars.resolution = 300;
%   
%   showfield_colormap(dfa2plot,sa.locs_2D,pars);
%   
%   print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_topo_a-p-cnt-rest_raw_f%d_v%d.jpg',ifoi,v_cnt))
%   
%   dfa2plot = dfa2plot.*stats.mask;
%   
%   showfield_colormap(dfa2plot,sa.locs_2D,pars);
%   
%   print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_topo_a-p-cnt-rest_f%d_v%d.jpg',ifoi,v_cnt))
%   
%   save(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_var_clusterstat_a-p-cnt-rest_f%d_v%d.mat',FOI(ifoi),v_cnt),'stats');
%   


%% BEHAVIOR

% 
% ifoi = 1;
% icontr = 1;
% 
% clear im tmp dfa var_cnt_all var_res_all
% 
% for isubj = SUBJLIST
%   for m = 1 : 3
%     
%     im = find(ord(isubj,:)==m);
%     
%     tmp_cnt = pcbi_cnt(isubj);
%     
%     cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
%     
%     load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
%     
%     dfa_cnt_all(:,isubj,im) = nanmean(par.dfa,2);
%     var_cnt_all(:,isubj,im) = nanmean(par.var,2);
%     
%     clear par
%     
%     load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
%     
%     dfa_res_all(:,isubj,im) = nanmean(par.dfa,2);
%     var_res_all(:,isubj,im) = nanmean(par.var,2);
%     
%   end
% end
% 
% % load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_var_stat_clusterstat_contr%d_f%d_v%d.mat',icontr,ifoi,v));
% 
% cnt = cnt(SUBJLIST,:);
% var_cnt_all = var_cnt_all(:,SUBJLIST,:);
% var_res_all = var_res_all(:,SUBJLIST,:);
% dfa_cnt_all = dfa_cnt_all(:,SUBJLIST,:);
% dfa_res_all = dfa_res_all(:,SUBJLIST,:);
% %%
% 
% load(sprintf('~/pconn_bttn/proc/pconn_bttn_mediandur_v%d.mat',4),'dur');
% 
% clc
% 
% % WHOLE BRAIN CORRELATION
% 
% a     = squeeze(nanmean(var_cnt_all,1));
% slp   = pconn_regress(a(:)',cnt(:));
% [r,p] = corrcoef(a(:),cnt(:));
% 
% figure; set(gcf,'color','white'); hold on
% line([min(a(:)) max(a(:))],[slp(2)*min(a(:))+slp(1) slp(2)*max(a(:))+slp(1)],'linewidth',7.5)
% scatter(a(:),cnt(:),200,'facebolor','r','markeredgecolor','w');
% set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025],'fontsize',24);
% title(sprintf('Whole-brain: r = %.2f | p = %.2f',r(1,2),p(1,2)))
% xlabel('Difference DFA'); ylabel('Difference switches');
% 
% axis([min(a(:))-0.1*min(a(:)) max(a(:))+0.1*max(a(:)) min(cnt(:))-2*min(cnt(:)) max(cnt(:))+0.2*max(cnt(:))]);
% 
% print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_behav_wb_f%d_v%d.eps',ifoi,v_cnt))
% 
% if sum(stats.mask)==0
%   error('No significant clusters found!');
% end
% 
% % CHANGE OF SWITCHES AFTER ATOMOX
% 
% % d_dfa = squeeze(nanmean(var_cnt_all(:,:,2),1)'-nanmean(var_cnt_all(:,:,1),1)');
% d_dfa = squeeze(nanmean(dfa_cnt_all(:,:,2),1)'-nanmean(dfa_cnt_all(:,:,1),1)');
% 
% % d_cnt = cnt(:,2)-cnt(:,1);
% d_cnt = [dur(2,:)-dur(1,:)]';
% slp   = pconn_regress(d_dfa',d_cnt);
% 
% figure; set(gcf,'color','white'); hold on
% line([min(d_dfa) max(d_dfa)],[slp(2)*min(d_dfa)+slp(1) slp(2)*max(d_dfa)+slp(1)],'linewidth',7.5)
% scatter(d_dfa,d_cnt,200,'facebolor','r','markeredgecolor','w')
% set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025]);
% axis([-0.05 0.15 -35 55]);
% 
% xlabel('Difference DFA');ylabel('Difference switchs');
% 
% [r,p]=corrcoef(d_dfa,d_cnt);
% title(sprintf('Cluster: r = %.2f | p = %.2f',r(1,2),p(1,2)))
% 
% print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_var_behav_pharmachange_f%d_c%d_v%d.eps',ifoi,icontr,v_cnt))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 



