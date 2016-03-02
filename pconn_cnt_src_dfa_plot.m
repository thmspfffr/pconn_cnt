%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 3;
v_rawdata = 2;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
foi       = [2 4; 4 8; 8 12; 12 24];
i_fit     = [3 28];
i_calc    = [.8 40];
gridsize  = 'coarse';
% --------------------------------------------------------




restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
addpath ~/pconn/matlab/


ft_defaults

indir     = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir    = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir   = '/home/tpfeffer/pconn_cnt/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

%% PLOT DFA
   
str = 'dfa'; 
ord   = pconn_randomization;

for ifoi = 1 : 4
  
cnt = 0;
for isubj = SUBJLIST
  
   d = dir(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m*_f%d_v%d.mat'],isubj,ifoi,v));
    if length(d) < 3
      warning('continue')
      continue
    end
    cnt = cnt + 1; disp(cnt);
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);

    load(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
   
    dfa_all(:,im,cnt,ifoi)  = nanmean(par.dfa,2);
    var_all(:,im,cnt,ifoi)  = nanmean(par.var,2);
   	cvar_all(:,im,cnt,ifoi) = nanmean(par.cvar,2); clear par
    
    load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
%    
    dfa_all_res(:,im,cnt,ifoi)  = nanmean(par.dfa,2);
    var_all_res(:,im,cnt,ifoi)  = nanmean(par.var,2);
   	cvar_all_res(:,im,cnt,ifoi) = nanmean(par.cvar,2); clear par
    

  end
end

end
clear d


%%

ifoi = 5;

if strcmp(str,'dfa')
  d_all = dfa_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  d_all = var_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  d_all = cvar_all(:,:,:,ifoi);
end


load sa_meg_template;

if strcmp(gridsize,'coarse')
  grid  = sa_meg_template.grid_coarse;
elseif strcmp(gridsize,'cortex')
  grid  = sa_meg_template.grid_cortex3000;
elseif strcmp(gridsize,'xcoarse')
  grid  = sa_meg_template.grid_xcoarse;
end

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
%  GENERATE SOURCE PLOT

para                  = [];
para.orientation      = 'axial';
para.colormaps        = {'jet'};
% para.colorlimits = [0.5 0.7];
h = figure; hold on
set(h,'color',' k');

% ------------------------
% PLACEBO - ATOMOXETINE
% ------------------------
[~,p] = ttest(d_all(:,2,:),d_all(:,1,:),'dim',3); m = p<0.001;
d = squeeze(nanmean(d_all(:,2,:),3))-squeeze(nanmean(d_all(:,1,:),3));
r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [-r r];
showmri_transp_v3(mri,para,[grid d]);



print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_a-p_f%d_v%d.jpg',ifoi,v)) 

% ------------------------
% PLACEBO - DONEPEZIL
% ------------------------
d = squeeze(nanmean(d_all(:,3,:),3))-squeeze(nanmean(d_all(:,1,:),3));
r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [-r r];
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_d-p_f%d_v%d.jpg',ifoi,v)) 

% ------------------------
% ATOMOX - DONEPEZIL
% ------------------------
d = squeeze(nanmean(d_all(:,2,:),3))-squeeze(nanmean(d_all(:,3,:),3));
r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [-r r];
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_a-d_f%d_v%d.jpg',ifoi,v)) 
% ------------------------

% ------------------------
% ALL CONDITIONS
% ------------------------
d = squeeze(nanmean(nanmean(d_all,3),2));
% r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [min(d) max(d)];
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_all_f%d_v%d.jpg',ifoi,v)) 
% ------------------------

% ------------------------
% PLACEBO
% ------------------------

d = squeeze(nanmean(d_all(:,1,:),3));
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_m1_f%d_v%d.jpg',ifoi,v)) 
% ------------------------

% ------------------------
% ATOMOXETINE
% ------------------------
d = squeeze(nanmean(d_all(:,2,:),3));
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_m2_f%d_v%d.jpg',ifoi,v)) 
% ------------------------

% ------------------------
% DONEPEZIL
% ------------------------
d = squeeze(nanmean(d_all(:,3,:),3));
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfa_m3_f%d_v%d.jpg',ifoi,v)) 
% ------------------------


%% TASK VS REST

str = 'dfa';
ifoi = 4;

if strcmp(str,'dfa')
  d1 = squeeze(nanmean(dfa_all(:,:,:,ifoi),2));
  d2 = squeeze(nanmean(dfa_all_res(:,:,:,ifoi),2));
elseif strcmp(str,'var')
  d1 = squeeze(nanmean(var_all(:,:,:,ifoi),2));
  d2 = squeeze(nanmean(var_all_res(:,:,:,ifoi),2));
elseif strcmp(str,'cvar')
  d1 = squeeze(nanmean(cvar_all(:,:,:,ifoi),2));
  d2 = squeeze(nanmean(cvar_all_res(:,:,:,ifoi),2));
end

load sa_meg_template;

if strcmp(gridsize,'coarse')
  grid  = sa_meg_template.grid_coarse;
elseif strcmp(gridsize,'cortex')
  grid  = sa_meg_template.grid_cortex3000;
elseif strcmp(gridsize,'xcoarse')
  grid  = sa_meg_template.grid_xcoarse;
end

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
%  GENERATE SOURCE PLOT

para                  = [];
para.orientation      = 'axial';
para.colormaps        = {'jet'};
% para.colorlimits = [0.5 0.7];
h = figure; hold on
set(h,'color',' k');

% ------------------------
% TASK - REST
% ------------------------
[~,p] = ttest(d1,d2,'dim',2); m = p<0.01;
d = squeeze(nanmean(d1,2))-squeeze(nanmean(d2,2));
r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [-r r];
showmri_transp_v3(mri,para,[grid d]);

print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_%s_task-rest_f%d_v%d.jpg',str,ifoi,v)) 

% ------------------------
% TASK
% ------------------------

d = squeeze(nanmean(d1,2));
para.colorlimits = [min(d) max(d)];
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_%s_task_f%d_v%d.jpg',str,ifoi,v)) 
% ------------------------

% ------------------------
% REST
% ------------------------
d = squeeze(nanmean(d2,2));
para.colorlimits = [min(d) max(d)];
pconn_showsrc(mri,para,[grid d]);
print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_cnt_src_%s_rest_f%d_v%d.jpg',str,ifoi,v)) 
% ------------------------


%% COMPARE ALL MEASURES

ifoi = 4;

p_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
p_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
p_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);

r(1) = corr(p_dfa,p_var);
r(2) = corr(p_dfa,p_cvar);
r(3) = corr(p_cvar,p_var);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Variance')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('Variance'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))

close 

p_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
p_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
p_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);

r(1) = corr(p_dfa,p_var);
r(2) = corr(p_dfa,p_cvar);
r(3) = corr(p_cvar,p_var);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Variance')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('Variance'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))

close 



%%  CORRELATION OF MEASURES ACROSS CONDITIONS

ifoi = 1;

p1_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
p1_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
p1_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);

p2_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
p2_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
p2_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);

r(1) = corr(p1_dfa,p2_dfa);
r(2) = corr(p1_var,p2_var);
r(3) = corr(p1_cvar,p2_cvar);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p1_dfa,p2_dfa,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA (rest)'); ylabel('DFA (task)')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(log10(p1_var),log10(p2_var),20,'facecolor','k','markeredgecolor','w')
xlabel('Variance (rest)'); ylabel('Variance (task)')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p1_cvar,p2_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('C. Variance (rest)'); ylabel('C. Variance (task)')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_resttaskcorr_f%d_v%d.jpg',ifoi,v))

close 

%% 

%% SHOW SURFACE

g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;
dd = .01;
m2 = spatfiltergauss(m,g1,dd,g2);

para.colorlimits = [-0.03 0.03];

figure;  pconn_showsurface(a,[],z2(idx))


