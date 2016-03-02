%% FIND CLUSTERS
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
   
clear dfa_all var_all cvar_all

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
   
    dfa_all_res(:,im,cnt,ifoi)  = nanmean(par.dfa,2);
    var_all_res(:,im,cnt,ifoi)  = nanmean(par.var,2);
   	cvar_all_res(:,im,cnt,ifoi) = nanmean(par.cvar,2); clear par
    

  end
end

end
clear d

edi
%% CLUSTER BASED PERMUTATION

str     = 'var';
for icontr = 1 : 3
% icontr  = 1 :;
for ifoi = 1:4
% ifoi    = 4;


load sa_meg_template;



if strcmp(str,'dfa')
  d_all = dfa_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  d_all = var_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  d_all = cvar_all(:,:,:,ifoi);
end

if v == 3 
  grid = sa_meg_template.grid_coarse;
elseif v == 2 || v == 1
  grid = sa_meg_template.grid_cortex3000;
end

n =  get_neighbours(grid);

contrasts = [2 1; 3 1; 2 3]; 

z(:,:,1)=d_all(:,contrasts(icontr,1),:);
z(:,:,2)=d_all(:,contrasts(icontr,2),:);

%
para.method     = 'dependentT';
para.neigh      = n;
para.minneigh   = 1;
para.clusteralpha = 0.025;
para.perm = 0;
clear ss clusters ss_num ss
s = find_cluster_sub(z,para);

% PERMUTATION
c = 'Permutation';
nperm = 2000;
para.perm = 1;
para.nperm = 2000;
for iperm = 1 : nperm
  
  fprintf('Permutation %d ...\n',iperm);

  perm = find_cluster_sub(z,para);
  
  if ~isempty(perm.stat_pos)
    
    for iclust = 1 : size(perm.stat_pos,2)  
      tmp(iclust) = perm.stat_pos{iclust}.stat;
    end
    
    pos_perm(iperm) = max(tmp); clear tmp
    
  else
    pos_perm(iperm) = 0;
  end
  
  if ~isempty(perm.stat_neg)
    
    for iclust = 1 : size(perm.stat_neg,2)
      tmp(iclust) = perm.stat_neg{iclust}.stat;
    end
    neg_perm(iperm) = max(tmp); clear tmp
  else
    neg_perm(iperm) = 0;
  end
   
end

for iclust = 1 : size(s.stat_pos,2)
  s.stat.p_pos(iclust) = sum(pos_perm>s.stat_pos{iclust}.stat)/nperm;
end
for iclust = 1 : size(s.stat_neg,2)
  s.stat.p_neg(iclust) = sum(neg_perm<s.stat_neg{iclust}.stat)/nperm;
end

s.mask = zeros(size(z,1),1);

if ~isempty(s.stat_pos)
  tmp = find(s.stat.p_pos<0.05);
  if ~isempty(tmp)
    for itmp = 1 : length(tmp)
      s.mask(s.stat_pos{itmp}.chan) = 1;
    end
  end
end
if ~isempty(s.stat_neg)
  tmp = find(s.stat.p_neg<0.025);
  if ~isempty(tmp)
    for itmp = 1 : length(tmp)
      s.mask(s.stat_neg{itmp}.chan) = 1;
    end
  end
end
 
%% PLOT SOURCE

% load sa_meg_template;

f = ifoi;
% f = [35];
a = squeeze(nanmean(nanmean(z(:,contrasts(icontr,1),:),3),4));
b = squeeze(nanmean(nanmean(z(:,contrasts(icontr,2),:),3),4));
d = a-b;
% d = 
  


m = zeros(size(z,1),1);
% k=logical(ttest(z(:,:,2),z(:,:,1),'dim',2))

if sum(s.mask)~=0
  m(find(s.mask)) = 1;
end


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
r = max(abs(min(d)),abs(max(d)));
para.colorlimits = [-r r];
showmri_transp_v3(mri,para,[grid d]);

print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_%s_raw_c%d_f%d_v%d.jpg',str,icontr,ifoi,v))

d = d.*m;
para.colorlimits = [-r r];
showmri_transp_v3(mri,para,[grid d]);

print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_%s_mask_c%d_f%d_v%d.jpg',str,icontr,ifoi,v))
end
end
%   end
% end
%%
% 
% g1 = sa_meg_template.grid_cortex3000;
% g2 = sa_meg_template.cortex10K.vc;
% dd = .01;
% m2 = spatfiltergauss(m,g1,dd,g2);
% 
% z2 = spatfiltergauss(d,g1,dd,g2);
% a.tri = sa_meg_template.cortex10K.tri;
% 
% idx = []; 
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,1)>size(a.vc,1))];
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,2)>size(a.vc,1))];
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,3)>size(a.vc,1))];
% 
% a.tri(unique(idx),:)=[];
% 
% clear idx;
% 
% idx = find(a.vc(:,1)>0);
% 
% a.vc = a.vc(idx,:);
% a.normals = a.normals(idx,:);
% a.tri = delaunay(a.vc);
% % a.tri = a.tri(:,2:4);
% % a.tri = get_neighbours(a.vc);
% % z2(z2<0.0001) =0;
% para.colorlimits = [-0.03 0.03];
% % a.vc = 
% figure;  pconn_showsurface(a,[],z2(idx))
% 
% 
% %%
% 
% load ~/pconn_bttn/proc/pconn_bttn_mediandur_v4.mat
% 
% d=nanmean(z(find(m),:,2),1)-nanmean(z(find(m),:,1),1)
% dd= dur(2,:)'-dur(1,:)''
% 
% scatter(nanmean(nanmean(z,3),1),nanmean(dur))
% 



