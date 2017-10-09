%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_cnt_src_pow

% requires execution of pconn_src_compute_common_filter.m first.
% this is where the common dipole orientation is determined.

clear

% -----------------------------------------------------
% VERSiON 1
% -----------------------------------------------------
v       = 1;
FOI     = [2 12; 8 12; 12 46; 54 100; 2 8];
smoo    = [4; 2; 17; 23; 1];
v_filt  = 1;
meth    = 'eloreta';
% -----------------------------------------------------
% VERSiON 2
% -----------------------------------------------------
% v       = 2;
% FOI     = [2 12; 8 12; 12 46; 54 100; 2 8];
% smoo    = [4; 2; 17; 23; 1];
% v_filt  = 2;
% meth    = 'lcmv';
% -----------------------------------------------------


SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

tp_addpaths

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/src/';

%%

for ifoi = 1: length(FOI)
  for m = 1 : 3
    for isubj = SUBJLIST
      
      if ~exist(sprintf([outdir 'pconn_cnt_src_pow_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_src_pow_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      
      clear d
      
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      d=dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
      
      if length(d)==1
        if v < 10
          blocks = str2num(d(1).name(end-7));
        else
          blocks = str2num(d(1).name(end-8));
        end
      elseif length(d) == 2
        blocks = [1 2];
      end
            
      for iblock = blocks
        
        if length(blocks) == 1
          if iblock == 1
            par.pow(:,2) = nan(3000,1); 
          else
            par.pow(:,1) = nan(3000,1); 
          end
        end
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        if strcmp(meth,'eloreta')
          
          if length(d) == 2
            load(['~/pconn_cnt/proc/preproc/' d(iblock).name]);
          else
            load(['~/pconn_cnt/proc/preproc/' d(1).name]);
          end

          cfg = [];
          cfg.length =  1/0.2;
          cfg.overlap = 0;

          data      = ft_redefinetrial(cfg,data);

          cfg               = [];
          cfg.method        = 'mtmfft';
          cfg.output        = 'powandcsd';
          cfg.taper         = 'dpss';
          cfg.pad           = 'nextpow2';
          cfg.foi           = mean(FOI(ifoi,:));
          cfg.keeptrials    = 'no';
          cfg.tapsmofrq     = smoo(ifoi);

          [~,  csd] = ft_freqanalysis(cfg, data); clear data

        end
        % watch out! L_coarse is correct even for cortex grid
        load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
        load(sprintf([outdir 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,v_filt));
  
        A = mkfilt_eloreta_v2(sa.L_coarse);

        if strcmp(meth,'eloreta')
        
          for ivox = 1 : size(A,2)

            Aloc = squeeze(A(:,ivox,:));
            [u s vv]=svd(real(cs_src_all(:,:,ivox)));
            A1(:,ivox)=Aloc*u(:,1);

          end

          par.pow(:,iblock) = diag(A1'*real(csd)*A1);

       	elseif strcmp(meth,'lcmv') || strcmp(meth,'dics')
                   
          % contains the 267 sensors to include
          load('~/pconn/matlab/idx_all.mat')
          
          cs_all_tmp = complex(zeros(274,274));
          cs_all     = complex(zeros(274,274));
          
          for im = 1 : 3
            
            if isubj >= 32
              load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
            elseif isubj < 4
              load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
            elseif isubj == 17
              load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
            else
              load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,im))
            end
            
            cfb = complex(zeros(274,274));
            cnt = 0;
            for ib = 1 : 2
              if ~isempty(csd{im}{ib})
                cnt = cnt + 1;
                cfb(find(idx),find(idx)) = cfb(find(idx),find(idx))+csd{im}{ib};
              end
            end 
            cs_all_tmp = cs_all_tmp + (cfb./cnt);
          end
          
          cs_all = cs_all_tmp./3;
          
          if isubj >= 32
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
          elseif isubj < 4
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
          elseif isubj == 17
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
          else
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
          end
          
          lf = zeros(274,3000,3);
          lf(find(idx),:,:) = sa.L_coarse;
          lf = lf(idx_all,:,:);
          lcmv = tp_lcmv(lf,real(cs_all(idx_all,idx_all)),0.05);
               
          cs(find(idx),find(idx))   = csd{m}{ib};
          cs                        = cs(idx_all,idx_all);
          
          par.pow(:,iblock)         = diag(lcmv.filt'*real(cs)*lcmv.filt);  
          
        end
        
      	clear csd A A1 
        
      end
      
      save(sprintf([outdir 'pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
      
      clear par
      
    end
  end
end

error('STOP')

%% CLEAN NON PROCESSED FILES
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap  = cbrewer('seq', 'YlOrRd', 500,'pchip');
% cmap  = cmap(end:-1:1,:);
% TEST PLOT
load sa_meg_template;

grid  = sa_meg_template.grid_cortex3000;
g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = sa_meg_template.grid_cortex3000;
par_interp = spatfiltergauss(pow,g1,dd,g2);

para = [] ;
para.colorlimits = [min(par_interp) max(par_interp)];

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);

%% PLOT ALL

ifoi = 5
v = 2

ord	= pconn_randomization;
cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    
      im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn_cnt/proc/src/' 'pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
      
      pow_all1(:,isubj,m) = nanmean(par.pow,2); clear par
      
    	load(sprintf(['~/pconn_cnt/proc/src/' 'pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,2));
      pow_all2(:,isubj,m) = nanmean(par.pow,2);

      
  end
end

pow1 = pow_all1(:,SUBJLIST,:);
pow2 = pow_all2(:,SUBJLIST,:);


% d=;
% for isubj = 1 : size(pow,2);
%   d = nanmean(nanmean(pow(:,isubj,:),2),3);


%%
% d = (nanmean(pow(:,:,3),2)-nanmean(pow(:,:,1),2))./(nanmean(pow(:,:,2),2)+nanmean(pow(:,:,1),2));
[t2,p2,~,s]=ttest(pow(:,:,3),pow(:,:,1),'dim',2);
d = s.tstat; d = d.*t2;

para.nperm = 10000;
para.tail = 0;
para.paired=1;
para.neigh = get_neighbours(g1);
para.method = 'dependentT'
para.alpha = 0.025;
para.clusteralpha = 0.05;
para.minneigh = 2;

s = tp_clusterperm(pow(:,:,[2 1]),para)

  par_interp = spatfiltergauss(d,g1,dd,g2);
  
  para = [] ;
  para.colorlimits = [-1.96 1.96];
%   para.colorlimits = [min(par_interp) max(par_interp)];

  % PLOT RESULTS
  tp_showsource(par_interp,cmap,sa_meg_template,para); drawnow
  
  
  %%
for isubj = 1 : 28
    p1 = pow1(:,isubj,2)-pow1(:,isubj,1);
    p2 = pow2(:,isubj,2)-pow2(:,isubj,1);
  
    r(isubj)=corr(p1,p2)
%     par_interp = spatfiltergauss(p,g1,dd,g2);
%     para.colorlimits = [log10(min(par_interp)) log10(max(par_interp))];
%     tp_showsource(log10(par_interp),hot,sa_meg_template,para); drawnow
% print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pow_s%d_f%d_v%d.mat.jpg',isubj,ifoi,v))
end
  
  
