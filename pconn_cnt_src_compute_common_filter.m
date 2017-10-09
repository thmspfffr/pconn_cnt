%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_cnt_src_compute_common_filter

clear 

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% -----------------------------------------------------
% VERSiON 1
% -----------------------------------------------------
v       = 1;
FOI     = [2 12; 8 12; 12 46; 54 100; 2 4; 1 41];
smoo    = [4; 2; 17; 23; 1; 20];
v_filt  = 1;
meth    = 'eloreta';
% -----------------------------------------------------
% VERSiON 2
% -----------------------------------------------------
% v       = 2;
% FOI     = [2 12; 8 12; 12 46; 54 100; 2 4; 1 41];
% smoo    = [4; 2; 17; 23; 1];
% v_filt  = 2;
% meth    = 'lcmv';
% -----------------------------------------------------

tp_addpaths

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/src/';

%%

for ifoi = 1 : length(FOI)
  for isubj = SUBJLIST
    
    if ~exist(sprintf([outdir 'pconn_cnt_src_compute_common_filter_s%d_f%d_v%d_processing.txt'],isubj,ifoi,v))
      system(['touch ' outdir sprintf('pconn_cnt_src_compute_common_filter_s%d_f%d_v%d_processing.txt',isubj,ifoi,v)]);
    else
      continue
    end
    
    clear csd
    
    for m = 1 : 3
      
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
      
      clear d
      
      for iblock = blocks
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
        
        cfg = [];
        cfg.length =  1/0.2;
        cfg.overlap = 0;
        
        data      = ft_redefinetrial(cfg,data);
        
        cfg               = [];
        cfg.method        = 'mtmfft';
        cfg.output        = 'powandcsd';
        cfg.taper         = 'dpss';
        cfg.channel       = 'MEG';
        cfg.pad           = 'nextpow2';
        cfg.foi           = mean(FOI(ifoi,:));
        cfg.keeptrials    = 'no';
        cfg.tapsmofrq     = smoo(ifoi);
        
        [~,  csd{m}{iblock}] = ft_freqanalysis(cfg, data); clear data
        
        load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
        
        % watch out! L_coarse is correct even for cortex grid
        if strcmp(meth,'eloreta')
          A{m}{iblock} = mkfilt_eloreta_v2(sa.L_coarse);
        elseif strcmp(meth,'lcmv')
          if ~(m == 3 && iblock == 2)
            continue
          else
            save(sprintf([outdir 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,v),'csd','-v7.3');
          end
        elseif strcmp(meth,'dics')
          para.iscs = 1;
          A{m}{iblock} = pconn_beamformer(csd,sa.L_coarse,para);
        end        
      end
      
      if strcmp(meth,'lcmv')
        
        continue
        
      end
      
      if length(blocks) > 1
        
        for ivox = 1 : size(A{m}{1},2)
          
          Aloc1               = squeeze(A{m}{1}(:,ivox,:));
          Aloc2               = squeeze(A{m}{2}(:,ivox,:));
          cs_src1             = Aloc1'*real(csd{m}{1})*Aloc1;
          cs_src2             = Aloc2'*real(csd{m}{2})*Aloc2;
          cs_src{m}(ivox,:,:) = (cs_src1+cs_src2)./2;
          
        end
        
      else
        
        for ivox = 1 : size(A{m}{blocks},2)
          
          Aloc1               = squeeze(A{m}{blocks}(:,ivox,:));
          cs_src1             = Aloc1'*real(csd{m}{blocks})*Aloc1;
          cs_src{m}(ivox,:,:) = cs_src1;
          
        end        
      end
    end
    
    if ~strcmp(meth,'lcmv')
      for ivox = 1 : 3000     
        cs_src_all(:,:,ivox)  = squeeze((cs_src{1}(ivox,:,:)+cs_src{2}(ivox,:,:)+cs_src{3}(ivox,:,:))./3);     
      end
      
      save(sprintf([outdir 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,v),'cs_src_all','-v7.3');
    else
      save(sprintf([outdir 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,v),'csd','-v7.3');
    end
  end
  
  clear cs_src1 cs_src2 Aloc1 Aloc2 csd A
  
end

error('Done!');

%% PLOT CERTAIN THINGS
m = 1;

load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);

if strcmp(meth,'eloreta')
  A = mkfilt_eloreta_v2(sa.L_coarse);
elseif strcmp(meth,'lcmv')
  A = mkfilt_lcmv(sa.L_coarse,csd);
elseif strcmp(meth,'dics')
  para.iscs = 1;
  A = pconn_beamformer(csd,sa.L_coarse,para);
end

for ivox = 1 : size(A,2)
  
  Aloc = squeeze(A(:,ivox,:));
  [u s vv]=svd(cs_src_all(:,:,ivox));
  A1(:,ivox)=Aloc*u(:,1);
  
end

load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
clear csd

cfg = [];
cfg.length =  1/0.2;
cfg.overlap = 0;

data      = ft_redefinetrial(cfg,data);

cfg               = [];
cfg.method        = 'mtmfft';
cfg.output        = 'powandcsd';
cfg.taper         = 'dpss';
cfg.channel       = 'MEG';
cfg.pad           = 'nextpow2';
cfg.foi           = mean(FOI(ifoi,:));
cfg.keeptrials    = 'no';
cfg.tapsmofrq     = smoo(ifoi);

[~,  csd] = ft_freqanalysis(cfg, data); clear data

pow               = diag(A1'*real(csd)*A1);
par_interp        = log10(spatfiltergauss(pow,g1,dd,g2));
para              = [] ;
para.colorlimits  = [min(par_interp) max(par_interp)];

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para); colormap(hot);

clear cs_src1 cs_src2 cs_src3 cs_src cs_avg cs1 cs2 cs3 f A2 A3 A csd

%           par.pow(:,iblock) = diag(A1'*real(cs)*A1);
%         end

% save(sprintf([outdir 'pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
%         clear expo tmp
%
%       end
%     end
%   end

% error('STOP')

% CLEAN NON PROCESSED FILES




