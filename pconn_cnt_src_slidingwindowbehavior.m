%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_cnt_src_dfa


clear

tp_addpaths

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';

for v = [2]
  
  if v==2
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v         = 2;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2 12; 8 12; 12 46; 54 100;2 8];
    i_fit     = [3 50];
    dfa_overlap = 0.5;
    gridsize  = 'cortex';
    filt      = 'eloreta';
  elseif  v==22
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v         = 2;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2 12; 8 12; 12 46; 54 100; 2 8];
    i_fit     = [3 50];
    dfa_overlap = 0.5;
    gridsize  = 'cortex';
    filt      = 'lcmv';
  end
 
  
  %%
  
  for ifoi = 1 : 5
    for m = 1 : 3
      for isubj = SUBJLIST
        %
        if ~exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          continue
        end
        
        disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
        
        d=dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
        
        if length(d)==1
          if v < 10
            blocks = str2num(d.name(end-7));
          else
            blocks = str2num(d.name(end-8));
          end
        elseif length(d) == 2
          blocks = [1 2];
        end
        
        for iblock = blocks
          
          disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
          
          if length(d) > 1
            load(['/home/tpfeffer/pconn_cnt/proc/preproc/' d(iblock).name])
          else
            load(['/home/tpfeffer/pconn_cnt/proc/preproc/' d.name])
          end
          
          clear data_hi
          [mydata,epleng] = megdata2mydata(data);
          clear data
          
          % ------------------------------------------
          % COMPUTE DFA IN SOURCE SPACE
          % ------------------------------------------
          
          clear L grid
          
          
          if strcmp(filt,'lcmv')
            if strcmp(gridsize,'xcoarse')
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
              [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
            elseif strcmp(gridsize,'coarse')
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
              [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
            elseif strcmp(gridsize,'cortex')
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
              [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
              A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
            end
            
          elseif strcmp(filt,'eloreta')
            if strcmp(gridsize,'xcoarse')
              % --------------------------------------------
              % XCOARSE GRID
              % --------------------------------------------
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
              A = mkfilt_eloreta_v2(sa.L_coarse);
              
              load(sprintf(['~/pconn_cnt/proc/src/' 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,1));
              
              for ivox = 1 : size(A,2)
                
                Aloc       = squeeze(A(:,ivox,:));
                [u s vv]   = svd(real(cs_src_all(:,:,ivox)));
                A1(:,ivox) = Aloc*u(:,1);
                
              end
              
              clear A
              
            elseif strcmp(gridsize,'coarse')
              % --------------------------------------------
              % COARSE GRID
              % --------------------------------------------
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
              A = mkfilt_eloreta_v2(sa.L_coarse);
              
              load(sprintf(['~/pconn_cnt/proc/src/' 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,1));
              
              for ivox = 1 : size(A,2)
                
                Aloc       = squeeze(A(:,ivox,:));
                [u s vv]   = svd(real(cs_src_all(:,:,ivox)));
                A1(:,ivox) = Aloc*u(:,1);
                
              end
              
              clear A
            elseif strcmp(gridsize,'cortex')
              % --------------------------------------------
              % CORTEX GRID
              % --------------------------------------------
              
              load([indir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
              A = mkfilt_eloreta_v2(sa.L_coarse);
              
              load(sprintf(['~/pconn_cnt/proc/src/' 'pconn_cnt_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,1));
              
              for ivox = 1 : size(A,2)
                
                Aloc       = squeeze(A(:,ivox,:));
                [u s vv]   = svd(real(cs_src_all(:,:,ivox)));
                A1(:,ivox) = Aloc*u(:,1);
                
              end
              
              clear A
            end
          end
          
          if length(d) == 1
            if iblock == 1
              par.dfa(1:size(A1,2),2)  = nan(3000,1);
              par.amp(1:size(A1,2),2)  = nan(3000,1);
              par.var(1:size(A1,2),2)  = nan(3000,1);
              par.cvar(1:size(A1,2),2) = nan(3000,1);
            else
              par.dfa(1:size(A1,2),1)  = nan(3000,1);
              par.amp(1:size(A1,2),1)  = nan(3000,1);
              par.var(1:size(A1,2),1)  = nan(3000,1);
              par.cvar(1:size(A1,2),1) = nan(3000,1);
            end
          end
          
          clear cs_src1 cs_src2 cs_src3 cs_src cs_avg cs1 cs2 cs3 f A2 A3 A
          
          siginfo = nbt_Info;
          siginfo.converted_sample_frequency = fsample;
          
          % compute bp-filtered signal
          ampenv = single(nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,2/foi(ifoi,1)));
          ampenv = hilbert(ampenv);
          
          clear mydata data data_hi
          % project bp-filtered signal into source space
          ampenv=ampenv*A1; clear A1
          
          for ivox = 1 : size(ampenv,2)
            fprintf('vox%d...\n',ivox)
            tmp=abs(ampenv(:,ivox));
            %           ampenv1(:,ivox) = resample(tmp,25,8);
            %             [par.acorr_acf(:,ivox),par.acorr_lags] = autocorr(tmp,10*400);
            tmp_dfa  = tp_dfa(tmp,i_fit,400,dfa_overlap,15);
            par.dfa(ivox,iblock) = tmp_dfa.exp;
            % FILTER RESAMPLED SIGNAL
            par.amp(ivox,iblock)  = nanmean(tmp);
            par.var(ivox,iblock)  = nanvar(tmp);
            par.cvar(ivox,iblock) = nanstd(tmp)./nanmean(tmp);
            
            clear tmp
          end
          
          
          clear ampenv ampenv1
          
        end
        
        save(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
        clear expo tmp
        
      end
    end
  end
end
error('STOP')

%% CLEAN NON PROCESSED FILES
outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';

cnt = 0;
v = 2;
cnt_exist = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1:5
      ifoi
      if exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v)) && exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        cnt_exist = cnt_exist + 1;
        
        continue
      elseif exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        
      elseif exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && ~exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v))
        warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
        delete(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        cnt = cnt + 1;
      elseif ~exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v)) && exist(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_src_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        warning('Nothing exists')
        cnt = cnt+1;
      end
      
    end
  end
end
cnt