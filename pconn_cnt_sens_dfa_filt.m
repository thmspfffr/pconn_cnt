%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_cnt_sens_dfa_filt

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.

clear all

for v = [1 2 9]
  if v == 1
    % --------------------------------------------------------
    % VERSION 1
    % --------------------------------------------------------
    v         = 1;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 100];
    i_calc    = [1.5 140];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 2
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v         = 2;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 50];
    i_calc    = [1.5 70];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 3
    % --------------------------------------------------------
    % VERSION 3
    % --------------------------------------------------------
    v         = 3;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 200];
    i_calc    = [1.5 250];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 4
    % --------------------------------------------------------
    % VERSION 4
    % --------------------------------------------------------
    v         = 4;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2:2:148; 4:2:150]';
    i_fit     = [2 100];
    i_calc    = [1.5 130];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 5
    % --------------------------------------------------------
    % VERSION 5
    % --------------------------------------------------------
    v         = 5;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 60];
    i_calc    = [1.5 80];
    dfa_overlap = 0.8;
    filt_ord = 2;
  elseif v == 6
    % --------------------------------------------------------
    % VERSION 6
    % --------------------------------------------------------
    v         = 6;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 70];
    i_calc    = [1.5 85];
    dfa_overlap = 0.8;
    filt_ord = 2;
  elseif v == 7
    % --------------------------------------------------------
    % VERSION 7
    % --------------------------------------------------------
    v         = 7;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 30];
    i_calc    = [1.5 50];
    dfa_overlap = 0.8;
    filt_ord = 2;
    % --------------------------------------------------------
  elseif v == 8
    % --------------------------------------------------------
    % VERSION 7
    % --------------------------------------------------------
    v         = 8;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [2 100];
    i_calc    = [1.5 140];
    dfa_overlap = 0.5;
    filt_ord = 2;
    % --------------------------------------------------------
   elseif v == 9
    % --------------------------------------------------------
    % VERSION 7
    % --------------------------------------------------------
    v         = 9;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
    foi       = [2 4; 4 8; 8 12; 12 36];
    i_fit     = [1 100];
    i_calc    = [0.5 140];
    dfa_overlap = 0.5;
    filt_ord = 2;
    % --------------------------------------------------------
  end
  
  addpath ~/Documents/MATLAB/
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
  addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
  
  indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
  outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
  % mridir = '/home/gnolte/neuconn/mri_data/';
  plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';
  
  addpath ~/pconn/matlab/
  
  run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
  %%
  for m = 1 : 3
    for isubj = SUBJLIST
      for ifoi = 1 : length(foi)
%         
        if ~exist(sprintf([outdir 'pconn_cnt_sens_dfa_filt_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pconn_cnt_sens_dfa_filt_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          continue
        end
        % %
        if isubj == 3
          load ~/pconn/matlab/pconn_sensorlabels.mat
        end
        
        if (isubj~=3 && isubj~=2 &&  isubj~=17) && isubj <= 24
          d = dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
        else
          d = dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
        end
        
        for iblock = 1 : length(d)
          
          disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
          
          load(['/home/tpfeffer/pconn_cnt/proc/preproc/' d(iblock).name])
          
          if isubj == 3 || isubj == 2
            [~,idx_lab]=intersect(data_low.label,lab);
            data_low.trial{1}=data_low.trial{1}(idx_lab,:);
            data_low.label = lab;
            save(['~/pconn/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
          end
          
          clear data_hi
          
          [mydata,epleng] = megdata2mydata(data_low);
          
          mydata = mydata(1:end-1000,:);
          
          siginfo = nbt_Info;
          siginfo.converted_sample_frequency = fsample;
          
          % compute bp-filtered signal
          tmp     = tp_filt(mydata,[foi(ifoi,1) foi(ifoi,2)],fsample,10,'ord',4);
          ampenv  = abs(hilbert(zscore(tmp))); clear tmp mydata
          
          % second filter
          ampenv    = zscore(tp_filt(ampenv,[0.01 2],fsample,10,'ord',1));
          
          [dfa,~,DFA_y_all] = nbt_doDFA(ampenv, siginfo, i_fit,i_calc,dfa_overlap,0,0,[]);

          clear data_low
          
          if isubj >= 32
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
          elseif isubj < 4
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
          elseif isubj == 17
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
          else
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
          end
          
          % INTERPOLATES THE RESULTS
          par.dfa(:,iblock)   = pconn_sens_interp274(idx,dfa.MarkerValues);
          par.amp(:,iblock)   = pconn_sens_interp274(idx,mean(ampenv));
          par.var(:,iblock)   = pconn_sens_interp274(idx,var(ampenv));
          %         par.pow(:,iblock)   = pconn_sens_interp274(idx,nanmean(p,2));
          par.cvar(:,iblock)  = sqrt(par.var(:,iblock))./(par.amp(:,iblock));
          
          clear dfa dat ampenv p
        end
        
        save(sprintf([outdir 'pconn_cnt_sens_dfa_filt_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
        clear par
        
      end
    end
  end
end

error('STOP')

%%

  outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
  % mridir = '/home/gnolte/neuconn/mri_data/';
  plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';
  
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];


v= 1;
    ord           = pconn_randomization;
% 
	for isubj = SUBJLIST
      for ifoi = 4;
        for m = 1 : 3
          	im = find(ord(isubj,:)==m);

            load(sprintf([outdir 'pconn_cnt_sens_dfa_filt_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));

            allpar(:,isubj,m) = nanmean(par.dfa,2);
            
        end
      end
  end
  
%   
%  
%   
%   
%   
  allpar = allpar(:,SUBJLIST,:);
%   
    %%
    

    
