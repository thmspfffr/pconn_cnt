%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_cnt_sens_planar_dfa

% --------------------------------------------------------
% VERSION 19
% --------------------------------------------------------
% v         = 19;
% v_rawdata = 1;
% is_src    = 0;
% fsample   = 400;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
% foi       = [2 4; 4 8; 8 12; 12 36];
% i_fit     = [1 100];
% i_calc    = [0.5 150];
% dfa_overlap = 0.5;
% filt_ord = 2;
% --------------------------------------------------------
% VERSION 20
% --------------------------------------------------------
v         = 20;
v_rawdata = 2;
is_src    = 0;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
foi       = [2 4; 4 8; 8 12; 12 36];
i_fit     = [1 100];
i_calc    = [0.5 150];
dfa_overlap = 0.5;
filt_ord = 2;
% --------------------------------------------------------



addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath ~/Documents/MATLAB/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
for m = 1 : 3
  for isubj =SUBJLIST
    for ifoi = 1 : length(foi)
      
      if ~exist(sprintf([outdir 'pconn_cnt_sens_planar_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_sens_planar_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      if isubj == 3
        load ~/pconn/matlab/pconn_sensorlabels.mat
      end
      
      for iblock = 1 : 2
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
        
        cfg = [];
        cfg.method      = 'template';
        cfg.template    = 'CTF275_neighb.mat';
        cfg.neighbours  = ft_prepare_neighbours(cfg);
        cfg.method      = 'sincos';
        data_pl         = ft_megplanar(cfg,data_low);
        
        if isubj == 3
          [~,idx_lab]=intersect(data_low.label,lab);
          data_low.trial{1}=data_low.trial{1}(idx_lab,:);
          data_low.label = lab;
          save(['~/pconn_cnt/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
        end
        
        clear data_hi clear data_low
        
        [mydata,epleng] = megdata2mydata(data_pl);
       
        mydata = mydata(1:end-1000,:);
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        
        % compute bp-filtered signal
        tmp    = single(nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1)));
                
        % LABELING OF PLANAR GRADIOMETERS
        planar    = planarchannelset(data_pl);
        [~, sel_dH]    = match_str(planar(:,1), data_pl.label);  % indices of the horizontal channels
        [~, sel_dV]    = match_str(planar(:,2), data_pl.label);  % indices of the vertical   channels
 
        ampenv = abs(hilbert(sqrt(tmp(:,sel_dH).^2+ tmp(:,sel_dV).^2)));
                 
      	pha = angle(hilbert(sqrt(tmp(:,sel_dH).^2+ tmp(:,sel_dV).^2)));
        
        r = abs(sum(exp(i*pha),2)/size(mydata,2)); clear tmp mydata  pha
        
        [dfa,~,DFA_y_all] = nbt_doDFA(ampenv, siginfo, i_fit,i_calc,dfa_overlap,0,0,[]);
        
        [r_dfa,~,~] = nbt_doDFA(r, siginfo, i_fit,i_calc,dfa_overlap,0,0,[]);

        win = hanning(400);
        p = zeros(size(data_pl.trial{1},1),length(foi(ifoi,1):foi(ifoi,2)));
        
        for ichan = 1 : size(data_pl.trial{1},1)
          [p(ichan,:),f]=pwelch(data_pl.trial{1}(ichan,:),win,50,foi(ifoi,1):foi(ifoi,2),fsample);
        end     
        
      
        par.dfa(:,iblock)   = dfa.MarkerValues;
        par.amp(:,iblock)   = mean(ampenv);
        par.var(:,iblock)   = nanvar(ampenv);
        par.cvar(:,iblock)  = sqrt(par.var(:,iblock))./(par.amp(:,iblock));    
        par.pow(:,iblock)   = nanmean(p(sel_dH,:)+p(sel_dV,:),2);
        par.r_dfa(iblock)   = r_dfa.MarkerValues;
        par.r_mean(iblock)  = mean(r);
        par.r_std(iblock)   = std(r);
                        
        clear dfa dat ampenv r
        
      end
      
      save(sprintf([outdir 'pconn_cnt_sens_planar_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
      clear par r

    end
  end
end


error('STOP')

