%% COMPUTE POWER SPECTRUM FOR ALL SENSORS
% fits a straight line to the powerspectrum in order
% to identify peaks in the spectrum for subsequent
% DFA analysis.

% pconn_cnt_powspec

% last update: 18-02-2015, tpfeffer

% to be implemented:
% ***nothing left to implement***

% --------------------------------------------------------
% VERSION 6 - NO PLANAR TRANSFORMATION
% --------------------------------------------------------
v         = 1;
d         = 'data_low.trial{1}+data_hi.trial{1};';
v_rawdata = 2;
fsample   = 400;
SUBJLIST 	= [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
FOI       = 2:0.2:200;
T         = 0.2;
planar    = 0;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/fieldtrip-20150318/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

ft_defaults

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/powspec/';
plotdir = '/home/tpfeffer/pconn_cnt/plots/';

%%

for m = 1 : 3
  for isubj = SUBJLIST
    
    if ~exist(sprintf([outdir 'pconn_cnt_powspec_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pconn_cnt_powspec_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    disp(sprintf('Processing s%d m%d ...', isubj,m))
        
    for iblock = 1 : 2
      
      disp(sprintf('Loading MEG data ...'));
      
        
      load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
      data_low.trial{1} = eval(d);
      
      cfg = [];
      cfg.length = 1/T;
      cfg.overlap = 0;
      
      data      = ft_redefinetrial(cfg,data_low); clear data_low

      cfg         = [];      
      cfg.trials  = 1:size(data.trial,2)-1;
      data        = ft_redefinetrial(cfg,data);

      cfg             = [];
      cfg.method      = 'mtmfft';
      cfg.output      = 'pow';
      cfg.taper       = 'hanning';
      cfg.channel     = {'MEG'};
      cfg.foi         = FOI;
      cfg.keeptrials  = 'yes';
      dat             = ft_freqanalysis(cfg, data); % run freqanalysis here
      
      save([outdir sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v)],'dat');
      
      tp(:,:,iblock) = nanmean(dat.powspctrm,1);
      
    end
  end
end

error('Stop here!')




