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
v_rawdata = 1;
fsample   = 400;
SUBJLIST 	= [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
FOI       = 2:0.2:200;
T         = 0.2;
planar    = 0;
% --------------------------------------------------------

restoredefaultpath

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
ft_defaults

outdir   = '/home/tpfeffer/pconn_cnt/proc/powspec/';

%%

for m = 3 : 3
  for isubj = 33
%     
%     if ~exist(sprintf([outdir 'pconn_cnt_powspec_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' outdir sprintf('pconn_cnt_powspec_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end
%     
    disp(sprintf('Processing s%d m%d ...', isubj,m))
        
    for iblock = 2 : 2
      
      disp(sprintf('Loading MEG data ...'));
      
      load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
      
      cfg = [];
      cfg.length = 1/T;
      cfg.overlap = 0;
      
      data      = ft_redefinetrial(cfg,data); clear data_low

      cfg         = [];      
      cfg.trials  = 1:size(data.trial,2)-1;
      data        = ft_redefinetrial(cfg,data);

      cfg             = [];
      cfg.method      = 'mtmfft';
      cfg.output      = 'pow';
      cfg.taper       = 'dpss';
      cfg.channel     = {'MEG'};
      cfg.foi         = FOI;
      cfg.keeptrials  = 'yes';
      cfg.tapsmofrq   = 2;
%       cfg.trials = 1;
      dat             = ft_freqanalysis(cfg, data); % run freqanalysis here
      
      save([outdir sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v)],'dat','-v7.3');
      
%       tp(:,:,iblock) = nanmean(dat.powspctrm,1);
      
    end
  end
end

error('Stop here!')

%% 
% 
% for m = 1 : 3
%   for isubj = SUBJLIST
%     isubj
%     for iblock = 1  : 2
%     
%     load([outdir sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v)]);
%     
%     pow(:,m,isubj,iblock) = squeeze(nanmean(nanmean(dat.powspctrm,1),2));
%     
%     end
%   end
% end
% 
% pow = squeeze(nanmean(nanmean(nanmean(pow(:,:,SUBJLIST,:),4),3),2));
% 
% 


