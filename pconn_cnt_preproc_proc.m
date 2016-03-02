%% pconn_cnt_preproc_proc.m
% pconn_cnt_preproc_proc

% last update: 01-03-2015, tpfeffer
% implemented: trigger resample

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 1
% -------------------------------------------------------------------------
v_in          = 1;
v_out         = 1;
pad           = 1200/2; % 500 ms
CFG.dftfreq   = [];
CFG.dftfilter = 'no';
HPFREQ        = 0.5;
planar        = 0;
% -------------------------------------------------------------------------

outdir = '/home/tpfeffer/pconn_cnt/proc/preproc/';

addpath /home/tpfeffer/pconn/matlab/
% addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')
addpath('/home/tpfeffer/fieldtrip-20150318/')

ft_defaults

%%
for m = 1 : 2
  for isubj =2 :3%24
    for ibl = 1 : 2
      
      a(1) = exist([outdir sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_in)]);
%       a(2) = ~exist(sprintf([outdir 'pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d.mat'],isubj,m,ibl,v_out));
      a(3) = ~exist(sprintf([outdir 'pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d_processing.txt'],isubj,m,ibl,v_out));
      
%       a = [1 1 1]
      
%       if a(1) && a(3)
        
        system(['touch ' outdir sprintf('pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d_processing.txt',isubj,m,ibl,v_out)]);
        
        disp(sprintf('Processing s%d b%d m%d',isubj,ibl,m));
        
        load([outdir sprintf('pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_in)])
        load([outdir sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_in)])

        if planar
          cfg             = [];
          cfg.method      = 'template';
          cfg.layout      = 'CTF275';
          neighbours       = ft_prepare_neighbours(cfg);

          cfg              = [];
          cfg.feedback     = 'no';
          cfg.method       = 'template';
          cfg.planarmethod = 'sincos';
          cfg.channel      = {'MEG'};
          cfg.neighbours   = neighbours;
          data             = ft_megplanar(cfg, data);

          cfg = [];
          dat = ft_combineplanar(cfg,data);
        end

      
        cfg = [];
       	cfg = cfgs;

        if v_out ~= 4 && v_out ~= 14
          cfg.channel     = {'MEG'};
        elseif v_out == 5
          cfg.channel     = {'EEG002';'EEG003';'EEG013';'EEG014';'EEG059'};
        end
        cfg.padding     = 10;
        cfg.continuous  = 'yes';
        
        cfg.bsfilter    = 'yes';
        cfg.hpfilter    = 'yes';
        cfg.lpfilter    = 'no';
        cfg.dftfilter   = CFG.dftfilter;
        cfg.dftfreq     = CFG.dftfreq;
        
        cfg.bsfreq      = [49 51; 99 101; 149 151];
        cfg.hpfreq      = 40;
        
%         data.Fs         = 1200;
        data_hi         = ft_preprocessing(cfg);
        
        cfg.hpfiltord   = 4;
        cfg.lpfiltord   = 4;
        
        cfg.hpfreq      = HPFREQ;
        data_low        = ft_preprocessing(cfg);
        
        % ADDED 01-03: TRIGGER CHANNEL NEEDS TO BE PROCESSED AS WELL!!!
      	cfg = [];
        cfg = cfgs;
        cfg.channel     = {'UPPT001'};
        cfg.continuous  = 'yes';
        trig            = ft_preprocessing(cfg);        
        
        save([outdir sprintf('pconn_cnt_trig_raw_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_out)],'trig','-v7.3')
        %
        %%
        % -------------------------------------------------------------------------
        % REJECT ARTIFACTS
        % -------------------------------------------------------------------------
%         cfg = [];
%         cfg.toilim = [1]
%         data_low = ft_redefinetrial(cfg,data_low);
%         
        art(art>data.sampleinfo(end))=data.sampleinfo(end);
        art(art<1) = 1;
        
        for iart = 1 : size(art,1)
          data_low.trial{1}(:,art(iart,1):art(iart,2)) = NaN;          
        end
        
        %% PROBLEM SOLUTION
        
        % CUT OUT LAST BIT OF IDX. it does not matter anyway
        % then data and idx are of same length
        % add idx as second sensor for ft_resample
        % test

        idx = isnan(data_low.trial{1}(1,:));

        data_low.trial{1}(:,idx)=[];
        data_low.time{1}(:,idx) =[];
        if v_out > 2
          data_low.trial{1}(:,1:20*1200) =[];
          data_low.time{1}(:,1:20*1200)  =[];
        end
%         
        data_low.sampleinfo = [1 length(data_low.time{1})];
%         
        for iart = 1 : size(art,1)
          data_hi.trial{1}(:,art(iart,1):art(iart,2)) = NaN;
        end
        
        data_hi.trial{1}(:,idx)=[];
        data_hi.time{1}(:,idx) =[];
        
        if v_out > 2
          data_hi.trial{1}(:,1:20*1200) =[];
          data_hi.time{1}(:,1:20*1200)  =[];
        end
        
        data_hi.sampleinfo = [1 length(data_hi.time{1})];
        
        %% RESAMPLE DATA
        
        cfg2            = [];
        cfg2.resamplefs = 400;
        cfg2.detrend    = 'yes'; % resampling with NaNs does not work if 'yes'
        cfg2.demean     = 'no'; % resampling with NaNs does not work if 'yes'
        data_hi         = ft_resampledata(cfg2, data_hi);
        data_low        = ft_resampledata(cfg2, data_low);
        trig            = ft_resampledata(cfg2, trig);
        data_low.idx    = resample(double(idx),400,1200);
        data_hi.idx     = resample(double(idx),400,1200);
        
        save([outdir sprintf('pconn_cnt_trig_resampled_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_out)],'trig','-v7.3')
        save([outdir sprintf('pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v_out)], 'data_hi','data_low','idx', '-v7.3')

        disp(sprintf('Saved.'));
        
%       else
% %         system(['touch ' outdir sprintf('pconn_preproc_artifacts_s%d_m%d_b%d_v%d_FILENOTFOUND.txt',isubj,m,ibl,v_out)]);
% %         continue
%       end
    end
  end
end