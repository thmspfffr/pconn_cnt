% =========================================================================
% INDEPENDENT COMPONENT ANALYSIS
% =========================================================================
% pconn_cnt_preproc_ica

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 1 - FASTICA
% -------------------------------------------------------------------------
v = 1;
v_out = 1;
CFG.method                = 'fastica';
split                     = 0;
CFG.fastica.approach      = 'defl';
numOfIC       = [128 32];
CFG.fastica.g             = 'pow3';
CFG.fastica.stabilization = 'on';
CFG.fastica.verbose       = 'on';
CFG.fastica.firstEig      = 1;
CFG.fastica.displaymode   = 'off';
% -------------------------------------------------------------------------

outdir = '/home/tpfeffer/pconn_cnt/proc/preproc/';

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

ft_defaults

%%
for m = 1 : 3
  for isubj =  34 : 34
    for iblock = 1 : 2
%       try
        
        if exist(sprintf([outdir 'pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d_processing.txt'],isubj,m,iblock,v_out)) 
          continue
        end
%         
        system(['touch ' outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d_processing.txt',isubj,m,iblock,v_out)]);

        disp(sprintf('Processing s%d, m%d, b%d',isubj,m,iblock))

        % READ DATA
        load(sprintf([outdir 'pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d.mat'],isubj,m,iblock,v));
        
        % Subtract high-freq data from low freq data
        % this way, data can be summed again later
        data_low.trial{1} = data_low.trial{1}-data_hi.trial{1};
        
        % -----------------------------------------------------
        % ICA SETTINGS
        % -----------------------------------------------------
        
        cfg = [];
        cfg = CFG;

        if v_out == 6
          cfg2          = [];
          cfg2.gradient = 'G3BR';
          data_low      =	ft_denoise_synthetic(cfg2, data_low);
         	data_hi       =	ft_denoise_synthetic(cfg2, data_hi);  
        end
        
         if strcmp(cfg.method,'fastica')
          cfg.fastica.lastEig = size(data_hi.trial{1},1);
        end
        
        cfg.channel = 'MEG';
        cfg.numcomponent = numOfIC(1);
        comp_low = ft_componentanalysis(cfg,data_low);
        
        cfg.numcomponent = numOfIC(2);
        comp_hi = ft_componentanalysis(cfg,data_hi);
        
        save([outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_out)],'comp_low','comp_hi','-v7.3')
        
%       catch me
%         system(['touch ' outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d_processing_ERR.txt',isubj,m,iblock,v_out)]);        
%       end
    end
  end
end
  
  
  
  
