%% PCONN SUBTRACT ICA COMPONENTS
% after pconn_viewcomps_v2.m

% pconn_preproc_rejcomp

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 1
% -------------------------------------------------------------------------
v     = 1; % ICA version
v_art = 1; % preproc version
split = 0;
FS    = 400;
CFG.method = 'fastica';
v_out = 1;
v_rej = 1;
% -------------------------------------------------------------------------
% VERSION 2
% -------------------------------------------------------------------------
% v     = 2; % ICA version
% v_art = 2; % preproc version
% split = 0;
% ses   = 1;
% hc    = 0;
% FS    = 400;
% CFG.method = 'runica';
% v_out = v;
% % -------------------------------------------------------------------------
% VERSION 2
% -------------------------------------------------------------------------
% v     = 1; % ICA version
% v_art = 1; % preproc version
% split = 0;
% FS    = 400;
% CFG.method = 'fastica';
% v_out = 2;
% v_rej = 5;
% % -------------------------------------------------------------------------


outdir  = sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/');

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/') 
ft_defaults

NSUBJ = 24;

%%
for im = 1 : 3
  for isubj = [3 17]
    for iblock = 1 : 2
      
      try
      ex(1) = exist([outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)]);
      
      if ex(1)
        
        load([outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)])
        
%         try
          for ifr = 1 : 1
            
            ex(2) = exist(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifr,v_rej));
            
            if ex(2)
              disp(sprintf('Processing s%d m%d b%d f%d ...',isubj,im,iblock,ifr))
              
              % load rejected components
              load(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifr,v_rej));
              
              % BE CAREFUL WITH VERSION HERE!!!
              % loads pre-ica data
              load(sprintf([outdir 'pconn_cnt_preproc_artifacts_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,v_art));
              data_low.trial{1} = data_low.trial{1}-data_hi.trial{1};
              
              cfg = [];
              cfg.component = find(~rej_comp);
              
              if ifr == 1
                data_art_low	= ft_rejectcomponent(cfg,comp_low); clear comp_low
                data_low.trial{1} = data_low.trial{1} - data_art_low.trial{1};
              elseif ifr == 2
                data_art_hi   = ft_rejectcomponent(cfg,comp_hi); clear comp_hi
                data_hi.trial{1} = data_hi.trial{1} - data_art_hi.trial{1};
              end
              
            else
              error('rej_comp missing for s%d m%d f%d b%d ...',isubj,im,ifr,iblock)
              continue
            end
          end
          if ~exist('data_hi','var')
            save([outdir sprintf('pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v_out)],'data_low','-v7.3')
          else
            save([outdir sprintf('pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v_out)],'data_low','data_hi','-v7.3')
          end
      else
        warning('rej_comp missing for s%d m%d f%d b%d ...',isubj,im,ifr,iblock)

        continue
      end
      catch me
      end
    end
  end
end








