%% COMPUTE VARIABLES FOR COMPONENT BROWSER

% pconn_cnt_viewcomp_compute_vars

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 1
% -------------------------------------------------------------------------
v     = 1; % ICA version
v_art = 1; % preproc version
split = 0;
ses   = 1;
FS    = 400;
CFG.method = 'fastica';
v_out = v;
% % -------------------------------------------------------------------------

indir   = sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/');
outdir  = sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/');

addpath /home/gnolte/neuconn/matlab/rest/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

NSUBJ = 24;

for im = 1 : 3
  for isubj = 2: 3
    for iblock = 1 : 2
      for ifreq = 1: 2
        
%         if ~exist(sprintf([outdir 'pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d_processing.txt'],isubj,im,iblock,ifreq,v))
%         	system(['touch ' outdir sprintf('pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d_processing.txt',isubj,im,iblock,ifreq,v)]);
%         else
%           continue
%         end

        try
          
          disp(sprintf('Processing s%d m%d b%d f%d ...',isubj,im,iblock,ifreq))
                 
          load([indir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v_out)])
          clear smax smoothed comp_var freq

          if ifreq == 1
            hilow = 'low'; clear comp_hi
          else
            hilow = 'hi'; clear comp_low
          end        

          % COMPUTE VARIANCE FOR 2s-WINDOWS
          siz = eval(sprintf('size(comp_%s.trial{1},1)',hilow));
          for i = 1 : siz
            
            smax =  eval(sprintf('floor(size(comp_%s.trial{1},2)/100);',hilow));
            for s = 1 : smax
              eval(sprintf('comp_var(i,s)=var(comp_%s.trial{1}(i,(s-1)*100+1:s*100));',hilow));
            end
            % -----------------
            if strcmp(hilow,'hi')
              [smoothed(:,i), freq]= pwelch(comp_hi.trial{1}(i,:),hann(1500),[],[],400);
            else
              [smoothed(:,i), freq]= pwelch(comp_low.trial{1}(i,:),hann(1500),[],[],400);
            end
          end

          save([outdir sprintf('pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifreq,v)], 'smoothed','freq','comp_var','smax','-v7.3')
          
        catch me
          disp(sprintf('Error processing s%d m%d b%d f%d ...',isubj,im,iblock,ifreq))
        	system(['touch ' outdir sprintf('pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d_processing_ERR.txt',isubj,im,iblock,ifreq,v)]);
          save([outdir sprintf('pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d_err.mat',isubj,im,iblock,ifreq,v)], 'me','-v7.3')
        end
        
      end
    end
  end
end
        
        
        
        