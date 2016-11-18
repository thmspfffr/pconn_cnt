%% COMPUTE AMPLITUDE ENVELOPE CORRELATIONS
% pconn_cnt_src_cs

% function nc_src_ana(rand_num)

% rng(rand_num,'twister')

clear all
% --------------------------------------------------------
% VERSION 1 - ELORETA
% --------------------------------------------------------
time      = 0;          % time resolved?
v         = 1;          % version of head model
v_cs      = 1;          % version of cross spectrum
v_filt    = 1;
v_out     = 1;          % version of output
SEGLENG   = 400;        % segment length
maxfreqbin = 200;
v_postproc = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
addpath /home/gnolte/neuconn/matlab/rest/
addpath ~/pconn/matlab

ft_defaults

outdir = '/home/tpfeffer/pconn_cnt/proc/src/';

freq = 1;


%%
for im = 1 : 3
  for isubj = 33:34
    for iblock = 1 : 2
      
      clear A L grid cs nave coh
      
      % if no output yet and file not currently being processed
      if ~exist([outdir sprintf('pconn_cnt_sens_cs_s%d_m%d_b%d_f%d_v%d_processing.txt',isubj,im,iblock,freq,v_cs)])
%         try
          system(['touch ' outdir sprintf('pconn_cnt_sens_cs_s%d_m%d_b%d_f%d_v%d_processing.txt',isubj,im,iblock,freq,v_out)]);
          
          disp(sprintf('Processing s%d m%d f%d ...',isubj,im,freq));
          disp(sprintf('Processing block %d ...',iblock));
          
          % -----------------------------------------------------------------
          % Load data
          % -----------------------------------------------------------------
          disp(sprintf('Loading MEG data ...'));
          load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1));
          % -----------------------------------------------------------------
          

          [mydata,epleng] = megdata2mydata(data);
          
          disp(sprintf('Computing cross spectrum ...'));
          
          segleng         = SEGLENG;
          segshift        = segleng/2;
          
          
          [cs, coh, nave] = pconn_data2cs_event(mydata,segleng,segshift,epleng,maxfreqbin);
          
          disp(sprintf('Saving cross spectrum ...'));
          save([outdir sprintf('pconn_cnt_sens_cs_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,freq,v_cs)],'cs');
          save([outdir sprintf('pconn_cnt_sens_coh_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,freq,v_cs)],'coh','nave');
          

      else
        continue
      end
      % -----------------------------------------------------------------
      % Catch error message and save
      % -----------------------------------------------------------------
      
    end
    clear data_hi data_low
  end
  
end



