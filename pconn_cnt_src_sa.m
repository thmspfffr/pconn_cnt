%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% pconn_cnt_src_sa

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v         = 1;
% grid_size = 'coarse';
% v_rawdat = 1;
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v         = 2;
% grid_size = 'xcoarse';
% v_rawdat = 1;
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
v         = 3;
grid_size = 'cortex';
v_rawdat  = 2;
% --------------------------------------------------------
% VERSION 4
% --------------------------------------------------------
% v         = 4;
% grid_size = 'aal';
% v_rawdat  = 2;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
addpath /home/gnolte/neuconn/OLD/matlab/rest/

ft_defaults

outdir = '/home/tpfeffer/pconn_cnt/proc/src/';

%%
for im = 1 : 3
  for isubj = 32:34
    
    indir  = sprintf('/home/tpfeffer/pconn_cnt/rawdata/meg/p%d/s%d/',isubj,im);
%     
      if ~exist(sprintf([outdir 'pconn_cnt_sa_s%d_m%d_v%d_processing.txt'],isubj,im,v))
        system(['touch ' outdir sprintf('pconn_cnt_sa_s%d_m%d_v%d_processing.txt',isubj,im,v)]);
      else
        continue
      end
%       
      disp(sprintf('Processing s%d m%d ...',isubj,im));
      
      % ----- ---------------------------------------------------
      % SELECT MRI DATA
      % --------------------------------------------------------
      
      load sa_meg_template;
      
      disp(sprintf('Looking for MRI ...'));
     
      mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/mri/',isubj);
      c_mrid = dir(mridir);
      if isempty(c_mrid)
        mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/MRI/',isubj);
        c_mrid = dir(mridir);
      end
      
      if ~isempty(c_mrid)
        if length(c_mrid) > 2
          for imri = 3 : length(c_mrid)
            if strcmp(c_mrid(imri).name(end-5:end),'V2.mri')
              mri_data = [mridir c_mrid(imri).name];
              disp(sprintf('Looking for MRI ... Found!'));
              br = 0;
            end
          end
        else
          br = 1;
        end
      else
        disp(sprintf('Looking for MRI ... Not Found!'));
%         continue
        br = 1;
      end
      addpath ~/pconn/matlab
      
      if br == 0
        sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_data);
      else
        mri_data = '/home/tpfeffer/pconn/rawdata/mri/p10/mri/aheitmann_V2.mri';
        disp(sprintf('Computing SA based on s10 ...'))
        sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_data);
      end
      
      % --------------------------------------------------------
      % SELECT SUBJECT & BLOCK
      % --------------------------------------------------------
      
      for iblock = 1 : 2
        clear meg_data
        
        disp(sprintf('Loading MEG-Data ...'));
        
        load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
        clear data
        meg_data = cfg1.headerfile;
        
        
        disp(sprintf('Processing block %d MEG-Data ...',iblock));
        
        sa          = mk_sa_meg_forward(sa_meg1, meg_data);
        
        if      strcmp(grid_size,'medium')
          
          L            = grid2L(sa.grid_medium_indi,sa.fp_indi);
          sa.L_medium  = L;
          sa.leadfield = 'medium';
          
        elseif  strcmp(grid_size,'coarse')
          
          L            = grid2L(sa.grid_coarse_indi,sa.fp_indi);
          sa.L_coarse  = L;
          sa.leadfield = 'coarse';
          
        elseif strcmp(grid_size,'cortex')
          
          L             = grid2L(sa.grid_cortex3000_indi,sa.fp_indi);
          sa.L_coarse   = L;
          sa.leadfield  = 'cortex';
          
        elseif strcmp(grid_size,'xcoarse')
          
          L             = grid2L(sa.grid_xcoarse_indi,sa.fp_indi);
          sa.L_xcoarse  = L;
          sa.leadfield  = 'xcoarse';
          
        elseif strcmp(grid_size,'aal')
          
          load aalmask_grid_medium
          load sa_meg_template;
          para.grid = 'medium';
          pos = tp_grid2aal(sa.grid_medium_indi',para);
          sa.grid_aal_indi = pos';

          L             = grid2L(sa.grid_aal_indi,sa.fp_indi);
          sa.L_aal  = L;
          sa.leadfield  = 'aal';
            
            
        end
        
        save([outdir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)],'sa');
        close all
      end
%     catch me
%       save([outdir sprintf('pconn_sa_s%d_m%d_v%d_ERROR.mat',isubj,im,v)],'me');
%       continue
%     end
  end
  
end

