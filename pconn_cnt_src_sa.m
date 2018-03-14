%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% pconn_cnt_src_sa

clear

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
% v_rawdat = 6;
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v         = 3;
% grid_size = 'cortex';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 4
% --------------------------------------------------------
% v         = 4;
% grid_size = 'aal';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 5
% --------------------------------------------------------
% v         = 5;
% grid_size = 'medium';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 6
% --------------------------------------------------------
% v         = 6;
% grid_size = 'aal_6mm';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
% v         = 7;
% grid_size = 'aal_4mm';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 8
% --------------------------------------------------------
% v         = 8;
% grid_size = 'm758_4mm';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 9
% --------------------------------------------------------
% v         = 9;
% grid_size = 'cortex_lowres';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 10
% --------------------------------------------------------
% v         = 10;
% grid_size = 'vtpm_4mm';
% v_rawdat  = 6;
% --------------------------------------------------------
% VERSION 11
% --------------------------------------------------------
v         = 11;
grid_size = 'vtpm_6mm';
v_rawdat  = 6;
% --------------------------------------------------------

restoredefaultpath

addpath ~/pconn/matlab
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
% addpath /home/gnolte/neuconn/OLD/matlab/rest/
addpath(genpath('/home/gnolte/meth'));

ft_defaults

outdir = '/home/tpfeffer/pconn_cnt/proc/src/';

%%
for im = 1 : 3
  for isubj = 4:34
        
      if ~exist(sprintf([outdir 'pconn_cnt_sa_s%d_m%d_v%d_processing.txt'],isubj,im,v))
        system(['touch ' outdir sprintf('pconn_cnt_sa_s%d_m%d_v%d_processing.txt',isubj,im,v)]);
      else
        continue
      end
%       
      fprintf('Processing s%d m%d ...\n',isubj,im);
      
      % ----- ---------------------------------------------------
      % SELECT MRI DATA
      % --------------------------------------------------------
      
      load sa_meg_template;
      
      fprintf('Looking for MRI ...\n');
     
      mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/mri/',isubj);
      c_mrid = dir(mridir);
      if isempty(c_mrid)
        mridir = sprintf('/home/tpfeffer/pconn/rawdata/mri/p%d/MRI/',isubj);
        c_mrid = dir(mridir);
      end
      
      br = 0;
      
      if ~isempty(c_mrid)
        if length(c_mrid) > 2
          for imri = 3 : length(c_mrid)
            if strcmp(c_mrid(imri).name(end-5:end),'V2.mri')
              mri_data = [mridir c_mrid(imri).name];
              fprintf('Looking for MRI ... Found!\n');
              br = 0;
            end
          end
        else
          fprintf('Looking for MRI ... Not Found!\n');
          br = 1;
        end
      else
        fprintf('Looking for MRI ... Not Found!\n');
%         continue
        br = 1;
      end
      
      if br == 0
        
          aal = tp_aalgrid();
          sa_meg_template.grid_aal4mm = aal.grid_4mm;
          sa_meg_template.grid_aal6mm = aal.grid_6mm;
                  
          m758 = tp_create_grid('m758');
          sa_meg_template.grid_m758_4mm = m758.grid_4mm;
          sa_meg_template.grid_m758_6mm = m758.grid_6mm;

          vtpm = tp_create_grid('vtpm');
          sa_meg_template.grid_vtpm_4mm = vtpm.grid_4mm;
          sa_meg_template.grid_vtpm_6mm = vtpm.grid_6mm;
          
          sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_data);
   
      end
        
      saa = 0;
      
      % --------------------------------------------------------
      % SELECT SUBJECT & BLOCK
      % --------------------------------------------------------
      
      for iblock = 1 : 2
        
        fprintf('Processing block %d MEG-Data ...\n',iblock);
        fprintf('Loading MEG-Data ...\n');

        if br == 1
          mri_data = '/home/tpfeffer/pconn/rawdata/mri/p10/mri/aheitmann_V2.mri';
%           disp(sprintf('Using template MRI ...'))
          load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
          load /home/gnolte/meth/templates/sa_template.mat
          
          aal = tp_aalgrid();
          sa_template.grid_aal4mm = aal.grid_4mm;
          sa_template.grid_aal6mm = aal.grid_6mm;
          
          m758 = tp_create_grid('m758');
          sa_template.grid_m758_4mm = m758.grid_4mm;
          sa_template.grid_m758_6mm = m758.grid_6mm;
          
          vtpm = tp_create_grid('vtpm');
          sa_meg_template.grid_vtpm_4mm = vtpm.grid_4mm;
          sa_meg_template.grid_vtpm_6mm = vtpm.grid_6mm;
        
          
          sa  = tp_mk_sa_meg_withoutmri(sa_template,cfg1.headerfile);
          saa = 1;
          
          
        end
        clear meg_data cfg1 data cfg2 
        
        
        load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
        
        clear data
        meg_data = cfg1.headerfile;
        
        if ~saa
          sa 	= mk_sa_meg_forward(sa_meg1, meg_data);
        end
        
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
          sa.L_aal      = L;
          sa.leadfield  = 'aal';
          
        elseif strcmp(grid_size,'aal_6mm')
           
          L             = grid2L(sa.grid_aal6mm_indi,sa.fp_indi);
          sa.L_aal_6mm  = L;
          sa.leadfield  = 'aal_6mm';
          sa.aal_centr  = aal.centroids(1:91,:);
          sa.aal_label  = aal.tissue_6mm(aal.tissue_6mm>0 & aal.tissue_6mm<92);
          
        elseif strcmp(grid_size,'aal_4mm')
           
          L             = grid2L(sa.grid_aal4mm_indi,sa.fp_indi);
          sa.L_aal_4mm  = L;
          sa.leadfield  = 'aal_4mm';
          sa.aal_centr  = aal.centroids(1:91,:);
          sa.aal_label  = aal.tissue_4mm(aal.tissue_4mm>0 & aal.tissue_4mm<92);
            
        elseif strcmp(grid_size,'m758_4mm')
           
          L              = grid2L(sa.grid_m758_4mm_indi,sa.fp_indi);
          sa.L_m758_4mm  = L;
          sa.leadfield   = 'm758_4mm';
            
        elseif strcmp(grid_size,'m758_6mm')
           
          L             = grid2L(sa.grid_m758_6mm_indi,sa.fp_indi);
          sa.L_m758_6mm = L;
          sa.leadfield  = 'm758_6mm';
            
        elseif strcmp(grid_size,'cortex_lowres')
           
          L             = grid2L(sa.grid_cortex_lowres_indi,sa.fp_indi);
          sa.L_coarse   = L;
          sa.leadfield  = 'cortex_lowres';
          
        elseif strcmp(grid_size,'vtpm_4mm')
        
            L              = grid2L(sa.grid_vtpm_4mm_indi,sa.fp_indi);
            sa.L_vtpm_4mm  = L;
            sa.leadfield   = 'vtpm_4mm';
            
        elseif strcmp(grid_size,'vtpm_6mm')
          
          L              = grid2L(sa.grid_vtpm_6mm_indi,sa.fp_indi);
          sa.L_vtpm_6mm  = L;
          sa.leadfield   = 'vtpm_6mm';
          
        end
        
        
        
        save([outdir sprintf('pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)],'sa');
        close all
      end
      br = 0;
%     catch me
%       save([outdir sprintf('pconn_sa_s%d_m%d_v%d_ERROR.mat',isubj,im,v)],'me');
%       continue
%     end
  end
  
end

