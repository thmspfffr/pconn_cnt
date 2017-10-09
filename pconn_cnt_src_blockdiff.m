%% pconn_cnt_src_blockdiff.m
% pconn_cnt_src_blockdiff


clear
restoredefaultpath
rng('shuffle','twister')

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 2;
v_stat = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

para.minneigh     = 2;
para.clusteralpha = 0.05;
para.alpha        = 0.025;
para.nperm        = 10000;
% --------------------------------------------------------
% VERSION 22
% --------------------------------------------------------
% v_stat         = 22;
% v = 2;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
%
% para.minneigh     = 2;
% para.clusteralpha = 0.05;
% para.alpha        = 0.025;
% para.nperm        = 10000;
% --------------------------------------------------------


outdir = '/home/tpfeffer/pconn/proc/dfa/';
indir  = '/home/tpfeffer/pconn/proc/preproc/';

addpath /home/tpfeffer/pconn/matlab/


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab/

load sa_meg_template;
ord   = pconn_randomization;

contr = [2 1; 3 1];

%% LOAD DATA
allstr       = {'dfa';'cvar';'pow'};
all_behav    = {'count';'numb_switches';'var_dur'};
all_cond     = {'rst';'tsk'};
all_meth     = {'correlation_pearson';'correlation_spearman'};
g1    = sa_meg_template.grid_cortex3000;

for ifoi = [2 5]
  for istr = 1 : 3
    str         = allstr{istr};
    for icond = 2 : 2
      cond        = all_cond{icond};
      for imeth = 1 : 1
        para.method = all_meth{imeth};
        % % %
        if ~exist(sprintf(['~/pconn_all/proc/' 'pconn_cnt_src_blockdiff_str%d_cond%d_meth%d_f%d_v%d_processing.txt'],istr,icond,imeth,ifoi,v_stat))
          system(['touch ' '~/pconn_all/proc/' sprintf('pconn_cnt_src_blockdiff_str%d_cond%d_meth%d_f%d_v%d_processing.txt',istr,icond,imeth,ifoi,v_stat)]);
        else
          continue
        end
        %
        fprintf('Computing f%d s%d cond%d ...\n',ifoi,istr,icond)
        
        str_behav   = all_behav{1};
        
        % READ IN BEHAVIORAL DATA
        if ~strcmp(str_behav,'avg')
          para.str_behav = str_behav;
          par_behav = pconn_read_behavioral_data(SUBJLIST,para);
        else
          for ib = 1 : 2
            para.str_behav   = all_behav{ib};
            par_behav_tmp{ib} = pconn_read_behavioral_data(SUBJLIST,para);
          end
          
          str_behav = 'avg';
          par_behav = cat(3,par_behav_tmp{1},par_behav_tmp{2});
          
        end
        
        par_diff = [par_behav(:,:,2)-par_behav(:,:,1)]';
        
        [nanidx,~]=find(isnan(par_diff));
        %% READ NEURAL DATA
        
        for isubj = SUBJLIST
          fprintf('Reading subj %d ...\n',isubj);
          for m = 1 : 3
            
            im = find(ord(isubj,:)==m);
            if strcmp(cond,'rst')
              if ~strcmp(str,'pow')
                load(sprintf('~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,v));
                eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
              else
                load(sprintf('~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,1));
                eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
              end
            elseif strcmp(cond,'tsk')
              if ~strcmp(str,'pow')
                load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,v));
                eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
              else
                load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,1));
                eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
              end
            end
          end
        end
        
        all_par   = all_par(:,:,:,SUBJLIST);
        
        %% COMPUTE STATSISTICS
        
        all_brain = squeeze(nanmean(all_par(:,2,:,:)-all_par(:,1,:,:),3));
        
        if ~isempty(nanidx)
          all_brain(:,nanidx) = [];
          par_diff(nanidx,:)  = [];
        end
        
        n                 = get_neighbours(g1);
        para.neigh        = n;
        para.extvar       = nanmean(par_diff,2);
        
        stats           = tp_clusterperm(all_brain,para);
        
        save(sprintf(['~/pconn_all/proc/' 'pconn_cnt_src_blockdiff_avg_str%d_cond%d_meth%d_f%d_v%d.mat'],istr,icond,imeth,ifoi,v_stat),'stats','-v7.3');
        
        
      end
    end
  end
end


error('!')
%%  PLOT STUFF

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

v_stat = 2;

for ifoi = [2 5]
  for istr = 1 : 3
    for icond = 2 : 2
      for imeth = 1 : 1
        
        load(sprintf(['~/pconn_all/proc/' 'pconn_cnt_src_blockdiff_avg_str%d_cond%d_meth%d_f%d_v%d.mat'],istr,icond,imeth,ifoi,v_stat));
        
        par_interp = spatfiltergauss(stats.stat,g1,dd,g2);
   
        para = [];         
        para.colorlimits = [-3 3];
        
        tp_showsource(par_interp,cmap,sa_meg_template,para);

        print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_blockdiff_mask0_f%d_str%d_cond%d_meth%d_v%d.jpg',ifoi,istr,icond,imeth,v))

        par_interp = spatfiltergauss(stats.mask.*stats.stat,g1,dd,g2);
        
        tp_showsource(par_interp,cmap,sa_meg_template,para);
        
        print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_blockdiff_mask1_f%d_str%d_cond%d_meth%d_v%d.jpg',ifoi,istr,icond,imeth,v))

      end
    end
  end
end

%% CLEAN CRASHED FILES


outdir   = '/home/tpfeffer/pconn_all/proc/';
cnt = 0;
v_stat = 2;
cnt_exist = 0;

for ifoi = 1 : 5
  for istr = 1 : 3
    for icontr = 1 : 2
      for icond = 1 : 2
        for imeth = 1 : 1
          for ibehav = 1 : 3
            
            ifoi
            
            txt_file  = exist(sprintf(['~/pconn_all/proc/all_behav_src_stats_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d_processing.txt'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat));
            mat_file1 = exist(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pharm_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat));
            mat_file2 = exist(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pbo_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat));
            
            if txt_file>0 & mat_file1>0 & mat_file2>0
              cnt_exist = cnt_exist + 1;
              continue
              
            elseif mat_file1>0 & mat_file2>0 & ~(txt_file>0)
              system(['touch ' sprintf(['~/pconn_all/proc/all_behav_src_stats_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d_processing.txt'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat)]);
              
            elseif ~(mat_file1>0) & ~(mat_file2>0) & txt_file>0
              warning(sprintf('Deleting stuff: s%d m%df %d',isubj,m,ifoi))
              delete(sprintf(['~/pconn_all/proc/all_behav_src_stats_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d_processing.txt'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat))
              cnt = cnt + 1;
            else
              warning('Nothing exists')
              %               cnt = cnt+1;
            end
          end
        end
      end
    end
  end
end
cnt



