%% nc_preproc_artvec.m
% pconn_preproc_artvec

clear all

restoredefaultpath



% -------------------------------------------------------------------------
% VERSION 1
% -------------------------------------------------------------------------
v   = 1;
pad = 1200/2; % 500 ms
m   = 3; % which session: 1/2/3?
% -------------------------------------------------------------------------

indir  = '/home/tpfeffer/pconn/rawdata/meg/';
outdir = '/home/tpfeffer/pconn_cnt/proc/preproc/';

addpath /home/tpfeffer/pconn_cnt/matlab/
addpath('/home/tpfeffer/fieldtrip-20150318/')

ft_defaults

%%

% isubj = 3;

  for im = 2 : 3
%   try
for isubj = 1


for ibl = 1 : 2
  
  clear cfgs
  
  if ~exist([outdir sprintf('pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)])
%     continue
  end
    
    if exist([outdir sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)])
%       continue
    end
% %       
      disp(sprintf('Processing s%d m%d ...',isubj,im))
      
      load([outdir sprintf('pconn_cnt_preproc_data_count_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)])
      
      cfg1.trl = cfg1.trl;
      cfgs.trl = cfg1.trl;
      
      % CUT OUT 10 MINUTES
      if size(data.trial,2)>800
%         for itrl = 1 : size(data.trial,2)
%           dat((itrl-1)*data.fsample+1:itrl*data.fsample)  = data.trial{itrl}(strcmp(data.label,'UPPT001'),:);
%         end
      else
%         for itrl = 1 : size(data.trial,2)
%           dat((itrl-1)*data.fsample+1:itrl*data.fsample)  = data.trial{itrl}(strcmp(data.label,'UPPT001'),:);
%         end
%       	strt = find(round(dat)==90,1,'first');
      end
      
      
      % VISUAL ARTIFACT REJECTION
      
      cfg = [];
      cfg.method    = 'summary';
      cfg.channel   = {'MEG'};
      tmp_data      = ft_rejectvisual(cfg,data);
      artifact_vis  = tmp_data.cfg.artfctdef.summary.artifact; tmp_data
%       artifact_vis = data.sampleinfo(end,:);
      
      cfg = [];
      cfg.channel   = {'MEG'};
      cfg.artfctdef.rejvis.artifact = artifact_vis;
      cfg.viewmode = 'vertical';
      cfg = ft_databrowser(cfg,data);
      artifact_vis = cfg.artfctdef.rejvis.artifact;

      % -------------------------------------------------------------------------
      % DETECT JUMP ARTIFACTS
      % -------------------------------------------------------------------------
      cfg = cfgs;
      cfg.trl=cfg.trl(1:end-1,:);
      cfg.continuous  = 'yes';
      cfg.artfctdef.zvalue.channel     = 'MEG';
      cfg.artfctdef.zvalue.cutoff      = 60;
      cfg.artfctdef.zvalue.trlpadding  = 0;
      cfg.artfctdef.zvalue.artpadding  = 0;
      cfg.artfctdef.zvalue.fltpadding  = 0;
      cfg.artfctdef.zvalue.boxcar      = 0;
      cfg.artfctdef.zvalue.interactive = 'yes';
      
      cfg.artfctdef.zvalue.cumulative     = 'yes';
      cfg.artfctdef.zvalue.medianfilter   = 'yes';
      cfg.artfctdef.zvalue.medianfiltord  = 9;
      cfg.artfctdef.zvalue.absdiff        = 'yes';

      [cfg, artifact_jump] = ft_artifact_zvalue(cfg);
%     
      % -------------------------------------------------------------------------
      % DETECT MUSCLE ARTIFACTS
      % -------------------------------------------------------------------------
      cfg = cfgs;
      cfg.continuous  = 'yes';
      cfg.artfctdef.zvalue.channel     = {'MRT','MLT'};
      cfg.artfctdef.zvalue.cutoff      = 8;
      cfg.artfctdef.zvalue.trlpadding  = 0;
      cfg.artfctdef.zvalue.artpadding  = 0.3;
      cfg.artfctdef.zvalue.fltpadding  = 0;
      cfg.artfctdef.zvalue.bpfilter    = 'yes';
      cfg.artfctdef.zvalue.bpfreq      = [110 140];
      cfg.artfctdef.zvalue.bpfiltord   = 7;
      cfg.artfctdef.zvalue.bpfilttype  = 'but';
      cfg.artfctdef.zvalue.hilbert     = 'yes';
      cfg.artfctdef.zvalue.boxcar      = 0.2;
      cfg.artfctdef.zvalue.interactive = 'no';
      
      [cfg, artifact_muscle] = ft_artifact_zvalue(cfg);
      
      % SAVE ALL ARTIFACTS
      if ~isempty(artifact_jump)
        art = [artifact_muscle; artifact_vis; artifact_jump];
      else
        art = [artifact_muscle; artifact_vis];
      end
      art(:,1) = art(:,1)-pad; art(:,2)=art(:,2)+pad;
      art(art>data.sampleinfo(end))=data.sampleinfo(end);
      art(art<0) = 1;
      
      if exist('strt','var')
        save([outdir sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)], 'art','strt')
      else
        save([outdir sprintf('pconn_cnt_preproc_artvec_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)], 'art')        
      end
    end
  end
end

%   end
% end
% end
% end