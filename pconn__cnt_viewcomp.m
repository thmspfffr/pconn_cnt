%% NEUCONN RESTING-STATE VIEW ICA COMPONENTS
% v2: difference to older version: ica was computed w/o the use of fieldtrip.
% therefore, data needs to be party restructured.

% views ica components in gui
% components can be marked for later rejection
% thomas pfeffer, february 2014

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 2
% -------------------------------------------------------------------------
v     = 1; % ICA version
split = 0;
ses   = 1;
hc    = 0;
FS    = 400;
CFG.method = 'runica';
v_out = v;
% -------------------------------------------------------------------------

indir   = sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/');
outdir  = sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/');

addpath /home/gnolte/neuconn/matlab/rest/
% addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')
addpath('/home/tpfeffer/fieldtrip-20150318/')
ft_defaults

NSUBJ = 33;

all_color = [linspace(1,0,13)' linspace(0,1,13)' zeros(13,1)];

% =========================================================================
% =========================================================================
% DO NOT MODIFY FOLLOWING SEGMENTS
% =========================================================================
% =========================================================================

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab
%% TEST

% 
% if isubj== 1 || isubj == 2 || isubj == 3
  load ~/pconn_cnt/proc/src/pconn_cnt_sa_s2_m1_b1_v3.mat
% else
%   load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
% end

total_err = 0;

for im = 1 : 3
  for iblock = 1 : 2
    for is = 1 : NSUBJ
      ex(im,iblock,is) = exist([outdir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',is,im,iblock,v)])>0;
    end
  end
end

subj_col = {'r';'g'};

while total_err == 0
  back = 1;
  while back == 1
    back = 0;
    err = 0;
    % INTERFACE
    figure('units','normalized','outerposition',[0.3 0.3 0.5 0.5]);
    
    for is = 1 : NSUBJ
      subj_cb{is} = '';
    end
    
    for is = 1 : 20
      cnt = 0;
      for im = 1 : 3;
        for ibl = 1 : 2;
          for iff = 1 : 2
            cnt = cnt + 1;
            a(cnt) = exist(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],is,im,ibl,iff,v));
          end
        end
      end
      
      subj_col = all_color(sum(a>1)+1,:);
      
%       helpbl
%       if ~any(a<1)
%         subj_col = {'r';'g'};
%       else
%         subj_col = {'r';'w'};
%       end
      
      subj{is} = eval(sprintf('uicontrol(''Units'',''normalized'',''Position'',[-0.05+is/20 0.90 0.05 0.05],''Style'',''pushbutton'',''Backgroundcolor'',subj_col,''String'',''S%d'',''Callback'',''isubj = %d'')',is,is));
    end
    
    for is = 21 : NSUBJ
      cnt = 0;
      for im = 1 : 3;
        for ibl = 1 : 2;
          for ifr = 1 : 2;
            cnt = cnt + 1;
            a(cnt) = exist(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],is,im,ibl,iff,v));
          end
        end
      end
      
     	subj_col = all_color(sum(a>1)+1,:);

      
      subj{is} = eval(sprintf('uicontrol(''Units'',''normalized'',''Position'',[-0.05+(is-20)/20 0.85 0.05 0.05],''Style'',''pushbutton'',''Backgroundcolor'',subj_col,''String'',''S%d'',''Callback'',''isubj = %d; clear comp_var'')',is,is));
    end
    
    clear im
    
    bl1 = uicontrol('Units','normalized','Position',[0.38 0.67 0.1 0.1],'Style','pushbutton','String','Block 1','Callback','iblock = 1;');
    bl2 = uicontrol('Units','normalized','Position',[0.52 0.67 0.1 0.1],'Style','pushbutton','String','Block 2','Callback','iblock = 2;');
    
    ses1 = uicontrol('Units','normalized','Position',[0.40 0.57 0.06 0.05],'Style','pushbutton','String','S1','Callback','im = 1;');
    ses2 = uicontrol('Units','normalized','Position',[0.47 0.57 0.06 0.05],'Style','pushbutton','String','S2','Callback','im = 2;');
    ses3 = uicontrol('Units','normalized','Position',[0.54 0.57 0.06 0.05],'Style','pushbutton','String','S3','Callback','im = 3;');
    
    seg1 = uicontrol('Units','normalized','Position',[0.1 0.3 0.3 0.2],'Style','pushbutton','String','Low Frequencies','Callback','which_freq = 1; close');
    seg2 = uicontrol('Units','normalized','Position',[0.6 0.3 0.3 0.2],'Style','pushbutton','String','High Frequencies','Callback','which_freq = 2; close');
    
    abort = uicontrol('Units','normalized','Position',[0.4 0.01 0.2 0.1],'Style','pushbutton','String','QUIT','Callback','close; err = 1;');
    
    uiwait
    
    if err == 1
      error('Aborted!');
    end
    
    clear comp_var smoothed
    % LOAD DATA
    load([indir sprintf('pconn_cnt_preproc_ica_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v_out)])
    
    if exist(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,v_out))
      col = 'r';
    else
      col = 'g';
    end
    
  end
  
  if which_freq == 1
    hilow = 'low';
  else
    hilow = 'hi';
  end
  
  subpl = 4;
  if strcmp(CFG.method,'fastica') && strcmp(hilow,'low')
    rej_comp = zeros(size(comp_low.trial{1},1),1);
  elseif strcmp(CFG.method,'fastica') && strcmp(hilow,'hi')
    rej_comp = zeros(size(comp_hi.trial{1},1),1);
  elseif strcmp(CFG.method,'runica')
    if strcmp(hilow,'hi')
      rej_comp = zeros(size(comp_hi.trial{1},1),1);
    else
      rej_comp = zeros(size(comp_low.trial{1},1),1);
    end
  end
  
  if exist(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,which_freq,v_out))
    load(sprintf([outdir 'pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,which_freq,v_out));
  end
  clear smax comp_var
  
  load([indir sprintf('pconn_cnt_preproc_viewcomp_vars_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,which_freq,v_out)])
  
  l = 0;
  
  cnt = 1;
  il = 0;
  
  while il < subpl
    
    il = il + 1;
    i = (cnt-1)*subpl+il;
    
    if mod(i-1,subpl)==0
      if exist('manpos')
        figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
      else
        figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
      end
      l = l + 1;
    end
    
    % PLOT POWER SPECTRUM
    if strcmp(hilow,'hi')
      subcomp{1}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-2);
      loglog(freq(find(freq>40,1,'first'):end),smoothed(find(freq>40,1,'first'):end,i)); grid on;
      set(gca,'TickDir','out','XTick',[freq(find(freq>50,1,'first')) freq(find(freq>100,1,'first')) freq(find(freq>150,1,'first'))],'XTickLabel',[50 100 150])
      xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');
      axis tight
    else
      subcomp{1}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-2);
      loglog(freq(find(freq>2,1,'first'):find(freq>40,1,'first')),smoothed(find(freq>2,1,'first'):find(freq>40,1,'first'),i)); grid on;
      set(gca,'TickDir','out','XTick',[freq(find(freq>10,1,'first')) freq(find(freq>30,1,'first'))],'XTickLabel',[10 30])
      xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');
      axis tight
    end
    
    % PLOT VARIANCE OVER TIME
    subcomp{2}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-1);
    
    eval(sprintf('scatter(1:smax,comp_var(i,:),''k.'');'))
    xlabel('Time'); ylabel('Variance');
    axis tight
    
    %     clear smoothed comp_var
    
    % PLOT COMPONENT TOPOGRAPHY
    subcomp{3}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3);
    cfg = [];
    cfg.component = [i];       % specify the component(s) that should be plotted
    cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.highlight = 'off';
    cfg.marker    = 'off';
    
%     [~,in] = max(abs(comp_low.topo(:,i)))
    
%     pars = [];
%     pars.markersize = 0;
%     pars.linewidth = 5;
%     pars.cbar =0;
%     range = eval(sprintf('max(abs(min(comp_%s.topo(:,i))),abs(max(comp_%s.topo(:,i))))',hilow,hilow));
%     pars.scale = [-range+(0.2*range) range-(0.2*range)];
%     
    eval(sprintf('ft_topoplotIC(cfg, comp_%s);',hilow));
%     eval(sprintf('showfield(comp_%s.topo(:,i),sa.locs_2D,pars);',hilow));

    if strcmp(hilow,'hi')
      maxnum = size(comp_hi.trial{1},1)-3;
    else
      maxnum = size(comp_low.trial{1},1)-3;
    end
    
    if mod(i,subpl)==0 || i == maxnum+3
      
      pos = [0.76 0.73 0.075 0.035; ...
        0.76 0.51 0.075 0.035; ...
        0.76 0.29 0.075 0.035; ...
        0.76 0.07 0.075 0.035];
      
      rej_callback1 = ['if (rej_comp(i-3) == 0), set(rej1,''Backgroundcolor'',''r''),rej_comp(i-3)=1;' ...
        'else set(rej1,''Backgroundcolor'',''g''), rej_comp(i-3)=0; end'];
      rej_callback2 = ['if (rej_comp(i-2) == 0), set(rej2,''Backgroundcolor'',''r''),rej_comp(i-2)=1;' ...
        'else set(rej2,''Backgroundcolor'',''g''), rej_comp(i-2)=0; end'];
      rej_callback3 = ['if (rej_comp(i-1) == 0), set(rej3,''Backgroundcolor'',''r''),rej_comp(i-1)=1;' ...
        'else set(rej3,''Backgroundcolor'',''g''), rej_comp(i-1)=0; end'];
      rej_callback4 = ['if (rej_comp(i-0) == 0), set(rej4,''Backgroundcolor'',''r''),rej_comp(i-0)=1;' ...
        'else set(rej4,''Backgroundcolor'',''g''), rej_comp(i-0)=0; end'];
      
      for ibgc = 1 : 4
        if rej_comp(i+(ibgc-4)) == 1
          bgc{ibgc} = 'r';
        else
          bgc{ibgc} = 'g';
        end
      end
      
      rej1 = uicontrol('Units','normalized','Position',pos(1,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{1},'Callback',rej_callback1);
      rej2 = uicontrol('Units','normalized','Position',pos(2,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{2},'Callback',rej_callback2);
      rej3 = uicontrol('Units','normalized','Position',pos(3,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{3},'Callback',rej_callback3);
      rej4 = uicontrol('Units','normalized','Position',pos(4,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{4},'Callback',rej_callback4);
      
      
      tc1_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-3]; ft_databrowser(cfg, comp_%s);',hilow);
      tc2_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-2]; ft_databrowser(cfg, comp_%s);',hilow);
      tc3_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-1]; ft_databrowser(cfg, comp_%s);',hilow);
      tc4_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-0]; ft_databrowser(cfg, comp_%s);',hilow);
      
      
      tc1 = uicontrol('Units','normalized','Position',[0.86 0.73 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc1_cb);
      tc2 = uicontrol('Units','normalized','Position',[0.86 0.51 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc2_cb);
      tc3 = uicontrol('Units','normalized','Position',[0.86 0.29 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc3_cb);
      tc4 = uicontrol('Units','normalized','Position',[0.86 0.07 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc4_cb);
      
      % SAVE COMPONENTS
      % ------------------------------------------------
      sc_cb1 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{1},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{1},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{1},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots_f%d_s%d_b%d_comp%d_v%d.pdf'']); close(h)',which_freq,isubj,iblock,1,v);
      sc_cb2 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{2},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{2},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{2},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots_f%d_s%d_b%d_comp%d_v%d.pdf'']); close(h)',which_freq,isubj,iblock,2,v);
      sc_cb3 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{3},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{3},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{3},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots_f%d_s%d_b%d_comp%d_v%d.pdf'']); close(h)',which_freq,isubj,iblock,3,v);
      sc_cb4 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{4},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{4},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{4},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots_f%d_s%d_b%d_comp%d_v%d.pdf'']); close(h)',which_freq,isubj,iblock,4,v);
      
      savecomp1 = uicontrol('Units','normalized','Position',[0.86 0.78 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb1);
      savecomp2 = uicontrol('Units','normalized','Position',[0.86 0.56 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb2);
      savecomp3 = uicontrol('Units','normalized','Position',[0.86 0.34 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb3);
      savecomp4 = uicontrol('Units','normalized','Position',[0.86 0.12 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb4);
      % ------------------------------------------------
      
      if i > 4
        prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','Prev','Callback','cnt = cnt - 1; il = 0; l = l - 2; manpos = get(gcf,''Position''); close');
      else
        prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','');
      end
      if i < maxnum
        next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','Next','Callback','cnt = cnt + 1; il = 0; close');
      else
        next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','');
      end
      
      save_callback = ['idx = find(rej_comp==1); save(sprintf([outdir ''pconn_cnt_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat''],isubj,im,iblock,which_freq,v_out),''idx'',''rej_comp''); close;'];
      
      s = uicontrol('Units','normalized','Position',[0.90 0.01 0.075 0.05],'Style','pushbutton','String','Save','Callback',save_callback);
      err = 0;
      quit = uicontrol('Units','normalized','Position',[0.80 0.01 0.075 0.05],'Style','pushbutton','String','Quit','Callback','close; err = 1;');
      orient landscape
      uiwait
      
      if err == 1
        error('Abort!')
      end
    end
  end
end


