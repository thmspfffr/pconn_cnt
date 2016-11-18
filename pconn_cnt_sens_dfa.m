%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_cnt_sens_dfa

clear all

for v = [2 8 24 84] 
  if v == 2
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v         = 2;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2 4; 4 8; 8 12; 12 36; 50 100];
    i_fit     = [3 50];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 8
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v         = 8;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2 4; 4 8; 8 12; 12 36; 50 100];
    i_fit     = [5 75];
    dfa_overlap = 0.5;
    filt_ord = 2;
    % --------------------------------------------------------
 	elseif v == 24
    % --------------------------------------------------------
    % VERSION 24 - COMPUTE DFA ACROSS FREQUENCIES
    % --------------------------------------------------------
    v         = 24;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2:2:148; 4:2:150]';
    i_fit     = [3 50];
    dfa_overlap = 0.5;
    filt_ord = 2;
  elseif v == 84
    % --------------------------------------------------------
    % VERSION 24 - COMPUTE DFA ACROSS FREQUENCIES
    % --------------------------------------------------------
    v         = 84;
    v_rawdata = 6;
    is_src    = 0;
    fsample   = 400;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    foi       = [2:2:148; 4:2:150]';
    i_fit     = [5 75];
    dfa_overlap = 0.5;
    filt_ord = 2;
  end
  
  addpath ~/pconn/matlab/
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
  
  outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
  plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';
  
  run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
  %%
  for m = 1 : 3
    for isubj = SUBJLIST
      for ifoi = 1 : length(foi)
        
        if ~exist(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
          system(['touch ' outdir sprintf('pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
        else
          continue
        end
        %
        if isubj == 3
          load ~/pconn/matlab/pconn_sensorlabels.mat
        end
        
        if (isubj~=3 && isubj~=2 &&  isubj~=17) && isubj <= 24
          d = dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b*_v%d.mat',isubj,m,v_rawdata));
        else
          d = dir(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b*_v%d.mat',isubj,m,v_rawdata));
        end
        
        if ~isempty(d)
        
        for iblock = 1 : length(d)
          
          disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
          
          load(['/home/tpfeffer/pconn_cnt/proc/preproc/' d(iblock).name])
          
          if isubj == 3
            [~,idx_lab]=intersect(data_low.label,lab);
            data_low.trial{1}=data_low.trial{1}(idx_lab,:);
            data_low.label = lab;
            save(['~/pconn_cnt/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
          end
          
          clear data_hi clear
          [mydata,epleng] = megdata2mydata(data_low);
          
          mydata = mydata(1:end-1000,:);
          
          siginfo = nbt_Info;
          siginfo.converted_sample_frequency = fsample;
          
          % compute bp-filtered signal
          tmp    = single(nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1)));
          ampenv = abs(hilbert(tmp));
        	tmp_dfa   = tp_dfa(ampenv,i_fit,400,dfa_overlap,15);

          clear data_low
          
          if isubj >= 32
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
          elseif isubj < 4
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
          elseif isubj == 17
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
          else
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
          end
          
          par.dfa(:,iblock)   = pconn_sens_interp274(idx,tmp_dfa);
          par.amp(:,iblock)   = pconn_sens_interp274(idx,mean(ampenv));
          par.var(:,iblock)   = pconn_sens_interp274(idx,var(ampenv));
          par.cvar(:,iblock)  = sqrt(par.var(:,iblock))./(par.amp(:,iblock));
          
          clear dfa dat ampenv r amp
          
        end
        
        save(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
        clear par r
        
        end
        
      end
    end
  end
end
  
  error('STOP')
  
  %% PLOT ACROSS FREQ
ord = pconn_randomization;

for isubj = SUBJLIST
  for m = 1 : 3
  	im = find(ord(isubj,:)==m);

    for ifoi = 1 : 74
      
      load(sprintf(['~/pconn_cnt/proc/dfa/' 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,4));
      
      alldfa(isubj,m,ifoi) = nanmean(nanmean(par.dfa,2));
      allvar(isubj,m,ifoi) = nanmean(nanmean(par.var,2));
      allcvar(isubj,m,ifoi) = nanmean(nanmean(par.cvar,2));
      allamp(isubj,m,ifoi) = nanmean(nanmean(par.amp,2));

   
    end
  end
end

alldfa = alldfa(SUBJLIST,:,:);
allvar = allvar(SUBJLIST,:,:);
allcvar = allcvar(SUBJLIST,:,:);
allamp = allamp(SUBJLIST,:,:);

  %% 
  
f = mean(foi,2); 
col = {[0.7 0.7 0.7];[1 0.4 0];[0 0.5 1]}

figure;set(gcf,'color','w'); hold on; set(gca,'tickdir','out')
alltitle = {'DFA';'VAR';'CVAR';'AMP'}

for iall = 1 : 4
  if iall == 1
    par = alldfa;
  elseif iall == 2
    par = allvar;
  elseif iall == 3
    par = allcvar;
  elseif iall == 4
    par = allamp;
  end
  
  subplot(2,2,iall); hold on; title(alltitle{iall})

  for i = 1 : 3
    plot(log2(f),log2(squeeze(nanmean(par(:,i,1:end),1))),'color',col{i},'linewidth',3)
  end

%   t = logical(squeeze(ttest(par(:,2,:),par(:,1,:),'dim',1)));
%   plot(log2(f(t)),log2(ones(sum(t),1)*min(par(:))),'*','color',col{2})
%   t = logical(squeeze(ttest(par(:,3,:),par(:,1,:),'dim',1)));
%   plot(log2(f(t)),log2(ones(sum(t),1)*min(par(:))),'*','color',col{3})
  
  set(gca,'xtick',[1  3  5  7],'xticklabel',[2  8  32  128],'tickdir','out')
  set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
end

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_para_across_freq_v%d.eps',v))

figure;set(gcf,'color','w'); hold on;
alltitle = {'DFA';'VAR';'CVAR';'AMP'}

for iall = 1 : 4
  if iall == 1
    par = alldfa;
  elseif iall == 2
    par = allvar;
  elseif iall == 3
    par = allcvar;
  elseif iall == 4
    par = allamp;
  end
  
  subplot(2,2,iall); hold on; title(alltitle{iall})

  [t1,p1] = ttest(par(:,2,:),par(:,1,:),'dim',1);
  [t2,p2] = ttest(par(:,3,:),par(:,1,:),'dim',1);

      line([1 7],[1.3 1.3],'color','k','linestyle','--')

    plot(log2(f),-log10(squeeze(p1)),'color',col{2},'linewidth',3)
   plot(log2(f),-log10(squeeze(p2)),'color',col{3},'linewidth',3)
    axis([0 8 -0.2 3.5])
set(gca,'xtick',[1 3  5  7],'xticklabel',[2 8  32  128],'tickdir','out')
     
     set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

end

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_para_across_freq_pval_v%d.eps',v))

  
