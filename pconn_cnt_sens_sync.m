%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_cnt_sens_sync

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.


% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_rawdata = 1;
is_src    = 0;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
foi       = [2 4; 4 8; 8 12; 12 36];
filt_ord  = 2;
i_fit     = [1 100];
i_calc    = [0.5 150];
dfa_overlap = 0.5;
% --------------------------------------------------------
% VERSION 2 (with DFA)
% --------------------------------------------------------


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1 : length(foi)
      
      if ~exist(sprintf([outdir 'pconn_cnt_sens_sync_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_sens_sync_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      if isubj == 3
        load ~/pconn/matlab/pconn_sensorlabels.mat
      end
      
      for iblock = 1 : 2
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
                
        if isubj == 3
          [~,idx_lab]=intersect(data_low.label,lab);
          data_low.trial{1}=data_low.trial{1}(idx_lab,:);
          data_low.label = lab;
          save(['~/pconn/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
        end
                
        clear data_hi
        
        [mydata,~] = megdata2mydata(data_low);
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        
        % compute bp-filtered signal
        
        tmp    = nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1));
        
        pha    = angle(hilbert(tmp));
        
        r      = abs(sum(exp(i*pha),2)/268);
        
        dfa = nbt_doDFA(r, siginfo, i_fit,i_calc,dfa_overlap,0,0,[]);        
        
        save(sprintf([outdir 'pconn_cnt_sens_sync_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'dfa','r','-v7.3');
        
        clear r tmp pha dfa mydata
        
      end  
    end
  end
end


error('STOP')

%%

clear all_sync_std all_sync_mean all_sync_dfa

v = 1;

addpath /home/tpfeffer/pconn/matlab/

ord       = pconn_randomization;

ifoi = 1;

for m = 1 : 3
  for isubj = SUBJLIST
    	
    im = find(ord(isubj,:)==m);

    for iblock = 1 : 2

      load(sprintf('/home/tpfeffer/pconn_cnt/proc/dfa/pconn_cnt_sens_sync_s%d_b%d_m%d_f%d_v%d.mat',isubj,iblock,im,ifoi,v))

      all_sync_std(iblock,m,isubj)  = std(r); 
    	all_sync_mean(iblock,m,isubj) = mean(r);
      all_dfa(iblock,m,isubj) = dfa.MarkerValues;
      
      clear dfa r
      
      load(sprintf('/home/tpfeffer/pconn/proc/dfa/pconn_sens_sync_s%d_b%d_m%d_f%d_v%d.mat',isubj,iblock,im,ifoi,v))

      all_sync_std_rest(iblock,m,isubj)  = std(r); 
    	all_sync_mean_rest(iblock,m,isubj) = mean(r);
      all_dfa_rest(iblock,m,isubj) = dfa.MarkerValues;
      
     	clear dfa 
      
    end
  end
end

all_sync_std  = all_sync_std(:,:,SUBJLIST);
all_sync_mean = all_sync_mean(:,:,SUBJLIST);
all_sync_dfa  = all_dfa(:,:,SUBJLIST);

all_sync_std_rest  = all_sync_std_rest(:,:,SUBJLIST);
all_sync_mean_rest = all_sync_mean_rest(:,:,SUBJLIST);
all_sync_dfa_rest  = all_dfa_rest(:,:,SUBJLIST);


%

par_std  = [squeeze(nanmean(all_sync_std(:,1,:),1)) squeeze(nanmean(all_sync_std(:,2,:),1)) squeeze(nanmean(all_sync_std(:,3,:),1))];
par_mean = [squeeze(nanmean(all_sync_mean(:,1,:),1)) squeeze(nanmean(all_sync_mean(:,2,:),1)) squeeze(nanmean(all_sync_mean(:,3,:),1))];
par_dfa  = [squeeze(nanmean(all_sync_dfa(:,1,:),1)) squeeze(nanmean(all_sync_dfa(:,2,:),1)) squeeze(nanmean(all_sync_dfa(:,3,:),1))];

par_std_rest  = [squeeze(nanmean(all_sync_std_rest(:,1,:),1)) squeeze(nanmean(all_sync_std_rest(:,2,:),1)) squeeze(nanmean(all_sync_std_rest(:,3,:),1))];
par_mean_rest = [squeeze(nanmean(all_sync_mean_rest(:,1,:),1)) squeeze(nanmean(all_sync_mean_rest(:,2,:),1)) squeeze(nanmean(all_sync_mean_rest(:,3,:),1))];
par_dfa_rest  = [squeeze(nanmean(all_sync_dfa_rest(:,1,:),1)) squeeze(nanmean(all_sync_dfa_rest(:,2,:),1)) squeeze(nanmean(all_sync_dfa_rest(:,3,:),1))];
%%

% STD OF R
[~,p_std(1)]=ttest(par_std(:,2),par_std(:,1));
[~,p_std(2)]=ttest(par_std(:,3),par_std(:,1));
[~,p_std(3)]=ttest(par_std(:,3),par_std(:,2));

% MEAN OF R
[~,p_mean(1)]=ttest(par_mean(:,2),par_mean(:,1));
[~,p_mean(2)]=ttest(par_mean(:,3),par_mean(:,1));
[~,p_mean(3)]=ttest(par_mean(:,3),par_mean(:,2));

% DFA OF R
[~,p_dfa(1)]=ttest(par_dfa(:,2),par_dfa(:,1));
[~,p_dfa(2)]=ttest(par_dfa(:,3),par_dfa(:,1));
[~,p_dfa(3)]=ttest(par_dfa(:,3),par_dfa(:,2));


% BAR PLOTS

m_std  = nanmean(par_std);
m_mean = nanmean(par_mean);
m_dfa  = nanmean(par_dfa);

sem_std  = nanstd(par_std)/sqrt(length(par_std));
sem_mean = nanstd(par_mean)/sqrt(length(par_std));
sem_dfa  = nanstd(par_dfa)/sqrt(length(par_std));

%

figure; set(gcf,'color','white'); hold on; title(sprintf('STD(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_std(1),p_std(2),p_std(3)))

lim = [min(m_std-sem_std)-0.1*min(m_std-sem_std) max(m_std+sem_std)+0.1*max(m_std+sem_std)];

bar([1 1.5 2],m_std,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_std(1)-sem_std(1) m_std(1)+sem_std(1)],'linewidth',5)
line([1.5 1.5],[m_std(2)-sem_std(2) m_std(2)+sem_std(2)],'linewidth',5)
line([2 2],[m_std(3)-sem_std(3) m_std(3)+sem_std(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_std_f%d_v%d.eps',ifoi,v))

figure; set(gcf,'color','white'); hold on; title(sprintf('MEAN(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_mean(1),p_mean(2),p_mean(3)))

lim = [min(m_mean-sem_mean)-0.1*min(m_mean-sem_mean) max(m_mean+sem_mean)+0.1*max(m_mean+sem_mean)];

bar([1 1.5 2],m_mean,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_mean(1)-sem_mean(1) m_mean(1)+sem_mean(1)],'linewidth',5)
line([1.5 1.5],[m_mean(2)-sem_mean(2) m_mean(2)+sem_mean(2)],'linewidth',5)
line([2 2],[m_mean(3)-sem_mean(3) m_mean(3)+sem_mean(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_mean_f%d_v%d.eps',ifoi,v))

figure; set(gcf,'color','white'); hold on; title(sprintf('DFA(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_dfa(1),p_dfa(2),p_dfa(3)))

lim = [min(m_dfa-sem_dfa)-0.1*min(m_dfa-sem_dfa) max(m_dfa+sem_dfa)+0.1*max(m_dfa+sem_dfa)];

bar([1 1.5 2],m_dfa,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_dfa(1)-sem_dfa(1) m_dfa(1)+sem_dfa(1)],'linewidth',5)
line([1.5 1.5],[m_dfa(2)-sem_dfa(2) m_dfa(2)+sem_dfa(2)],'linewidth',5)
line([2 2],[m_dfa(3)-sem_dfa(3) m_dfa(3)+sem_dfa(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_dfa_f%d_v%d.eps',ifoi,v))

%% COUNT AND SYNC
addpath ~/pcbi/
cnt = pcbi_cnt(SUBJLIST);
cnt = nanmean(cnt,2);
[r,p]=corr(nanmean(par_mean,2),cnt)




