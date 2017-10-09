%% PCONN SUBTRACT ICA COMPONENTS
%  pconn_cnt_postpostproc


clear all
restoredefaultpath
d         = 'data_low.trial{1}+data_hi.trial{1};';

v = 1;


%%
for im = 3 : 3
  for isubj = 33
    for iblock = 1 : 1
%       
%       if exist(['~/pconn_cnt/proc/preproc/' sprintf('pconn_cnt_postpostproc_s%d_m%d_b%d_v%d_processing.txt',isubj,im,iblock,v)])
%         continue
%       else
%         system(['touch ' '~/pconn_cnt/proc/preproc/' sprintf('pconn_cnt_postpostproc_s%d_m%d_b%d_v%d_processing.txt',isubj,im,iblock,v)]);
%       end
      
      disp(sprintf('Processing s%d m%d b%d ...',isubj,im,iblock))
           
      load(sprintf(['~/pconn_cnt/proc/preproc/' 'pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,v));
      
      % this is necessary as *preproc_proc.m does not create low freq
        data_low.trial{1} = eval(d);
      
         
      clear data_hi
      
       data = data_low; clear data_low
      
      save(sprintf(['~/pconn_cnt/proc/preproc/' 'pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,1),'data');
            
    end
  end
end
        
            
           
    
error('!')
  
  
  
  
  %%
%   win = hanning(400);
%   
%   [p_pup,f_pup]=pwelch(nanmean(data_low.trial{1},1)+nanmean(data_hi.trial{1},1),win,50,2:1:150,400);
%   
%   p_pup(find(f_pup==49):find(f_pup==51)) = NaN;
%   p_pup(find(f_pup==49):find(f_pup==51)) = NaN;
