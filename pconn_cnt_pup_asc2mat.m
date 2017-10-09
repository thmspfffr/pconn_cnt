%% EYE TRACKING DATA
% pconn_cnt_pup_asc2mat

clear

indir   = '~/pconn/rawdata/edf/';
outdir  = '~/pconn_cnt/proc/';

trl = [];
samples=[];
dat=[];
v = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34]
%%
for isubj = SUBJLIST
  disp(sprintf('Processing subject %d ...',isubj))
  % -------------------------------------------------------------------------
  % SAMPLES
  % -------------------------------------------------------------------------
  
  for m = 1:3
%     
    if ~exist(sprintf([outdir 'pconn_cnt_pup_asc2mat_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pconn_cnt_pup_asc2mat_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
%     try
      clear d
      sampledir = [indir sprintf('p%d/s%d/',isubj,m) 'edf_samples/'];
      nblocks = 0;
      
      while 1
        if isubj < 25
          
          if isubj == 6 && m == 1
            for idir = [2 5]
              d_samples = dir([sampledir sprintf('*-s%d_*b%d*.asc',m,idir)]);
              if ~isempty(d_samples)
                nblocks=nblocks+1;
                d{nblocks} = d_samples;
              end
            end
          else
            
            for idir = [3 6]
              d_samples = dir([sampledir sprintf('*-s%d_*b%d*.asc',m,idir)]);
              if ~isempty(d_samples)
                nblocks=nblocks+1;
                d{nblocks} = d_samples;
              end
            end
          end
          break
          
        elseif isubj >= 25
          for idir = [1 4]
            d_samples = dir([sampledir sprintf('*-s%d_*b%d*.asc',m,idir)]);
            if ~isempty(d_samples)
              nblocks=nblocks+1;
              d{nblocks} = d_samples;
            end
          end
          break
        end
      end
           
          
      for ifile = 1 : length(d)
        
        fid = fopen([sampledir d{ifile}.name]);
        
        str     = d{ifile}.name;
        pat     = 'b(\w*)_20';
        t       = regexp(str, pat, 'tokens');
        
        iblockkk  = str2num(cell2mat(t{1}));
        
        if isubj < 25
          if isubj == 6 && m == 1
            if iblockkk == 2
              iblock = 1;
            elseif iblockkk == 5
              iblock = 2;
            end
          else
            if iblockkk == 3
              iblock = 1;
            elseif iblockkk == 6
              iblock = 2;
            end
          end
        else
          if iblockkk == 1
            iblock = 1;
          elseif iblockkk == 4
            iblock = 2;
          end
        end
          
        
        disp(sprintf('Processing samples file %d / %d...',ifile,length(d)))
        data=textscan(fid,'%q%q%q%q%*q','headerlines',0,'delimiter','\t');
        
        fclose(fid);
        
        % Speed improvements possible here!
        for icell = 1 : numel(data)
%           disp(sprintf('Processing column %d ...',icell));
          if icell == 1
            C = data{icell};
            S = sprintf('%s*', C{:});
            samples(:,icell) = sscanf(S, '%f*')';
          else
            samples(:,icell)=zeros(length(data{icell}),1);
            samples(:,icell)=str2double(data{icell});
          end
        end
        
        clear data
        
        save(sprintf([outdir 'pconn_cnt_pup_samples_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'samples','-v7.3');
        
                
        clear temp_samples samples
      end
      % -------------------------------------------------------------------------
      % Events
      % -------------------------------------------------------------------------
      
      eventdir  = [indir sprintf('p%d/s%d/',isubj,m) 'edf_events/'];
%       d_events  = dir([eventdir sprintf('*-s%d_*rest*.asc',m)]);
      
      clear d
      nblocks =0;
      while 1
        if isubj < 25
          
          if isubj == 6 && m == 1
            for idir = [2 5]
              d_events = dir([eventdir sprintf('*-s%d_*b%d*.asc',m,idir)]);
              if ~isempty(d_events)
                nblocks=nblocks+1;
                d{nblocks} = d_events;
              end
            end
          else
            
            for idir = [3 6]
              d_events = dir([eventdir sprintf('*-s%d_*b%d*.asc',m,idir)]);
              if ~isempty(d_events)
                nblocks=nblocks+1;
                d{nblocks} = d_events;
              end
            end
          end
          break
          
        elseif isubj >= 25
          for idir = [1 4]
            d_events = dir([eventdir sprintf('*-s%d_*b%d*.asc',m,idir)]);
            if ~isempty(d_events)
              nblocks=nblocks+1;
              d{nblocks} = d_events;
            end
          end
          break
        end
      end
      
      for ifile = 1 : length(d)
        
        clear fid
        
        fid = fopen([eventdir d{ifile}.name]);
        
      	disp(sprintf('Processing events file %d ...',ifile))

        str     = d{ifile}.name;
        pat     = 'b(\w*)_20';
        t       = regexp(str, pat, 'tokens');
       
        iblockkk  = str2num(cell2mat(t{1}));
        
        if isubj < 25
          if isubj == 6 && m == 1
            if iblockkk == 2
              iblock = 1;
            elseif iblockkk == 5
              iblock = 2;
            end
          else
            if iblockkk == 3
              iblock = 1;
            elseif iblockkk == 6
              iblock = 2;
            end
          end
        else
          if iblockkk == 1
            iblock = 1;
          elseif iblockkk == 4
            iblock = 2;
          end
        end
        
        i=0;j=0;l=0;k=0;trlNum=0;
        
        while ~feof(fid)
          
          tmp = fgetl(fid);
          
          if ~isempty(strfind(tmp,'pretrial_sync'))
            i = i + 1;
            tmp_trl(i,1)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
          elseif ~isempty(strfind(tmp,'posttrial_sync'))
            j = j + 1;
            tmp_trl(j,2)=cell2mat(textscan(tmp, '%*s%d%*[\n]'));
          elseif ~isempty(strfind(tmp, 'EBLINK'))
            k = k + 1;
            buff = textscan(tmp, '%s%s%d%d%*[\n]');
            blinks(k,1:2)=cell2mat(buff(3:4));
            clear buff
          elseif ~isempty(strfind(tmp, 'ESACC'))
            l = l + 1;
            buff = textscan(tmp, '%s%s%d%d%*[\n]');
            saccs(l,1:2)=cell2mat(buff(3:4));
            clear buff
            
          end
        end
        
          if exist('blinks','var') && exist('saccs','var')
            save(sprintf([outdir 'pconn_cnt_pup_events_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'saccs', 'blinks');
          elseif ~exist('blinks','var') && exist('saccs','var')
            save(sprintf([outdir 'pconn_cnt_pup_events_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'saccs');
          elseif exist('blinks','var') && ~exist('saccs','var')
            save(sprintf([outdir 'pconn_cnt_pup_events_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'blinks');
          end
          
          clear saccs blinks
        fclose(fid);     
      end
%     catch me
%         warning('Something went wrong...')
%         disp(me.message);
%     end
  end
  
  clear tmp saccs blinks
end


