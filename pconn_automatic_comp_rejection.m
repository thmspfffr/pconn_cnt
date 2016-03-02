%% CORRECTS STUFF

cnt = 1;


if cnt
  ss = 'pconn_cnt';

  indir   = '~/pconn_cnt/proc/preproc/';
  outdir  =  '~/pconn_cnt/proc/preproc/';
  
  v_ica = 1;
  v_in = [1 4];
  v_out = 5;
  iff = 1;
  
else
  
  ss = 'pconn';
  
 	indir   = '~/pconn/proc/preproc/';
  outdir  =  '~/pconn/proc/preproc/';

  v_ica = 3;
  v_in = [3];
  v_out = 4;
  iff = 1;
  
end




for isubj = 21
  for im = 1 : 3
    for iblock = 1 : 2
      
    fprintf('Processing subj%d m%d b%d ...\n',isubj,im,iblock)

    load([indir sprintf('%s_preproc_ica_s%d_m%d_b%d_v%d.mat',ss,isubj,im,iblock,v_ica)])
    
    for icomp = 1 : 128
      rej(icomp) = std(abs(fft(comp_low.topo(:,icomp))));
      if isubj == 21
        [~,in] = max(abs(comp_low.topo(:,icomp)));
        if any(in==[123 124 125 126]) || any(in==[237 244])
          rej(icomp) = 1;
          icomp
        end
      end
    end
    
    

    rej = find(rej<mean(rej)/2);

    if ~exist(sprintf([outdir '%s_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],ss,isubj,im,iblock,iff,v_in(1)))
      
      rej_comp = zeros(128,1);
      rej_comp(rej) = 1;
      
    else
      load(sprintf([outdir '%s_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],ss,isubj,im,iblock,iff,v_in(1)));

      rej_comp(rej) = 1;
    end
    
    save(sprintf([outdir '%s_rejected_comps_s%d_m%d_b%d_f%d_v%d.mat'],ss,isubj,im,iblock,iff,v_out),'rej_comp');
 
    end
  end
end