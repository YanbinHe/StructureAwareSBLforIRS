% =========================================================================
% Title       : Consolidate simulation results
% File        : consolidate.m
% -------------------------------------------------------------------------
% Description :
%   Collects specified simulation results and averages bit-error rate
%   (BERs) and frame-error rate (FER) for all files. 
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   11-dec-11  1.1      studer  cleanup for reproducible research
%   10-mar-08  1.0      studer  file created
% -------------------------------------------------------------------------
%   (C) 2006-2011 Communication Theory Group                      
%   ETH Zurich, 8092 Zurich, Switzerland                               
%   Author: Dr. Christoph Studer (e-mail: studer@rice.edu)     
% =========================================================================

function Out = consolidate(filename,maxfiles)

  % -- read in specified files
  disp(sprintf('### consolidate %s',filename));
  k = 0;
  for i=1:maxfiles
    fullname = ['results/',filename,'_',num2str(i-1),'.mat'];    
    if exist(fullname)>0
      tmp = load(fullname);
      if k==0
        ptype = tmp;
        BER = ptype.Results.BER;
        FER = ptype.Results.FER;
        Trials = ptype.Results.Trials;
      else
        BER = BER + tmp.Results.BER;
        FER = FER + tmp.Results.FER;
        Trials = Trials + tmp.Results.Trials;
      end
      k=k+1;
    end
  end    
  
  % -- compute output  
  ptype.Results.BER = BER / k;
  ptype.Results.FER = FER / k;
  ptype.Results.Trials = Trials;
  disp(sprintf('### %i trials of %i completed',Trials,k*ptype.Results.TxRx.Sim.nr_of_channels)); 
  disp(sprintf('### %i file(s) consolidated\n',k)); 
  Out = ptype.Results; 
  
return
  