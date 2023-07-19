function [RandResult] = orig_kRandTest(DataA,DataB, alphaThresh,  NPerms, UseMean)
%RandResult = randomisation_test2dABS(DataA,DataB, alpha,  NPerms)
%
% Randomisation t-test (twosample unpaired or onesample), omnibus corrected
% for multiple comparisons.
% If you supply DataA and DataB this will run a 2-sample unpaired test
% If you set DataB to be empty i.e. [] it runs a 1-sample test
%
% If you want to run a paired-t-test, generate the differences socores
% DataA-DataB for each subject and then run this as a 1-sample t-test on the
% differences - this is equivalent to a paired-t test
%
%Inputs - DataA = Subjects x Variables (e.g. Timepoins)
%
%       - DataB = Subjects x Variables (e.g. Timepoints) or [] for a 1-sample test on DataA
%
%         alpha - your desired significance level   
%         NPerms = NumbJUNEer of permutations to resample
%                  I suggest 5000 or more;
%
%
%Outputs 
%        RandResult.t_series = Raw t test values (no thresholding performed)
%
%        RandResult.tseries_corr = t test values (omnibus corrected)
%
%        RandResult.t_crit = the critical t test value score
%
%        RandResult.pseries_uncorr - the pvalues, uncorrected per test
%        RandResult.pseries_corr - the p-values, omnibus corrected for multiple
%        comparisons
%        RandResult.t_hist_max  - the omnibus max-T distribution for each
%        randomisation, as a histogram.

%        Note this script searches for the ABS maximum t-statistic on each
%        randomisation stage to define the omnibus test null distribution.
%        Then abs(T-scores) are tested against this null.
%
%
% SM Jul 2010
% Krish August 2016

Mean=0;
try
    if (UseMean~=0),
        fprintf(1,'Forcing mean instead of TTEST\n');
        Mean=1;
    end
end

N1 = size(DataA, 1);

if isempty(DataB) == 1
   Data_Cat = DataA;
   NTot = N1;
   OneSample=1;
else
    N2 = size(DataB, 1);
    NTot = N1 + N2;
    Data_Cat = [DataA ; DataB];
    OneSample=0;
end

NVars = size(DataA, 2);

TMax = zeros([1 NPerms]);

%Calculate the Obtained t value
if (OneSample==0),
    if (Mean==0),
        [H,P,CI,STATS] = ttest2(Data_Cat(1:N1,:),Data_Cat(N1+1:end,:));
    end
     if (Mean~=0),
        TrueMean = nanmean(Data_Cat(1:N1,:))-nanmean(Data_Cat(N1+1:end,:)); %fixed! Should be nanmean...
     end
end

if (OneSample~=0),
    if (Mean==0),
        [H,P,CI,STATS] = ttest(Data_Cat);
    end
    if (Mean~=0),
        TrueMean = nanmean(Data_Cat);
    end

end


if (Mean==0),
    tseries = STATS.tstat;
end

if (Mean~=0),
    tseries=TrueMean;
end


Pcount = zeros([1 NVars]);

%Now run the permutations
fprintf(1,'\n');

for i = 1 : NPerms
   
   if(OneSample==0),
        permuted_data = shuffle(Data_Cat, 1);
   
        if (Mean==0),
            [H,P,CI,STATS] = ttest2(permuted_data(1:N1,:),permuted_data(N1+1:end,:));
        end
        
        if (Mean~=0),
            PermMean = nanmean(permuted_data(1:N1,:))-nanmean(permuted_data(N1+1:end,:));%Fixed! Should be nanmean!
        end
        

   end
   
   if(OneSample~=0),
        RandSigns=sign(randn(NTot,1));
        RandSignsRep=repmat(RandSigns,[1,NVars]);
        permuted_data = RandSignsRep.*Data_Cat;
        if(Mean==0),
            [H,P,CI,STATS] = ttest(permuted_data);
        end
        if (Mean~=0),
            PermMean=nanmean(permuted_data);
        end
   end
   
   if (Mean==0),
        Tmax(i) = max(abs(STATS.tstat));
        TSignedmax(i) = max((STATS.tstat));
        TSignedmin(i) = min((STATS.tstat));
   end
   if (Mean~=0),
        Tmax(i) = max(abs(PermMean));
        TSignedmax(i) = max((PermMean));
        TSignedmin(i) = min((PermMean));
  end
   
  
   if  mod(i, floor(NPerms/20)) == 0
      fprintf('Calculated permutation %d of %d, abs maximum=%f\r', i, NPerms,Tmax(i));
   end
    
    for j = 1 : length(tseries),
       if(Mean==0),
           Val=STATS.tstat(j);
       end
       if(Mean~=0),
           Val=PermMean(j);
       end
       if abs(Val) > abs(tseries(j))
          Pcount(j) = Pcount(j) + 1; 
       end
    end
end

%The histogram of obtained t values
t_hist_max = sort(Tmax);
t_hist_signedmax = sort(TSignedmax);
t_hist_signedmin = sort(TSignedmin,'descend');


%The critical value
  thresh_index = round((1 - (alphaThresh)) * NPerms);
  ActualPThresh=1-(thresh_index/NPerms);
  fprintf(1,'alpha=%f NPerms=%d ThreshIndex=%d ActualPThresh=%f\n',alphaThresh,NPerms,thresh_index,ActualPThresh);
  thresh_signed = round((1 - (alphaThresh/2)) * NPerms);
  
  tcrit =  t_hist_max(thresh_index);
  tcrit_pos=t_hist_signedmax(thresh_signed);
  tcrit_neg=t_hist_signedmin(thresh_signed);
  
fprintf(1,'\ntcrit(abs)=%f (actual pthrsh=%f) tcrit_neg=%f tcrit_pos=%f\n',tcrit,ActualPThresh,tcrit_neg,tcrit_pos);

%Now create the Omnibus thresholded time-series

tseries_corr = tseries;

for i = 1 : length(tseries_corr(:)),
    pseries_corr(i)=1;
    %NumBigger=length(find(t_hist_max <= abs(tseries(i))));
%     Test=abs(tseries(i));
%     NumBigger=0;
%     for j=1:length(t_hist_max)
%         if(t_hist_max(j)<Test),
%             NumBigger=NumBigger+1;
%         end
%     end
    NumBigger=findthenearest(abs(tseries(i)),t_hist_max);
    if (isempty(NumBigger)),
        NumBigger=NPerms;
    end
    if length(NumBigger)>1
        NumBigger=NumBigger(1);
    end
    
    pseries_corr(i)=1-(NumBigger/NPerms);
    
     if abs(tseries_corr(i)) < tcrit       
         tseries_corr(i) = 0;
     end
%      if tseries_corr(i) < tcrit_pos && tseries_corr(i) > tcrit_neg,       
%           tseries_corr(i) = 0;
%      end
end
pseries_corr(isnan(tseries))=1;
tseries_corr(isnan(tseries))=0;
tseries(isnan(tseries))=0;

pseries_uncorr = Pcount / NPerms;

if (Mean==0)
    fprintf(1,'\nFIXED Max abs(T-stat)=%f, T-critical=%f, Number of t-statistics that make omnibus significance=%d\n',max(abs(tseries(:))),tcrit,length(find(tseries_corr~=0)));
end

if (Mean~=0)
    fprintf(1,'\nFIXED Max abs(TrueMean)=%f, Mean-critical=%f, Number of Values that make omnibus significance=%d\n',max(abs(tseries(:))),tcrit,length(find(tseries_corr~=0)));
end

RandResult.tseries=tseries; 
RandResult.tseries_corr=tseries_corr; 
RandResult.tcrit=tcrit; 
RandResult.pseries_uncorr=pseries_uncorr;
RandResult.pseries_corr=pseries_corr;
RandResult.t_hist_max=t_hist_max;
RandResult.t_hist_signedmax=t_hist_signedmax;
RandResult.t_hist_signedmin=t_hist_signedmin;
RandResult.UseMeanInsteadOfT=Mean;




