function [RandResult] = kRandCorr(DataA,DataB, alphaThresh,  NPerms,Type)
%RandResult = randomisation_test2dABS(DataA,DataB, alpha,  NPerms)
%
% Linear correlation  omnibus corrected
% for multiple comparisons between DataA and DataB
%
%Inputs - DataA = Subjects x Variables (e.g. Timepoins)
%
%       - DataB = A row vector to correlate each variable with i.e 1 x Subjects
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
% 
% Krish August 2017


N1 = size(DataA, 1);

NVars = max([ size(DataA, 2) size(DataB, 2)]);

RMax = zeros([1 NPerms]);

%Calculate the Obtained r value
[rseries,Ptrue]=corr(DataA,DataB,'rows','pairwise','type',Type);


Pcount = zeros([1 NVars]);

%Now run the permutations
fprintf(1,'\n');

[~,puc]=corr(DataA,DataB,'rows','pairwise','type',Type);

for i = 1 : NPerms
   
        permuted_data = aconnectivity.shuffle(DataB, 1);
   
        [R,P]=corr(DataA,permuted_data,'rows','pairwise','type',Type);
        
        
        Rmax(i) = max(abs(R(:)));
        RSignedmax(i) = max((R(:)));
        RSignedmin(i) = min((R(:)));

        if  mod(i, floor(NPerms/20)) == 0
            fprintf('Calculated permutation %d of %d, abs maximum=%f\r', i, NPerms,Rmax(i));
        end
    
    for j = 1 : length(rseries),
        Val=R(j);
        if abs(Val) > abs(rseries(j))
          Pcount(j) = Pcount(j) + 1; 
       end
    end
end

%The histogram of obtained t values
r_hist_max = sort(Rmax);
r_hist_signedmax = sort(RSignedmax);
r_hist_signedmin = sort(RSignedmin,'descend');


%The critical value
  thresh_index = round((1 - (alphaThresh)) * NPerms);
  ActualPThresh=1-(thresh_index/NPerms);
  fprintf(1,'alpha=%f NPerms=%d ThreshIndex=%d ActualPThresh=%f\n',alphaThresh,NPerms,thresh_index,ActualPThresh);
  thresh_signed = round((1 - (alphaThresh/2)) * NPerms);
  
  rcrit =  r_hist_max(thresh_index);
  rcrit_pos=r_hist_signedmax(thresh_signed);
  rcrit_neg=r_hist_signedmin(thresh_signed);
  
fprintf(1,'\nrcrit(abs)=%f (actual pthrsh=%f) rcrit_neg=%f rcrit_pos=%f\n',rcrit,ActualPThresh,rcrit_neg,rcrit_pos);

%Now create the Omnibus thresholded time-series

rseries_corr = rseries;

for i = 1 : length(rseries_corr(:)),
    pseries_corr(i)=1;

    NumBigger=aconnectivity.findthenearest(abs(rseries(i)),r_hist_max);
    if (isempty(NumBigger)),
        NumBigger=NPerms;
    end
    if length(NumBigger)>1
        NumBigger=NumBigger(1);
    end
    
    pseries_corr(i)=1-(NumBigger/NPerms);
    
     if abs(rseries_corr(i)) < rcrit       
         rseries_corr(i) = 0;
     end
end
pseries_corr(isnan(rseries))=1;
rseries_corr(isnan(rseries))=0;
rseries(isnan(rseries))=0;

%pseries_uncorr = Pcount / NPerms;
pseries_uncorr = puc;

    fprintf(1,'\nFIXED Max abs(R)=%f, R-critical=%f, Number of R-statistics that make omnibus significance=%d\n',max(abs(rseries(:))),rcrit,length(find(rseries_corr~=0)));

% reshape p series corr
pseries_corr = reshape(pseries_corr,[size(DataA,2),size(DataB,2)]);
 
RandResult.rseries=rseries; 
RandResult.rseries_corr=rseries_corr; 
RandResult.rcrit=rcrit; 
RandResult.pseries_uncorr=pseries_uncorr;
RandResult.pseries_corr=pseries_corr;
RandResult.r_hist_max=r_hist_max;
RandResult.r_hist_signedmax=r_hist_signedmax;
RandResult.r_hist_signedmin=r_hist_signedmin;
RandResult.Ptrue=Ptrue;
RandResult.Type=Type;





