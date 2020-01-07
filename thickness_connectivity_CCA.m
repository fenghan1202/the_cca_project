%% 
% Run CCA between connectivity and thickness data
% focus on the 818 subjects with NETMAT data but 32K resolution
% Copyright by Feng Han & Xiao Liu.

%%
clear all;

load('/path/to/thickness/data.mat');%%%%thickness data
load('/path/to/connectivity/data.mat');%%%%connectivity data
load('/path/to/sub_label/data.mat');%%%% subj label
%
%%
%%% input data
var_key = importdata('/path/to/column_headers.txt');        % needs to be a subjects X subjectmeasures text file - see www.fmrib.ox.ac.uk/analysis/HCP-CCA
clear varidx
for li = 1:length(var_key)
      temp = find(ismember(keyname_s,var_key{li}));  
      if isempty(temp)
          varidx(li) = nan;
      else
          varidx(li) = temp;
      end
end
varidx(isnan(varidx)) = [];
%
%
sbjsidx3 = find(ismember(behdat_c(:,1),sbjs));  % index in behavior data
sbjsidx4 = find(ismember(sbjs900,sbjs));            % index in motion data
%
varsQconf= ~strcmp(behdat_c(sbjsidx3,23),'r177');
meanMot = cat(1,RL1.meanMot,RL2.meanMot,LR1.meanMot,LR2.meanMot);
meanMot = meanMot(:,sbjsidx4);
rfmrimotion = nanmean(meanMot);
%
vars = zeros(length(sbjsidx3),length(var_key));
vars(:,[1 3:6 8:end]) = behdat_d(sbjsidx3,varidx);
vars(:,2) = varsQconf;
vars(:,7) = rfmrimotion';
%
sex_all=double(strcmp(behdat_c(:,4),'M'));
vars(:,3)=sex_all(sbjsidx3);

vars(:,4) = behdat_d(sbjsidx3,551);
vars(:,sum(isnan(vars)==0)<60)=NaN;  % pre-delete vars with LOADS of missing data

%% load netmats from HCP PTN release
% set NET to be the subjects X unwrapped network-matrices (partial correlations)
%
netpath = 'HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d200_ts2/';
sbjs1k = textread('HCP_PTN1200/subjectIDs.txt');
sbjs1k = cellstr(num2str(sbjs1k));
sbjsidx = find(ismember(sbjs1k,sbjs)); 

tmp = dlmread([netpath,'netmats2.txt'],' ');
tmp = reshape(tmp,length(sbjs1k),200,200);
tmp = tmp(sbjsidx,:,:);

clear NET;
for li = 1:size(tmp,1)
    NET(li,:) = squareform(squeeze(tmp(li,:,:)));
    li
end
%%
sbjsidx2 = find(ismember(sbjthk,sbjs));            % index in motion data
%
Thk=Thk(:,sbjsidx2)';


%% setup confounds matrix
conf=palm_inormal([ varsQconf vars(:,[3 4 7 14 15 22 23 25]) vars(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,2:end).^2]);  % add on squared terms and renormalise

%% number of components to feed into CCA
Nkeep=100;

%% prepare permutation scheme using PALM - for more details see:
%%% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks#Generating_the_set_of_permutations_only
Nperm=1000;                                                                       % in the paper we used 100000 but 10000 should be enough
EB=hcp2blocks_xl('/path/to/restricted/file.csv',[ ], false, vars(:,1)); % change the filename to your version of the restricted file
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
amNET=abs(mean(NET));  NET3=nets_demean(NET./repmat(amNET,size(NET,1),1));  NET3(:,amNET<0.1)=[];
NET3=NET3/std(NET3(:)); % norm by mean of columns, removing badly conditioned ones
grot=[NET1 NET3]; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu1,ss1,vv1]=nets_svds(double(NETd),Nkeep); % SVD reduction

%% prepare main netmat matrix - we have a range of normalisation possibilities
Thk1=nets_demean(Thk);  Thk1=Thk1/std(Thk1(:)); % no norm
amThk=abs(mean(Thk));  Thk3=nets_demean(Thk./repmat(amThk,size(Thk,1),1));  Thk3(:,amThk<0.1)=[];
Thk3=Thk3/std(Thk3(:)); % norm by mean of columns, removing badly conditioned ones
grot=[Thk1 Thk3]; % concat horizontally
Thkd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu3,ss3,vv3]=nets_svds(double(Thkd),Nkeep); % SVD reduction

%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values
badvars=[];
for i=1:size(vars,2)
  Y=vars(:,i); grotKEEP=~isnan(Y);  
  grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?
  if (sum(grotKEEP)>250) & (std(Y(grotKEEP))>0) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95) & (grot<100)
    i=i; % do nothing
  else
    [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
    badvars=[badvars i];
  end
end

%% get list of which SMs to feed into CCA
varskeep=setdiff([1:size(vars,2)],[1 6 267:457 ...                            % SMs we generally ignore (ID, race, FreeSurfer)
 2 7 14 15 22 23 25 265 266  ...                                              % confound SMs
 11 12 13 17 19 27 29 31 34 40 204 205 212:223 229:233 236 238 242 477 ...    % some more SMs to ignore for the CCA
 3 4 5 8 9 10 16 18 20 21 24 26 28 30 32 33 35:39 458 459 460 463 464 ...     % some more SMs to ignore for the CCA
 badvars]);                                                                   % the "bad" vars auto-detected above

%% "impute" missing vars data - actually this avoids any imputation
varkeep_noalc=varskeep([1:2 4:50 72:end]);

%varsd=palm_inormal(vars(:,varskeep)); % Gaussianise
varsd=palm_inormal(vars(:,varkeep_noalc)); % Gaussianise

for i=1:size(varsd,2) % deconfound ignoring missing data
  grot=(isnan(varsd(:,i))==0); 
  grotconf=nets_demean(conf(grot,:)); 
  varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsdCOV=zeros(size(varsd,1));
for i=1:size(varsd,1) % estimate "pairwise" covariance, ignoring missing data
  for j=1:size(varsd,1)
    grot=varsd([i j],:); 
    grot=cov(grot(:,sum(isnan(grot))==0)'); 
    varsdCOV(i,j)=grot(1,2);
  end
end
varsdCOV2=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[uu,dd]=eigs(varsdCOV2,Nkeep);  % SVD (eigs actually)
%%

uu2=uu3-conf*(pinv(conf)*uu3);    % deconfound again just to be safe

%% CCA
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,uu2);

%% CCA permutation testing
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,uu2(PAPset(:,j),:));
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%% netmat weights for CCA mode 1
grotAA = corr(grotU(:,1),NET)';
 % or
grotAAd = corr(grotU(:,1),NETd(:,1:size(NET,2)))'; % weights after deconfounding

%% Thickness weights for CCA mode 1
grotBB = corr(grotV(:,1),Thk)';
 % or
grotBBd = corr(grotV(:,1),Thkd(:,1:size(Thk,2)))'; % weights after deconfounding

%%
save('mats/test.mat','grotA','grotB','grotU','grotV','grotBB','grotBBd','grotAAd','grotAA','sbjs','var_key','grotR','grotRpval');
