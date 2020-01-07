% Run CCA between connectivity and behavior data
% focus on the 818 subjects with NETMAT data but 32K resolution
% Copyright by Feng Han & Xiao Liu.

%%
clear all;
load('/path/to/behavioral/data.mat');%%%%behavioral data
load('/path/to/connectivity/data.mat');%%%%connectivity data
load('/path/to/sub_label/data.mat');%%%% subj label
%
cd /the/path/to/HCP_PTN1200/
netpath = 'HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d200_ts2/';
sbjs1k = textread('HCP_PTN1200/subjectIDs.txt');
sbjs1k = cellstr(num2str(sbjs1k));
sbjsidx = find(ismember(sbjs1k, sbjthk));
sbjs = sbjs1k(sbjsidx);
sbjs(257) = []; % remove the subj 168240 who has a missing information 
sbjsidx = find(ismember(sbjs1k, sbjs)); % index in netmat data
sbjsidx2 = find(ismember(sbjthk, sbjs)); % index in thickness data

%%
%%% input data
var_key = importdata('/path/to/your/column_headers.txt');        % needs to be a subjects X subjectmeasures text file - see www.fmrib.ox.ac.uk/analysis/HCP-CCA
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
vars = zeros(length(sbjsidx),length(var_key));
vars(:,[1 3:6 8:end]) = behdat_d(sbjsidx3,varidx);
vars(:,2) = varsQconf;
vars(:,7) = rfmrimotion';

sex_all=double(strcmp(behdat_c(:,4),'M'));
vars(:,3)=sex_all(sbjsidx3);

vars(:,4) = behdat_d(sbjsidx3,551);

%
vars(:,sum(isnan(vars)==0)<60)=NaN;  % pre-delete vars with LOADS of missing data

%% load netmats from HCP PTN release
% set NET to be the subjects X unwrapped network-matrices (partial correlations)
netpath = 'HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d200_ts2/';

sbjs1k = textread('HCP_PTN1200/subjectIDs.txt');
sbjs1k = cellstr(num2str(sbjs1k));

tmp = dlmread([netpath,'netmats2.txt'],' ');
tmp = reshape(tmp,length(sbjs1k),200,200);
tmp = tmp(sbjsidx,:,:);

clear NET;
for li = 1:size(tmp,1)
    NET(li,:) = squareform(squeeze(tmp(li,:,:)));
    li
end

%% setup confounds matrix
conf=palm_inormal([ varsQconf vars(:,[3 4 7 14 15 22 23 25]) vars(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,3:end).^2]);  % add on squared terms and renormalise

%% number of components to feed into CCA

Nkeep=100;

%% prepare permutation scheme using PALM - for more details see:
%%% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks#Generating_the_set_of_permutations_only
Nperm=1000;       % permutation #
EB=hcp2blocks_xl('/path/to/restricted/file.csv',[ ], false, vars(:,1)); % change the filename to your version of the restricted file
%%%%            Defines whether dizygotic twins should be
%               treated as ordinary siblings (true), or be a category
%               on its own (false). Default = false.%%%%permu for the var(:,1)
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
amNET=abs(mean(NET));  NET3=nets_demean(NET./repmat(amNET,size(NET,1),1));  NET3(:,amNET<0.1)=[];
NET3=NET3/std(NET3(:)); % norm by mean of columns, removing badly conditioned ones %%%%A(:) reshapes all elements of A into a single column vector. This has no effect if A is already a column vector.

grot=[NET1 NET3]; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean%%%%pinv:Moore-Penrose pseudoinverse
[uu1,ss1,vv1]=nets_svds(double(NETd),Nkeep); % SVD reduction  %%%%SVD: RAM/time-efficient wrapper for eig/eigs

%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values
badvars=[];
for i=1:size(vars,2)%grotB
  Y=vars(:,i); grotKEEP=~isnan(Y);  %%%% each column was gotten,and stat # of nan in each column
  grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers? %%%% extract the non-nan ones in a column and de-median,maybe the way to detect outliers
  if (sum(grotKEEP)>250) & (std(Y(grotKEEP))>0) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95) & (grot<100)
     
    i=i; % do nothing
  else
    [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]%%%% the last bad column is as grot is too large; 
    badvars=[badvars i];
  end
end

%% get list of which SMs to feed into CCA
varskeep=setdiff([1:size(vars,2)],[1 6 267:457 ...                            % SMs we generally ignore (ID, race, FreeSurfer)
 2 7 14 15 22 23 25 265 266  ...                                              % confound SMs
 11 12 13 17 19 27 29 31 34 40 204 205 212:223 229:233 236 238 242 477 ...    % some more SMs to
3 4 5 8 9 10 16 18 20 21 24 26 28 30 32 33 35:39 458 459 460 463 464 ...      % some more SMs to ignore for the CCA
 badvars]);                                                                   % the "bad" vars auto-detected above
varkeep_noalc=varskeep([1:2 4:50 72:end]);

varsd=palm_inormal(vars(:,varkeep_noalc));% Gaussianise

%%
for i=1:size(varsd,2) % deconfound ignoring missing data
  grot=(isnan(varsd(:,i))==0); 
  grotconf=nets_demean(conf(grot,:)); 
  varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsdCOV=zeros(size(varsd,1));
for i=1:size(varsd,1) % estimate "pairwise" covariance, ignoring missing data%%%%*******
  for j=1:size(varsd,1)
    grot=varsd([i j],:); 
    grot=cov(grot(:,sum(isnan(grot))==0)'); %%%%sum(isnan(grot))==0) means the num of grot is not nan.
    varsdCOV(i,j)=grot(1,2);
  end
end
varsdCOV2=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix%%%% converted to the nearest Symmetric Positive Definite Matrix.
[uu,dd]=eigs(varsdCOV2,Nkeep);  % SVD (eigs actually)
uu2=uu-conf*(pinv(conf)*uu);    % deconfound again just to be safe

%% CCA
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,uu2);


%% CCA permutation testing
clear tmpU tmpV;
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,uu2(PAPset(:,j),:));
  tmpU(:,j) = grotUr(:,1);
  tmpV(:,j) = grotVr(:,1);
end
for i=1:Nkeep  % get FWE-corrected pvalues

  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end


grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%% netmat weights for CCA mode 1
grotAA = corr(grotU(:,1),NET)';
 % or
grotAAd = corr(grotU(:,1),NETd(:,1:size(NET,2)))'; % weights after deconfounding

%% SM weights for CCA mode 1

grotBB = corr(grotV(:,1),palm_inormal(vars),'rows','pairwise');

 % or 
varsgrot=palm_inormal(vars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grotBBd = corr(grotV(:,1),varsgrot,'rows','pairwise')'; % weights after deconfounding
%% making cifti for the connectivity strength modulation 







%% 
netpath = 'HCP_PTN1200/netmats/3T_HCP1200_MSMAll_d200_ts2/';

icapath = 'HCP_PTN1200/groupICA/groupICA_3T_HCP1200_MSMAll_d200.ica/';




%%
sbjs1k = textread('HCP_PTN1200/subjectIDs.txt');
sbjs1k = cellstr(num2str(sbjs1k));
sbjsidx = find(ismember(sbjs1k, sbjs)); % index in netmat data
tmp = dlmread([netpath,'netmats.txt'],' ');
tmp = reshape(tmp,length(sbjs1k),200,200);
tmp = tmp(sbjsidx,:,:);

clear NET;
for li = 1:size(tmp,1)
    NET(li,:) = squareform(squeeze(tmp(li,:,:)));
    li
end


%%
conn = mean(NET)';
connmtx = squareform(conn);

sgn = sign(mean(NET)');
sgnmtx = squareform(sgn);
% sgnmtx(sgnmtx==1) = 0;
grotAAmat = squareform(grotAAd.*sgn);
grotAAmat2 = squareform(grotAAd);
%
[sgAA, sind] = sort(grotAAmat);

% dcwgt = mean(sgAA(1:50,:));
% icwgt = mean(sgAA(151:end,:));
% 
% diwgt = abs(dcwgt)-icwgt;

dcwgt2 = sum(sgAA<-0.068);% z=-atanh(0.068)*(sqrt(818-3))
icwgt2 = sum(sgAA>0.068);
diwgt2 = icwgt2-dcwgt2;



%%

icsurf = ciftiopen([icapath, 'melodic_IC.dscalar.nii'],'/path/to/your/WorkBench/wb_command',1);
[~, labels] = max(icsurf.cdata,[],2);
icmask = zeros(size(icsurf.cdata));
for li = 1:size(icmask,2)
    icmask(:,li) = (labels==li);  
end
icsurf.cdata(:,1) = icmask*(-dcwgt2');
 
ciftisave(icsurf,'test.dscalar.nii','/path/to/your/WorkBench/wb_command');

icsurf = ciftiopen([icapath, 'melodic_IC.dscalar.nii'],'/path/to/your/WorkBench/wb_command',1);
[~, labels] = max(icsurf.cdata,[],2);
icmask = zeros(size(icsurf.cdata));
for li = 1:size(icmask,2)
    icmask(:,li) = (labels==li);  
end
icsurf.cdata(:,1) = icmask*(icwgt2');
ciftisave(icsurf,'test.dscalar.nii','/path/to/your/WorkBench/wb_command');

icsurf = ciftiopen([icapath, 'melodic_IC.dscalar.nii'],'/path/to/your/WorkBench/wb_command',1);
[~, labels] = max(icsurf.cdata,[],2);
icmask = zeros(size(icsurf.cdata));
for li = 1:size(icmask,2)
    icmask(:,li) = (labels==li);  
end
icsurf.cdata(:,1) = (icmask)*(diwgt2');
ciftisave(icsurf,'test.dscalar.nii','/path/to/your/WorkBench/wb_command');

