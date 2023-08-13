clear;clc;
atlas=246;
p_threshold=0.05;
slope=load('fnfa246_6-13slope.mat').slope;
slope(:,3)=slope(:,3)*atlas;
slope(:,12)=slope(:,12)*atlas;
for i=1:3
    sig_regions_model1{i}=find(slope(1+(i-1)*atlas:i*atlas,3)<p_threshold & ((slope(1+(i-1)*atlas:i*atlas,4)-slope(1+(i-1)*atlas:i*atlas,5))<0 | slope(1+(i-1)*atlas:i*atlas,12)>p_threshold));
    sig_regions_model2{i}=find(slope(1+(i-1)*atlas:i*atlas,12)<p_threshold & (slope(1+(i-1)*atlas:i*atlas,4)-slope(1+(i-1)*atlas:i*atlas,5))>0);
end
outindex=[30,31,39,40,51,117,134,147,165,169,208,211:246];
bootnum=10000;
perm=100;
load('100DS246scaledRobustSigmoidNSGRNAseqQC1wholeBrain_ROI_NOdistCorrEuclidean.mat');

%% spatial autocorrelation corrected permutation test to assess the significance of PLS component variance explained ratios
outindex=[30,31,39,40,51,117,134,147,165,169,208,211:246];
S=slope(1:atlas,1);
S(outindex,:)=[];
dNeParcelExpression=parcelExpression;
tempindex1=sum(isnan(dNeParcelExpression(:,2:end)),2)==10027;
dNeParcelExpression(tempindex1,:)=[];
index=1:246;
index(tempindex1)=[];
X=dNeParcelExpression(:,2:end);
Y=S;
X(isnan(X))=0;
Y(isnan(Y))=0;
X=zscore(X);
Y=zscore(Y);
genesSymbol = probeInformation.GeneSymbol;
geneID = probeInformation.EntrezID;
load(['ge_slope_surrogates.mat']);
dim=15;
for i=1:dim
    [R(i),~]=corr(XS(:,i),Y);
end

for counter=1:size(surrogates,1)
    Yp=surrogates(counter,:)';
    [XL,YL,XS,YS,BETA,PCTVAR2,MSE,stats]=plsregress(X,Yp,dim);
    [R_s(counter,:),~]=corr(XS,Y);
end
for i=1:dim
    p(i)=sum(abs(R_s(:,i))>R(i))/size(surrogates,1);
end

%% get gene weights and p values
% Do PLS in 1 dimensions (with 1 components):
dim = 1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
PLS1_Group_D = [XS(:,1),Y];
R=corr(XS(:,1),YS(:,1));
%align PLS components with desired direction for interpretability 
if R<0  
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end

%Order the genes accroding to the initial PLS weights
PLS1w=stats.W(:,1);
PLS1weights=[];

bootnum=10000;
for counter=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1); 
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 

    [XL2,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim);

    newWei=stats.W(:,1);
%     order the newly obtained weights the same way as initial PLS 

%     As the sign of PLS components is arbitrary, make sure this aligns between runs
    if corr(PLS1w,newWei)<0 
        newWei=-1*newWei;
    end
    PLS1weights=[PLS1weights,newWei];
end

%standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

%z-score weights
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
Z1_p=ones(size(Z1))*0.5;
Z1_p(Z1>0)=1-normcdf(Z1(Z1>0));
Z1_p(Z1<0)=normcdf(Z1(Z1<0));

