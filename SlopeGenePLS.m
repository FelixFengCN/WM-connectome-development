clc
clear
load('AHBA_processed.mat');
slope=load(['slope.mat']);
load(['surrogates.mat']);
outindex=[30,31,39,40,51,117,134,147,165,169,208,211:246];
slope(outindex,:)=[];
dNeParcelExpression=parcelExpression;
X=dNeParcelExpression(:,2:end);
Y=S;
X(isnan(X))=0;
Y(isnan(Y))=0;
X=zscore(X);
Y=zscore(Y);
genesSymbol = probeInformation.GeneSymbol;
geneID = probeInformation.EntrezID;

dim=15;
[XL,YL,XS,YS,BETA,PCTVAR2,MSE,stats]=plsregress(X,Y,dim);
for i=1:dim
    [R(i),~]=corr(XS(:,i),Y);
end
for counter=1:size(surrogates,1)
    Yp=surrogates(counter,:)';
    [XL,YL,XS,YS,BETA,PCTVAR2,MSE,stats]=plsregress(X,Yp,dim);
    [R_s(counter,:),~]=corr(XS,Y);
end
for i=1:dim
    p(i)=sum(abs(R_s(:,i))>abs(R(i)))/size(surrogates,1);
end



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
parfor counter=1:bootnum
%     for counter=1:bootnum
    counter
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
OrderedGeneSymbol=genesSymbol(ind1);
OrderedGeneID=geneID(ind1);
Z1_p=ones(size(Z1))*0.5;
Z1_p(Z1>0)=1-normcdf(Z1(Z1>0));
Z1_p(Z1<0)=normcdf(Z1(Z1<0));













