function [ SugiCorr , SugiR , LM , SugiY , SugiX , origY , origX ] = SugiLM( X , Y , tau , E , LMN )
% Calculating Sugihara's CMM, L and M causality. 
%
% References:
% - for Sugihara's CCM method > Sugihara, George, et al., Detecting Causality in Complex Ecosystems, Science 26 October 2012, Vol. 338, no. 6106, pp. 496-500.   
% - for L method > Chicharro, Daniel, and Ralph G. Andrzejak, Reliable detection of directional couplings using rank statistics, Physical Review E 80.2 (2009): 026217.
% - for M method > Andrzejak, Ralph G., et al., Bivariate surrogate techniques: necessity, strengths, and caveats, Physical review E 68.6 (2003): 066202.	
% 
% Inputs:
% X,Y - time series with the same length
% tau - time step for the reconstruction
% E   - dimension of the reconstruction  
% LMN - number of neighborhoods for L and M methods 
%       the number of neighborhoods for Sugihara's CCM method is E+1  
%
% Outputs:
% SugiCorr - correlation between the CCM estimation of original data and original data
% SugiR    - sqrt((sum((originaldata-CCMestimaleddata).^2)/numel(originaldata)))/std(origaldata)
% LM - results for L and M methods
% SugiY, SugiX - the CCM estimate of original data 
% origY, origX - original data
switch nargin
    case 5
    case 4
        LMN = E+1;
    otherwise
        error('Bad input')
end
L=length(X);
T=1+(E-1)*tau;
Xm=zeros((L-T+1),E);
Ym=zeros((L-T+1),E);
SugiN=E+1;
N = L-T+1;
%% RECONTRUCTIONS OF ORIGINAL SYSTEMS
for t=1:(L-T+1)
    Xm(t,:)=X((T+t-1):-tau:(T+t-1-(E-1)*tau));
    Ym(t,:)=Y((T+t-1):-tau:(T+t-1-(E-1)*tau));
end
%%
LMj= zeros(2,2,N);
SugiX=zeros(N,1);
SugiY=zeros(N,1);
origY=Y(T:end);
origX=X(T:end);
parfor j=1:N
%% neighborhood search 
[n1,d1]=knnsearch(Xm,Xm(j,:),'k',N);
[n2,d2]=knnsearch(Ym,Ym(j,:),'k',N);
susY=origY(n1(2:end)); 
susX=origX(n2(2:end));
%% CMM
SugsusY=susY(1:SugiN);
SugsusX=susX(1:SugiN);
Sugid1=d1(:,2:SugiN+1);
Sugid2=d2(:,2:SugiN+1);
u1=exp(-Sugid1./(Sugid1(:,1)*ones(1,SugiN)));
u2=exp(-Sugid2./(Sugid2(:,1)*ones(1,SugiN)));
w1=u1./(sum(u1,2)*ones(1,SugiN));
w2=u2./(sum(u2,2)*ones(1,SugiN));
SugiY(j)= w1*SugsusY;
SugiX(j)= w2*SugsusX;
%% LM
LMn1=n1(n1~=j);
LMn2=n2(n2~=j);
LMd1=d1(n1~=j);
LMd2=d2(n2~=j);
   
susXY=arrayfun(@(x) find(LMn1(:) == x,1,'first'), LMn2(1:LMN) );
susYX=arrayfun(@(x) find(LMn2(:) == x,1,'first'), LMn1(1:LMN) );
    
sum1=sum(LMd1(:))/(N-1);
sum2=sum(LMd2(:))/(N-1);
        
LMj(:,:,j) = [(N/2-sum(susXY)/LMN)/(N/2-(LMN+1)/2) , (sum1-sum(LMd1(susXY))/LMN)/(sum1-sum(LMd1(1:LMN))/LMN) ; (N/2-sum(susYX)/LMN)/(N/2-(LMN+1)/2) , (sum2-sum(LMd2(susYX))/LMN)/(sum2-sum(LMd2(1:LMN))/LMN)]; 
    
   % (GN-G(Y|X))/(GN-Gk) => L(Y|X)
   % (GN-G(X|Y))/(GN-Gk) => L(X|Y)
  
   %(RNY-Rcond(Y|X))/(RNY-RkY) => M(Y|X)
   %(RNX-Rcond(X|Y))/(RNX-RkX) => M(X|Y)
end
SugiCorr1=corrcoef(origY,SugiY);
SugiCorr(2,1)=SugiCorr1(1,2);
SugiCorr2=corrcoef(origX,SugiX);
SugiCorr(1,1)=SugiCorr2(1,2);
SugiR(2,1)=sqrt((sum((origY-SugiY).^2)/numel(origY)))/std(origY);
SugiR(1,1)=sqrt((sum((origX-SugiX).^2)/numel(origX)))/std(origX);
LM = squeeze(mean(LMj,3));
end