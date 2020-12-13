% Author: Rong Wang
% Date: 2020.Oct.29
tic
clear;

% Date Info
dd2020=313; % number of days in 2020 (to be updated: 313 for 8 Nov; 330 for 25 Nov)
nt=floor((dd2020-1)/7)+1;

% Prediction for the future starting from 25 Nov to 21 Jan
dd2020case=330; 

% Country Info
cninfo2=load('input\countryinfo1031.txt');

% Total GDP
gdpcn=mean(cninfo2(:,11:14),2);

% dNO2 + weekly cases
load('input_tvcdbyregion\tvcdcn.dat','-mat'); % tvcdcn=ones(338,nt,9); tvcdcn=tvcdcn.*(-999);
load('input_dailycasesbyregion\daily_covid19.dat','-mat'); % input1207xiong_unconfirmed cases=zeros(ddtotal,338,8); 1 cumul case; 2 cumul death; 3 case; 4 death; 5 log cases; 6 velocity; 7 p value; 8 number of cases>=1
load('input_tvcdbyregion\mete2020.dat','-mat'); % mete2020=zeros(338,nt,8); 3 pbl; 4 temp; 5 press; 6 preci; 7 rh; 8 uwind; 9 vwind; 10 wind

% Missing data
missdata=zeros(256*5+3,1);
vss=[1 2 4 5 8];
for cn=1:256
    for v=1:5
        a1=mete2020(cn,1:end,vss(v)); idxa=find(a1>-990);
        if size(idxa,2)>1
            missdata(cn+256*(v-1),1)=mean(a1(idxa),2);
        else
            missdata(cn+256*(v-1),1)=-999;
        end
    end
end
for v=1:5
    a2=missdata((v*256-255):(v*256),1); idxa=find(a2>-990); a2mean=mean(a2(idxa),1);
    for i=(v*256-255):(v*256)
        if missdata(i,1)<-990
            missdata(i,1)=a2mean;
        end
    end
end
a3=cninfo2(:,5); idxa=find(a3>-990 & a3<120); missdata(256*5+1,1)=mean(a3(idxa),1);
a3=cninfo2(:,19); idxa=find(a3>-990 & a3<12); missdata(256*5+2,1)=mean(a3(idxa),1);
a3=cninfo2(:,10); idxa=find(a3>-990 & a3<520); missdata(256*5+3,1)=mean(a3(idxa),1);

% Preparing data for ML...
xy1=ones(nt*2,338); xy1=xy1.*(-999);
xy2=zeros(338+7,4*3); xy2gdp=zeros(7,1);
xy3=ones(338,15); xy3=xy3.*(-999);
xy4=zeros(10000,13); e1=0;
xy5=ones(338,18,nt); xy5=xy5.*(-999);
for cn=1:338
    
    region=cninfo2(cn,1);
    if cninfo2(cn,1)==-999 || tvcdcn(cn,1,1)<-100 || gdpcn(cn,1)<0
        continue;
    end    
    
    % Moving average of tvcdcn 
    tv=zeros(nt,13);
    for v=1:13
        tv(1:nt,v)=tvcdcn(cn,1:nt,v); % 1-4 for 360x720; 1 OMI-RF; 2 OMI-PLS; 3 TRO-RF; 4 TRO-PLS; 5-8 for 90x180
        moving=zeros(nt-2,2); moving(1:(nt-2),1)=tvcdcn(cn,1:(nt-2),v); moving(1:(nt-2),2)=tvcdcn(cn,3:nt,v);
        tv(2:(nt-1),v)=(tv(2:(nt-1),v)+moving(:,1)+moving(:,2))/3;
    end    
    
    % Cases by week
    weekcase=zeros(50,7); 
    for i=1:nt
        weekcase(i,1)=mean(cases((i*7-6):min(i*7,dd2020case),cn,3),1); % N, day-1
        weekcase(i,2)=mean(cases((i*7-6):min(i*7,dd2020case),cn,6),1); % velocity = slope of logN
        weekcase(i,3)=mean(cases((i*7-6):min(i*7,dd2020case),cn,8),1); % number of data points (should be 15)
        weekcase(i,4)=sqrt(mean(cases((i*7-6):min(i*7,dd2020case),cn,7).*cases((i*7-6):min(i*7,dd2020case),cn,7),1)); % se of velocity
        % log rate
        i1=(i*7-13); i2=min(i*7+7,dd2020case); tt=([i1:i2])';
        if i>1 && weekcase(i,3)>=10
            tt2=cases(i1:i2,cn,6); idx=find(tt2>-990); [b4,bint4,r4,rint4,stats4]=regress(tt2(idx),[ones(size(idx,1),1) tt(idx,1)]);
            weekcase(i,5)=b4(2,1); % acceleration velocity, day-2
            weekcase(i,6)=(bint4(2,2)-bint4(2,1))/2/1.96; % se of acceleration velocity
            tt2=cases(i1:i2,cn,5); idx=find(tt2>-990); [b4,bint4,r4,rint4,stats4]=regress(tt2(idx),[ones(size(idx,1),1) tt(idx,1)]);
            weekcase(i,7)=b4(2,1); % velocity, day-1
            weekcase(i,8)=(bint4(2,2)-bint4(2,1))/2/1.96; % se of velocity
        end
    end
    
    % Weeks after the outbreak of covid-19
    idx=find(cases(:,cn,9)==1);
    if size(idx,1)==0
        continue;
    end
    w1=floor((idx(1)-1)/7)+1;
    w2=floor((idx(end)-1)/7)+1;    
    
    % Loading the continuous data of tvcd and covid-19
    continuity=0;
    qi=1;
    q=zeros(nt,2);
    q(1,1)=w1;
    xy=ones(nt,12);
    xy=xy.*(-999);
    
    % Loading data by week...
    for i=(w1+1):w2
        if weekcase(i,3)>=10 && weekcase(i-1,3)>=10
            % User-defined variable: 1 OMI-RF; 2 OMI-PLS; 3 TRO-RF; 4 TRO-PLS;
            v=1; 
            
            % Loading data ...
            xy(i,1)=tv(i,v); % tvcd
            xy(i,2)=tv(i,v+9); % SE of tvcd or min(tv(i,v+9),0.4)
            xy(i,3)=weekcase(i,7); % velocity = slope of logN, % day-1
            xy(i,4)=weekcase(i,8); % se of velocity
            xy(i,5)=weekcase(i,5); % acceleration velocity, % day-1
            xy(i,6)=weekcase(i,6); % se of acceleration velocity
            continuity=continuity+1;
        else
            qi=qi+1;
            q(qi,1)=i;
            q(qi,2)=continuity;
            continuity=0;
        end
    end
    
    % Loading data...
    qi=qi+1; q(qi,1)=i+1; q(qi,2)=continuity;
    if qi>1
        [B, IX] = sort(-q(1:qi,2), 1);
        for i=(w1+1):w2
            if i>q(IX(1)-1,1) && i<q(IX(1),1)
            else
                xy(i,1:end)=-999;
            end
        end
    end
   
    % User-defined variable: 8 for velocity; 9 for acceleration velocity
    vv=9; 
    
    % Determine the continuity variable (k) defined in our apper
    for i=(w1+1):w2
        if xy(i,1)>-990 && xy(i,5)>-990
            xy1(i,cn)=xy(i,1);
            xy1(i+nt,cn)=xy(i,3);
            % v = 1 tvcd, 2 velocity and 3 acceleration velocity
            for v=1:3
                if (xy(i,v*2-1)+xy(i,v*2))<0
                    xy(i,v+6)=4;
                    if (xy(i-1,v*2-1)+xy(i-1,v*2))<0
                        xy(i,v+6)=5;
                        if (xy(i-2,v*2-1)+xy(i-2,v*2))<0
                            xy(i,v+6)=6;
                        end
                    end
                else
                    xy(i,v+6)=3;
                    if (xy(i-1,v*2-1)+xy(i-1,v*2))>=0 && xy(i-1,1)>-990
                        xy(i,v+6)=2;
                        if (xy(i-2,v*2-1)+xy(i-2,v*2))>=0 && xy(i-2,1)>-990
                            xy(i,v+6)=1;
                        end
                    end
                end
            end
            % coupling or decoupling
            for level=1:3
                if xy(i,7)<1 || xy(i,vv)<1
                    continue;
                end
                if xy(i,7)>=(level+3)
                    if xy(i,vv)>=(level+3)
                        xy(i,level+9)=1; % low tvcd & low risk
                    else
                        xy(i,level+9)=4; % low tvcd & high risk
                    end
                else
                    if xy(i,vv)>=(level+3)
                        xy(i,level+9)=2; % high tvcd & low risk
                    else
                        xy(i,level+9)=3; % high tvcd & high risk
                    end
                end
            end
        end
    end
    
    % Loading socioeconomic and environmental indicators
    qs=zeros(100,14);
    q1=0;
    q2=0;
    if cn>=257 && cn<=287
        cn2=44;
    elseif cn>=288
        cn2=235;
    else
        cn2=cn;
    end
    for i=(w1+1):w2
        if xy(i,7)>=1 && xy(i,7)<=6
            e1=e1+1;
            xy4(e1,1)=cn;
            xy4(e1,2)=xy(i,vv*2-13); % velocity or acceleration velocity
            xy4(e1,3)=xy(i,7); % tvcd level 1-6
            xy4(e1,4)=xy(i,1); % tvcd
            if cninfo2(cn2,5)>-990
                xy4(e1,5)=min(120,cninfo2(cn2,5)); % time travel
            else
                xy4(e1,5)=missdata(256*5+1,1);
            end
            if cninfo2(cn,19)>-990
                xy4(e1,6)=min(12,cninfo2(cn,19)); % gdp per cap
            else
                xy4(e1,6)=missdata(256*5+2,1);
            end
            if cninfo2(cn2,10)>-990
                xy4(e1,7)=min(520,cninfo2(cn2,10)); % pop density
            else
                xy4(e1,7)=missdata(256*5+3,1);
            end
            for v=1:5
                if mete2020(cn,i,vss(v))>-990
                    xy4(e1,v+7)=mete2020(cn,i,vss(v)); % pbl,temp,precip,rh,wind
                else
                    xy4(e1,v+7)=missdata(cn2+(v-1)*256,1);
                end
            end
        end
        if xy(i,7)==1 || xy(i,7)==6
            if xy(i,9)==1 || xy(i,9)==6
                if xy(i,7)==1
                    q1=q1+1; % no control for 3 consecutive weeks
                    if xy(i-1,1)>-990 && xy(i-2,1)>-990
                        qs(q1,1)=(xy(i,1)+xy(i-1,1)+xy(i-2,1))/3;
                    else
                        qs(q1,1)=xy(i,1);
                    end
                    if xy(i-1,vv*2-13)>-990 && xy(i-2,vv*2-13)>-990
                        qs(q1,2)=(xy(i,vv*2-13)+xy(i-1,vv*2-13)+xy(i-2,vv*2-13))/3;
                    else
                        qs(q1,2)=xy(i,vv*2-13);
                    end
                    for v=1:5
                        qs(q1,v+2)=mete2020(cn,i,vss(v));
                    end
                else
                    q2=q2+1; % control for 3 consecutive weeks
                    if xy(i-1,1)>-990 && xy(i-2,1)>-990
                        qs(q2,8)=(xy(i,1)+xy(i-1,1)+xy(i-2,1))/3;
                    else
                        qs(q2,8)=xy(i,1);
                    end
                    if xy(i-1,vv*2-13)>-990 && xy(i-2,vv*2-13)>-990
                        qs(q2,9)=(xy(i,vv*2-13)+xy(i-1,vv*2-13)+xy(i-2,vv*2-13))/3;
                    else
                        qs(q2,9)=xy(i,vv*2-13);
                    end
                    for v=1:5
                        qs(q2,v+9)=mete2020(cn,i,vss(v));
                    end
                    if xy3(cn,15)<0
                        xy3(cn,15)=min(15,i-w1); % first week with control
                    end
                end
            end
        end
        for level=1:3
            if xy(i,level+9)>0
                xy2(cn+7,xy(i,level+9)+(level-1)*4)=xy2(cn+7,xy(i,level+9)+(level-1)*4)+1;
            end
        end
    end
    
    for i=(w1+1):nt
            xy5(cn,1,i)=cn;
            xy5(cn,2,i)=weekcase(i,5); % velocity or acceleration velocity
            v1=1;
            tvcdup=tv((i-2):i,v1)+tv((i-2):i,v1+9); % 1 OMI-ML; 2 OMI-PLS; 3 TRO-ML; 4 TRO-PLS;
            if tvcdup(3)<0
                xy5(cn,3,i)=4;
                if tvcdup(2)<0
                    xy5(cn,3,i)=5;
                    if tvcdup(1)<0
                        xy5(cn,3,i)=6;
                    end
                end
            else
                xy5(cn,3,i)=3;
                if tvcdup(2)>=0
                    xy5(cn,3,i)=2;
                    if tvcdup(1)>=0
                        xy5(cn,3,i)=1;
                    end
                end
            end
            xy5(cn,4,i)=tv(i,v1); % 1 OMI-ML; 2 OMI-PLS; 3 TRO-ML; 4 TRO-PLS;
            xy5(cn,16,i)=tv(i,v1+9); % SE of tvcd
            xy5(cn,17,i)=tv(i-1,v1+9); % i-1
            xy5(cn,18,i)=tv(i-2,v1+9); % i-2
            if cninfo2(cn2,5)>-990
                xy5(cn,5,i)=min(120,cninfo2(cn2,5)); % time travel
            else
                xy5(cn,5,i)=missdata(256*5+1,1);
            end
            if cninfo2(cn,19)>-990
                xy5(cn,6,i)=min(12,cninfo2(cn,19)); % gdp per cap
            else
                xy5(cn,6,i)=missdata(256*5+2,1);
            end
            if cninfo2(cn2,10)>-990
                xy5(cn,7,i)=min(520,cninfo2(cn2,10)); % pop density
            else
                xy5(cn,7,i)=missdata(256*5+3,1);
            end
            for v=1:5
                if mete2020(cn,i,vss(v))>-990
                    xy5(cn,v+7,i)=mete2020(cn,i,vss(v)); % pbl,temp,precip,rh,wind
                else
                    xy5(cn,v+7,i)=missdata(cn2+(v-1)*256,1);
                end
            end
            if xy3(cn,15)>0
                xy5(cn,13,i)=xy3(cn,15); % first week with control
            else
                xy5(cn,13,i)=15;
            end
            xy5(cn,14,i)=floor(weekcase(i,1)); % number of cases day-1
            xy5(cn,15,i)=weekcase(i,7); % velocity = slope of logN, % day-1
    end
    
    % no control for 3 consecutive weeks
    if q1>3
        for vq=1:7
            xy3(cn,vq)=mean(qs(1:q1,vq),1); % tvcd
        end
    end
    
    % control for 3 consecutive weeks
    if q2>0
        for vq=8:14
            xy3(cn,vq)=mean(qs(1:q2,vq),1); % tvcd
        end
    end
    
    % TVCD weighted by GDP
    if cn~=44 && cn~=235
        for level=1:3
            typesum=xy2(cn+7,level*4-3)+xy2(cn+7,level*4-2)+xy2(cn+7,level*4-1)+xy2(cn+7,level*4);
            for type=1:4
                if typesum>0
                    xy2(region+1,type+(level-1)*4)=xy2(region+1,type+(level-1)*4)+gdpcn(cn,1)*xy2(cn+7,type+(level-1)*4)/typesum;
                    xy2(1,type+(level-1)*4)=xy2(1,type+(level-1)*4)+gdpcn(cn,1)*xy2(cn+7,type+(level-1)*4)/typesum;
                    % weighted average
                    if level==1 && type==1
                        xy2gdp(region+1,1)=xy2gdp(region+1,1)+gdpcn(cn,1);
                        xy2gdp(1,1)=xy2gdp(1,1)+gdpcn(cn,1);
                    end
                end
            end
        end
    end
end

% Average TVCD weighted by GDP
for i=1:7
    xy2(i,1:12)=xy2(i,1:12)./xy2gdp(i);
end

% Number of weeks with d_TVCD>0 after the outbreak
xy4=xy4(1:e1,1:end);
for i=1:e1
    if xy3(xy4(i,1),15)>-990
        xy4(i,13)=xy3(xy4(i,1),15);
    else
        xy4(i,13)=max(xy3(:,15),[],1);
    end
end

% MLR regression
[b4,bint4,r4,rint4,stats4]=regress(xy4(:,2),[ones(e1,1) xy4(:,3:13)]);
[R P]=corrcoef(xy4(:,2),xy4(:,4));
xy4pre=xy4(:,2)-r4;

% Finding the optimal model parameters in gradient-boosting-decision-tree regression
learnRate = [0.1 0.25 0.5 1];
numLR = numel(learnRate);
mnsPlot = [1 8 64];
numTrees = 150;
Mdl = cell(3,numLR);

for k = 1:numLR
    for j = 1:3
        t = templateTree('MaxNumSplits',mnsPlot(j),'Surrogate','on');
        Mdl{j,k} = fitrensemble(xy4(:,3:13),xy4(:,2),'NumLearningCycles',numTrees, ...
            'Learners',t,'KFold',4,'LearnRate',learnRate(k));
    end
end

% Estiamte the mean squared error in the four-fold cross validation
kflAll = @(x)kfoldLoss(x,'Mode','cumulative');
errorCell = cellfun(kflAll,Mdl,'Uniform',false);
error = reshape(cell2mat(errorCell),[numTrees 3 numel(learnRate)]);
for k = 1:3
    if k<=2
        subplot(2,4,k);
    else
        subplot(2,4,5);
    end
    plot(squeeze(error(:,k,:)),'LineWidth',2)
    axis tight
    xlabel('Number of trees')
    ylabel('Cross-validated MSE')
    title(sprintf('MaxNumSplits = %0.3g', mnsPlot(k)))
end

% Using the optimal model parameters in gradient-boosting-decision-tree regression
idxNumTrees=7;
idxlearnRate=0.25;
idxk=64;
tFinal = templateTree('MaxNumSplits',idxk,'Surrogate','on');
Mdl=fitrensemble(xy4(:,3:13),xy4(:,2),'NumLearningCycles',idxNumTrees,'Learners',tFinal,'LearnRate',idxlearnRate);

% Machine Learning without cross-validation
xy4preml=predict(Mdl,xy4(:,3:13));

% hold-one-out cross validation by SLR for simple linear regression, MLR for multiple linear regression and ML for machine learning
xy4preml_holdone=zeros(e1,3);
for i=1:e1
    % SLR
    b4=regress(xy4one(1:(e1-1),2),[ones(e1-1,1) xy4one(1:(e1-1),4)]);
    xy4preml_holdone(i,2)=b4(1)+b4(2)*xy4one(e1,4);
    
    % MLR
    b5=regress(xy4one(1:(e1-1),2),[ones(e1-1,1) xy4one(1:(e1-1),3:13)]);
    xy4preml_holdone(i,3)=b5(1);
    for k=1:11
        xy4preml_holdone(i,3)=xy4preml_holdone(i,3)+b5(1+k)*xy4one(e1,2+k);
    end
    
    % ML
    xy4one=zeros(e1,13+4);
    xy4one(1:end,1:13)=xy4;
    xy4one(i:(e1-1),1:13)=xy4((i+1):e1,1:13);
    xy4one(e1,1:13)=xy4(i,1:13);
    MdlFinal = fitrensemble(xy4one(1:(e1-1),3:13),xy4one(1:(e1-1),2),'NumLearningCycles',idxNumTrees,'Learners',tFinal,'LearnRate',idxlearnRate);
    Mdl=fitrensemble(xy4one(1:(e1-1),3:13),xy4one(1:(e1-1),2),'Method','LSBoost','NumLearningCycles',10); % 2019-2019
    xy4preml_holdone(i,1)=predict(MdlFinal,xy4one(e1,3:13)); % 2019-2020
end

xy5pre=ones(338,11*2+1); xy5pre=xy5pre.*(-999);
xy5prediscrete=ones(338,11*2+1); xy5prediscrete=xy5prediscrete.*(-999);
xy5tvcd=zeros(338,2);
for cn=1:338
    % Region ID
    region2=cninfo2(cn,1);
    
    % Excluding the national data of China and US, where provincial and
    % state data are used for ML
    if cninfo2(cn,1)==-999 || tvcdcn(cn,1,1)<-100 || gdpcn(cn,1)<0 || cn==44 || cn==235
        continue;
    end
    
    % Days with covid-19
    idx=find(cases(:,cn,9)==1);
    if size(idx,1)==0
        continue;
    end
    w1=floor((idx(1)-1)/7)+1;
    if cn==76
        w1=25; % low peak in the first turn of pandemic
    elseif cn==177
        w1=28; % suspicious data before the 26th week
    elseif cn==8 || cn==205 || cn==232
        w1=w1+1; % flutuation in the first 10 days
    end
    
    % Simulation starting at the end of the first 3 weeks
    istart=max(3,w1)+3;
    if (istart+8)>nt
        continue;
    end
    for ww=1:2
        if ww==2
            istart=nt;
        end
        a=xy5(cn,:,istart); % a=xy5(3,1:13,nt); a(1,4)=-0.3; a(1,3)=4; a(1,13)=1; predict(Mdl,a(1,3:13));
        if min(a,[],2)<-990
            continue;
        end
        if ww==1
            a(1,13)=1;
            for wk=1:8
                xy5tvcd(cn,1)=xy5tvcd(cn,1)+min(0,xy5(cn,4,(istart+wk)))/8; % exclude positive tvcd
            end
        else
            xy5tvcd(cn,2)=min(0,xy5(cn,4,nt)); % exclude positive tvcd
        end
        a(1,3)=1; a(1,4)=0;
        A0=predict(Mdl,a(1,3:13))*100;
        xy5pre(cn,1+(ww-1)*11)=A0;
        
        % Consider that the control is discrete (level=1)
        xy5prediscrete(cn,1+(ww-1)*11)=A0;
        for k=1:10
            totalc=-k/10*3;
            
            % Scenario 1 - control in the last week
            A1=A0;
            A1discrete=A0;
            for k2=1:min(3*k,10)
                a(1,4)=-k2/10;
                if (a(1,4)+a(1,16))<0
                    a(1,3)=4;
                else
                    a(1,3)=1;
                end
                Aa=predict(Mdl,a(1,3:13))*100;
                if Aa<A1
                    A1=Aa; % find the minimum
                end
                a2=a;a2(1,3)=1; Aadiscrete=predict(Mdl,a2(1,3:13))*100;
                if Aadiscrete<A1discrete
                    A1discrete=Aadiscrete; % find the minimum
                end
            end
            
            % Scenario 2 - control in the last 2 weeks
            A2=A0;
            A2discrete=A0;
            if (totalc+a(1,16)+a(1,17))<0
                a(1,3)=5;
                kn1=floor(a(1,16)*10)+1;
                kn2=floor((-totalc-a(1,17))*10);
                for k2=kn1:kn2
                    a(1,4)=-k2/10;
                    Aa=predict(Mdl,a(1,3:13))*100;
                    if Aa<A2
                        A2=Aa; % find the minimum
                    end
                    a2=a;a2(1,3)=1; Aadiscrete=predict(Mdl,a2(1,3:13))*100;
                    if Aadiscrete<A2discrete
                        A2discrete=Aadiscrete; % find the minimum
                    end
                end
            end
            
            % Scenario 3 - control in the last 3 weeks
            A3=A0;
            A3discrete=A0;
            if (totalc+a(1,16)+a(1,17)+a(1,18))<0
                a(1,3)=6;
                kn1=floor(a(1,16)*10)+1;
                kn2=floor((-totalc-a(1,17)-a(1,18))*10);
                for k2=kn1:kn2
                    a(1,4)=-k2/10;
                    Aa=predict(Mdl,a(1,3:13))*100;
                    if Aa<A3
                        A3=Aa; % find the minimum
                    end
                    a2=a;a2(1,3)=1; Aadiscrete=predict(Mdl,a2(1,3:13))*100;
                    if Aadiscrete<A3discrete
                        A3discrete=Aadiscrete; % find the minimum
                    end
                end
            end
            xy5pre(cn,1+k+(ww-1)*11)=min([A1,A2,A3],[],2);
            
            % We consider that the control is discrete
            xy5prediscrete(cn,1+k+(ww-1)*11)=min([A1discrete,A2discrete,A3discrete],[],2);
        end
    end
    xy5pre(cn,23)=1;
end

subplot(2,4,3);
for cn=1:338
   if xy5pre(cn,1)>-990
       plot(xy5pre(cn,12:22)); hold on;
   end
end

subplot(2,4,4);
for cn=1:338
   if xy5pre(cn,1)>-990
       plot(xy5pre(cn,12:22)); hold on;
   end
end

% Sensitivity of cases to tvcd
xy9case0=zeros(338,25);
xy90=zeros(338,24,2);
xy9case=zeros(338,25);
xy9=zeros(338,24,2);
for cn=1:338
    % Excluding the national data of China and US, where provincial and
    % state data are used for ML
    if xy5pre(cn,23)~=1 || gdpcn(cn,1)<=0 || cn==44 || cn==235
        continue;
    end
    
    % Days with covid-19
    idx=find(cases(:,cn,9)==1);
    if size(idx,1)==0
        continue;
    end
    w1=floor((idx(1)-1)/7)+1;
    if cn==76
        w1=25; % low peak in the first turn of pandemic
    elseif cn==177
        w1=28; % suspicious data before the 26th week
    elseif cn==8 || cn==205 || cn==232
        w1=w1+1; % flutuation in the first 10 days
    end
    
    % Simulation starting at the end of the first 3 weeks
    istart=max(3,w1)+3;
    if (istart+8)>nt
        continue;
    end
    astart=xy5(cn,:,istart);
    if min(astart,[],2)<-990
        continue;
    end
    
    % Starting at the end of first 3 weeks
    wss=1; % 1 starting at the end of first 3 weeks; 2 starting from Nov 9
    for sce=0:24
        d0=1;
        s1=floor((sce-1)/3);
        s2=sce-s1*3;
        simu=zeros(57,3); % 1 Nt; 2 Vt; 3 At
        if cases(istart*7,cn,5)>-0.7
            simu(d0,1)=cases(istart*7,cn,5); % lnN
        else
            simu(d0,1)=log(0.5);
        end
        simu(d0,2)=min(0.1,cases(istart*7,cn,6)); % velocity
        simu(d0,3)=astart(1,2); % A
        for w=1:8
            if sce==0
                Aa=xy5pre(cn,1+(wss-1)*11)/100; % no control
            elseif sce<=24
                if w<(s1+1)
                    Aa=xy5pre(cn,4+(wss-1)*11)/100; % 30% control in early week
                elseif w==(s1+1)
                    Aa=xy5pre(cn,1+s2+(wss-1)*11)/100; % s2*10% in this week
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100; % no control
                end
            end
            for d=1:7
                d0=d0+1;
                if simu(d0-1,1)>=0
                    simu(d0,1)=simu(d0-1,1)+simu(d0-1,2); % N
                else
                    simu(d0,1)=log(0.5);
                end
                simu(d0,2)=simu(d0-1,2)+simu(d0-1,3); % velocity
                simu(d0,3)=Aa; % A velocity
            end
        end
        xy9case0(cn,sce+1)=floor(exp(simu(end,1)));
        if sce>0
            xy90(cn,sce,1)=(xy9case0(cn,sce)-xy9case0(cn,sce+1))/1000; % reduced cases, thousand
            xy90(cn,sce,2)=gdpcn(cn,1)/365*7*0.1/1e9; % gdp * tvcd, billion USD
        end
    end
    
    % starting from Nov 26
    wss=2; % 1 starting at the end of first 4 weeks; 2 starting from Nov 9
    for sce=0:24
        s1=floor((sce-1)/3);
        s2=sce-s1*3;
        d0=1;
        simu=zeros(57,3); % 1 Nt; 2 Vt; 3 At
        if cases(dd2020case,cn,5)>-0.7
            simu(d0,1)=cases(dd2020case,cn,5); % lnN
        else
            simu(d0,1)=log(0.5);
        end
        simu(d0,2)=max(-0.05,min(0.05,cases(dd2020case,cn,6))); % velocity
        simu(d0,3)=0; % A
        if cn==265
            simu(d0,1)=log(2); % http://www.bjfsh.gov.cn/ztzx/2020/yqfk/zxyq/202011/t20201123_40009834_fs.shtml
            simu(d0,2)=(log(2)-log(0.5))/10; % V 0.05 or log(1)-log(0.5)
            simu(d0,3)=0; % A
        end
        for w=1:8
            if sce==0
                Aa=xy5pre(cn,1+(wss-1)*11)/100; % no control
            elseif sce<=24
                if w<(s1+1)
                    Aa=xy5pre(cn,4+(wss-1)*11)/100; % 30% control in early week
                elseif w==(s1+1)
                    Aa=xy5pre(cn,1+s2+(wss-1)*11)/100; % s2*10% in this week
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100; % no control
                end
            end
            for d=1:7
                d0=d0+1;
                if simu(d0-1,1)>=0
                    simu(d0,1)=simu(d0-1,1)+simu(d0-1,2); % N
                else
                    simu(d0,1)=log(0.5);
                end
                simu(d0,2)=simu(d0-1,2)+simu(d0-1,3); % velocity
                simu(d0,3)=Aa; % A velocity
            end
        end
        xy9case(cn,sce+1)=floor(exp(simu(end,1)));
        if sce>0
            xy9(cn,sce,1)=(xy9case(cn,sce)-xy9case(cn,sce+1))/1000; % reduced cases, thousand
            xy9(cn,sce,2)=gdpcn(cn,1)/365*7*0.1/1e9; % gdp * tvcd, billion USD
        end
    end
end

% Finding the 20% or 10% regions where confinement is the most effective
xy9casereduced=zeros(338,2);
xy9casereduced(:,1)=sum(xy90(:,1:24,1),2);
xy9casereduced(:,2)=sum(xy9(:,1:24,1),2);
m9=zeros(6,4);
for region=1:6
    for ws=1:2
        a9=zeros(200,1); i9=0;
        for cn=1:338
            if cninfo2(cn,1)==region && xy9casereduced(cn,ws)>0
                i9=i9+1;
                a9(i9,1)=xy9casereduced(cn,ws);
            end
        end
        m9(region,(ws*2-1):(ws*2)) = prctile(a9(1:i9,1),[80 90]);
    end
end

xy6case=zeros(57,4*7);
xy7case=zeros(57,9*7);
xy8case=zeros(57,8*7);
for cn=1:338
    region2=cninfo2(cn,1);
    % Excluding the national data of China and US, where provincial and
    % state data are used for ML
    if cninfo2(cn,1)==-999 || tvcdcn(cn,1,1)<-100 || gdpcn(cn,1)<0 || xy5pre(cn,1)<-990 || xy5pre(cn,12)<-990 || cn==44 || cn==235
        continue;
    end
    
    % Days with covid-19
    idx=find(cases(:,cn,9)==1);
    if size(idx,1)==0
        continue;
    end
    w1=floor((idx(1)-1)/7)+1;
    
%     We consider the first day in week w0 as the outbreak of COVID-19 in this region. 
%     It is noticed that the first peak of Nw (4263 day-1 on 1–-7 Apr.) in France 
%     reported by ECDPC is far lower than the second one (48720 day-1 on 4–-11 Nov.), 
%     and we derive w1000 from the second peak. Because there are two kinds of Covid-19 
%     tests by June in Philippines, we derive w0 using daily cases reported after July 
%     to reduce the impact of different tests 
%     (https://www.philstar.com/headlines/2020/05/21/2015542/there-are-two-kinds-covid-19-tests-used-philippines-how-are-they-different).
%     In addition, we notice a remarkable fluctuation in the reported daily cases in 
%     United Arab Emirates (0, 16, 0, 14, 15, 0 and 11 cases on 7–13 Mar., respectively), 
%     Serbia (13, 6, 17, 5, 9, 2 and 15 cases on 12–18 Mar., respectively) 
%     and Ukraine (2, 9, 5, 7, 0, 15 and 6 cases on 17–23 Mar., respectively), 
%     where the data in the first week is not used.
    if cn==76
        w1=25; % low peak in the first turn of pandemic
    elseif cn==177
        w1=28; % suspicious data before the 26th week
    elseif cn==8 || cn==205 || cn==232
        w1=w1+1; % flutuation in the first 10 days
    end
    
    % Simulation starting at the end of the first 3 weeks
    istart=max(3,w1)+3;
    if (istart+8)>nt
        continue;
    end
    astart=xy5(cn,:,istart);
    if min(astart,[],2)<-990
        continue;
    end
    
    % Starting at the end of first 3 weeks
    wss=1; % 1 starting at the end of first 3 weeks; 2 starting from Nov 9
    for sce=1:10
        d0=1;
        simu=zeros(57,6); % 1 Nt; 2 Vt; 3 At
        if cases(istart*7,cn,5)>-0.7
            simu(d0,1)=cases(istart*7,cn,5); % lnN: log(xy5(cn,14,istart)) or cases(istart*7,cn,5)
        else
            simu(d0,1)=log(0.5);
        end
        simu(d0,2)=min(0.1,cases(istart*7,cn,6)); % velocity: xy5(cn,15,istart) or cases(istart*7,cn,6)
        simu(d0,3)=astart(1,2); % A
        simu(d0,4)=simu(d0,1); % real N
        simu(d0,5)=cases(istart*7,cn,6); % real V
        simu(d0,6)=astart(1,2); % real A
        for w=1:8
            a=xy5(cn,:,istart+w);
            if min(a,[],2)<-990
                a=astart;
            end
            alast=xy5(cn,:,istart+w-1);
            if min(alast,[],2)<-990
                alast=a;
            end
            if (istart+w+1)<=nt
                anext=xy5(cn,:,istart+w+1);
                if min(anext,[],2)<-990
                    anext=a;
                end
            else
                anext=a;
            end
            
            if w==1
                a(1,3)=4; alast(1,3)=4; anext(1,3)=4;
            end
            if sce==1
                % Scenario 1: lower bounds
                Aa=predict(Mdl,a(1,3:13));
                Aanext=predict(Mdl,anext(1,3:13));
                Aalast=predict(Mdl,alast(1,3:13));
                Aa=min([Aa,Aanext,Aalast],[],2);
            elseif sce==2
                % Scenario 2: higher bounds
                Aa=predict(Mdl,a(1,3:13));
                Aanext=predict(Mdl,anext(1,3:13));
                Aalast=predict(Mdl,alast(1,3:13));
                Aa=max([Aa,Aanext,Aalast],[],2);
            elseif sce==3
                % Scenario 3: real tvcd
                Aa=predict(Mdl,a(1,3:13)); % predict(Mdl,a(1,3:13)); Aa=xy5pre(cn,1+min(10,floor(-min(0,a(1,4))*10)+1)+(wss-1)*11)/100;
            elseif sce==4
                % Scenario 4: 0% tvcd
                Aa=xy5pre(cn,1+(wss-1)*11)/100;
            elseif sce==5
                % Scenario 5: 10% tvcd
                Aa=xy5pre(cn,2+(wss-1)*11)/100;
            elseif sce==6
                % Scenario 6: 30% tvcd
                Aa=xy5pre(cn,4+(wss-1)*11)/100;
            elseif sce==7
                % Scenario 7: 50% tvcd
                Aa=xy5pre(cn,6+(wss-1)*11)/100;
            elseif sce==8
                % Scenario 8: 50% tvcd in a discrete scenario
                Aa=xy5prediscrete(cn,6+(wss-1)*11)/100;
            elseif sce==9
                % Scenario 9: tvcd toward top 10% sensitive countries
                if xy9casereduced(cn,wss)>=m9(region,wss*2)
                    Aa=xy5pre(cn,6+(wss-1)*11)/100;
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100;
                end
            elseif sce==10
                % Scenario 10: tvcd toward top 20% sensitive countries
                if xy9casereduced(cn,wss)>=m9(region,wss*2-1)
                    Aa=xy5pre(cn,6+(wss-1)*11)/100;
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100;
                end
            end
            for d=1:7
                d0=d0+1;
                if simu(d0-1,1)>=0
                    simu(d0,1)=simu(d0-1,1)+simu(d0-1,2); % N
                else
                    simu(d0,1)=log(0.5);
                end
                simu(d0,2)=simu(d0-1,2)+simu(d0-1,3); % velocity
                simu(d0,3)=Aa; % A velocity
                simu(d0,4)=max(-0.6931,cases(min(dd2020,istart*7+d0),cn,5)); % real N
                simu(d0,5)=cases(min(dd2020,istart*7+d0),cn,6); % real V
                simu(d0,6)=a(1,2); % real A
            end
        end
        if sce==1
            xy6case(:,4)=xy6case(:,4)+exp(simu(:,4));
            xy6case(:,4+region*4)=xy6case(:,4+region*4)+exp(simu(:,4));
            xy7case(:,9)=xy7case(:,9)+exp(simu(:,4)); % obs
            xy7case(:,9+region*9)=xy7case(:,9+region*9)+exp(simu(:,4)); % obs
        end
        if sce<=3
            xy6case(:,sce)=xy6case(:,sce)+exp(simu(:,1));
            xy6case(:,sce+region*4)=xy6case(:,sce+region*4)+exp(simu(:,1));
        end
        if sce>=3
            xy7case(:,sce-2)=xy7case(:,sce-2)+exp(simu(:,1));
            xy7case(:,sce-2+region*9)=xy7case(:,sce-2+region*9)+exp(simu(:,1));
        end
    end
    
    % Future Simulation Starting on Nov 26
    wss=2; % 1 starting at the end of first 4 weeks; 2 starting from Nov 9
    for sce=1:8
        d0=1;
        simu=zeros(57,3); % 1 Nt; 2 Vt; 3 At
        if cases(dd2020case,cn,5)>-0.7
            simu(d0,1)=cases(dd2020case,cn,5); % lnN: log(xy5(cn,14,istart)) or cases(istart*7,cn,5)
        else
            simu(d0,1)=log(0.5);
        end
        simu(d0,2)=max(-0.05,min(0.05,cases(dd2020case,cn,6))); % velocity: xy5(cn,15,istart) or cases(istart*7,cn,6) max(-0.05,min(0.01,cases(dd2020case,cn,6)))
        simu(d0,3)=0; % A
        if cn==265
            simu(d0,1)=log(2); % http://www.bjfsh.gov.cn/ztzx/2020/yqfk/zxyq/202011/t20201123_40009834_fs.shtml
            simu(d0,2)=(log(2)-log(0.5))/10; % V 0.05 or log(1)-log(0.5)
            simu(d0,3)=0; % A
        end
        for w=1:8
            a=xy5(cn,:,nt);
            if sce==1
                % Scenario 1: 0% tvcd
                Aa=xy5pre(cn,1+(wss-1)*11)/100; % Aa=predict(Mdl,a(1,3:13)); % [Aa,AaSTE] = Predict(Mdl,a(1,3:13),'Quantile',[0.025,0.975]);
            elseif sce==2
                % Scenario 2: 10% tvcd
                Aa=xy5pre(cn,2+(wss-1)*11)/100;
            elseif sce==3
                % Scenario 3: 20% tvcd
                Aa=xy5pre(cn,3+(wss-1)*11)/100;
            elseif sce==4
                % Scenario 4: 30% tvcd
                Aa=xy5pre(cn,4+(wss-1)*11)/100;
            elseif sce==5
                % Scenario 5: 30% tvcd in a discrete scenario
                Aa=xy5prediscrete(cn,4+(wss-1)*11)/100;
            elseif sce==6
                % Scenario 6: tvcd toward top 10% sensitive countries
                if xy9casereduced(cn,wss)>=m9(region,wss*2)
                    Aa=xy5pre(cn,4+(wss-1)*11)/100;
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100;
                end
            elseif sce==7
                % Scenario 7: tvcd toward top 20% sensitive countries
                if xy9casereduced(cn,wss)>=m9(region,wss*2-1)
                    Aa=xy5pre(cn,4+(wss-1)*11)/100;
                else
                    Aa=xy5pre(cn,1+(wss-1)*11)/100;
                end
            elseif sce==8
                % Scenario 8: tvcd in Nov19-25
                Aa=xy5pre(cn,1+min(10,floor(-min(0,a(1,4))*10)+1)+(wss-1)*11)/100;
            end
            for d=1:7
                d0=d0+1;
                if simu(d0-1,1)>=0
                    simu(d0,1)=simu(d0-1,1)+simu(d0-1,2); % N
                else
                    simu(d0,1)=log(0.5);
                end
                simu(d0,2)=simu(d0-1,2)+simu(d0-1,3); % velocity
                simu(d0,3)=Aa; % A velocity
            end
        end
        xy8case(:,sce)=xy8case(:,sce)+floor(exp(simu(:,1)));
        xy8case(:,sce+region*8)=xy8case(:,sce+region*8)+floor(exp(simu(:,1)));
    end
end
xy6caselog=log(max(0.5,xy6case));
xy7caselog=log(max(0.5,xy7case));
xy8caselog=log(max(0.5,xy8case));

% Optimization procedure (see descriptions in our paper)
% Threshold (H) for the reduced number on 21 Jan. 2021: 250000,249000, ..., 2000, 1000, 900, ..., 100,90, ..., 10
ts=zeros(1,318); 
ts(1,1:300)=[300:-1:1];
ts(1,301:309)=[0.9:-0.1:0.1];
ts(1,310:318)=[0.09:-0.01:0.01];
xy10=zeros(338,318);
xy10cost=zeros(319,4);
for k=1:318
    for cn=1:338
        if sum(xy9(cn,1:24,1),2)<=ts(k) || cn==44 || cn==235
            continue;
        end
        j=24;
        casesave=xy9(cn,j,1);
        while casesave<ts(k)
            j=j-1;
            casesave=casesave+xy9(cn,j,1);
        end
        xy10(cn,k)=j;
        for j2=1:j
            xy10cost(k,1)=xy10cost(k,1)+xy9(cn,j2,1); % reduced cases, thousand
            xy10cost(k,2)=xy10cost(k,2)+xy9(cn,j2,2); % gdp * tvcd, billion USD
        end
    end
end

% Total reduced cases after k steps
xy10cost(k+1,1)=xy10cost(k+1,1)+sum(sum(xy9(:,:,1),2),1); % cases

% GDP total
gdpsum=(sum(gdpcn,1)-gdpcn(44,1)-gdpcn(235,1))/365*56/1e9; % gdp * tvcd, billion USD
xy10cost(k+1,2)=gdpsum*0.3; % gdp * tvcd, billion USD
xy10cost(:,3)=xy10cost(:,2)./gdpsum;
xy10cost(:,4)=xy8case(end,1)-xy10cost(:,1)*1000;

% Optimization by region
xy11cost=zeros(5,4);
tvcdreal=zeros(338,1);
for cn=1:338
    % Excluding the national data of China and US, where provincial and
    % state data are used for ML
    if cn==44 || cn==235
        continue;
    end
    
    % Real-time scenario
    if xy5pre(cn,23)==1 && xy5(cn,4,nt)<0
        tvcdreal(cn,1)=tvcdreal(cn,1)+gdpcn(cn,1)/365*56/1e9*(-xy5(cn,4,nt));
    end
end

% Scenarios of T0, T1, T2 and T3
xy11cost(1,1)=(xy8case(end,1)-xy8case(end,8))/1000; % reduced cases, thousand
xy11cost(1,2)=sum(tvcdreal,1); % gdp * tvcd, billion USD
for sce=1:3
    xy11cost(2+sce,1)=(xy8case(end,1)-xy8case(end,sce+1))/1000; % reduced cases, thousand
    xy11cost(2+sce,2)=xy10cost(end,2)/3*sce; % gdp * tvcd, billion USD
end
xy11cost(:,3)=xy11cost(:,2)./gdpsum;
xy11cost(:,4)=xy8case(end,1)-xy11cost(:,1)*1000;

subplot(2,4,8);
plot(xy10cost(:,2),xy10cost(:,4),'LineStyle','-','LineWidth',1,'Color',[1 0 0]); hold on;
plot(xy11cost(2:5,2),xy11cost(2:5,4),'s','MarkerEdgeColor','none','MarkerFaceColor',[1 0.8 0],'MarkerSize',5); hold on;
plot(xy11cost(2:5,2),xy11cost(2:5,4),'LineStyle','-','LineWidth',1,'Color',[1 0.8 0]); hold on;
plot(xy11cost(1,2),xy11cost(1,4),'p','MarkerEdgeColor',[0.8 0 0.8],'MarkerFaceColor',[1 1 1],'MarkerSize',10);

