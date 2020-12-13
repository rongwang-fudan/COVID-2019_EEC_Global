% Author: Rong Wang
% Date: 2020.Oct.21
tic
clear;

% Date Info
dd2020=313; % number of days in 2020 (to be updated: 313 for 8 Nov; 330 for 25 Nov)
ddyear=[0,366,731,1096,1461,1461+dd2020];
d1=ddyear(4)+1;
d2=ddyear(6);
d3=ddyear(5);
nt=floor((dd2020-1)/7)+1;

% Pop Data
pop=load('input\pop.txt'); % 256 countries (cap)
load('input\popfrac05new.dat','-mat'); % popfrac05new(360,720,1);

% Dummy variables
yt=zeros(d2,2); % 1 for year; 2 for season
for y=1:5
    for d=(ddyear(y)+1):ddyear(y+1)
        yt(d,1)=y;
        s=floor((d-ddyear(y)+30)/90)+1; % 1 for winter; 2 for spring / summer; 3 for summer
        if s==1 || s==5
            yt(d,2)=1;
        elseif s==2 || s==4
            yt(d,2)=2;
        elseif s==3
            yt(d,2)=3;
        end
    end
end
v_dum=zeros(d2,nt);
for d=1:dd2020
    v_dum(ddyear(5)+d,floor((d-1)/7)+1)=1; % give a value of 1, elsewhere 0
end

% NO2 + Mete Info
load('input_meteorologydata\mydata.dat','-mat'); % mydata=zeros(34310,ddyear(5)+dd2020,9); input1207xiong_unconfirmed input1206confirmed
% 1 omi no2; 2 tropomi no2; 3 pbl; 4 temp; 5 press; 6 preci; 7 rh; 8 uwind; 9 vwind

% Fix-effects model started here
diff=ones(360,720,nt,20); diff=diff.*(-999);
regres=zeros(360,720,4);
gi=0;
for i=1:360
    display(i);
    for j=1:720
        if popfrac05new(i,j,1)<=2500 % > 1 capita per km2
            continue;
        end
        gi=gi+1;
        xy=zeros(d2,10); vi=0; vi89=0;
        for v=1:9
            aa=mydata(gi,1:d2,v);
            if var(aa)~=0
                vi=vi+1;
                xy(1:d2,vi)=aa;
            else
                if v==8 || v==9
                    vi89=1;
                end
            end
        end
        if vi89==0
            vi=vi+1; % ensure that data is available for v=8,9 
            xy(:,vi)=sqrt(xy(:,vi-2).*xy(:,vi-2)+xy(:,vi-1).*xy(:,vi-1));
        end
        
        % Fixed-effects model - using OMI
        X=xy(:,3:vi);
        
        % Machine learning for NO2-meteorology
        Mdl=fitrensemble(X(1:d3,:),xy(1:d3,1),'Method','LSBoost','NumLearningCycles',10);
        no2ml=predict(Mdl,X);
        no2=xy(:,1)-no2ml;
        [b2,bint2,r2,rint2,stats2]=regress(no2,[ones(d2,1) yt(:,1:2) v_dum]);
        [b4,bint4,r4,rint4,stats4]=regress(xy(:,1),[ones(d2,1) X yt(:,1:2) v_dum]);
        
        % Statistics
        regres(i,j,1)=1-resubLoss(Mdl)/var(xy(1:d3,1));
        regres(i,j,2)=stats4(1);
        
        % Change of NO2 by week
        for w=1:nt
            no2ref=zeros(100,8); c1=0; c2=0;
            for y=1:5
                for d=(w*7-6):(w*7)
                    di2=max(1,min(d2,ddyear(y)+d));
                    if y==5
                        c1=c1+1;
                        no2ref(c1,1)=no2(di2,1)-r2(di2,1)+no2ml(di2); % prediction machine learning with covid-19 2020
                        no2ref(c1,2)=no2(di2,1)-r2(di2,1)+no2ml(di2)-b2(w+3); % prediction machine learning without covid-19 2020
                        no2ref(c1,3)=xy(di2,1)-r4(di2,1); % prediction linear regression with covid-19 2020
                        no2ref(c1,4)=xy(di2,1)-r4(di2,1)-b4(w+vi+1); % prediction linear regression without covid-19 2020
                    else
                        c2=c2+1;
                        no2ref(c2,5)=no2(di2,1)-r2(di2,1)+no2ml(di2)+b2(w+3); % prediction machine learning with covid-19 2016-2019
                        no2ref(c2,6)=no2(di2,1)-r2(di2,1)+no2ml(di2); % prediction machine learning without covid-19 2016-2019
                        no2ref(c2,7)=xy(di2,1)-r4(di2,1)+b4(w+vi+1); % prediction linear regression with covid-19 2016-2019
                        no2ref(c2,8)=xy(di2,1)-r4(di2,1); % prediction linear regression without covid-19 2016-2019
                    end
                end
            end
            for y=1:4
                diff(i,j,w,y)=mean(no2ref(1:c1,y),1); % 2020
                diff(i,j,w,y+4)=mean(no2ref(1:c2,y+4),1); % 2016-2019
            end
            diff(i,j,w,9)=(bint2(w+3,2)-bint2(w+3,1))/b2(w+3)/2/1.96;
            diff(i,j,w,10)=(bint4(w+vi+1,2)-bint4(w+vi+1,1))/b4(w+vi+1)/2/1.96;
        end
        
        % Fixed-effects model - using tropomi
        X=xy(:,3:vi);
        
        % Machine learning for NO2-meteorology
        Mdl=fitrensemble(X(d1:d3,:),xy(d1:d3,2),'Method','LSBoost','NumLearningCycles',10); % 2019-2019
        no2ml=predict(Mdl,X(d1:d2,:)); % 2019-2020
        no2=xy(d1:d2,2)-no2ml; % 2019-2020
        [b2,bint2,r2,rint2,stats2]=regress(no2,[ones(d2-d1+1,1) v_dum(d1:d2,:)]); % 2019-2020
        [b4,bint4,r4,rint4,stats4]=regress(xy(d1:d2,2),[ones(d2-d1+1,1) X(d1:d2,:) v_dum(d1:d2,:)]); % 2019-2020
        
        % Statistics
        regres(i,j,3)=1-resubLoss(Mdl)/var(xy(d1:d3,1));
        regres(i,j,4)=stats4(1);
        
        % Change of NO2 by week
        for w=1:nt
            no2ref=zeros(100,8); c1=0; c2=0;
            for y=4:5
                for d=(w*7-6):(w*7)
                    di2=max(1,min(d2-d1+1,ddyear(y)+d-d1+1));
                    di3=max(1,min(d2-d1+1,ddyear(y)+d));
                    if y==5
                        c1=c1+1;
                        no2ref(c1,1)=no2(di2,1)-r2(di2,1)+no2ml(di2); % prediction machine learning with covid-19 2020
                        no2ref(c1,2)=no2(di2,1)-r2(di2,1)+no2ml(di2)-b2(w+1); % prediction machine learning without covid-19 2020
                        no2ref(c1,3)=xy(di3,1)-r4(di2,1); % prediction linear regression with covid-19 2020
                        no2ref(c1,4)=xy(di3,1)-r4(di2,1)-b4(w+vi-1); % prediction linear regression without covid-19 2020
                    else
                        c2=c2+1;
                        no2ref(c2,5)=no2(di2,1)-r2(di2,1)+no2ml(di2)+b2(w+1); % prediction machine learning with covid-19 2016-2019
                        no2ref(c2,6)=no2(di2,1)-r2(di2,1)+no2ml(di2); % prediction machine learning without covid-19 2016-2019
                        no2ref(c2,7)=xy(di3,1)-r4(di2,1)+b4(w+vi-1); % prediction linear regression with covid-19 2016-2019
                        no2ref(c2,8)=xy(di3,1)-r4(di2,1); % prediction linear regression without covid-19 2016-2019
                    end
                end
            end
            for y=1:4
                diff(i,j,w,y+10)=mean(no2ref(1:c1,y),1); % 2020
                diff(i,j,w,y+14)=mean(no2ref(1:c2,y+4),1); % 2016-2019
            end
            diff(i,j,w,19)=(bint2(w+1,2)-bint2(w+1,1))/b2(w+1)/2/1.96;
            diff(i,j,w,20)=(bint4(w+vi-1,2)-bint4(w+vi-1,1))/b4(w+vi-1)/2/1.96;
        end
    end
end
save('diff_weeks_omi360x720_machinelearning_bg.dat','diff');
save('regres360x720_machinelearning_bg.dat','regres');
clear
toc

