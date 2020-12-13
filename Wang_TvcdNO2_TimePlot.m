% Author: Rong Wang
% Date: 2020.Oct.21
tic
clear;

% Country Info
cninfo2=load('input\countryinfo1031.txt');

% Date Info
dd2020=313; % number of days in 2020 (to be updated: 313 for 8 Nov; 330 for 25 Nov)
nt=floor((dd2020-1)/7)+1;

% dNO2 + weekly cases
load('input_tvcdbyregion\tvcdcn.dat','-mat'); % tvcdcn=ones(338,nt,9); tvcdcn=tvcdcn.*(-999);
load('input_dailycasesbyregion\daily_covid19.dat','-mat'); % cases=zeros(ddtotal,338,8); 1 cumul case; 2 cumul death; 3 case; 4 death; 5 log cases; 6 velocity; 7 p value; 8 number of cases>=1


% Colorbar from -0.2 to 0.2
nbar=21;
mycolor=zeros(nbar*2+1,3);
mycolor(1,1:3)=[0.8 1 1];
mycolor(2:(nbar+1),1:3)=summer(nbar); % green to yellow
mycolor((nbar*2+1):-1:(nbar+2),1:3)=autumn(nbar); % red to yellow

% Time series of NO2 TVCD change and Velocity by week
aa1=ones(338,2); aa1=aa1.*(-999);
aa2=ones(nt,338,3); aa2=aa2.*(-999);
for cn=1:338
    if cninfo2(cn,1)==-999 || tvcdcn(cn,1,1)<-100
        continue;
    end
    
    % moving average of tvcd change
    aa2(1:nt,cn,1)=tvcdcn(cn,1:nt,1); % 1-4 for 360x720; 1 OMI-RF; 2 OMI-PLS; 3 TRO-RF; 4 TRO-PLS; 5-8 for 90x180
    moving=zeros(nt-2,2);moving(1:(nt-2),1)=tvcdcn(cn,1:(nt-2),1); moving(1:(nt-2),2)=tvcdcn(cn,3:nt,1);
    aa2(2:(nt-1),cn,1)=(aa2(2:(nt-1),cn,1)+moving(:,1)+moving(:,2))/3;
    
    % uncertainty (se) in tvcd
    aa2(1:nt,cn,3)=tvcdcn(cn,1:nt,10); % 1-4 for 360x720; 1 OMI-RF; 2 OMI-PLS; 3 TRO-RF; 4 TRO-PLS; 5-8 for 90x180
    moving=zeros(nt-2,2);moving(1:(nt-2),1)=tvcdcn(cn,1:(nt-2),10); moving(1:(nt-2),2)=tvcdcn(cn,3:nt,10);
    aa2(2:(nt-1),cn,3)=(aa2(2:(nt-1),cn,1)+moving(:,1)+moving(:,2))/3;
    
    % find the first week by region
    weekcase=zeros(nt,3);
    for i=1:nt
        weekcase(i,1)=sum(cases((i*7-6):min(i*7,dd2020),cn,3),1);   % number of cases in a week
        weekcase(i,2)=mean(cases((i*7-6):min(i*7,dd2020),cn,6),1);  % velocity
        weekcase(i,3)=mean(cases((i*7-6):min(i*7,dd2020),cn,8),1);
    end
    startcase=100;
    idx=find(weekcase(:,1)>startcase & weekcase(:,2)>-100);
    if size(idx,1)<1
        continue;
    end
    clear wi;
    if idx(1)==1
        wi=1;
    else
        w2=idx(1)-1;
        for i=w2:-1:1
            if weekcase(i,2)<-100 || weekcase(i,1)<=20
                wi=i+1; break;
            end
            if weekcase(i,1)>=weekcase(i+1,1)
                if (i-3)>0
                    if weekcase(i,1)>weekcase(i-1,1) && weekcase(i-1,1)>weekcase(i-2,1) && weekcase(i-2,1)>weekcase(i-3,1)
                    else
                        wi=i+1; break;
                    end
                else
                    wi=i+1; break;
                end
            end
        end
        if i==1
            wi=1;
        end
    end
    
    % number of weeks with covid-19
    idx2=find(weekcase(:,2)>-100);
    aa1(cn,1)=wi;
    aa1(cn,2)=idx2(end)-wi+1; 
    
    % velocity of covid-19
    for i=wi:idx2(end)
        jx=zeros(7,1);jn=0;
        for j=(i*7-6):(i*7)
            if cases(min(j,dd2020),cn,6)>-100
                jn=jn+1;
                jx(jn,1)=cases(min(j,dd2020),cn,6); % 1 cumul case; 2 cumul death; 3 case; 4 death; 5 log cases; 6 velocity; 7 p value; 8 number of cases>=1
            end
        end
        if weekcase(i,3)>=10
            if jn<1
                aa2(i,cn,2)=0;
            else
                aa2(i,cn,2)=mean(jx(1:jn,1),1);
            end
        end
    end
end

% Plot NO2 TVCD change and the velocity of daily cases by week
for region=1:6
    subplot(2,3,region);
    xylow=zeros(1000,2); k1=0;
    for cn=1:338
        for i=1:nt
            if cninfo2(cn,1)==region && aa1(cn,1)<20 && aa1(cn,1)>0 && cn~=44 && cn~=235 && aa2(i,cn,2)<-100
                k1=k1+1;
                xylow(k1,1)=i; xylow(k1,2)=aa2(i,cn,1); % no velocity
            end
        end
    end
    for cn=1:338
        if cninfo2(cn,1)==region && aa1(cn,1)<20 && aa1(cn,1)>0 && cn~=44 && cn~=235
            plot(aa2(1:nt,cn,1),'LineStyle','-','LineWidth',0.1,'Color',[0.8 0.8 0.8]); hold on;
            startweek=aa1(cn,1);
            for nwe=1:3
                if (aa2(startweek+nwe*2-1,cn,1)+min(0.4,aa2(startweek+nwe*2-1,cn,3)))<0
                end
            end
        else
            if region==1 && cn>=257 && cn<=287
                plot(aa2(1:nt,cn,1),'LineStyle','-','LineWidth',0.1,'Color',[0.8 0.8 0.8]); hold on;
                for i=1:nt
                    k1=k1+1;
                    xylow(k1,1)=i; xylow(k1,2)=aa2(i,cn,1); % no velocity
                end
            end
        end
    end
    plot(xylow(1:k1,1),xylow(1:k1,2),'o','MarkerEdgeColor','none','MarkerFaceColor',mycolor(1,1:3),'MarkerSize',4); hold on;
    for cn=1:338
        for i=1:nt
            if cninfo2(cn,1)==region && aa1(cn,1)<20 && aa1(cn,1)>0 && cn~=44 && cn~=235 && aa2(i,cn,2)>=-100
                colorn=min(nbar*2,max(1,floor(aa2(i,cn,2)*100)+nbar+1))+1; % from -0.2 to 0.2
                plot(i,aa2(i,cn,1),'o','MarkerEdgeColor','none','MarkerFaceColor',mycolor(colorn,1:3),'MarkerSize',4); hold on;
            elseif region==1 && cn>=257 && cn<=287 && aa2(i,cn,2)>=-100
                colorn=min(nbar*2,max(1,floor(aa2(i,cn,2)*100)+nbar+1))+1; % from -0.2 to 0.2
                plot(i,aa2(i,cn,1),'o','MarkerEdgeColor','none','MarkerFaceColor',mycolor(colorn,1:3),'MarkerSize',4); hold on;
            end
        end
    end
    axis([0 nt -1 1]);
%     set(gca,'xtick',-5:55:50);
%     set(gca,'ytick',-2:4:2);
%     set(gca,'color',[0.965,0.898,1]);
end
