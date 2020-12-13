% Author: Rong Wang
% Date: 2020.Oct.21
tic
clear;

% Date Info
dd2020=330; % number of days in 2020 (to be updated: 313 for 8 Nov; 330 for 25 Nov)
nt=floor((dd2020-1)/7)+1;
mday=[31,29,31,30,31,30,31,31,30,31,30,31;31,28,31,30,31,30,31,31,30,31,30,31;31,28,31,30,31,30,31,31,30,31,30,31;31,28,31,30,31,30,31,31,30,31,30,31;31,29,31,30,31,30,31,31,30,31,30,31;31,28,31,30,31,30,31,31,30,31,30,31;];
mdaysum=zeros(6,13);
dateinfo=zeros(366,2); i=0;
for m=1:12
    mdaysum(1:6,m+1)=mdaysum(1:6,max(1,m-1)+1)+mday(1:6,m);
    for d=1:mday(5,m)
        i=i+1;
        dateinfo(i,1)=m;
        dateinfo(i,2)=d;
    end
end

% Loading data of COVID-19 cases from different sources
ddtotal=dd2020;
cases=zeros(ddtotal,338,9); % 0 Aruba; 43 China; 234 US; 255 Zimbabwe; 256 Beijing; 286 Xinjiang; 287 Hawaii; 337 Alaska / last 256 for WHO
for cn=1:338
    if cn<=256
        filename=strcat('input\COVID_cases\World_ECDPC\ECDPCdata-',num2str(cn-1),'.dat');
    elseif cn<=287
        filename=strcat('input\COVID_cases\China\Chinadata-',num2str(cn-1),'.dat');
    elseif cn<=338
        filename=strcat('input\COVID_cases\USA\USAdata-',num2str(cn-1),'.dat');
    elseif cn<=594
        filename=strcat('input\COVID_cases\World_WHO\WHOdata-',num2str(cn-339),'.dat');
    end
    if exist(filename,'file')==0
        continue;
    else
        load(filename,'-mat');
    end
    
    % take away countries will few data
    if size(data,1)<5
        continue;
    end
    
    for i=1:size(data,1)
        if data(i,2)==0
            continue;
        end
        j=mdaysum(5,data(i,2))+data(i,1);
        cases(j,cn,1:2)=data(i,3:4);
        if i>1
            cases(j,cn,3:4)=data(i,3:4)-data(i-1,3:4);
        end
    end
    clear data;
end

% On Feb 20, shangdong reports the infections for a prison on one day,
% which is excluded: https://baijiahao.baidu.com/s?id=1659113674465060786&wfr=spider&for=pc
cases(51:ddtotal,271,1) = cases(51:ddtotal,271,1) - 200;
cases(51,271,3) = cases(51,271,3) - 200;
cases(51:ddtotal,44,1) = cases(51:ddtotal,44,1) - 200;
cases(51,44,3) = cases(51,44,3) - 200;

% Consider that in China "the ascertainment rate during the  
% early outbreak in Wuhan was estimated to be 0.23" Hao,2020,Nature
cases(:,44,1) = cases(:,44,1) + round(cases(:,273,1).* ((1-0.23)/0.23-1));
cases(:,44,3) = cases(:,44,3) + round(cases(:,273,3).* ((1-0.23)/0.23-1));
cases(:,273,1) = cases(:,273,1) + round(cases(:,273,1).* ((1-0.23)/0.23-1));
cases(:,273,3) = cases(:,273,3) + round(cases(:,273,3).* ((1-0.23)/0.23-1));

for cn=232:338
    % velocities
    [vv, xout]=velocity(cases(:,cn,3),3,7); % velocity(dailycase, half-time of movingaverage, half-time of velocity)
    cases(:,cn,5)=vv(:,1); % log of daily new cases
    cases(:,cn,6)=vv(:,2); % velocity, day-1
    cases(:,cn,7)=vv(:,3); % se of velocity
    cases(:,cn,8)=vv(:,4); % number of cases>=1
    
    % find the first week
    weekcase=zeros(nt,2);
    for i=1:nt
        weekcase(i,1)=mean(cases((i*7-6):min(i*7,dd2020),cn,3),1)*7;   % number of cases in a week
        weekcase(i,2)=mean(cases((i*7-6):min(i*7,dd2020),cn,6),1);  % average velocity
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
    idx2=find(weekcase(:,2)>-100);
    wi2=idx2(end); % number of weeks with covid-19
    dayend=min(wi2*7,dd2020);
    cases((wi*7-6):dayend,cn,9)=1;
end

save('input\daily_covid19.dat','cases');

clear
toc
clear
toc


