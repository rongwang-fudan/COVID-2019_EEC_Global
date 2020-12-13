% Author: Rong Wang
% Date: 2020.Oct.21
tic
clear;

% Date Info
dd2020=313; % number of days in 2020 (to be updated: 313 for 8 Nov; 330 for 25 Nov)
ddyear=[0,366,731,1096,1461,1461+dd2020];
d1=ddyear(4)+1;
d2=ddyear(6);
mydata=zeros(34310,ddyear(5)+dd2020,9);

% Loading data
load('input\popfrac05.dat','-mat');
for yy=2016:2020

% NO2 Info
load(strcat('no2_me\NO2_Metro_2016_2020\NO2_',num2str(yy),'.dat'),'-mat'); % OMI?molecular/cm^2 (360,720,366)
if yy==2020
    no2grid_day_omi = NO2grid_day ./ 6.02214e16 * 46; % mg/m2
    clear NO2grid_day;
else
    no2grid_day_omi = DATAgrid_day ./ 6.02214e16 * 46; % mg/m2
    clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\Temperature_',num2str(yy),'.dat'),'-mat'); % Celsium degree
if yy==2020
    temperature=Temperaturegrid_day; clear Temperaturegrid_day;
else
    temperature=DATAgrid_day; clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\Pressure_',num2str(yy),'.dat'),'-mat'); % Pa
if yy==2020
    pressure=Pressuregrid_day; clear Pressuregrid_day;
else
    pressure=DATAgrid_day; clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\precipitation_',num2str(yy),'.dat'),'-mat'); % 0-50
if yy==2020
    precipitation=precipitationgrid_day; clear precipitationgrid_day;
else
    precipitation=DATAgrid_day; clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\RH_',num2str(yy),'.dat'),'-mat'); % 0-100
if yy==2020
    RH=RHgrid_day; clear RHgrid_day;
else
    RH=DATAgrid_day; clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\Uwind_',num2str(yy),'.dat'),'-mat'); % m/s
if yy==2020
    uwind=Uwindgrid_day; clear Uwindgrid_day;
else
    uwind=DATAgrid_day; clear DATAgrid_day;
end

load(strcat('no2_me\NO2_Metro_2016_2020\Vwind_',num2str(yy),'.dat'),'-mat'); % m/s
if yy==2020
    vwind=Vwindgrid_day; clear Vwindgrid_day;
else
    vwind=DATAgrid_day; clear DATAgrid_day;
end

if yy>=2019
    load(strcat('no2_me\NO2_Metro_2016_2020\S5_NO2_',num2str(yy),'.dat'),'-mat'); % Tropomi?mol/m^2
    no2grid_day_tropomi = S5_NO2grid_day * 460000; % mg/m2
    clear S5_NO2grid_day;
    
    load(strcat('no2_me\NO2_Metro_2016_2020\PBL_',num2str(max(2019,yy)),'.dat'),'-mat'); % m
    pbl=PBLgrid_day; clear PBLgrid_day;
end

% compiling data in our format
gi=0;
yy2=yy-2015;
for i=1:360
    for j=1:720
        % country and urban masks
        if popfrac05(i,j,1)>2500 % > 1 capita per km2
            gi=gi+1;
            for d=(ddyear(yy2)+1):ddyear(yy2+1)
                mydata(gi,d,1)=no2grid_day_omi(i,j,d-ddyear(yy2));
                if yy>=2019
                    mydata(gi,d,2)=no2grid_day_tropomi(i,j,d-ddyear(yy2));
                    mydata(gi,d,3)=pbl(i,j,d-ddyear(yy2));
                end
                mydata(gi,d,4)=temperature(i,j,d-ddyear(yy2));
                mydata(gi,d,5)=pressure(i,j,d-ddyear(yy2));
                mydata(gi,d,6)=precipitation(i,j,d-ddyear(yy2));
                mydata(gi,d,7)=RH(i,j,d-ddyear(yy2));
                mydata(gi,d,8)=uwind(i,j,d-ddyear(yy2));
                mydata(gi,d,9)=vwind(i,j,d-ddyear(yy2));
            end
        end
    end
end

end

gi=0;
gi2=0;
mydata2=mydata;
popfrac05new=popfrac05;
for i=1:360
    for j=1:720
        % country and urban masks
        if popfrac05(i,j,1)>2500 % > 1 capita per km2
            gi=gi+1;
            gi2=gi2+1;
            
            % OMI NO2
            tt1=(mydata(gi,1:d2,1))';
            idx=find(tt1>0);
            if size(idx,1)<600
                gi2=gi2-1;
                popfrac05new(i,j,1:2)=0;
                continue;
            end
            mydata2(gi2,1:d2,1)=fillmissdata_no2(tt1);
            
            % TROPOMI NO2
            tt1=(mydata(gi,d1:d2,2))';
            idx=find(tt1>0);
            if size(idx,1)<200
                mydata2(gi2,d1:d2,2)=mydata2(gi2,d1:d2,1);
            else
                mydata2(gi2,d1:d2,2)=fillmissdata_no2(tt1);
            end
            mydata2(gi2,1:(d1-1),2)=mean(mydata2(gi2,d1:d2,2),2);
            
            % PBL
            tt1=(mydata(gi,d1:d2,3))';
            idx=find(tt1>0);
            if size(idx,1)<200
                mydata2(gi2,d1:d2,3)=0;
            else
                mydata2(gi2,d1:d2,3)=fillmissdata_mete(tt1);
            end
            mydata2(gi2,1:(d1-1),3)=mean(mydata2(gi2,d1:d2,3),2);
            
            for v=1:6
                tt1=(mydata(gi,1:d2,v+3))';
                idx=find(tt1>0);
                if size(idx,1)<600
                    mydata2(gi2,1:d2,v+3)=0;
                else
                    mydata2(gi2,1:d2,v+3)=fillmissdata_mete(tt1);
                end
            end
        end
    end
end

clear mydata;
mydata=mydata2(1:gi2,:,:);
save('input\mydata.dat','mydata');

clear
toc


