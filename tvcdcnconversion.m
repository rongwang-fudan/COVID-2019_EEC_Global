tic
clear;

% Country Info
cninfo2=load('input\countryinfo1031.txt');
gdpcn=mean(cninfo2(:,11:14),2);

% Date Info
dd2020=313; % 313 for Nov 8; 330 for Nov 25; 339 for Dec 4
dd2020case=330;
nt=floor((dd2020-1)/7)+1;

% loading coupled  data (Acase is coupled on Dno2)
load('input\tvcdcn.dat','-mat')
load('input\daily_covid19.dat','-mat');
load('input\grid_id_global.dat','-mat');

tvcd_case=zeros(338,45,2);
for cn=1:338
    region=cninfo2(cn,1);
    if cninfo2(cn,1)==-999 || tvcdcn(cn,1,1)<-100 || gdpcn(cn,1)<0
        continue;
    end
    % tvcdcn moving average
    tv=zeros(nt,13);
    for v=1:13
        tv(1:nt,v)=tvcdcn(cn,1:nt,v); % 1-4 for 360x720; 1 OMI-RF; 2 OMI-PLS; 3 TRO-RF; 4 TRO-PLS; 5-8 for 90x180
        moving=zeros(nt-2,2); 
        moving(1:(nt-2),1)=tvcdcn(cn,1:(nt-2),v); moving(1:(nt-2),2)=tvcdcn(cn,3:nt,v);
        tv(2:(nt-1),v)=(tv(2:(nt-1),v)+moving(:,1)+moving(:,2))/3;
    end
    % week cases
    weekcase=zeros(nt,8);
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
    % acceleration velocity moving average
    mweekcase=zeros(nt,1);
    mweekcase(end,1)=weekcase(end,5);
    casemoving=zeros(nt-2,2);casemoving(1:(nt-2),1)=weekcase(1:(nt-2),5); casemoving(1:(nt-2),2)=weekcase(3:nt,5);
    mweekcase(2:(nt-1),1)=(mweekcase(2:(nt-1),1)+casemoving(:,1)+casemoving(:,2))/3;
    idx=find(mweekcase~=0);
    if size(idx,1)<10
        continue;
    end
    %filling missed data (KNN)
    if size(idx,1)~=idx(end)-idx(1)+1
        idx1=find(mweekcase==0);
        mweekcase(idx1)=nan;
        mweekcase(idx(1):idx(end),1)=knnimpute(mweekcase(idx(1):idx(end)));
        idx2=find(isnan(mweekcase));
        mweekcase(idx2)=0;
    end  
    % days with covid-19
    idx=find(cases(:,cn,9)==1);
    if size(idx,1)==0
        continue;
    end
    tvcd_case(cn,:,1)=tv(:,1);
    tvcd_case(cn,:,2)=mweekcase;
end
A=zeros(338,5);
for cn=1:338
    A(cn,1)=cn;
    A(cn,2)=cninfo2(cn,1);
    if cn>=257 && cn<=287
        A(cn,2)=1;
    end
end
% The parameters of the CCM analysis
ccmp=zeros(338,45,2);Lt=zeros(338,45);
for cn=1:338
    t1=tvcd_case(cn,:,1)';% tvcd: OMI-RF;
    t2=tvcd_case(cn,:,2)';% acceleration velocity , % day-2
    idx11=find(t1~=0 & t2~=0);
    x=t1(idx11);y=t2(idx11);
    if size(x,1)<10
        continue;
    end
    % References:
    % - Clark A. T, et al. Spatial convergent cross mapping to detect causal relationships from short time series. Ecology 96, 1174-1181(2015).
    % - Ye H, et al. Distinguishing time-delayed causal interactions using convergent cross mapping. Sci. Rep. 5, 14750 (2015).
    [XR,eLag,eDim] = phaseSpaceReconstruction(x);
    tau = 1;
    E = eDim;
    LMN = E+1;
    r=corrcoef(x,y);A(cn,3)=r(1,2);A(cn,4)=tau;A(cn,5)=E;
    L=zeros(size(x,1),1);
    for i = E+2+(E-1)*tau:size(x,1)
        L(i-(E-1)*tau-E-1,1)=i;
    end
    idx1=find(L==0); L(idx1)=[];
    Lt(cn,1:size(L,1))=L;
    len=length(L);SC=zeros(2,len);
    for i = 1:len
        [SC(:,i)] = SugiLM(x(1:L(i)),y(1:L(i)),tau,E,LMN );
     end
    ccmp(cn,1:len,1)=SC(1,:);ccmp(cn,1:len,2)=SC(2,:);
    clear L;
end

tt=A(:,3);idx=find(tt==0);A(idx,:)=[];
idxnan=find(isnan(ccmp));ccmp(idxnan)=0;
vali=zeros(size(A,1),5);
for i=1:size(A,1)
    if A(i,3)>0
        vali(i,:)=A(i,:);
    end
end
idx=find(vali(:,1)==0);vali(idx,:)=[];
vaccmp=zeros(size(vali,1),45,2);
Lt1=zeros(size(vali,1),45);
cnvar1=zeros(size(vali,1),45,2);
for i=1:size(vali,1)
    vaccmp(i,:,:)=ccmp(vali(i,1),:,:);
    Lt1(i,:)=Lt(vali(i,1),:);
    cnvar1(i,:,:)=tvcd_case(vali(i,1),:,:);
end
Cnvar=zeros(size(vali,1),45,2);
for i=1:size(vali,1)
    t1=cnvar1(i,:,1);
    t2=cnvar1(i,:,2);
    idx11=find(t1~=0 & t2~=0);
    t11=t1(idx11);t22=t2(idx11);
    Cnvar(i,1:size(t11,2),1)=t11;
    Cnvar(i,1:size(t22,2),2)=t22;
end
% The random noises in the time series of ΔNO2 and A were reduced according to the upward and downward trends in the correlation series
merge1=zeros(size(vali,1),50);merge1(:,1:5)=vali;merge1(:,6:50)=vaccmp(:,:,1);
merge2=zeros(size(vali,1),50);merge2(:,1:5)=vali;merge2(:,6:50)=vaccmp(:,:,2);
Lt2=zeros(size(vali,1),50);Lt2(:,1:5)=vali;Lt2(:,6:50)=Lt1;
% diff if
diffcn_no2_acase=zeros(size(merge1,1),45,3);
w1=zeros(size(merge1,1),45,3);w2=zeros(size(merge1,1),45,3);w3=zeros(size(merge1,1),45,2);
for cn=1:size(merge1,1)
    tt=merge1(cn,6:50);idx=find(tt~=0);tt=tt(idx);
    
    % step 1 diff
    for i=1:size(tt,2)-1
        diffcn_no2_acase(cn,i,1)=tt(i+1)-tt(i);
        if tt(i+1)-tt(i)>-0.05
            w1(cn,i,1)=w1(cn,i,1)+3;
        end
    end
    % step 2 diff
    for i=1:size(tt,2)-2
        diffcn_no2_acase(cn,i,2)=tt(i+2)-tt(i);
        if tt(i+2)-tt(i)>-0.1
            w1(cn,i,2)=w1(cn,i,2)+2;
        end
    end
    % step 3 diff
    for i=1:size(tt,2)-3
        diffcn_no2_acase(cn,i,3)=tt(i+3)-tt(i);
        if tt(i+3)-tt(i)>-0.15
            w1(cn,i,3)=w1(cn,i,3)+1;
        end
    end
    
    % max diff location
    for j=1:3
        [v,index]=max(diffcn_no2_acase(cn,:,j));
        w2(cn,index,j)=w2(cn,index,j)+1;
    end
    w3(cn,:,1)=w1(cn,:,1)+w1(cn,:,2)+w1(cn,:,3);
    w3(cn,:,2)=w2(cn,:,1)+w2(cn,:,2)+w2(cn,:,3);
end
A=w3(:,:,1);
A1=w3(:,:,2);
no2_acase=zeros(size(merge1,1),50);
acase_no2=zeros(size(merge1,1),50);
Ltvar=zeros(size(merge1,1),50);
for cn=1:size(merge1,1)
    no2_acase(cn,1:5)=merge1(cn,1:5);
    acase_no2(cn,1:5)=merge2(cn,1:5);
    Ltvar(cn,1:5)=Lt2(cn,1:5);
    
    t1=w3(cn,:,1);idx=find(t1>=3);
    t2=merge1(cn,6:50); t21=merge2(cn,6:50);lt1=Lt2(cn,6:50);
    
    no2_acase(80,6:14+5)=merge1(80,6:14+5);
    acase_no2(80,6:14+5)=merge2(80,6:14+5);
    Ltvar(80,6:14+5)=Lt2(80,6:14+5);
    
    idiff=diff(idx); idx1=find(idiff~=1);
    if size(idiff,2)>1 && size(idx1,2)==0
        no2_acase(cn,6:size(idx,2)+5)=t2(idx);
        no2_acase(cn,size(idx,2)+5+1)=t2(idx(end)+1);
         
        acase_no2(cn,6:size(idx,2)+5)=t21(idx);
        acase_no2(cn,size(idx,2)+5+1)=t21(idx(end)+1);
        
        Ltvar(cn,6:size(idx,2)+5)=lt1(idx);
        Ltvar(cn,size(idx,2)+5+1)=lt1(idx(end)+1);
        
    else
        idx2=find(idx1>1);
        if size(idx2,2)>0
            idx3=find(idx1==idx1(idx2(1)));
            if idx3~=1
                t3=t2(idx1(find(idx1==idx1(idx2(1)))-1)+1:idx1(find(idx1==idx1(idx2(1))))+1);
                t4=t1(idx1(find(idx1==idx1(idx2(1)))-1)+1:idx1(find(idx1==idx1(idx2(1))))+1);
                
                t5=t21(idx1(find(idx1==idx1(idx2(1)))-1)+1:idx1(find(idx1==idx1(idx2(1))))+1);
                
                t6=lt1(idx1(find(idx1==idx1(idx2(1)))-1)+1:idx1(find(idx1==idx1(idx2(1))))+1);
            else
                t3=t2(1:idx1(find(idx1==idx1(idx2(1))))+1);
                t4=t1(1:idx1(find(idx1==idx1(idx2(1))))+1);
                
                t5=t21(1:idx1(find(idx1==idx1(idx2(1))))+1);
                
                t6=lt1(1:idx1(find(idx1==idx1(idx2(1))))+1);
            end
            
            if size(t3,2)>=4
                idx4=find(t4==6);
                if  size(idx4,2)>=4
                    t3=t3(idx4(1):end);
                    no2_acase(cn,6:size(t3,2)+5)=t3;
                    
                    t5=t5(idx4(1):end);
                    acase_no2(cn,6:size(t5,2)+5)=t5;
                    
                    t6=t6(idx4(1):end);
                    Ltvar(cn,6:size(t5,2)+5)=t6;
                    
                end
            end
        end
    end
end

tt=no2_acase(:,6);
idx=find(tt==0);
no2_acase(idx,:)=[];
acase_no2(idx,:)=[];
Ltvar(idx,:)=[];
Cnvar(idx,:,:)=[];
load('ccm\noadj_no2.dat','-mat');
num_valid_countryoregion=size(Cnvar,1);
no2_acase=zeros(num_valid_countryoregion,45,3);
Levar=zeros(num_valid_countryoregion,48);Levar(:,1)=Ltvar(:,1);
Levar(:,2:48)=Ltvar(:,4:end);
for i=1:size(Ltvar,1)
    t1=Cnvar(i,:,:);idx1=find(t1(:,:,1)~=0);t1=t1(1,idx1,:);
    t2=Ltvar(i,:);idx2=find(t2~=0);t2=t2(idx2);
    x1=t1(:,1:t2(end),1);x2=noadj_no2(i,1:t2(end));y=t1(:,1:t2(end),2);
    no2_acase(i,1:size(x1,2),1)=x1;
    no2_acase(i,1:size(x1,2),2)=x2;
    no2_acase(i,1:size(x1,2),3)=y;
end

idx=find(no2_acase==0);no2_acase(idx)=-999;

save('ccm\no2_acase.dat','no2_acase');
save('ccm\Levar.dat','Levar');

toc;


