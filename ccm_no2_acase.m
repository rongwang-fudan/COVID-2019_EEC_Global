tic
clear;

% Loading data
load('ccm\no2_acase.dat','-mat');% no2_acase=zeros(78,45,3): 1 ΔNO2; 2 ΔNO2*; 3 The acceleration of COVID-19 cases (A, in % day-2).
load('ccm\Levar.dat','-mat');% Levar=zeros(78,48): 1 territories code; 2 reconstruction time step; 3 The optimal embedding dimension (E);4:48 Time step.

% The convergent cross-mapping method (CCM)
% The convergent cross mapping algorithm package of MATLAB using the SugiLM function
% (https://www.mathworks.com/matlabcentral/fileexchange/52964-convergent-cross-mapping).
ccm=zeros(78,45,4); % 1 CCM analysis of ΔNO2 and A;2 CCM analysis of A andΔNO2;3 CCM analysis of ΔNO2* and A;4 CCM analysis of A andΔNO2*.
for i=1:78
    % The ΔNO2,ΔNO2* and A information for each territories.
    x1=no2_acase(i,:,1);x2=no2_acase(i,:,2);y=no2_acase(i,:,3);
    idx=find(x1~=-999 & x2~=-999 & y~=-999);
    x1=x1(idx);x2=x2(idx);y=y(idx);
    % Parameters of the SugiLM function.
    arg=Levar(i,:);tau = arg(2);E = arg(3);LMN = E+1;
    % Time step series for CCM analysis in each territories.
    ts=arg(4:end);L=ts(idx);
    % The CCM analysis for ΔNO2,ΔNO2* and A in each territories.
    len=length(L);
    adj_no2_acase=zeros(2,len);% 1 CCM analysis of ΔNO2 and A; 2 CCM analysis of A andΔNO2.
    noadj_no2_acase=zeros(2,len);% 1 CCM analysis of ΔNO2* and A; 2 CCM analysis of A andΔNO2*.
    for m = 1:len
        [adj_no2_acase(:,m)] = SugiLM(x1(1:L(m))',y(1:L(m))',tau,E,LMN );
        [noadj_no2_acase(:,m)] = SugiLM(x2(1:L(m))',y(1:L(m))',tau,E,LMN );
    end
    ccm(i,1:len,1)=adj_no2_acase(1,:);ccm(i,1:len,2)=adj_no2_acase(2,:);
    ccm(i,1:len,3)=noadj_no2_acase(1,:);ccm(i,1:len,4)=noadj_no2_acase(2,:);
    clear L;
end

% plot fig.2 A B D E
colormap1=[107 49 9; 166 76 14;223 102 19;241 148 85; 248 199 166; 90 90 90 ;207 205 205;];colormap1=colormap1./255;% fig.2 A B
colormap2=[20 49 76;41 107 167;64 139 208;105 164 217; 191 216 239; 90 90 90; 207 205 205;];colormap2=colormap2./255;% fig.2 D E
cncode=[33 22 52 61 21 38];% 33 Namibia; 22 Israel; 52 Hubei; 61 California; 21 Iran; 38 Paraguay
for k=1:4
    subplot(2,2,k);t1=ccm(:,:,k);
    % CCM analysis of ΔNO2 and A
    if k<=2
        % Each territories
        for cn=1:size(ccm,1)
            if ismember(cn,cncode)==0
                t2=t1(cn,:);idx=find(t2~=0 & ~isnan(t2));t2=t2(idx);
                plot(t2,'LineStyle','-','LineWidth',0.2,'Color',colormap1(7,1:3)); hold on;
            end
        end
        % Namibia; Israel; Hubei; California; Iran; Paraguay
        for cn=1:6
            t2=t1(cncode(cn),:);idx=find(t2~=0 & ~isnan(t2));t2=t2(idx);
            plot(t2,'LineStyle','--','LineWidth',1.5,'Color',colormap1(cn,1:3)); hold on;
        end
    % CCM analysis of ΔNO2* and A;
    else
        % Each territories
        for cn=1:size(ccm,1)
            if ismember(cn,cncode)==0
                t2=t1(cn,:);idx=find(t2~=0 & ~isnan(t2));t2=t2(idx);
                plot(t2,'LineStyle','-','LineWidth',0.2,'Color',colormap2(7,1:3)); hold on;
            end
        end
        % Namibia; Israel; Hubei; California; Iran; Paraguay
        for cn=1:6
            t2=t1(cncode(cn),:);idx=find(t2~=0 & ~isnan(t2));t2=t2(idx);
            plot(t2,'LineStyle','--','LineWidth',1.5,'Color',colormap2(cn,1:3)); hold on;
        end
    end
    axis([1 30 -1 1]);
end