% PhD Thesis- Case 1- variable Power 
% Calculation of the transmitted power s.t. both users try to maximize the
% rate/R_min ratio for every subchannel
% Calculation of MUD rate and determination of MUD mode for Sx12
% subchannels (S resource blocks)

% PROJECTED SUBGRADIENT ALGORITHM_N=3_08/05/2017
% FILE IS BASED ON METHODS= FOR ALL STEPS

% Subgradient Projection Method composed with Bisection Method for lambda
% 1st BS has U(1) users\ 2nd BS has U(2) users 
% 1st step: SC Assignment wrt rate/Rmin with assuming equal power
% SCA== proportional to SIR ratio
% SCA step with MUD rate
% 2nd step: Power Allocation with the SC assignment performed in the 1st step
%           PA proportional to Rmin requirements.
% lambda_opt12 discarded. lambda_opt12=lambda_opt11

% clc
clear
tic
% load test_channel.mat
% load test_channel_S100.mat
% load N2U2_500sim_channelmatrix.mat;
% load N2U2_500sim_radius2000m_channelmatrix.mat;
% load N2U2_500sim_channelmatrix_radius2000m_randomusers.mat;
% load N2U2_500sim_channelmatrix_radius2000m_inner_outer_randomusers.mat; %Scenario IIb-inner-outer
% load N2U2_1sim_S6_Wyner_channelmatrix_BESTcase.mat
% load N2U2_1sim_S6_Wyner_channelmatrix_WORSTcase.mat

%load N3U2_500sim_radius2000_channelmatrix.mat
% load N3U2_500sim_channelmatrix_radius2000m_inner_outer_randomusers.mat
load N3U2_500sim_channelmatrix_rb2000m_NONOVERLAP_inner_outer_randomusers.mat
% load N3U2_500sim_channelmatrix_rb5000m_NONOVERLAP_inner_outer_randomusers.mat

% load S5_N3U2_500sim_channelmatrix_rb2000m_NONOVERLAP_inner_outer_randomusers.mat

simulate=250;
Nonconvergence=zeros(1,simulate);
Rmin_vector= 1;%:5:21;
% Rmin_vector1= [20 10];
% Rmin_vector2= [20 10];

% for Pt_max=20 W
% epsilon= 100;
% unf_delta= 50;

Pt_max_vector= [0.1 1 3 5 7.5 10 15 20 30 40 50 100]*1e3; % General use
% Pt_max_vector= [0.1 1 5 10 15 20 50 100]*1e3; % Small set
% Pt_max_vector= [0.1 1 3 5]*1e3; %for diagnosis
% Pt_max_vector= 20e3;

% epsilon= 30e-3;
% epsilon= 0.0059;
epsilon= 15e-3;
iter_max= 25;

Used= zeros(K,S); %Initialize the used matrix
P_MAT= zeros(K,U(1),iter_max);
R_MAT= zeros(K,U(1),iter_max);

for bbb=1:size(Pt_max_vector,2)

Pt_max= Pt_max_vector(bbb)*ones(K,1);  %mW
Pt_max= Pt_max/12; % 1 RB= 12 sub-channels
Antenna_Gain= 15 ; %dBi
Pt_max= Pt_max*10^(Antenna_Gain/10);

for aaa= 1:simulate

H_dB_temp(:,:,:,:)= H_dB(aaa,:,:,:,:);

H= 10.^(H_dB_temp./10);
p= mean(Pt_max)/S*ones(K,S); % initial powers
rate= 0.1*ones(K,S); % initial rate for each SC in each cell
MUD_way= {}; 

% Rmin= [Rmin_vector*ones(1,max(U));   % min rate cons for each user in each cell
%        Rmin_vector*ones(1,max(U));];
    
Rmin= ones(3,2);
% Rmin= [2 1;
%        2 1;
%        2 1;];

Used(1,:)= zeros(1,S); % the matrix of the usage state of the subchannels in each cell
                  % the value of (k,s) is the user index to which s is assigned in cell k    
       
pow_convergence= 0;

% pow_vector= zeros(2,max(U),S);
pow_count=1;
%user_rate_matrix= zeros(2,max(U));

p_uniform1= Pt_max(1)/S; % assume equal power distribution for SC assignment purpose
p_uniform2= Pt_max(2)/S; % assume equal power distribution for SC assignment purpose
p_uniform3= Pt_max(3)/S; % assume equal power distribution for SC assignment purpose

%Initial Power Distribution
p(1,:)=p_uniform1;
p(2,:)=p_uniform2;
p(3,:)=p_uniform3;

%% power convergence for that run
while pow_convergence== 0 && pow_count< iter_max

p_temp= p;

%% Calculate SCA and PA+MUD for every BS iteratively

for BS=1:K

%% Heuristic SCA w/ OPT PA
% [Used(BS,:),rate(BS,:)] = SCA(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]=SubgradProj(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% MUD-SCA w/ OPT PA
[Used(BS,:),rate(BS,:)] = SCA_wrtRATE(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
[p(BS,:),rate(BS,:)]=SubgradProj(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% NoINTERFERENCE case
% [Used(BS,:),rate(BS,:)] = SCA_NoINT(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [Used(BS,:),rate(BS,:)] = SCA_wrtRATE_NoINT(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]=SubgradProj_NoINT(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% Shen's performance
% [Used(BS,:),rate(BS,:)] = SCA_Shen(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]=SubgradProj_Shen(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% Shen +WF performance
% [Used(BS,:),rate(BS,:)] = SCA_Shen(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]=SubgradProj_SD(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% Case 2 performance
% [Used(BS,:),rate(BS,:)] = SCA_case2(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]=SubgradProj_case2(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);

%% Heuristic SCA+PA
% [Used(BS,:),rate(BS,:)] = SCA(BS,H,p,rate,Pt_max(BS),Z,Rmin(BS,:),U,S);
% [p(BS,:),rate(BS,:)]= PAMarMin(BS,H,p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S);
end

%% Determine Convergence of the System
if norm(p_temp-p)/sum(sum(p(:,:)))< epsilon %|| pow_count > 50
   pow_convergence= 1;
end
% % another convergence criteria
% if norm(max(p'-p_temp')./sum(p'))< 0.01
%    pow_convergence= 1;
% end 

% P_MAT(:,:,)= zeros;
% R_MAT(pow_count,:,:)= zeros;
for i=1:K
   for s=1:S
       P_MAT(i,Used(i,s),pow_count)= P_MAT(i,Used(i,s),pow_count) + p(i,s);
       R_MAT(i,Used(i,s),pow_count)= R_MAT(i,Used(i,s),pow_count) + rate(i,s);
   end
end

pow_count= pow_count+1;
fark1(pow_count-1)= norm(p_temp(1,:)-p(1,:));
fark2(pow_count-1)= norm(p_temp(2,:)-p(2,:));


% % norm(p_temp-p)/sum(sum(p(:,:)))
end

R_MAT= 12*Subcarrier_bw.*R_MAT;

if pow_count> 28
   Nonconvergence(bbb,aaa)=1;
end

final_user_rate_matrix= zeros(K,max(U));
% rate= 12*Subcarrier_bw.*rate; % 1 RB= 12 subcarriers of 15 kHz bw.

for i=1:K
   for s=1:S
       final_user_rate_matrix(i,Used(i,s))= final_user_rate_matrix(i,Used(i,s)) + rate(i,s); 
   end    
end

final_user_rate_matrix= 12*Subcarrier_bw.*final_user_rate_matrix;

BIG_user_rate_matrix(bbb,aaa,:,:)= final_user_rate_matrix;

aaa=aaa+1;

[bbb,aaa]
end
%%
bbb=bbb+1;


end

final_user_rate_matrix
toc
% p
% rate
% 
% Used
% sum_rate
% count
% MUD_way
% MaxMin_Capacity= [Capacity(H(1,1,s)*p(1,s)/(1+ H(2,1,s)*p(2,s))) Capacity(H(1,1,s)*p(1,s));
%                  Capacity(H(2,2,s)*p(2,s)/(1+ H(1,2,s)*p(1,s))) Capacity(H(2,2,s)*p(2,s));]