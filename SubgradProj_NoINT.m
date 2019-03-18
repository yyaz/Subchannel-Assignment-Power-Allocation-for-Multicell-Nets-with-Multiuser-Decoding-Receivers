function [p1,rate1] = SubgradProj_NoINT(BS,H,p,rate,Pt_max,Pt_max_vector,Used,Z,Rmin,S)

% = SubgradProj(BS,H(:,BS,:,:),p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S)
%                1   3x1xUx1       3xS 3xS       1            1               1xS      1 NxU  1

% Subgradient Projection Method Initialization for Power Allocation
subgrad_say1_limit= 2000;
subgrad_say1= 1;
subgrad_convergence1= 0;
Lagrangiandual1= zeros;

subgrad_delta= 5e-6;
step_size= 1e-3;
if Pt_max_vector<=1e3;
  step_size_lambda= step_size*2e0/Pt_max_vector;
else
  step_size_lambda= step_size*2e-2/Pt_max_vector;
end

x1= [3.5/5*20/Pt_max_vector 1 1]'; 
% x1_preprevious= rand(3,1);
x1_previous= rand(3,1);

X1MAT= zeros(3,1);

while subgrad_convergence1==0 && subgrad_say1 < subgrad_say1_limit

lambda_opt11= x1(1);
% q_11= x1(2); 
% q_12= x1(3); 

check_user_rate_matrix1= zeros(1,2);
for s= 1:S
    q_BS1=x1(Used(1,s)+1);
    [p(BS,s),rate(BS,s)] = calc_power_rateSC_NoINT(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z,q_BS1,lambda_opt11);
    check_user_rate_matrix1(1,Used(1,s))= check_user_rate_matrix1(1,Used(1,s)) + rate(BS,s); 
end

%% keep track of the best solution found
Lagrangiandual1(subgrad_say1)= [Pt_max-sum(p(BS,:)) check_user_rate_matrix1]*x1;

d_l1= Pt_max-sum(p(BS,:));    % Subgradient
d_q1= (check_user_rate_matrix1(1,:))';  % Subgradient

A1= [Rmin];
b1= 1;

x1_preprevious= x1_previous;
x1_previous= x1;
x1(1)= x1(1)- step_size_lambda*d_l1;
x1(2:3)= x1(2:3)- step_size*(eye(2)-A1'*1/(A1*A1')*A1)*d_q1; % - A1'*1/(A1*A1')*b1

if norm(x1(2:3)- x1_previous(2:3))/sum(sum(x1(2:3,:)))< subgrad_delta || norm(x1- x1_preprevious)/sum(sum(x1(:,:))) < subgrad_delta 
   subgrad_convergence1= 1;
end

subgrad_say1= subgrad_say1+ 1;

X1MAT(:, subgrad_say1)=x1; 
end

lambda_opt11= x1(1);

% Calculate rate of BS1 users with opt parameters
for s= 1:S
    q_BS1=x1(Used(1,s)+1);
    [p(BS,s),rate(BS,s)] = calc_power_rateSC_NoINT(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z,q_BS1,lambda_opt11);
end

p1=p(BS,:);
rate1= rate(BS,:);