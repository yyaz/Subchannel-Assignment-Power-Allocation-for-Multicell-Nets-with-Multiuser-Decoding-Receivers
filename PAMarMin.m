function [p1,rate1] = PAMarMin(BS,H,p,rate,Pt_max,Pt_max_vector,Used,Z,Rmin,S)

% = PAMarMin(BS,H(:,BS,:,:),p,rate,Pt_max(BS),Pt_max_vector(bbb),Used(BS,:),Z,Rmin(BS,:),S)
%                1   3x1xUx1       3xS 3xS       1            1               1xS      1 NxU  1

% Bisection Method Initialization for Power Allocation
delta= 1e-5;
Pt_epsilon= 30*Pt_max_vector*1e-4;
% unf_delta= 10+ (Pt_max_vector-100)*1.0040e-004;
unf_delta= 0.2+ (Pt_max_vector-100)*1.9e-004;

Unfairness1= 100;

Pt_max11= Pt_max/2;
Pt_max12= Pt_max/2;

while abs(Unfairness1) > unf_delta

lambda_opt11_max= 1;
lambda_opt11_min= 0;
lambda_opt12_max= 1;
lambda_opt12_min= 0;

while lambda_opt11_max - lambda_opt11_min> delta
 
lambda_opt11= (lambda_opt11_max + lambda_opt11_min)/2;
 
for s= find(Used(1,:)==1)
    [p(BS,s),rate(BS,s)] = calc_power_rateSC(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z,1,lambda_opt11);
end
 
if sum(p(BS,find(Used(1,:)==1)))- Pt_max11 < 0;
    lambda_opt11_max= lambda_opt11; 
else
   lambda_opt11_min= lambda_opt11; 
end
 
end
% End of Calculate SCA and PA for BS1- u11
%% PA for BS1- u12
while lambda_opt12_max - lambda_opt12_min> delta
 
lambda_opt12= (lambda_opt12_max + lambda_opt12_min)/2;

for s= find(Used(1,:)==2)
    [p(BS,s),rate(BS,s)] = calc_power_rateSC(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z,1,lambda_opt12);
end

if sum(p(BS,find(Used(1,:)==2)))- Pt_max12 < 0;
   lambda_opt12_max= lambda_opt12; 
else
   lambda_opt12_min= lambda_opt12; 
end
 
end

check_user_rate_matrix1= zeros(size(Rmin));
for s= 1:S
    check_user_rate_matrix1(1,Used(1,s))= check_user_rate_matrix1(1,Used(1,s)) + rate(BS,s); 
end

Unfairness1= (check_user_rate_matrix1(1,1)./Rmin(1,1))-(check_user_rate_matrix1(1,2)./Rmin(1,2)); %% max(user_rate_matrix(1,:)./Rmin(1,:))-min(user_rate_matrix(1,:)./Rmin(1,:)); for U>2 users

if Unfairness1 > 0
   Pt_max11= Pt_max11 - Pt_epsilon;
   Pt_max12= Pt_max12 + Pt_epsilon;
else
   Pt_max11= Pt_max11 + Pt_epsilon;
   Pt_max12= Pt_max12 - Pt_epsilon; 
end

end
%%
p1=p(BS,:);
rate1= rate(BS,:);