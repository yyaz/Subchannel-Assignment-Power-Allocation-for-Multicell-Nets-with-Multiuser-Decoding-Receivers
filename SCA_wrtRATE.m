function [Used,rate1] = SCA_wrtRATE(BS,H,p,rate,Pt_max,Z,Rmin,U,S)

% First SC assignment is based on the min (SNR/INR) ratio
% SCA= SCA(BS, H(:,BS,Used(1,s),s), p(:,s),rate(:,s),Z,Rmin(BS,:), U, S)
%          1             3x1x1x1     3x1      3x1    1    1xU     Nx1 1


assignment1=zeros(U(BS),S);  % temp assignment matrix
user_rate_SC_matrix= zeros(U(BS),S); % temp rate for each SC matrix

Used(1,:)= zeros(1,S);       % final assignment vector    
user_rate_matrix= zeros(1,U(BS)); % temp user rate matrix

p(BS,:)= Pt_max/S; % Assume uniform power allocation

if BS==1;
   otherBS1=2;
   otherBS2=3;
elseif BS==2;
   otherBS1=1;
   otherBS2=3; 
elseif BS==3;
   otherBS1=1;
   otherBS2=2;
end

% SxU complexity
for s=1:S
    for ee=1:U(BS)
        user_rate_SC_matrix(ee,s)= calcrateSC_rateregion(BS,H(:,BS,ee,s),p(:,s),rate(:,s),Z);
    end
end

assignment1_convergence= 0;
while assignment1_convergence== 0

prop_rate= user_rate_matrix./Rmin;
min_prop_rate_user= find(prop_rate== min(prop_rate));
% a measure if the proportional rates are same, i.e. there are more than 1 minimum. randomly selects the user
% min_prop_rate_user=  min_prop_rate_user(floor(rand*size(min_prop_rate_user,2))+1);    
min_prop_rate_user=  min_prop_rate_user(1); 

strongest_channels1= min(find(user_rate_SC_matrix(min_prop_rate_user,:)==max(user_rate_SC_matrix(min_prop_rate_user,:)))); % find best rate subchannel of the min rate user
assignment1(min_prop_rate_user,strongest_channels1)= 1;   % update the assignment matrix
Used(1,strongest_channels1)= min_prop_rate_user; % update the SC usage matrix
rate(BS,strongest_channels1) = user_rate_SC_matrix(min_prop_rate_user,strongest_channels1);
user_rate_matrix(1,min_prop_rate_user)= user_rate_matrix(1,min_prop_rate_user)+ rate(BS,strongest_channels1);

user_rate_SC_matrix(:,strongest_channels1)= zeros; % delete the assigned subchannel from possible assignment list
%H11_temp(1,1,:,strongest_channels1)= zeros; % delete the assigned subchannels from possible assignment list

%%
Assigned_SC1= find(Used(1,:)> 0);  % update the # of assigned SCs

if length(Assigned_SC1)== S
   assignment1_convergence= 1;    
end

end % assignment1 convergence

rate1=rate(BS,:);