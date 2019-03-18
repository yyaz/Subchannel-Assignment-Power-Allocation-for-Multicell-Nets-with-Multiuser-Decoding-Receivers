function [Used,rate1] = SCA_Shen(BS,H,p,rate,Pt_max,Z,Rmin,U,S)

% SCA= SCA(BS, H(:,BS,Used(1,s),s), p(:,s),rate(:,s),Z,Rmin(BS,:), U, S)
%          1             3x1x1x1     3x1      3x1   1  NxU Nx1 1

assignment1=zeros(U(BS),S);
Used(1,:)= zeros(1,S);
user_rate_matrix= zeros(1,U(BS));

% sub-carrier allocation
%% assignment of sub-channels proportional to rate/Rmin of each user

% % assignment1_temp= assignment1;
% % rate_temp(1,:)= rate(1,:);

p(BS,:)= Pt_max/S; % Assume uniform power allocation

% assign for each user the strongest channel of that user, first
% RP_temp1(1,1,:,:)= (2.^[rate(2,:); rate(2,:)]-1);
% H11_temp= H(1,1,:,:)+ H(1,1,:,:)./(H(2,1,:,:));
% H11_temp= H(1,1,:,:)./(H(2,1,:,:));

% % if BS==1;
% %    otherBS1=2;
% %    otherBS2=3;
% % elseif BS==2;
% %    otherBS1=1;
% %    otherBS2=3; 
% % elseif BS==3;
% %    otherBS1=1;
% %    otherBS2=2;
% % end

H11_temp= H(BS,BS,:,:);

for dd=1:U(BS)
    strongest_channels1= find(H11_temp(1,1,dd,:)== max(H11_temp(1,1,dd,:))); % find strongest channel
    assignment1(dd,strongest_channels1)= 1;   % update the assignment matrix
    Used(1,strongest_channels1)= dd; % update the SC usage matrix
    H11_temp(1,1,:,strongest_channels1)= zeros; % delete the assigned subchannels from possible assignment list
end

Assigned_SC1= find(Used(1,:)> 0);

for s= Assigned_SC1 % for already assigned channels, calculate the respective rate

rate(BS,s) = calcrateSC_rateregion_Shen(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z);
user_rate_matrix(Used(1,s))= user_rate_matrix(Used(1,s))+ rate(BS,s);
end

assignment1_convergence= 0;
Not_Assigned_SC1= size(setdiff(1:S, Assigned_SC1),2);

while assignment1_convergence== 0

prop_rate= user_rate_matrix./Rmin;
min_prop_rate_user= find(prop_rate== min(prop_rate));
% a measure if the proportional rates are same, i.e. there are more than 1 minimum. randomly selects the user
% min_prop_rate_user=  min_prop_rate_user(floor(rand*size(min_prop_rate_user,2))+1);    
min_prop_rate_user=  min_prop_rate_user(1); 

strongest_channels1= find(H11_temp(1,1,min_prop_rate_user,:)==max(H11_temp(1,1,min_prop_rate_user,:))); % find strongest channel of the min rate user
assignment1(min_prop_rate_user,strongest_channels1)= 1;   % update the assignment matrix
Used(1,strongest_channels1)= min_prop_rate_user; % update the SC usage matrix
H11_temp(1,1,:,strongest_channels1)= zeros; % delete the assigned subchannels from possible assignment list

% % HH(:,1,strongest_channels1)= H(:,1,min_prop_rate_user,strongest_channels1);
% % pow_vector(1,min_prop_rate_user,strongest_channels1)=p(1,strongest_channels1);

% rate(1,strongest_channels1)= Capacity(HH(1,1,strongest_channels1)*p(1,strongest_channels1)/(Z+HH(2,1,strongest_channels1)*p(2,strongest_channels1)));
%% new trial
s= strongest_channels1;

rate(BS,s) = calcrateSC_rateregion_Shen(BS,H(:,BS,Used(1,s),s),p(:,s),rate(:,s),Z);
%%
user_rate_matrix(1,min_prop_rate_user)= user_rate_matrix(1,min_prop_rate_user)+ rate(BS,s);
                                                                                %eski-yanlis hali rate(1,s)
Assigned_SC1= find(Used(1,:)> 0);  % update the # of assigned SCs

if length(Assigned_SC1)== S
   assignment1_convergence= 1;    
end

end % assignment1 convergence

rate1=rate(BS,:);