function [rate] = ComputeRateThomasType1user(P,Lambda,m,sigma)
% This function computes average ergodic rate of two tier HetNet described in for Type-1
% users when Tier 1 is modelled as a Thomas cluster process.
% Input ---
% P: 1x2 vector specifying the power levels
% Lambda: 1x2 vector specifying BS intensities
% m: average number of BSs of the TCP
% sigma: cluster standard deviation of TCP
% Output: average ergodic rate
P_1  = P(1);
P_2 =  P(2);


l_p_2 = Lambda(2);
l_p_1 = Lambda(1);
        
alpha = 4;

Noise = 0;

f_d1 = @(x,r)    ricepdf(x,r,sigma);

%% Compute laplace transform of power %%
 
 %%% Tier 1 term %%%
  
  %%%% P-s %%%%
  P1_P_11 = 1;
  P1_P_21 = (P_2/P_1)^(1/alpha);
  
  %%%% D-s %%%% (For Sum-Product)
  P1_D1 = @(r,z,s) exp(-m*(1+(integral(@(y)f_d1(y,z),0,P1_P_11*r,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)) ...
      -(integral(@(y)f_d1(y,z)./(1+P_1.*s./(y.^alpha)),0,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));
  P1_D_arr1 = @(r,z,s) arrayfun(@(r,z,s) P1_D1(r,z,s),r,z,s);
  
  %%%% T-s %%%% (For Sum-Product)
  P1_T1 = @(r,s) 2*pi*l_p_1.*integral(@(z) f_d1(r,z).*P1_D_arr1(r,z,s).*z,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
  
  %%%% Q-s: combine all PGFLs %%%%
  P1_Q1 = @(r,s) exp(-2*pi*l_p_1*integral(@(z) (1-P1_D_arr1(r,z,s)).*z,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
  P1_Q2 = @(r,s) exp(-pi*l_p_2*(P1_P_21^2.*r.^2+(P_2.*s).^(2/alpha).*((2/alpha).*beta(2/alpha,1-2/alpha))));
  P1_f1 = @(r,s) P1_T1(r,s).*P1_Q1(r,s).*P1_Q2(r,s);
  P1_f1_arr = @(r,s) arrayfun(@(r,s) P1_f1(r,s),r,s);
  N_1 = @(s) m*integral(@(r) P1_f1_arr(r,s),0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
 
 %%% Tier 2 term %%%
  
  %%%% P-s %%%%
  P2_P_12 = (P_1/P_2)^(1/alpha);
  P2_P_22 = 1;
  
  %%%% D-s %%%% (For Sum-Product)
  P2_D1 = @(r,z,s) exp(-m*(1+(integral(@(y)f_d1(y,z),0,P2_P_12*r,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)) ...
      -(integral(@(y)f_d1(y,z)./(1+P_1.*s./(y.^alpha)),0,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));  
  P2_D_arr1 = @(r,z,s) arrayfun(@(r,z,s) P2_D1(r,z,s),r,z,s);
  
  %%%% Q-s: combine all PGFLs %%%%
  P2_Q1 = @(r,s) exp(-2*pi*l_p_1*integral(@(z) (1-P2_D_arr1(r,z,s)).*z,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
  P2_Q2 = @(r,s) 2*pi*l_p_2.*r.*exp(-pi*l_p_2*(P2_P_22^2.*r.^2+(P_2.*s).^(2/alpha).*((2/alpha).*beta(2/alpha,1-2/alpha))));
  P2_f2 = @(r,s) P2_Q1(r,s).*P2_Q2(r,s);
  P2_f2_arr = @(r,s) arrayfun(@(r,s) P2_f2(r,s),r,s);
  N_2 = @(s) integral(@(r) P2_f2_arr(r,s),0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
  
%% Compute laplace transform of interference %% 

 %%% Tier 1 term %%%

  %%%% P-s %%%%
  I1_P_11 = 1;
  I1_P_21 = (P_2/P_1)^(1/alpha);
  
  %%%% C-s %%%% (For Sum-Product)
  I1_C1= @(r,s,z) exp(-m*(1-(integral(@(y)f_d1(y,z)./(1+s.*P_1./(y.^alpha)),I1_P_11*r,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));
  I1_C_arr1 = @(r,s,z) arrayfun(@(r,s,z) I1_C1(r,s,z),r,s,z);
  
  %%%% T-s %%%%
  I1_T1=@(r,s) 2*pi*l_p_1.*integral(@(z)f_d1(r,z).*I1_C_arr1(r,s,z).*z,0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
   
  %%%% M-s: combine all PGFLs %%%%
  I1_M1= @(r,s) exp(-2*pi*l_p_1*integral(@(z)(1-I1_C_arr1(r,s,z)).*z,0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
  I1_M2=@(r,s) exp(-pi*l_p_2*I1_P_21.^2.*r.^2.*(1+(2.*s.*P_1./((alpha-2).*r.^alpha)).*hypergeom([1,1-2/alpha],2-2/alpha,-s.*P_1./(r.^alpha))));
  I1_f1 = @(r,s) I1_T1(r,s).*I1_M1(r,s).*I1_M2(r,s);
  I1_f1_arr=@(r,s) arrayfun(@(r,s)I1_f1(r,s),r,s);
  J_1= @(s) m.*integral(@(r) I1_f1_arr(r,s),0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3); 
   
 %%% Tier 2 term %%%

  %%%% P-s %%%%
  I2_P_12 = (P_1/P_2)^(1/alpha);
  I2_P_22 = 1;
  
  %%%% C-s %%%% (For Sum-Product)
  I2_C1= @(r,s,z) exp(-m*(1-(integral(@(y)f_d1(y,z)./(1+s.*P_1./(y.^alpha)),I2_P_12*r,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));
  I2_C_arr1 = @(r,s,z) arrayfun(@(r,s,z) I2_C1(r,s,z),r,s,z);
 
  %%%% M-s: combine all PGFLs %%%%
  I2_M1= @(r,s) exp(-2*pi*l_p_1*integral(@(z)(1-I2_C_arr1(r,s,z)).*z,0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
  I2_M2=@(r,s) 2*pi*l_p_2.*r.*exp(-pi*l_p_2*I2_P_22.^2.*r.^2.*(1+(2.*s.*P_2./((alpha-2).*r.^alpha)).*hypergeom([1,1-2/alpha],2-2/alpha,-s.*P_2./(r.^alpha))));
  I2_f2 = @(r,s) I2_M1(r,s).*I2_M2(r,s);
  I2_f2_arr=@(r,s) arrayfun(@(r,s)I2_f2(r,s),r,s);
  J_2= @(s) integral(@(r) I2_f2_arr(r,s),0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3); 
    
%% Compute Average Ergodic Rate %%
J_arr1= @(s) arrayfun(@(s)J_1(s),s);
J_arr2= @(s) arrayfun(@(s)J_2(s),s);
N_arr1= @(s) arrayfun(@(s)N_1(s),s);
N_arr2= @(s) arrayfun(@(s)N_2(s),s);
rate = integral(@(s) ((exp(-s.*Noise))./s).*(J_arr1(s)+J_arr2(s)-N_arr1(s)-N_arr2(s)),0,20*sigma, 'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
fprintf('%f',rate);