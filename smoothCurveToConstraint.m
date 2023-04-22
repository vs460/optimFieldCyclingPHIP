function [s_sm] = smoothCurveToConstraint(s0,dt,alpha,beta,minB0,maxB0)
N = length(s0);
C_kine=set_MRI_constraints_RIV(alpha,beta,dt); 

C_linear.n_linear_constraints = 2;
C_linear.A = zeros(2,N);
C_linear.A(1,N) = 1;
C_linear.A(2,1) = 1;
C_linear.AT = C_linear.A';
C_linear.v = [maxB0;minB0];
C_linear.PI = C_linear.AT;

% Algorithm parameters
Algo_param.nit = 30;  % number of iterations
Algo_param.L=36;    % Lipschitz constant of the gradient
Algo_param.discretization_step=dt;

% optional parameters (default 0)
Algo_param.show_progression = 0; % 0 = no progression , can be really slow
Algo_param.display_results  = 0;

%% Project curve with Rotation-Invariant Constraints
s_sm=Project_Curve_Affine_Constraints(s0,C_kine,C_linear,Algo_param);


end

