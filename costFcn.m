function[fitness] = costFcn(population,pars)
% tic
[npop,nvar] = size(population);
gamma_C       = pars.gX;                 % Xnuc gyromagnetic ratio
gamma_H       = pars.gH;                 % 1H gyromagnetic ratio
tau = pars.tau;
alpha = pars.alpha;
beta  = pars.beta;
minB0 = pars.minB0;
maxB0 = pars.maxB0;
H0_perB0 = pars.H0_perB0;
H_JHs = pars.H_JHs;
H_JCH = pars.H_JCH;
rho_init = pars.rho_init;
coil = pars.Iz;
traceNorm = pars.traceNorm;
%% fitness evaluation
fitness  = zeros(npop,1);
parfor i = 1:npop
    parameters = (population(i,1:nvar));
    parameters(1) = pars.minB0;
    parameters(end) = pars.maxB0;
    parameters = min(pars.maxB0,max(pars.minB0,parameters));
    parameters = 1e-6*smoothCurveToConstraint(1e6*parameters',tau,alpha,beta,1e6*minB0,1e6*maxB0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 4: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Remagnetization from 50nT to 2uT %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho = rho_init;
    for j = 1:nvar
        B0 = parameters(j);
        %%% 1H-Xnuc pulse %%%
        H = B0*H0_perB0+H_JHs+H_JCH;                                                          %coupling
        PH  = expm(-1i*H*tau);
        PHt = PH';
        rho = PH*rho*PHt;
    end
    signal_tmp  = trace(coil'*rho)/traceNorm;
    fitness(i)	= (-abs(real(signal_tmp)));
end
% toc
end















