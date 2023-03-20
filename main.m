% Genetic algorithm optimized pulse design
% (2019) Vencel Somai -> vs460@cam.ac.uk

clear all
close all
rng default
%% parameters
gX    = 10.71e6;    % Xnuc gyromagnetic ratio [Hz/T]
gH    = 42.57e6;    % 1H gyromagnetic ratio   [Hz/T]
T     = 2;          % sequence time length [s]
minB0 = 50e-9;
maxB0 = 2e-6;
pars.prepPhase = true;
%% plot options
plotOpt    = {@gaplotbestf};

%% population options
nvar       = 100;
pop.Type   = 'doubleVector';
pop.Size   = 32;
pop.Create = @gacreationuniform;
pop.Range = [minB0;maxB0];
population = rand(pop.Size,nvar).*repmat(pop.Range(2)-pop.Range(1),[pop.Size,1]) + repmat(pop.Range(1),[pop.Size,1]);

%% selection options
select.Type = @selectionroulette;
ec          = 2;
crossFrac   = 0.8;
%% mutation options
mutFcn      = @mutationgaussian;
Scale       = 0.5;
Shrink      = 0.5;
%% crossover options
crssover    = 'crossoverheuristic';
%% number of iterations
NofGens     = 3e4;

%% create options input
if size(population,2) < nvar
    population = [population,zeros(size(population,1),nvar - size(population,2))];
end
population(1,:) = linspace(minB0,maxB0,nvar);
pars.gX  = gX; pars.gH = gH; pars.tau = T/nvar; 
pars.minB0 = minB0; pars.maxB0 = maxB0; 
%%
[pars.rho_init,pars.H0_perB0,pars.H_JHs,pars.H_JCH,pars.Iz,pars.traceNorm] = calcOps(pars);

%%
options = optimoptions('ga','InitialPopulationMatrix',population,'PlotFcn',plotOpt,'PopulationType',pop.Type,'PopulationSize',pop.Size,...
    'CreationFcn',pop.Create,'InitialPopulationRange',pop.Range,'SelectionFcn',{select.Type},'EliteCount',ec,...
    'CrossoverFraction',crossFrac,'MutationFcn',{mutFcn,Scale,Shrink},'CrossoverFcn',{crssover},'MaxGenerations',NofGens,...
    'UseVectorized',true,'MaxStallGenerations',NofGens);
B0s = ga(@(x) costFcn(x,pars),nvar,[],[],[],[],[],[],[],options);
%%
figure
B0s = min(pars.maxB0,max(pars.minB0,B0s));
t = linspace(0,T,nvar);
subplot(1,2,1); plot(t,B0s*1e6)
xlabel('Time [ms]'); ylabel('B_0 [uT]')
rho = pars.rho_init;
for j = 1:nvar
    B0 = min(pars.maxB0,max(pars.minB0,B0s(j)));
    H = B0*pars.H0_perB0+pars.H_JHs+pars.H_JCH;                                                          %coupling
    PH  = expm(-1i*H*pars.tau);
    PHt = PH';
    rho = PH*rho*PHt;
    signal_tmp(j)  = trace(pars.Iz'*rho)/pars.traceNorm;
end
subplot(1,2,2); plot(t,real(signal_tmp))
xlabel('Time [ms]'); ylabel('I_z operator amplitude')