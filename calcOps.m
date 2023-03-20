function[rho_init,H0_perB0,H_JHs,H_JCH,Iz,traceNorm] = calcOps(pars)

gamma_C       = pars.gX;                 % Xnuc gyromagnetic ratio
gamma_H       = pars.gH;                 % 1H gyromagnetic ratio
tau = pars.tau;
minB0 = pars.minB0;
maxB0 = pars.maxB0;
B0_earth = 50e-6;
% Scalar-coupling values
JH1H2  = 10.5;
JH1H3  = 17.26;
JH2H3  = -1.58;
JH1H45 = 5.56;
JH2H45 = 1.27;
JH3H45 = 1.61;
JH45C  = 3.15;

%%%
n_of_I_spins = 1;
n_of_S_spins = 5;
nspins = n_of_I_spins + n_of_S_spins;

% Define Pauli matrices
sigma_x = ([0 1/2; 1/2 0]);
sigma_y = ([0 -1i/2; 1i/2 0]);
sigma_z = ([1/2 0; 0 -1/2]);
unit    = ([1 0; 0 1]);

% Calculate spin operators
L = zeros(2^nspins,2^nspins,nspins,4);
for l = 1:nspins
    Lx_current = 1;Ly_current = 1;Lz_current = 1;
    for s = 1:nspins
        if s == l
            Lx_current = kron(Lx_current,sigma_x);
            Ly_current = kron(Ly_current,sigma_y);
            Lz_current = kron(Lz_current,sigma_z);
        else
            Lx_current = kron(Lx_current,unit);
            Ly_current = kron(Ly_current,unit);
            Lz_current = kron(Lz_current,unit);
        end
    end
    L(:,:,l,1) = eye(2^nspins); L(:,:,l,2) = Lx_current; L(:,:,l,3) = Ly_current; L(:,:,l,4) = Lz_current;
end

% Carbon spin operators
Ix = sum(L(:,:,1,2),3);
Iy = sum(L(:,:,1,3),3);
Iz = sum(L(:,:,1,4),3);
I = [Ix;Iy;Iz];
% H1 spin operators
Sx1 = sum(L(:,:,2,2),3);
Sy1 = sum(L(:,:,2,3),3);
Sz1 = sum(L(:,:,2,4),3);
S1 = [Sx1;Sy1;Sz1];
% H2 spin operators
Sx2 = sum(L(:,:,3,2),3);
Sy2 = sum(L(:,:,3,3),3);
Sz2 = sum(L(:,:,3,4),3);
S2 = [Sx2;Sy2;Sz2];
% H3 spin operators
Sx3 = sum(L(:,:,4,2),3);
Sy3 = sum(L(:,:,4,3),3);
Sz3 = sum(L(:,:,4,4),3);
S3 = [Sx3;Sy3;Sz3];
% H4 spin operators
Sx4 = sum(L(:,:,5,2),3);
Sy4 = sum(L(:,:,5,3),3);
Sz4 = sum(L(:,:,5,4),3);
S4 = [Sx4;Sy4;Sz4];
% H4' spin operators
Sx5 = sum(L(:,:,6,2),3);
Sy5 = sum(L(:,:,6,3),3);
Sz5 = sum(L(:,:,6,4),3);
S5 = [Sx5;Sy5;Sz5];
% sum H operators
Sx = sum(L(:,:,2:6,2),3);
Sy = sum(L(:,:,2:6,3),3);
Sz = sum(L(:,:,2:6,4),3);
S = [Sx;Sy;Sz];

% construct cartesian product operator basis
% basis elements are represented as a 1x4 row vector:
% entry 1 = E, 2 = Lx, 3 = Ly, 4 = Lz
% and the basis operators are the product of the permutations (with repetitions) of these
for bb = 1:nspins
    permutation_mtx_tmp(:,bb) = reshape(repmat(mod((1:4^(nspins-bb+1))-1,4)'+1,[1,4^(bb-1)])',[1,4^nspins])';
end
permutation_mtx = permutation_mtx_tmp;

% constructing cell with the basis operators and the relaxation supermatrix
B             = cell((2^nspins)^2,1);
R             = zeros((2^nspins)^2,1);
RelMtx        = zeros(size(Ix));
RelMtx_struct = zeros(size(Ix));

% define Hamiltonians
H_JHs   = 2*pi*(JH1H2*S1'*S2+...
              JH1H3*S1'*S3+...
              JH2H3*S2'*S3+...
              JH1H45*S1'*(S4+S5)+...
              JH2H45*S2'*(S4+S5)+...
              JH3H45*S3'*(S4+S5));
H_JCH   = 2*pi*(JH45C*I'*(S4+S5));
H_earth = 2*pi*(-gamma_H*B0_earth*Sz-gamma_C*B0_earth*Iz) + H_JHs + JH45C*Iz'*(Sz4+Sz5);
H_zf    = 2*pi*(-gamma_H*minB0*Sz-gamma_C*minB0*Iz) + H_JHs + H_JCH;
H0_perB0 = 2*pi*(-gamma_H*Sz-gamma_C*Iz);
%% initialization
rho_init  = S1'*S2;
traceNorm = trace(rho_init'*rho_init);
coil      = Iz;
if pars.prepPhase
    rho = rho_init;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% average density matrix after 1s time evolution at earth’s field %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = 1000;rho_avg = rho/N; tau_earth = 1/N;
    for n = 2:N
        PH_earth  = expm(-1i*H_earth*tau_earth);
        PHt_earth = PH_earth';
        rho = PH_earth*rho*PHt_earth;
        rho_avg = rho_avg + rho/N;
    end
    rho = rho_avg;%/trace(rho_avg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% diabatic transport from earth’s field to zero field %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nH=100; tau_dia = 1e-3;
    for n=1:19
        nH=nH-10; % H Larmor frequency: the speed of this passage is changed by varying
        % nH (-1, -5, -10…) and the delay for each step
        nH1=nH;
        nX=nH/3.9772;
        H_dia  = 2*pi*(-nH*Sz-nX*Iz) + H_JHs + JH45C*Iz'*(Sz4+Sz5);
        U_dia = expm(-1i*2*pi*H_dia*tau_dia); U_diat = expm(1i*2*pi*H_dia*tau_dia);
        for m=1:5 % time delay for each passage: 1ms, it can be increased
            rho= U_dia*rho*U_diat;
            T(m,n)=(real(trace(Iz'*rho))); % 13C net magnetization
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 3: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% evolution at zero field 1s %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho_avg = rho/N; tau_zf = 1/N;
    U_zf = expm(-1i*2*pi*H_zf*tau_zf); U_zft =expm(1i*2*pi*H_zf*tau_zf);
    for m=2:N
        rho = U_zf*rho*U_zft;
        rho_avg = rho_avg+rho/N;
    end
    rho_init = rho_avg;%/trace(rho_avg);
end

end