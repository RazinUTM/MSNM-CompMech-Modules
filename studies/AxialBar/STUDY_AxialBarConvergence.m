% Title: Code for an Axial Bar by 1D Finite Elements: the errors against the mesh
% Attempt: to solve a uniform & a tapered bar by 2- & 3-noded elements on
%          five meshes, & to measure the rates of convergence against the
%          closed form (P2 of the ladder)
% Created on: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's answers of 23/11/2026 on the choices of P2 (2 changed):
%   these are the lines of problems/AxialBar/MAIN_AxialBar as they stood on
%   24/11/2026, on the kernels KM_Line1D, SV_Partition, SF_Line1D &
%   UT_AxialBar, with its two loops, over the bars (Bar, one row [a p]
%   each) & over the meshes (Ne), its data, its guards & its rule for the
%   rates that are not given. They were that main file until this day,
%   with these numbers; the main file is now one bar on one mesh, & as it
%   stands it is the entry [a = 0.5, p = 1, Ne = 8] of this study.
%   The names are those of the main file & not the Loop1, Error1, Rate1 of
%   the other studies, so that an entry reads against the main file:
%     eN, eU, eE, eR   one row per bar, one column per mesh: the errors of
%                      the main file (at the nodes, L2 norm, energy norm; the
%                      reaction with the loads)
%     UL               u at x = L on every mesh, & the exact one (last column)
%     RateN, RateU, RateE   the observed order of eN, eU & eE between one
%                      mesh & the next (RateN is new with this file)
%   What it shows: with a constant EA the elements are exact at the nodes &
%   the 3-noded ones everywhere; on the tapered bar the L2 norm falls at the
%   order p+1, the energy norm at p & the error at the nodes at 2p.
%   With P = q*L*(1/a - 1) the strain of the tapered bar is uniform & the
%   elements are exact: such a bar has no rate to show.
%   The load is uniform: that the nodes of a uniform bar are exact "for any
%   load" (the catalogue) is shown for this load only.
%   differences from the lines of 24/11/2026: the shape functions of the
%   element by SF_FEMLG1D (the lines wrote [s 0] to call SF_FEMLG), the rate
%   at the nodes, & what eR checks, in its comment
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 25/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 25/11/2026 (his answers on the twelve choices of
%   P4: choice 12, the rates of P2 on the kernel):
%   (1) the rates are PL_Rate (core/), with the rule of before (no rate where
%       the error of the coarser mesh is below 1e-10); the same numbers in
%       every bit
%   run in GNU Octave 8.4.0 on 25/11/2026; not yet run in MATLAB
% Double checked on: -

%% Project initialisation
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'core','assembly'),fullfile(projectRoot,'core','basis'), ...
    fullfile(projectRoot,'core','quadrature'),fullfile(projectRoot,'core','solver'), ...
    fullfile(projectRoot,'core','post'),fullfile(projectRoot,'problems','AxialBar'));

%% DATA INPUT

% the bar: x along the bar from the held end (N & m); q & P are positive along
% +x (away from the support) & N = E*A*du/dx is positive in tension
L = 2;      % length (m)
E = 200e9;  % Young's modulus (N/m^2)
A0 = 1e-3;  % area at x = 0 (m^2)
q = 50e3;   % distributed axial load, uniform (N/m)
P = 50e3;   % load at the free end x = L (N)
ubar = 0;   % displacement held at x = 0 (m): 0 as in (P2.4); with another
            % value every u moves by it through the term -K(f,c)*U(c), which
            % TEST_AxialBar holds on the main file, with ubar = 2e-3

% the bars solved, one row each [taper a , order p of the elements]:
% A(x) = A0*(1 - a*x/L) with 0 <= a < 1, a = 0 the uniform bar of the
% catalogue; p = 1 the linear element (2 nodes), p = 2 the quadratic (3 nodes).
% The closed form of the tapered bar loses about eps/a of its figures: a
% between 0 & about 1e-4 is not for the rates of the quadratic elements
Bar = [0   1;
       0   2;
       0.5 1;
       0.5 2];

% number of elements, one solution per entry (each twice the one before; two
% entries at the least: a rate is taken between one mesh & the next)
Ne = [2 4 8 16 32];

ng = 2;  % gauss points of an element for [k] & {r}: exact for these bars (degree 3); fewer than p make [K] singular
nge = 8; % gauss points of an element for the error norms (a higher rule: the exact u is not a polynomial)

% the solution drawn at the end [row of Bar , entry of Ne] (both must exist:
% change it with Bar & Ne)
Draw = [3 1];

% the ranges of the input (the closed form & the elements hold only there)
if any(Bar(:,1) < 0 | Bar(:,1) >= 1)
    error('STUDY_AxialBarConvergence: the taper a must be 0 <= a < 1 (row %d of Bar)',find(Bar(:,1) < 0 | Bar(:,1) >= 1,1))
end
if ng < max(Bar(:,2))
    error('STUDY_AxialBarConvergence: ng = %d is fewer than the order p = %d: [K] would be singular',ng,max(Bar(:,2)))
end
if any(Bar(:,1) > 0 & Bar(:,1) < 1e-4)
    warning('STUDY_AxialBarConvergence: a taper below 1e-4 leaves the closed form about eps/a of its figures: not for the rates')
end


%% SOLUTION OF EVERY BAR ON EVERY MESH
eN = zeros(size(Bar,1),numel(Ne)); % error at the nodes (largest, of the largest u)
eU = zeros(size(Bar,1),numel(Ne)); % error of the displacement over the bar (L2 norm)
eE = zeros(size(Bar,1),numel(Ne)); % error in the energy norm
eR = zeros(size(Bar,1),numel(Ne)); % reaction & loads, of the loads
UL = zeros(size(Bar,1),numel(Ne)+1); % u at x = L on every mesh, & exact
for b = 1:size(Bar,1)
    % loop over the bars
    a = Bar(b,1);
    p = Bar(b,2);

    % AREA & AXIAL STIFFNESS: A(x) = A0*(1 - a*x/L), k = E*A for KM_Line1D;
    % the exact u, N & du/dx are those of UT_AxialBar
    Aex = @(x) A0*(1 - a*x/L);
    EAx = @(x) E*Aex(x);
    qx = @(x) q + 0*x;

    % SHAPE FUNCTIONS OF THE ELEMENT in s from -1 to 1, the xi of the catalogue
    % (the same for every element): its nodes from end to end, its gauss
    % points & weights for [k] & {r}, & those for the error norms
    sn = linspace(-1,1,p+1)';
    [G,W] = GQ_Points(ng);
    [Ns,dNs] = SF_FEMLG1D(G,sn);
    [Ge,We] = GQ_Points(nge);
    [Nse,dNse] = SF_FEMLG1D(Ge,sn);

    for m = 1:numel(Ne)
        % loop over the meshes

        % NODES & ELEMENTS: equal elements, the nodes of each from end to end
        nN = p*Ne(m)+1;
        xi = linspace(0,L,nN)';
        Elem = p*((1:Ne(m))'-1) + (1:p+1);

        % STIFFNESS MATRIX & FORCE VECTOR: the integrals of E*A*dNi/dx*dNj/dx &
        % Ni*q over the elements (KM_Line1D), & the load at x = L
        [K,F] = KM_Line1D(xi,Elem,Ns,dNs,W,EAx,qx);
        F(nN) = F(nN) + P;

        % SUPPORT, SOLUTION & REACTION: the node at x = 0 held at ubar, by
        % partition, K(f,f)*U(f) = F(f) - K(f,c)*U(c) & R = K(c,:)*U - F(c)
        [U,R] = SV_Prescribed(K,F,1,ubar);

        % COMPARE WITH THE EXACT SOLUTION
        % u & du/dx of the elements at the gauss points of the error norms
        % [xQ], with their weights [WxQ] (one column per element)
        [xQ,WxQ,Uq,dUq] = SF_Line1D(xi,Elem,U,Nse,dNse,We);
        [Uex,~,dUex] = UT_AxialBar(xQ,L,E,A0,a,q,P,ubar);
        Uxi = UT_AxialBar(xi,L,E,A0,a,q,P,ubar);
        % largest error at the nodes, of the largest displacement
        eN(b,m) = max(abs(U-Uxi))/max(abs(Uxi));
        % error of the displacement over the bar (L2 norm)
        eU(b,m) = sqrt(sum(sum(WxQ.*(Uq-Uex).^2))/sum(sum(WxQ.*Uex.^2)));
        % error in the energy norm, the integral of E*A*(du/dx - duex/dx)^2
        eE(b,m) = sqrt(sum(sum(WxQ.*E.*Aex(xQ).*(dUq-dUex).^2))/sum(sum(WxQ.*E.*Aex(xQ).*dUex.^2)));
        % equilibrium: the reaction & the applied loads add up to zero. This
        % is a check of the assembly & of the partition, not of the accuracy:
        % it holds for any U that solves the free equations, on any mesh
        eR(b,m) = abs(R + P + q*L)/(abs(P) + abs(q*L));
        UL(b,m) = U(nN);

        % kept for the figure: u at the nodes, & the axial force E*A*du/dx
        % along every element from end to end (one column per element)
        if b == Draw(1) && m == Draw(2)
            sD = linspace(-1,1,11)';
            [NsD,dNsD] = SF_FEMLG1D(sD,sn);
            [xD,~,~,dUD] = SF_Line1D(xi,Elem,U,NsD,dNsD,[]);
            ND = E*Aex(xD).*dUD;
            xDn = xi; UDn = U;
            xDe = linspace(0,L,201)';
            [UDe,NDe] = UT_AxialBar(xDe,L,E,A0,a,q,P,ubar);
        end
    end
    UL(b,end) = UT_AxialBar(L,L,E,A0,a,q,P,ubar);
end


%% RATES OF CONVERGENCE

% observed order between one mesh & the next, log(e1/e2)/log(l1/l2) with the
% element length l = L/Ne, & no rate where the error of the coarser mesh of a
% pair is already at the level of rounding (below 1e-10): the elements are
% exact there (the rounding grows with the mesh) (PL_Rate, the meshes along
% the 2nd dimension)
RateN = PL_Rate(eN,L./Ne,1e-10);
RateU = PL_Rate(eU,L./Ne,1e-10);
RateE = PL_Rate(eE,L./Ne,1e-10);


%% DISPLAY RESULTS
fprintf('Reaction & loads: largest |R + P + q*L| = %.1e of the loads\n',max(eR(:)))
for b = 1:size(Bar,1)
    fprintf('Bar %d: A(L) = %g*A0, elements of order %d, u(L) exact = %.8e m\n',b,1-Bar(b,1),Bar(b,2),UL(b,end))
    fprintf('   Ne   u(L) (m)         at the nodes (rate)   L2 norm (rate)     energy norm (rate)\n')
    for m = 1:numel(Ne)
        if m == 1
            fprintf('  %3d   %.8e   %.3e             %.3e          %.3e\n',Ne(m),UL(b,m),eN(b,m),eU(b,m),eE(b,m))
        else
            fprintf('  %3d   %.8e   %.3e (%.2f)      %.3e (%.2f)   %.3e (%.2f)\n',Ne(m),UL(b,m),eN(b,m),RateN(b,m-1),eU(b,m),RateU(b,m-1),eE(b,m),RateE(b,m-1))
        end
    end
end
fprintf('=> Rates between the last two meshes, one bar per row [at the nodes , L2 norm , energy norm]\n')
fprintf('   %.3f   %.3f   %.3f\n',[RateN(:,end) RateU(:,end) RateE(:,end)]')

% the errors against the element length, every bar that has a rate
l = L./Ne;
figure; Leg = {};
for b = 1:size(Bar,1)
    if ~isnan(RateN(b,end))
        loglog(l,eN(b,:),'^:'); hold on
        Leg{end+1} = sprintf('a = %g, p = %d: at the nodes (%.2f)',Bar(b,1),Bar(b,2),RateN(b,end));
    end
    if ~isnan(RateU(b,end))
        loglog(l,eU(b,:),'o-'); hold on
        Leg{end+1} = sprintf('a = %g, p = %d: L2 norm (%.2f)',Bar(b,1),Bar(b,2),RateU(b,end));
    end
    if ~isnan(RateE(b,end))
        loglog(l,eE(b,:),'s--'); hold on
        Leg{end+1} = sprintf('a = %g, p = %d: energy norm (%.2f)',Bar(b,1),Bar(b,2),RateE(b,end));
    end
end
xlabel('element length l (m)'); ylabel('relative error'); grid on
legend(Leg,'Location','southeast')
title('Axial bar: error against the element length (rate between the last two meshes)')

% the displacement & the axial force of the solution Draw: the force is not
% continuous from one element to the next
PL_AxialBar(xDn,UDn,xD,ND,xDe,UDe,NDe, ...
    sprintf('Axial bar: a = %g, %d elements of order %d',Bar(Draw(1),1),Ne(Draw(2)),Bar(Draw(1),2)));
