% Title: Code for an Axial Bar by 1D Finite Elements (linear or quadratic)
% Attempt: to solve a tapered bar under a distributed load & an end load by
%          2-noded (or 3-noded) elements on one mesh, & to compare with the
%          closed form (P2 of the ladder)
% Created on: 23/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 23/11/2026 (the next problem of his ladder): P2
%   of "A Graded Catalogue of Computational Mechanics Problems for Civil
%   Engineering" (v0.1, 22/11/2026). The functions of an element are those
%   of SF_FEMLG on its nodes in s (one row of nodes, so they are the 1D
%   Lagrange functions) & its gauss points are those of GQ_Points; the
%   element matrices are assembled as in MAIN_SpringSystem & the node at
%   x = 0 is held by partition. Nothing is printed or drawn inside the loops.
%   With P = q*L*(1/a - 1) the strain of the tapered bar is uniform & the
%   elements are exact: such a bar has no rate to show.
%   The load is uniform: that the nodes of a uniform bar are exact "for any
%   load" (the catalogue) is shown for this load only
%   run in GNU Octave 8.4.0 on 23/11/2026
%   run in MATLAB by Razin, reported on 23/11/2026 (runValidation: 22 of 22 test files passed in 174.6 s)
% Last updated: 23/11/2026 by Claude
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (the codes written his way, by his
%   skill "my way of writing codes"), after an AI review against it:
%   (1) no rate where the error of the coarser mesh of a pair is already at
%       the level of rounding, as in MAIN_SteadyConduction1D (the rule of the
%       first mesh printed rates of rounding once the meshes were fine)
%   (2) the file stops where a taper is not 0 <= a < 1 or ng is fewer than
%       p, & warns where a taper is below 1e-4
%   (3) a comment: the sign convention
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 23 of 23 test files passed in 114.4 s)
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (use the kernels of core/; where there
%   is none, write one, modular, & keep it in core/):
%   (4) the exact solution is UT_AxialBar, the assembly KM_Line1D, the
%       partition SV_Partition & u, du/dx at the gauss points of the error
%       norms & of the figure SF_Line1D (core/), in place of their lines here;
%       no recorded number moved
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's answers of 23/11/2026 on the choices of P2 (2 changed):
%   (5) one bar on one mesh (a = 0.5, p = 1, Ne = 8). The loops over the bars
%       & the meshes, the rates & the figure of the errors are now
%       studies/AxialBar/STUDY_AxialBarConvergence, which keeps the names of
%       this file: this file as it stands is its entry [a = 0.5, p = 1, Ne = 8]
%   (6) the shape functions by SF_FEMLG1D (the file wrote the points & the
%       nodes as [s 0] to call SF_FEMLG), so sn is a vector
%   (7) the comment at eR says what it checks
%   (8) the guards & the warning are written for the one bar (a, p); the
%       rule of the rates went with the rates to the study
%   no number moved: eN, eU, eE & eR are those of the file of before for this
%   bar & mesh in every figure (Python & GNU Octave)
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
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
a = 0.5;    % taper: A(x) = A0*(1 - a*x/L) with 0 <= a < 1; a = 0 is the uniform bar of the catalogue
q = 50e3;   % distributed axial load, uniform (N/m)
P = 50e3;   % load at the free end x = L (N)
ubar = 0;   % displacement held at x = 0 (m): 0 as in (P2.4); with another
            % value every u moves by it through the term -K(f,c)*U(c) of
            % SV_Prescribed, which TEST_AxialBar holds with ubar = 2e-3

% the elements
p = 1;      % order of the elements: p = 1 the linear element (2 nodes); p = 2 gives the 3-noded (quadratic) element
Ne = 8;     % number of elements

ng = 2;  % gauss points of an element for [k] & {r}: exact for this bar (degree 3); fewer than p make [K] singular
nge = 8; % gauss points of an element for the error norms (a higher rule: the exact u is not a polynomial)

% the ranges of the input (the closed form & the elements hold only there)
if a < 0 || a >= 1
    error('MAIN_AxialBar: the taper a must be 0 <= a < 1 (a = %g)',a)
end
if ng < p
    error('MAIN_AxialBar: ng = %d is fewer than the order p = %d: [K] would be singular',ng,p)
end
if a > 0 && a < 1e-4
    warning('MAIN_AxialBar: a taper below 1e-4 leaves the closed form about eps/a of its figures')
end


%% AREA & AXIAL STIFFNESS

% A(x) = A0*(1 - a*x/L), k = E*A for KM_Line1D; the exact u, N & du/dx are
% those of UT_AxialBar
Aex = @(x) A0*(1 - a*x/L);
EAx = @(x) E*Aex(x);
qx = @(x) q + 0*x;


%% SHAPE FUNCTIONS OF THE ELEMENT

% in s from -1 to 1, the xi of the catalogue (the same for every element): its
% nodes from end to end, its gauss points & weights for [k] & {r}, & those for
% the error norms
sn = linspace(-1,1,p+1)';
[G,W] = GQ_Points(ng);
[Ns,dNs] = SF_FEMLG1D(G,sn);
[Ge,We] = GQ_Points(nge);
[Nse,dNse] = SF_FEMLG1D(Ge,sn);


%% NODES & ELEMENTS

% equal elements, the nodes of each from end to end
nN = p*Ne+1;
xi = linspace(0,L,nN)';
Elem = p*((1:Ne)'-1) + (1:p+1);


%% STIFFNESS MATRIX & FORCE VECTOR

% the integrals of E*A*dNi/dx*dNj/dx & Ni*q over the elements (KM_Line1D), &
% the load at x = L
[K,F] = KM_Line1D(xi,Elem,Ns,dNs,W,EAx,qx);
F(nN) = F(nN) + P;


%% SUPPORT, SOLUTION & REACTION

% the node at x = 0 held at ubar, by partition, K(f,f)*U(f) = F(f) - K(f,c)*U(c)
% & R = K(c,:)*U - F(c)
[U,R] = SV_Prescribed(K,F,1,ubar);


%% COMPARE WITH THE EXACT SOLUTION

% u & du/dx of the elements at the gauss points of the error norms [xQ], with
% their weights [WxQ] (one column per element)
[xQ,WxQ,Uq,dUq] = SF_Line1D(xi,Elem,U,Nse,dNse,We);
[Uex,~,dUex] = UT_AxialBar(xQ,L,E,A0,a,q,P,ubar);
Uxi = UT_AxialBar(xi,L,E,A0,a,q,P,ubar);
% largest error at the nodes, of the largest displacement
eN = max(abs(U-Uxi))/max(abs(Uxi));
% error of the displacement over the bar (L2 norm)
eU = sqrt(sum(sum(WxQ.*(Uq-Uex).^2))/sum(sum(WxQ.*Uex.^2)));
% error in the energy norm, the integral of E*A*(du/dx - duex/dx)^2
eE = sqrt(sum(sum(WxQ.*E.*Aex(xQ).*(dUq-dUex).^2))/sum(sum(WxQ.*E.*Aex(xQ).*dUex.^2)));
% equilibrium: the reaction & the applied loads add up to zero. This is a
% check of the assembly & of the partition, not of the accuracy: it holds for
% any U that solves the free equations, on any mesh
eR = abs(R + P + q*L)/(abs(P) + abs(q*L));


%% DISPLAY RESULTS
fprintf('Axial bar: A(L) = %g*A0, %d elements of order %d (%d nodes)\n',1-a,Ne,p,nN)
fprintf('  u(L) = %.8e m (exact %.8e m)\n',U(nN),Uxi(nN))
fprintf('  R    = %.8g N; R + P + q*L = %.1e of the loads\n',R,eR)
fprintf('=> Error: at the nodes %.3e, L2 norm %.3e, energy norm %.3e\n',eN,eU,eE)

% the displacement & the axial force against the exact ones: u at the nodes,
% & the axial force E*A*du/dx along every element from end to end (one column
% per element), which is not continuous from one element to the next
sD = linspace(-1,1,11)';
[NsD,dNsD] = SF_FEMLG1D(sD,sn);
[xD,~,~,dUD] = SF_Line1D(xi,Elem,U,NsD,dNsD,[]);
ND = E*Aex(xD).*dUD;
xDe = linspace(0,L,201)';
[UDe,NDe] = UT_AxialBar(xDe,L,E,A0,a,q,P,ubar);

fig = PL_AxialBar(xi,U,xD,ND,xDe,UDe,NDe,sprintf('Axial bar: a = %g, %d elements of order %d',a,Ne,p));

%% OPTIONAL EXPORT
outputDir = getenv('P2_OUTPUT_DIR');
if ~isempty(outputDir)
    [ok,message] = mkdir(outputDir);
    assert(ok,'%s',message)
    fid = fopen(fullfile(outputDir,'p2_displacement_matlab.csv'),'w');
    assert(fid>=0,'Could not open the displacement CSV.')
    fprintf(fid,'x_m,fem_m,exact_m\n');
    fprintf(fid,'%.16g,%.16g,%.16g\n',[xi U Uxi].');
    fclose(fid);
    elements = repelem((1:Ne)',size(xD,1));
    forceRows = [elements, xD(:), ND(:), P + q*(L - xD(:))];
    fid = fopen(fullfile(outputDir,'p2_force_matlab.csv'),'w');
    assert(fid>=0,'Could not open the axial force CSV.')
    fprintf(fid,'element,x_m,fem_N,exact_N\n');
    fprintf(fid,'%.16g,%.16g,%.16g,%.16g\n',forceRows.');
    fclose(fid);
    print(fig,fullfile(outputDir,'p2_results_matlab.png'),'-dpng','-r160');
end
