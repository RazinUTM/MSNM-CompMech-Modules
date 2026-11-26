% Title: Test of the axial bar by 1D finite elements
% Attempt: to check the element matrices against the closed forms of the
%          catalogue, & the main file against the exact solution of the
%          tapered bar on its one mesh (P2 of the ladder)
% Created on: 23/11/2026 by Claude (not yet double checked by Razin)
%   the main file is run as it stands, with ubar = 0: the term -K(f,c)*U(c)
%   of its partition is zero & not held here; TEST_SpringSystem holds the
%   same line of MAIN_SpringSystem, not this copy (a partition routine is on
%   the list for another session)
%   run in GNU Octave 8.4.0 on 23/11/2026
%   run in MATLAB by Razin, reported on 23/11/2026 (runValidation: 22 of 22 test files passed in 174.6 s)
% Last updated: 23/11/2026 by Claude
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   after an AI review against the author's way of writing codes:
%   (1) each term of the weak form alone, with the main file run on other
%       data in memory (q, P & L, which coincide in its data): against
%       (P2.4) & an independent integral of N/(EA)
%   (2) the uniform strain of the tapered bar with P = q*L*(1/a - 1), which
%       the header states
%   (3) the guards of the main file held
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 23 of 23 test files passed in 114.4 s)
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (use the kernels of core/; where there
%   is none, write one, modular, & keep it in core/), with the main file now
%   on KM_Line1D, SV_Partition, SF_Line1D & UT_AxialBar:
%   (4) the weights of the norms are a matrix (one column per element):
%       summed whole
%   (5) each term alone: the limit of the uniform bar is 1e-12 (it was
%       1e-13); [K] moves by an ulp & u(L) of 32 quadratic elements by
%       1.3e-13 in GNU Octave, rounding
%   (6) the values observed brought up to date; the partition of the header
%       is SV_Partition now, held by TEST_Line1D with prescribed values
%       away from zero
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's answers of 23/11/2026 on the choices of P2: the main file
%   is one bar on one mesh.
%   (7) the element matrices through SF_FEMLG1D
%   (8) the main file as it stands: the tapered bar on 8 linear elements
%   (9) what this file held of the four bars & the five meshes (the uniform
%       bar, the rates, the values by hand on 2 elements) is in
%       TEST_AxialBarStudy, each with the limit it had here
%   (10) each term alone, the uniform strain & the guards are held here on
%       the one bar of the main file, with a = 0 for (P2.4), & in
%       TEST_AxialBarStudy on the four bars
%   (11) the main file with ubar = 2e-3: the value held at x = 0 reaches
%       SV_Partition & every u moves by it
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Double checked on: -
% The same checks as tests/Validation/test_axial_bar.py

function TEST_AxialBar()
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'core','quadrature'),fullfile(projectRoot,'core','basis'),fullfile(projectRoot,'core','assembly'));
previous_visibility = get(0,'DefaultFigureVisible');
restore_visibility = onCleanup(@() set(0,'DefaultFigureVisible',previous_visibility));
previous_output = getenv('P2_OUTPUT_DIR');
restore_output = onCleanup(@() setenv('P2_OUTPUT_DIR',previous_output));
set(0,'DefaultFigureVisible','off');
setenv('P2_OUTPUT_DIR','');
example = fullfile(projectRoot,'problems','AxialBar');

%% THE ELEMENT MATRICES AGAINST THE CLOSED FORMS (P2.3)

% one element of length l with constant EA & q, by the lines of the main
% file: SF_FEMLG1D gives the 1D Lagrange functions on the nodes of the
% element in s
l = 0.4; EA = 3; q = 7;
[G,W] = GQ_Points(2);
kex = {EA/l*[1 -1; -1 1], EA/(3*l)*[7 -8 1; -8 16 -8; 1 -8 7]};
rex = {q*l/2*[1; 1], q*l/6*[1; 4; 1]};
for p = 1:2
    sn = linspace(-1,1,p+1)';
    [Ns,dNs] = SF_FEMLG1D(G,sn);
    Xe = 1.3 + linspace(0,l,p+1)';
    [ke,re] = KM_Line1D(Xe,1:p+1,Ns,dNs,W,@(x) EA + 0*x,@(x) q + 0*x);
    fprintf('AxialBar: element of order %d, [k] against (P2.3) %.1e, {r} %.1e\n',p, ...
        max(max(abs(ke-kex{p})))/max(max(abs(kex{p}))),max(abs(re-rex{p}))/max(abs(rex{p})));
    % observed: 3.6e-16 & 1.6e-16 (linear), 5.3e-16 & 4.8e-16 (quadratic;
    % 3.6e-16 & 3.3e-16 in GNU Octave)
    % (in MATLAB on 23/11/2026, through SF_FEMLG: 3.6e-16 & 1.6e-16, 3.6e-16 & 3.3e-16)
    assert(max(max(abs(ke-kex{p})))<1e-12*max(max(abs(kex{p}))),'AxialBar: the element matrix')
    assert(max(abs(re-rex{p}))<1e-12*max(abs(rex{p})),'AxialBar: the element load')
end

%% THE EXAMPLE AGAINST THE EXACT SOLUTION

r = resultsOfExample(fullfile(example,'MAIN_AxialBar.m'));
fprintf('AxialBar: reaction & loads %.1e, R + 150000 = %.1e of it; the weights of the norms %.1e from 2\n', ...
    r.eR,abs(r.R+150000)/150000,abs(sum(r.WxQ(:))-2));
fprintf('          at the nodes %.3e; L2 %.9e, energy %.9e; the energy norm by F.U %.1e of it\n', ...
    r.eN,r.eU,r.eE,abs(r.eE-sqrt(1-(r.F'*r.U)/(100+50*log(2))))/r.eE);

% the main file as it stands: the tapered bar (a = 0.5) on 8 linear
% elements, 9 nodes
assert(r.a==0.5 && r.p==1 && r.Ne==8 && r.nN==9 && isequal(size(r.U),[9 1]),'AxialBar: the bar & the mesh of the main file')
% the closed form by hand: the exact u at the nodes is 0 at the held end &
% at the free end (a = 0.5) u(L) = 2e-3 - 1e-3*log(2) m. observed: 0 & 0
assert(r.Uxi(1)==0 && abs(r.Uxi(end)-(2e-3-1e-3*log(2)))<1e-15,'AxialBar: the closed form')
% the held end stays where it is held (ubar = 0)
assert(r.U(1)==0,'AxialBar: the held end')
% the reaction & the loads add up to zero: R = -(P + q*L) = -150000 N.
% observed: 1.6e-15 of the loads (1.4e-15 in GNU Octave)
assert(r.eR<1e-10,'AxialBar: the reaction')
assert(abs(r.R+150000)<1e-10*150000,'AxialBar: the reaction by hand')
% the weights of the error norms add up to the length of the bar.
% observed: 1.3e-15 (0 in GNU Octave)
assert(abs(sum(r.WxQ(:))-2)<1e-12,'AxialBar: the weights of the norms')
% tapered bar: not exact at the nodes. observed: 3.7e-4
assert(r.eN>1e-9,'AxialBar: the nodes of the tapered bar')
% the energy norm by another route: the error of a Galerkin solution has the
% energy a(e,e) = a(u,u) - F.U (ubar = 0), with a(u,u) = 100 + 50*log(2) J by
% hand & F.U of the main file (134.633 J). observed: 4.6e-12 of it (3.0e-12
% in GNU Octave)
assert(abs(r.eE-sqrt(1-(r.F'*r.U)/(100+50*log(2))))<1e-8*r.eE,'AxialBar: the energy norm by the energy of the solution')
% recorded on 23/11/2026 (Python & GNU Octave agree to 3e-11 of them): the
% L2 & the energy norms of the error. They are those the file of before gave
% for this bar & mesh, in every figure; TEST_AxialBarStudy holds the same two
% as an entry of the study
% (in MATLAB on 23/11/2026, the file of before: the two in every figure printed)
% (in MATLAB on 24/11/2026, this file: the two in every figure printed; the
% reaction 1.9e-16, the energy norm by F.U 1.3e-12, each term alone 1.6e-15 at
% most, ubar = 2e-3 1.4e-16, the uniform strain 6.2e-16 at most)
rec = [7.503743828e-04 1.344212226e-02];
assert(abs(r.eU-rec(1))<1e-7*rec(1) && abs(r.eE-rec(2))<1e-7*rec(2),'AxialBar: the recorded errors')
% the figure: the axial force along every element, 11 points on each of the 8
assert(isequal(size(r.xD),[11 8]) && isequal(size(r.ND),[11 8]),'AxialBar: the axial force along the elements')

%% EACH TERM ALONE, WITH OTHER DATA (the main file run with lines of its input changed in memory)

% in the data of the main file q = P = 50e3 & q*L/2 = P, so a load put for
% another would not show: the distributed load alone, the end load alone & the
% two, on a uniform bar of 3 m (a = 0), whose nodes the elements give exactly:
% u(L) against (P2.4), u(L) = q*L^2/(2*E*A0) + P*L/(E*A0). observed: 3.2e-16 &
% 3.2e-16 at most (4.8e-16 & 4.8e-16 in GNU Octave)
main = fullfile(example,'MAIN_AxialBar.m');
for qP = [0 70e3; 30e3 0; 30e3 70e3]'
    rc = resultsOfChanged(main,{'L = 2;      %','L = 3;      %'; 'a = 0.5;    %','a = 0;      %'; 'q = 50e3;   %',sprintf('q = %g;   %%',qP(1)); 'P = 50e3;   %',sprintf('P = %g;   %%',qP(2))});
    uni = qP(1)*3^2/(2*2e8) + qP(2)*3/2e8;
    fprintf('AxialBar, uniform bar of 3 m with q = %g & P = %g: u(L) %.1e from (P2.4), the nodes %.1e\n',qP,abs(rc.U(end)-uni)/uni,rc.eN);
    assert(abs(rc.U(end)-uni)<1e-13*uni && rc.eN<1e-10 && rc.eR<1e-10,'AxialBar: each term alone')
end
% & the support away from zero: with ubar = 2e-3 m every u of the uniform bar
% moves by it, through the term -K(f,c)*U(c) of the partition, which is zero
% in the data of the main file: u(L) = 2e-3 + 1e-3 m. observed: 1.2e-15
% (1.4e-15 in GNU Octave)
rc = resultsOfChanged(main,{'a = 0.5;    %','a = 0;      %'; 'ubar = 0;   %','ubar = 2e-3; %'});
fprintf('AxialBar, uniform bar held at ubar = 2e-3: u(L) %.1e from 3e-3, the nodes %.1e\n',abs(rc.U(end)-3e-3)/3e-3,rc.eN);
assert(rc.U(1)==2e-3 && abs(rc.U(end)-3e-3)<1e-13*3e-3 && rc.eN<1e-10 && rc.eR<1e-10,'AxialBar: the support away from zero')
% P = q*L*(1/a - 1) = 100e3 N: N/A is constant, the strain uniform, & the
% elements are exact whatever their EA (the header of UT_AxialBar). observed:
% 1.5e-15, 1.7e-15 & 1.7e-15 (1.2e-15, 1.5e-15 & 1.7e-15 in GNU Octave)
rc = resultsOfChanged(main,{'P = 50e3;   %','P = 100e3;   %'});
fprintf('AxialBar, P = 100e3 (uniform strain): at the nodes %.1e, L2 %.1e, energy %.1e\n',rc.eN,rc.eU,rc.eE);
assert(rc.eN<1e-10 && rc.eU<1e-10 && rc.eE<1e-10,'AxialBar: the uniform strain of the tapered bar')

%% WHAT IS GUARDED

% a taper of 1, & fewer gauss points than the order, stop the main file with a message
m1 = messageOfChanged(main,{'a = 0.5;    %','a = 1;      %'});
m2 = messageOfChanged(main,{'p = 1;      %','p = 2;      %'; 'ng = 2;  %','ng = 1;  %'});
fprintf('AxialBar guards: %s | %s\n',m1,m2);
assert(~isempty(strfind(m1,'MAIN_AxialBar: the taper a must be 0 <= a < 1')) && ~isempty(strfind(m2,'[K] would be singular')),'AxialBar: the guards')

fprintf('TEST_AxialBar passed.\n');
end

%======================================================

%to run the main file & take the results it leaves in the workspace
%(the main file starts with clear, so nothing else is kept here)
function r = resultsOfExample(file)

evalc('run(file)');
names = {'a','p','Ne','nN','eN','eU','eE','eR','F','U','R','Uxi','WxQ','xD','ND'};
for k = 1:numel(names)
    r.(names{k}) = eval(names{k});
end
close all

end

%======================================================

%to run the main file with lines of its input changed in memory, each
%{old , new} found once, & take some of its results (the file starts
%with clear, so the names are given after it returns)
function r = resultsOfChanged(file,changes)

txt = strrep(fileread(file),"mfilename('fullpath')","file");
for k = 1:size(changes,1)
    assert(numel(strfind(txt,changes{k,1}))==1,'resultsOfChanged: a line to change is not found once')
    txt = strrep(txt,changes{k,1},changes{k,2});
end
evalc(txt);
names = {'U','eN','eU','eE','eR'};
for k = 1:numel(names)
    r.(names{k}) = eval(names{k});
end
close all

end

%======================================================

%to run the main file with lines of its input changed in memory & take what
%it stopped with ('' if it ran to the end)
function msg = messageOfChanged(file,changes)

txt = strrep(fileread(file),"mfilename('fullpath')","file");
for k = 1:size(changes,1)
    assert(numel(strfind(txt,changes{k,1}))==1,'messageOfChanged: a line to change is not found once')
    txt = strrep(txt,changes{k,1},changes{k,2});
end
try
    evalc(txt);
    msg = '';
catch err
    msg = err.message;
end
close all

end

%======================================================
