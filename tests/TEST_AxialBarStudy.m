% Title: Test of the study of the axial bar
% Attempt: to check that the study runs to the end & gives the four bars on
%          the five meshes against the exact solution, with the rates, & that
%          the main file is an entry of it (P2 of the ladder)
% Created on: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's answers of 23/11/2026 on the choices of P2
%   the checks of 24/11/2026 of TEST_AxialBar (each term alone on other
%   data, the uniform strain, the guards) are here for the four bars, with
%   its limits & its observed values of that day (the lines on KM_Line1D,
%   SV_Partition, SF_Line1D & UT_AxialBar)
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 24/11/2026 by Claude
% Double checked on: -
% The same checks as tests/Validation/test_axial_bar_study.py. The study is
% run as it stands, in a workspace of its own, & so is the main file of the
% problem. Each limit is stated with the value observed in Python and in GNU
% Octave 8.4.0 on 24/11/2026. What is held: the four bars on the five meshes
% against the exact solution (the uniform bar exact at the nodes, the
% quadratic elements exact everywhere; the rates; the values worked by hand;
% the four recorded errors), which TEST_AxialBar held until the main file
% became one bar on one mesh, each with the limit it had there; and, new with
% the study, the rate of the error at the nodes, and the main file as it
% stands as the entry [a = 0.5, p = 1, Ne = 8] of the study

function TEST_AxialBarStudy()
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'core','quadrature'),fullfile(projectRoot,'core','basis'));
previous_visibility = get(0,'DefaultFigureVisible');
restore_visibility = onCleanup(@() set(0,'DefaultFigureVisible',previous_visibility));
previous_output = getenv('P2_OUTPUT_DIR');
restore_output = onCleanup(@() setenv('P2_OUTPUT_DIR',previous_output));
set(0,'DefaultFigureVisible','off');
setenv('P2_OUTPUT_DIR','');
studyFolder = fullfile(projectRoot,'studies','AxialBar');
example = fullfile(projectRoot,'problems','AxialBar');

%% STUDY_AxialBarConvergence: THE FOUR BARS ON THE FIVE MESHES

r = resultsOfStudy(fullfile(studyFolder,'STUDY_AxialBarConvergence.m'));
fprintf('STUDY_AxialBarConvergence: reaction & loads %.1e; uniform bar at the nodes %.1e, quadratic everywhere %.1e & %.1e\n', ...
    max(r.eR(:)),max(max(r.eN(1:2,:))),max(r.eU(2,:)),max(r.eE(2,:)));
fprintf('          rates between the last two meshes [nodes L2 energy]: uniform linear %.4f %.4f %.4f, tapered linear %.4f %.4f %.4f, tapered quadratic %.4f %.4f %.4f\n', ...
    r.RateN(1,end),r.RateU(1,end),r.RateE(1,end),r.RateN(3,end),r.RateU(3,end),r.RateE(3,end),r.RateN(4,end),r.RateU(4,end),r.RateE(4,end));
fprintf('          rates at the nodes of the tapered bar: linear %.4f %.4f %.4f %.4f, quadratic %.4f %.4f %.4f %.4f\n',r.RateN(3:4,:)');
fprintf('          tapered bar on 8 elements: L2 %.9e & %.9e, energy %.9e & %.9e\n',r.eU(3:4,3),r.eE(3:4,3));

% the study as it stands: rows of Bar = uniform linear, uniform quadratic,
% tapered linear, tapered quadratic; five meshes; the solution drawn is the
% tapered bar on 2 linear elements
assert(isequal(r.Bar,[0 1; 0 2; 0.5 1; 0.5 2]) && isequal(r.Ne,[2 4 8 16 32]),'STUDY_AxialBarConvergence: the bars & the meshes of the study')
assert(isequal(r.Draw,[3 1]),'STUDY_AxialBarConvergence: the solution drawn')
assert(isequal(size(r.eN),[4 5]) && isequal(size(r.UL),[4 6]) && isequal(size(r.RateN),[4 4]),'STUDY_AxialBarConvergence: sizes')
% the closed forms by hand: uniform, u(L) = q*L^2/(2*E*A) + P*L/(E*A) =
% 5e-4 + 5e-4 m; tapered (a = 0.5), u(L) = 2e-3 - 1e-3*log(2) m
assert(max(abs(r.UL(1:2,end)-1e-3))<1e-15 && max(abs(r.UL(3:4,end)-(2e-3-1e-3*log(2))))<1e-15,'STUDY_AxialBarConvergence: the closed form')
% u(L) of the elements: at the closed form on every mesh of the uniform bar;
% by hand for the tapered bar on 2 linear elements, 125e3/175e6 = 1/1400 m at
% the middle node & 1/1400 + 75e3/125e6 = 23/17500 m at the end. observed:
% 9.6e-17 & 2.2e-19 (1.5e-16 & 0 in GNU Octave)
assert(max(max(abs(r.UL(1:2,1:end-1)-1e-3)))<1e-13 && abs(r.UL(3,1)-23/17500)<1e-13,'STUDY_AxialBarConvergence: u(L) of the elements')
% the reaction & the loads add up to zero on every mesh. observed: 9.2e-14
% (1.4e-13 in GNU Octave)
assert(all(r.eR(:)<1e-10),'STUDY_AxialBarConvergence: the reaction')
% the weights of the error norms (of the last bar & mesh) add up to the
% length of the bar. observed: 0 (8.9e-16 in GNU Octave)
assert(abs(sum(r.WxQ(:))-2)<1e-12,'STUDY_AxialBarConvergence: the weights of the norms')
% uniform bar: the elements are exact at the nodes, & the quadratic ones
% everywhere, so that bar has no rate. observed: 9.6e-14, 1.0e-13 & 1.0e-13
% (1.5e-13 in GNU Octave)
% (in MATLAB on 23/11/2026, with these lines in the main file: 1.5e-13; the reaction 1.3e-13)
assert(all(all(r.eN(1:2,:)<1e-10)),'STUDY_AxialBarConvergence: the nodes of the uniform bar')
assert(all(r.eU(2,:)<1e-10) && all(r.eE(2,:)<1e-10),'STUDY_AxialBarConvergence: the quadratic elements on the uniform bar')
assert(all(isnan(r.RateU(2,:))) && all(isnan(r.RateE(2,:))),'STUDY_AxialBarConvergence: no rate for a bar solved exactly')
% & the nodes of the uniform bar have no rate, with either element
assert(all(all(isnan(r.RateN(1:2,:)))),'STUDY_AxialBarConvergence: no rate for nodes that are exact')
% uniform bar, linear elements: the error between the nodes is the same
% curve on every element, q/(2EA)*s*(l-s), so the rates are 2 & 1 to
% rounding. observed: 3.2e-11 & 1.4e-15 from them (5.9e-11 in GNU Octave)
assert(all(abs(r.RateU(1,:)-2)<1e-8) && all(abs(r.RateE(1,:)-1)<1e-8),'STUDY_AxialBarConvergence: the rates of the uniform bar')
% & the errors themselves by hand: the L2 norm is l^2/(2*sqrt(204)) & the
% energy norm l/sqrt(52), with l = L/Ne. observed: 2.6e-11 & 1.2e-15 of them
% (5.4e-11 & 4.0e-16 in GNU Octave)
l = 2./r.Ne;
assert(all(abs(r.eU(1,:)-l.^2/(2*sqrt(204)))<1e-8*r.eU(1,:)),'STUDY_AxialBarConvergence: the L2 norm of the uniform bar by hand')
assert(all(abs(r.eE(1,:)-l/sqrt(52))<1e-8*r.eE(1,:)),'STUDY_AxialBarConvergence: the energy norm of the uniform bar by hand')
% tapered bar: no longer exact at the nodes. observed: 2.3e-5 at the least
% (linear), 2.7e-9 (quadratic)
assert(all(all(r.eN(3:4,:)>1e-9)),'STUDY_AxialBarConvergence: the nodes of the tapered bar')
% by hand on 2 linear elements: the largest error at the nodes is at x = L,
% 23/17500 against 2e-3 - 1e-3*log(2). observed: 1.7e-16 from it (0 in GNU Octave)
assert(abs(r.eN(3,1)-(23/17500-(2e-3-1e-3*log(2)))/(2e-3-1e-3*log(2)))<1e-10,'STUDY_AxialBarConvergence: the error at the nodes by hand')
% & the energy norm there by hand: a(e,e) = a(u,u) - F.U, with a(u,u) = 100 +
% 50*log(2) J & F.U = 50e3/1400 + 75e3*23/17500 = 940/7 J. observed: 2.4e-14
% of it
assert(abs(r.eE(3,1)-sqrt(1-(940/7)/(100+50*log(2))))<1e-8*r.eE(3,1),'STUDY_AxialBarConvergence: the energy norm of the tapered bar by hand')
% tapered bar: the rates between the last two meshes are p+1 in the L2 norm
% & p in the energy norm. observed: 1.9989 & 0.9995, 2.9978 & 1.9986
assert(abs(r.RateU(3,end)-2)<0.01 && abs(r.RateE(3,end)-1)<0.01,'STUDY_AxialBarConvergence: the rates of the linear elements')
assert(abs(r.RateU(4,end)-3)<0.01 && abs(r.RateE(4,end)-2)<0.01,'STUDY_AxialBarConvergence: the rates of the quadratic elements')
% & every error falls from one mesh to the next
assert(all(all(diff(r.eU([1 3 4],:),1,2)<0)) && all(all(diff(r.eE([1 3 4],:),1,2)<0)),'STUDY_AxialBarConvergence: the errors against the mesh')
% recorded on 23/11/2026 (Python & GNU Octave agree to 3e-11 of them): the
% tapered bar on 8 elements, linear & quadratic
% (in MATLAB on 23/11/2026, with these lines in the main file: the four in
% every figure printed, & the rates)
% (in MATLAB on 24/11/2026, the study: the four & every rate in every figure
% printed; the main file against its entry 0; each term alone 1.7e-13 & 1.0e-15
% at most)
rec = [7.503743828e-04 1.240235526e-05; 1.344212226e-02 3.418021231e-04];
assert(all(abs(r.eU(3:4,3)'-rec(1,:))<1e-7*rec(1,:)) && all(abs(r.eE(3:4,3)'-rec(2,:))<1e-7*rec(2,:)),'STUDY_AxialBarConvergence: the recorded errors')
% the solution drawn (tapered, 2 linear elements): the axial force of the
% elements is not continuous at the node between them; by hand it jumps
% from 750000/7 to 90000 N there (the exact force is 100000 N)
assert(abs(r.ND(end,1)-750000/7)<1e-4 && abs(r.ND(1,2)-90000)<1e-4,'STUDY_AxialBarConvergence: the axial force between the elements')

% new with the study: the rate of the error at the nodes.
% the rates that the study prints & draws are those of its errors. observed: 0
dN = log(r.Ne(2:end)./r.Ne(1:end-1));
assert(max(max(abs(r.RateN(3:4,:) - log(r.eN(3:4,1:end-1)./r.eN(3:4,2:end))./dN)))<1e-12,'STUDY_AxialBarConvergence: RateN')
assert(max(max(abs(r.RateU([1 3 4],:) - log(r.eU([1 3 4],1:end-1)./r.eU([1 3 4],2:end))./dN)))<1e-12,'STUDY_AxialBarConvergence: RateU')
assert(max(max(abs(r.RateE([1 3 4],:) - log(r.eE([1 3 4],1:end-1)./r.eE([1 3 4],2:end))./dN)))<1e-12,'STUDY_AxialBarConvergence: RateE')
% tapered bar: the error at the nodes falls on every mesh, at the order 2p:
% 2 with the linear elements & towards 4 with the quadratic ones, whose rate
% still rises on these meshes. observed: 1.9474, 1.9857, 1.9963, 1.9991
% (linear); 3.3885, 3.6698, 3.8291, 3.9132 (quadratic; 3.9131 in GNU Octave,
% where the error is 2.7e-9 of u)
assert(all(all(diff(r.eN(3:4,:),1,2)<0)),'STUDY_AxialBarConvergence: the error at the nodes against the mesh')
assert(abs(r.RateN(3,end)-2)<0.01,'STUDY_AxialBarConvergence: the rate at the nodes of the linear elements')
assert(all(diff(r.RateN(4,:))>0) && abs(r.RateN(4,end)-4)<0.15,'STUDY_AxialBarConvergence: the rate at the nodes of the quadratic elements')
% so the nodes are better than the bar between them: the error at the nodes
% over the L2 norm falls with the mesh for the quadratic elements. observed:
% 0.107 to 0.0137
assert(all(diff(r.eN(4,:)./r.eU(4,:))<0),'STUDY_AxialBarConvergence: the nodes against the bar between them')

%% MAIN_AxialBar AS AN ENTRY OF THE STUDY

% the main file as it stands (one bar on one mesh) is the entry of the study
% with its taper a, its order p & its number of elements Ne:
% [a = 0.5, p = 1, Ne = 8], row 3 of Bar & the third mesh
g = resultsOfStudy(fullfile(example,'MAIN_AxialBar.m'));
b = find(r.Bar(:,1)==g.a & r.Bar(:,2)==g.p);
m = find(r.Ne==g.Ne);
assert(isequal(b,3) && isequal(m,3),'STUDY_AxialBarConvergence: the entry of the main file')
fprintf('          the main file against the entry [%g %d %d]: eN %.1e, eU %.1e, eE %.1e of them, eR %.1e, u(L) %.1e & exact %.1e\n',g.a,g.p,g.Ne, ...
    abs(g.eN-r.eN(b,m))/g.eN,abs(g.eU-r.eU(b,m))/g.eU,abs(g.eE-r.eE(b,m))/g.eE,abs(g.eR-r.eR(b,m)),abs(g.U(end)-r.UL(b,m)),abs(g.Uxi(end)-r.UL(b,end)));
% the three errors & u(L), of the elements & exact. observed: 0, in every
% figure, in Python & in GNU Octave (the same lines on the same kernels)
assert(abs(g.eN-r.eN(b,m))<1e-12*g.eN && abs(g.eU-r.eU(b,m))<1e-12*g.eU && abs(g.eE-r.eE(b,m))<1e-12*g.eE,'STUDY_AxialBarConvergence: the errors of the main file')
assert(abs(g.U(end)-r.UL(b,m))<1e-12*abs(g.U(end)) && abs(g.Uxi(end)-r.UL(b,end))<1e-12*abs(g.Uxi(end)),'STUDY_AxialBarConvergence: u(L) of the main file')
% the reaction with the loads is the same rounding in both. observed: no
% difference (each is 1.6e-15; 1.4e-15 in GNU Octave)
assert(abs(g.eR-r.eR(b,m))<1e-10,'STUDY_AxialBarConvergence: the reaction of the main file')

%% EACH TERM ALONE, WITH OTHER DATA (the study run with lines of its input changed in memory)

% in the data of the study q = P = 50e3 & q*L/2 = P, so a load put for another
% would not show: the distributed load alone, the end load alone & the two, on
% a bar of 3 m. The uniform bar against (P2.4), u(L) = q*L^2/(2*E*A0) +
% P*L/(E*A0); the tapered bar against an independent integral of N/(E*A) by 40
% gauss points. observed: 1.9e-14 & 7.9e-16 at most (1.7e-13 & 4.5e-16 in GNU
% Octave; 2.4e-14 there until the assembly was KM_Line1D, which takes E*A as
% one function: [K] moves by an ulp & u(L) of 32 quadratic elements by
% 1.3e-13, rounding, so the limit is 1e-12 since 24/11/2026)
study = fullfile(studyFolder,'STUDY_AxialBarConvergence.m');
[t40,w40] = GQ_Points(40);
for qP = [0 70e3; 30e3 0; 30e3 70e3]'
    rc = resultsOfChanged(study,{'L = 2;      %','L = 3;      %'; 'q = 50e3;   %',sprintf('q = %g;   %%',qP(1)); 'P = 50e3;   %',sprintf('P = %g;   %%',qP(2))});
    uni = qP(1)*3^2/(2*2e8) + qP(2)*3/2e8;
    s40 = 3/2*(t40+1);
    tap = 3/2*sum(w40.*(qP(2) + qP(1)*(3 - s40))./(2e8*(1 - 0.5*s40/3)));
    fprintf('STUDY_AxialBarConvergence, q = %g & P = %g on 3 m: u(L) %.1e (uniform) & %.1e (tapered) from the references\n',qP,max(max(abs(rc.UL(1:2,:)-uni)))/uni,max(abs(rc.UL(3:4,end)-tap))/tap);
    assert(all(all(abs(rc.UL(1:2,:)-uni)<1e-12*uni)) && all(abs(rc.UL(3:4,end)-tap)<1e-13*tap) && all(rc.eR(:)<1e-10),'STUDY_AxialBarConvergence: each term alone')
end
% P = q*L*(1/a - 1) = 100e3 N: N/A is constant, the strain uniform, & the elements
% are exact whatever their EA (the header). observed: 9.9e-15 & 2.0e-14
rc = resultsOfChanged(study,{'P = 50e3;   %','P = 100e3;   %'});
assert(all(all(rc.eU(3:4,:)<1e-10)) && all(all(rc.eE(3:4,:)<1e-10)) && all(all(isnan(rc.RateU(3:4,:)))),'STUDY_AxialBarConvergence: the uniform strain of the tapered bar')

%% WHAT IS GUARDED

% a taper of 1 & fewer gauss points than the order stop the study, with a message
m1 = messageOfChanged(study,{'       0.5 1;','       1   1;'});
m2 = messageOfChanged(study,{'ng = 2;  %','ng = 1;  %'});
fprintf('STUDY_AxialBarConvergence guards: %s | %s\n',m1,m2);
assert(~isempty(strfind(m1,'taper a must be 0 <= a < 1 (row 3 of Bar)')) && ~isempty(strfind(m2,'[K] would be singular')),'STUDY_AxialBarConvergence: the guards')

fprintf('TEST_AxialBarStudy passed.\n');
end

%======================================================

%to run a study or a main file & take what it leaves in the workspace
%(it starts with clear, which also clears what this function holds; so the
%names are taken from the workspace after the run)
function g = resultsOfStudy(file)

evalc('run(file)');
vars = who;
for k = 1:numel(vars)
    g.(vars{k}) = eval(vars{k});
end
close all

end

%======================================================

%to run the study with lines of its input changed in memory, each
%{old , new} found once, & take some of its results (the file starts
%with clear, so the names are given after it returns)
function r = resultsOfChanged(file,changes)

txt = strrep(fileread(file),"mfilename('fullpath')","file");
for k = 1:size(changes,1)
    assert(numel(strfind(txt,changes{k,1}))==1,'resultsOfChanged: a line to change is not found once')
    txt = strrep(txt,changes{k,1},changes{k,2});
end
evalc(txt);
names = {'UL','eR','eU','eE','RateU'};
for k = 1:numel(names)
    r.(names{k}) = eval(names{k});
end
close all

end

%======================================================

%to run the study with lines of its input changed in memory & take what
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
