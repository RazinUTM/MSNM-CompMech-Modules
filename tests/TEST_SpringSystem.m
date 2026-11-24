% Title: Test of the system of springs by the direct stiffness method
% Attempt: to check the main file against the closed form of springs in
%          series, for the support at zero & away from zero (P1 of the ladder)
% Created on: 22/11/2026 by Claude (not yet double checked by Razin)
%   run in GNU Octave 8.4.0 on 22/11/2026
%   run in MATLAB by Razin, reported on 23/11/2026 (runValidation: 21 of 21 test files passed in 167.7 s)
% Last updated: 22/11/2026 by Claude
% Double checked on: -
% The same checks as tests/Validation/test_spring_system.py

function TEST_SpringSystem()
projectRoot = fileparts(fileparts(mfilename('fullpath')));
example = fullfile(projectRoot,'problems','SpringSystem');

%% THE EXAMPLE AGAINST THE CLOSED FORM

r = resultsOfExample(fullfile(example,'MAIN_SpringSystem.m'));
fprintf('SpringSystem: row sums of [K] %.1e; u3 = %.8g & %.8g m (exact 0.25 & 0.26)\n',r.Krow,r.U(3,:));
fprintf('              errors: partition %.1e & %.1e, without -K(f,c)U(c) %.1e & %.1e, penalty %.1e & %.1e\n',r.eU,r.eC,r.eP);
fprintf('              reactions %s N; condition numbers %.2f & %.2e\n',mat2str(r.R,8),r.condKff,r.condKp);

% before the support every row of [K] adds up to zero. observed: 0
assert(r.Krow<1e-15,'SpringSystem: [K] before the support')
% the closed form by hand for u1 = 0: u2 = 100/2000 = 0.05, u3 = 0.05 + 100/500 = 0.25
assert(max(abs(r.Uex(:,1)-[0; 0.05; 0.25]))<1e-15,'SpringSystem: the closed form')
% by partition, for u1 = 0 & u1 = 0.01 m. observed: 0 & 2.7e-17
assert(all(r.eU<1e-14),'SpringSystem: the displacements')
% the reaction is -F, & with the loads it adds up to zero. observed: 0
assert(max(abs(r.R-r.Rex))<1e-12 && all(r.eR<1e-14),'SpringSystem: the reaction')
% without the term -K(f,c)U(c) the answer is right only for u1 = 0: the check
% catches the mistake, which leaves u off by u1. observed: 0 & 3.8e-2
assert(r.eC(1)<1e-14 && abs(r.eC(2)-r.ubar(2)/max(abs(r.Uex(:,2))))<1e-12,'SpringSystem: the forgotten term')
% the penalty is off by about F/(Alpha max(diag(K))) & conditions [K] worse.
% observed: 1.6e-7 & 1.5e-7; 6.85 & 6.55e6
assert(all(r.eP>1e-8 & r.eP<1e-6),'SpringSystem: the penalty')
assert(abs(r.condKff-6.85)<0.01 && r.condKp>1e5*r.condKff,'SpringSystem: the condition numbers')

assert(max(abs(r.U(:,2)-[0.01; 0.06; 0.26]))<1e-14,'SpringSystem: displaced support')
assert(max(abs(r.extension(:)-[0.05; 0.20; 0.05; 0.20]))<1e-14,'SpringSystem: spring extensions')
assert(max(abs(r.spring_force(:)-100))<1e-12,'SpringSystem: spring forces')
assert(max(abs(r.equilibrium_residual(:)))<1e-12,'SpringSystem: nodal equilibrium')
fprintf('TEST_SpringSystem passed.\n');
end

function r = resultsOfExample(file)
% Suppress display and export while keeping the plotting path under test.
previous_visibility = get(0,'DefaultFigureVisible');
previous_output = getenv('P1_OUTPUT_DIR');
restore_visibility = onCleanup(@() set(0,'DefaultFigureVisible',previous_visibility));
restore_output = onCleanup(@() setenv('P1_OUTPUT_DIR',previous_output));
set(0,'DefaultFigureVisible','off');
setenv('P1_OUTPUT_DIR','');
evalc('run(file)');
close_figure = onCleanup(@() close(fig));
axes_handles = findall(fig,'Type','axes');
% Octave represents each legend as a separate axes with the tag legend.
assert(sum(~strcmp(get(axes_handles,'Tag'),'legend'))==4,'SpringSystem: plot panels')
names = {'U','Uex','R','Rex','eU','eC','eP','eR','Krow','condKff','condKp','ubar', ...
         'extension','spring_force','equilibrium_residual'};
for k = 1:numel(names)
    r.(names{k}) = eval(names{k});
end
end
