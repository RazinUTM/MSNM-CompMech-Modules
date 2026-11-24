% Title: Code for a System of Springs by the Direct Stiffness Method
% Attempt: to assemble the element matrices into [K], hold a node (at zero
%          or not), solve & recover the reaction, against the closed form
%          of springs in series (P1 of the ladder)
% Created on: 22/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 22/11/2026 (the problems of his ladder one by
%   one, from the simplest): P1 of "A Graded Catalogue of Computational
%   Mechanics Problems for Civil Engineering" (v0.1, 22/11/2026). The
%   sections are those of MAIN_PlateWithHole.m; no routine of core/ is
%   needed. The support is held by partition, & by penalty for comparison
%   run in GNU Octave 8.4.0 on 22/11/2026
%   run in MATLAB by Razin, reported on 23/11/2026 (runValidation: 21 of 21 test files passed in 167.7 s)
% Last updated: 22/11/2026 by Claude
% Double checked on: -

%% DATA INPUT

% springs [node i , node j , stiffness ke (N/m)], one row per spring
Spr = [1 2 2000;
       2 3  500];
nN = 3; % number of nodes

% Applied Forces [node , force (N)]
Fn = [3 100];

% Support conditions: node Su is held at the displacement ubar (m); the
% figure of the catalogue (0), & a value that is not 0, which shows the term
% -K(f,c)*U(c) that is most often forgotten
Su = 1;
ubar = [0 0.01];

% penalty factor, times max(diag(K)) (as the meshless problems hold their supports)
Alpha = 1e6;


%% STIFFNESS MATRIX (assembly of the element matrices)
K = zeros(nN);
for e = 1:size(Spr,1)
    % loop over springs
    ij = Spr(e,1:2);
    K(ij,ij) = K(ij,ij) + Spr(e,3)*[1 -1; -1 1];
end
% before the support every row adds up to zero: the rigid motion, [K] is singular
Krow = max(abs(sum(K,2)))/max(abs(K(:)));

% Force vector
F = zeros(nN,1);
F(Fn(:,1)) = F(Fn(:,1)) + Fn(:,2);


%% SUPPORTS, SOLUTION & REACTIONS (one column for each value of ubar)
c = Su;                         % constrained unknowns
f = setdiff((1:nN)',c);         % free unknowns
U = zeros(nN,numel(ubar));      % displacements by partition
R = zeros(numel(c),numel(ubar)); % reactions
Uw = zeros(nN,numel(ubar));     % the same without -K(f,c)*U(c): wrong unless ubar = 0
Up = zeros(nN,numel(ubar));     % the same by penalty
for s = 1:numel(ubar)
    % partition: K(f,f)*U(f) = F(f) - K(f,c)*U(c)
    U(c,s) = ubar(s);
    U(f,s) = K(f,f)\(F(f) - K(f,c)*U(c,s));
    R(:,s) = K(c,:)*U(:,s) - F(c);
    % the commonest mistake: the row & column of c deleted, the right-hand side not corrected
    Uw(c,s) = ubar(s);
    Uw(f,s) = K(f,f)\F(f);
    % penalty: Alpha*max(diag(K)) added to K(c,c), & times ubar to F(c)
    Kp = K; Fp = F;
    Kp(c,c) = Kp(c,c) + Alpha*max(diag(K))*eye(numel(c));
    Fp(c) = Fp(c) + Alpha*max(diag(K))*ubar(s);
    Up(:,s) = Kp\Fp;
end
condKff = cond(K(f,f));
condKp = cond(Kp);


%% COMPARE WITH THE EXACT SOLUTION

% springs in series (spring e joins nodes e & e+1), node 1 held at ubar, the
% force at the last node: u(I) = ubar + F*sum(1/ke) over the springs before
% node I, & R1 = -F
Uex = ubar + Fn(1,2)*[0; cumsum(1./Spr(:,3))];
Rex = -Fn(1,2);

% largest error of the displacements, of the largest displacement (one per ubar)
eU = max(abs(U-Uex),[],1)./max(abs(Uex),[],1);
eC = max(abs(Uw-Uex),[],1)./max(abs(Uex),[],1);
eP = max(abs(Up-Uex),[],1)./max(abs(Uex),[],1);

% equilibrium: the reactions & the applied forces add up to zero
eR = abs(sum(R,1) + sum(F))/abs(Fn(1,2));


%% POST-PROCESSING
% Positive extension and force mean tension, with each spring directed i -> j.
extension = U(Spr(:,2),:) - U(Spr(:,1),:);
spring_force = Spr(:,3).*extension;
equilibrium_residual = K*U - F;
equilibrium_residual(c,:) = equilibrium_residual(c,:) - R;
reactions = zeros(size(U));
reactions(c,:) = R;

%% DISPLAY RESULTS
fprintf('Row sums of [K] before the support = %.1e of its largest entry\n',Krow)
for s = 1:numel(ubar)
    fprintf('u1 = %g m:\n',ubar(s))
    fprintf('  u  = %s m (exact %s)\n',mat2str(U(:,s)',8),mat2str(Uex(:,s)',8))
    fprintf('  R1 = %.8g N (exact %g); R1 + sum(F) = %.1e of F\n',R(1,s),Rex,eR(s))
    fprintf('  Spring extensions (m): %s\n',mat2str(extension(:,s)',8))
    fprintf('  Spring forces (N, tension positive): %s\n',mat2str(spring_force(:,s)',8))
    fprintf('  Largest nodal equilibrium residual: %.1e N\n',max(abs(equilibrium_residual(:,s))))
    fprintf('  error: partition %.1e, penalty %.1e, without -K(f,c)*U(c) %.1e\n',eU(s),eP(s),eC(s))
end
fprintf('Condition number: K(f,f) %.2f, with the penalty %.2e\n',condKff,condKp)
fprintf('=> Displacement of node %d = %.8g m (exact %.8g)\n',nN,U(nN,1),Uex(nN,1))

%% PLOTS
nodes = (1:nN)';
elements = (1:size(Spr,1))';
fig = figure('Name','P1 Spring System','Position',[100 100 1000 700]);
support_labels = arrayfun(@(value) sprintf('Support = %g m',value),ubar,'UniformOutput',false);

for s = 1:numel(ubar)
    subplot(2,2,s)
    plot(nodes,Uex(:,s),'k-',nodes,U(:,s),'ro',nodes,Uw(:,s),'bx',nodes,Up(:,s),'g+')
    set(gca,'XTick',nodes,'XLim',[nodes(1)-0.4 nodes(end)+0.4])
    xlabel('Node'); ylabel('Displacement (m)'); grid on
    title(sprintf('Support displacement = %g m',ubar(s)))
    if s == 1
        legend('Closed form','Partition','Omitted support term','Penalty','Location','northwest')
    end
end

subplot(2,2,3)
plot(elements,spring_force(:,1),'bo',elements,spring_force(:,2),'rx')
set(gca,'XTick',elements,'XLim',[elements(1)-0.4 elements(end)+0.4])
xlabel('Spring'); ylabel('Axial force (N)')
title('Spring forces (tension positive)')
legend(support_labels{:},'Location','northwest')

subplot(2,2,4)
bar(nodes - 0.18,F,0.36)
hold on
plot(nodes,reactions(:,1),'bo',nodes,reactions(:,2),'rx')
hold off
set(gca,'XTick',nodes,'XLim',[nodes(1)-0.4 nodes(end)+0.4])
xlabel('Node'); ylabel('Force (N)')
force_limit = 1.1*max([1; abs(F); abs(R(:))]);
ylim([-force_limit force_limit])
title('Applied forces and reactions')
legend('Applied force',support_labels{:},'Location','northwest')

%% OPTIONAL EXPORT
output_dir = getenv('P1_OUTPUT_DIR');
if ~isempty(output_dir)
    [ok,message] = mkdir(output_dir);
    assert(ok,'%s',message)
    nodal_rows = zeros(nN*numel(ubar),9);
    spring_rows = zeros(numel(elements)*numel(ubar),7);
    for s = 1:numel(ubar)
        nodal_rows((s-1)*nN + nodes,:) = ...
            [repmat(ubar(s),nN,1), nodes, U(:,s), Uex(:,s), Up(:,s), Uw(:,s), F, reactions(:,s), equilibrium_residual(:,s)];
        spring_rows((s-1)*numel(elements) + elements,:) = ...
            [repmat(ubar(s),numel(elements),1), elements, Spr, extension(:,s), spring_force(:,s)];
    end

    fid = fopen(fullfile(output_dir,'p1_nodal_matlab.csv'),'w');
    assert(fid >= 0,'Could not open the nodal CSV for writing.')
    fprintf(fid,'support_m,node,partition_m,exact_m,penalty_m,omitted_term_m,applied_N,reaction_N,residual_N\n');
    fprintf(fid,'%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n',nodal_rows.');
    fclose(fid);

    fid = fopen(fullfile(output_dir,'p1_springs_matlab.csv'),'w');
    assert(fid >= 0,'Could not open the spring CSV for writing.')
    fprintf(fid,'support_m,spring,node_i,node_j,stiffness_N_per_m,extension_m,axial_force_N\n');
    fprintf(fid,'%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n',spring_rows.');
    fclose(fid);
    print(fig,fullfile(output_dir,'p1_results_matlab.png'),'-dpng','-r160');
end
