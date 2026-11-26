% Title: 2D SF Code considering Lagrange interpolation
% Attempt: to test Lagrange SF & 1st derivatives (fitting)
% Created on: 05/07/2021 by Razin
% Last updated: 12/11/2026 by Claude (not yet double checked by Razin)
%   (1) 2nd derivatives dNNxx, dNNxy, dNNyy added; NN, dNNx, dNNy unchanged
% Last updated: 18/11/2026 by Claude (not yet double checked by Razin)
%   (2) the 1D functions are worked out for all the points together (A holds
%       one row for each point) & every output is sized before it is filled;
%       the terms of each product & of each sum are those of before, in the
%       same order (the six outputs are the same in every figure in GNU Octave
%       8.4.0; 31 s became 0.09 s there for 11 x 11 nodes & 900 points).
%       Run in MATLAB by Razin, reported on 18/11/2026 (runValidation:
%       12 of 12 test files passed)
% Last updated: 22/11/2026 by Claude (not yet double checked by Razin)
%   (3) a note, no line of code changed. This is ONE Lagrange element through
%       all the nodes of a direction. Through equally spaced nodes its
%       functions grow at the edges (the largest |N| at the gauss points of a
%       plate is 28, 201, 1660 & 16000 for 11, 13, 15 & 17 nodes a side; along
%       an edge 11, 18, 99 & 557 for 12, 13, 16 & 19 nodes), so figures are lost
%       beyond about 12 nodes in a direction: in the plate buckling 13 nodes
%       still give k to 4 figures,
%       & 14 x 14 nodes give NaN (README, Findings 43e). No number of nodes
%       is an exact limit, so there is no warning: keep to 12 or fewer
% Last updated: 23/11/2026 by Claude (not yet double checked by Razin)
%   (4) an error where the nodes are not a full grid, each node once: the
%       columns are the nodes of the grid unique(x) x unique(y), so on other
%       nodes the columns were not the nodes & nothing said so; a note under
%       the function line says in which order the columns come
%   run in MATLAB R2026a on 23/11/2026 by Claude, on the author's computer (runValidation: 19 of 19 test files passed in 127.0 s)
%   run in MATLAB by Razin, reported on 23/11/2026 (runValidation: 21 of 21 test files passed in 167.7 s)
% Double checked on: 05/07/2021 by Razin (before the 12/11/2026 update)

function [NN,dNNx,dNNy,dNNxx,dNNxy,dNNyy] = SF_FEMLG(points,nodes)
% rows are points; columns are the nodes of the grid unique(x) x unique(y),
% each sorted, x running first: column k is node k only when the nodes are
% in that order (sortrows(nodes,[2 1]))

if numel(unique(nodes(:,1)))*numel(unique(nodes(:,2))) ~= size(nodes,1) || ...
        size(unique(nodes,'rows'),1) ~= size(nodes,1)
    error('SF_FEMLG:nodes','SF_FEMLG: the nodes must be a full grid: %d different x & %d different y make %d nodes, these are %d', ...
        numel(unique(nodes(:,1))),numel(unique(nodes(:,2))),numel(unique(nodes(:,1)))*numel(unique(nodes(:,2))),size(nodes,1));
end

% calculate SFunc for 1 element
% in 1D (x-direction)
Xp = points(:,1);
Xm = unique(nodes(:,1))';
Xj = Xm';
A = Xp-Xm;    %numerator (one row for each point of interest)
B = Xj-Xm;   %denominator
B(logical(eye(size(B)))) = 1;
Nx = zeros(length(Xp),length(Xm));
dNx = zeros(length(Xp),length(Xm));
ddNx = zeros(length(Xp),length(Xm));
for j = 1:length(Xm) %loops over the nodes of the direction
    dA = A;
    dA(:,j) = 1;
    Nx(:,j) = prod(dA./B(j,:),2);

    %derivative over-x
    for k = 1:length(Xm) %1st term taken out of the product
        if k~=j
            dA = A;
            dA(:,[j k]) = 1;
            dNx(:,j) = dNx(:,j) + prod(dA./B(j,:),2);
        end
    end

    %2nd derivative over-x
    for k = 1:length(Xm) %2nd term taken out of the product
        if k~=j
            ddNxjk = zeros(length(Xp),1);
            for l = 1:length(Xm)
                if l~=j && l~=k
                    dA = A;
                    dA(:,[j k l]) = 1;
                    ddNxjk = ddNxjk + prod(dA./B(j,:),2);
                end
            end
            ddNx(:,j) = ddNx(:,j) + ddNxjk;
        end
    end
end

% in 1D (y-direction)
Yp = points(:,2);
Ym = unique(nodes(:,2))';
Yj = Ym';
A = Yp-Ym;    %numerator (one row for each point of interest)
B = Yj-Ym;   %denominator
B(logical(eye(size(B)))) = 1;
Ny = zeros(length(Yp),length(Ym));
dNy = zeros(length(Yp),length(Ym));
ddNy = zeros(length(Yp),length(Ym));
for j = 1:length(Ym) %loops over the nodes of the direction
    dA = A;
    dA(:,j) = 1;
    Ny(:,j) = prod(dA./B(j,:),2);

    %derivative over-y
    for k = 1:length(Ym) %1st term taken out of the product
        if k~=j
            dA = A;
            dA(:,[j k]) = 1;
            dNy(:,j) = dNy(:,j) + prod(dA./B(j,:),2);
        end
    end

    %2nd derivative over-y
    for k = 1:length(Ym) %2nd term taken out of the product
        if k~=j
            ddNyjk = zeros(length(Yp),1);
            for l = 1:length(Ym)
                if l~=j && l~=k
                    dA = A;
                    dA(:,[j k l]) = 1;
                    ddNyjk = ddNyjk + prod(dA./B(j,:),2);
                end
            end
            ddNy(:,j) = ddNy(:,j) + ddNyjk;
        end
    end
end

% in 2D (x-y directions)
k = 0;
NN = zeros(size(points,1),length(Xm)*length(Ym));
dNNx = zeros(size(points,1),length(Xm)*length(Ym));
dNNy = zeros(size(points,1),length(Xm)*length(Ym));
dNNxx = zeros(size(points,1),length(Xm)*length(Ym));
dNNxy = zeros(size(points,1),length(Xm)*length(Ym));
dNNyy = zeros(size(points,1),length(Xm)*length(Ym));
for j = 1:length(Ym)
    for i = 1:length(Xm)
        k = k + 1;
        NN(:,k) = Nx(:,i).*Ny(:,j);
        dNNx(:,k) = dNx(:,i).*Ny(:,j);
        dNNy(:,k) = Nx(:,i).*dNy(:,j);
        dNNxx(:,k) = ddNx(:,i).*Ny(:,j);
        dNNxy(:,k) = dNx(:,i).*dNy(:,j);
        dNNyy(:,k) = Nx(:,i).*ddNy(:,j);
    end
end

end