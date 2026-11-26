% Title: Code for the matrix & the vector of a line of 1D elements
% Attempt: to assemble [K] & {F} of -d/dx(k*du/dx) = Q on 1D elements of
%          any order, as a kernel
% Created on: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (the kernels that core/ has not are
%   written, modular, & kept in core/): the loop over the elements that
%   MAIN_AxialBar (k = E*A, Q = q) & MAIN_SteadyConduction1D (k = k*A,
%   Q = Q*A) each wrote (README, open item 22 (b)), which now call this file;
%   the functions of the element are those of SF_FEMLG on one row of nodes,
%   its points those of GQ_Points
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 24/11/2026 by Claude
% Double checked on: -

function [K,F] = KM_Line1D(xi,Elem,Ns,dNs,W,kx,Qx)
% UDF for the matrix [K] & the vector {F} of -d/dx(k(x)*du/dx) = Q(x) on a
% line of 1D elements, by gauss quadrature on each element
%   K(i,j) = sum over the elements of the integral of k*dNi/dx*dNj/dx dx
%   F(i)   = sum over the elements of the integral of Ni*Q dx
% with x = N(s)*Xe, dx/ds = dN/ds*Xe & dN/dx = (dN/ds)/(dx/ds) on an element
% xi: the nodes (x), nN of them; Elem: the nodes of each element from end to
% end, one row per element; Ns, dNs: the functions of the element & their
% derivatives in s (from -1 to 1) at its gauss points, one row per point, with
% their weights W (a column); kx, Qx: k(x) & Q(x), functions of a column of x.
% Nothing of a face is added: a load at a node or a film coefficient is added
% by the caller

xi = xi(:);
W = W(:);
nN = numel(xi);
K = zeros(nN);
F = zeros(nN,1);
for el = 1:size(Elem,1)
    % loop over elements
    en = Elem(el,:);
    Xe = xi(en);       % nodes of the element
    J = dNs*Xe;        % dx/ds at the gauss points (half the length)
    if any(J <= 0)
        error('KM_Line1D: dx/ds <= 0 in element %d (its nodes must run from end to end, x increasing)',el)
    end
    x = Ns*Xe;         % the gauss points on the line
    B = dNs./J;        % dN/dx
    ke = B'*((W.*J.*kx(x)).*B); % integral of k*dNi/dx*dNj/dx
    re = Ns'*(W.*J.*Qx(x));     % integral of Ni*Q
    K = AS_Global(K,ke,en);
    F(en) = F(en) + re;
end

end
