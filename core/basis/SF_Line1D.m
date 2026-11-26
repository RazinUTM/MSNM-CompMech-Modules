% Title: Code for a field of 1D elements at their gauss points
% Attempt: to give the points, their weights, u & du/dx of every element of
%          a line, for the error norms & the figures, as a kernel
% Created on: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (the kernels that core/ has not are
%   written, modular, & kept in core/): the loops over the elements of the
%   error norms & the lines of the figure that MAIN_AxialBar &
%   MAIN_SteadyConduction1D each wrote (README, open item 22 (b)), which now
%   call this file; written for all the elements at once (one column each),
%   where those loops took one element at a time
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 24/11/2026 (runValidation: 25 of 25 test files passed in 134.0 s)
% Last updated: 24/11/2026 by Claude
% Double checked on: -

function [xQ,WxQ,Uq,dUq] = SF_Line1D(xi,Elem,U,Ns,dNs,W)
% UDF for the points [xQ], their weights [WxQ], & a field u of the elements
% with its derivative du/dx ([Uq], [dUq]) at the points s of the element on
% every element of a line
%   x = N(s)*Xe, dx/ds = dN/ds*Xe, u = N(s)*Ue, du/dx = (dN/ds*Ue)/(dx/ds)
% xi: the nodes (x); Elem: the nodes of each element from end to end, one row
% per element; U: the field at the nodes; Ns, dNs: the functions of the
% element & their derivatives in s (from -1 to 1) at the points, one row per
% point, with their weights W (a column; [] for points with no weight, as for
% a figure: WxQ is then []). One column per element, one row per point

xi = xi(:);
U = U(:);
Xe = xi(Elem');    % nodes of the elements, one column each
Ue = U(Elem');     % the field at those nodes
J = dNs*Xe;        % dx/ds at the points
xQ = Ns*Xe;
if isempty(W)
    WxQ = [];
else
    WxQ = W(:).*J;
end
Uq = Ns*Ue;
dUq = (dNs*Ue)./J;

end
