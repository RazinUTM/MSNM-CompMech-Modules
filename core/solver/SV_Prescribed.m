% SPDX-License-Identifier: MIT
function [U,R] = SV_Prescribed(K,F,c,values)
% Solve one load case with prescribed degrees of freedom by partitioning.
% K is square. F is the load vector. c contains unique one-based DOFs.
% values and the returned reactions R follow the order in c.

n = size(K,1);
assert(ismatrix(K) && size(K,2)==n && (isvector(F) || isempty(F)) && numel(F)==n, ...
    'K must be square and F must match its size.')
assert((isvector(c) || isempty(c)) && (isvector(values) || isempty(values)) && numel(c)==numel(values), ...
    'Provide one prescribed value per constrained DOF.')
c = c(:);
F = F(:);
assert(isreal(c) && all(c>=1 & c<=n & c==fix(c)) && numel(unique(c))==numel(c), ...
    'Constrained DOFs must be unique, valid integer indices.')

f = setdiff((1:n)',c);
U = zeros(n,1);
U(c) = values(:);
if ~isempty(f)
    U(f) = K(f,f)\(F(f) - K(f,c)*U(c));
end
R = K(c,:)*U - F(c);
end
