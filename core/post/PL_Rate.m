% Title: Code for the observed order of a sequence of errors
% Attempt: to give the rate between one mesh & the next, or the order of an
%          iteration, with no rate where the error is already at the level of
%          rounding, as a kernel
% Created on: 24/11/2026 by Claude (not yet double checked by Razin)
%   on the author's word of 24/11/2026 (use the kernels of core/; where there
%   is none, write one, modular, & keep it in core/; "make modular as much
%   as possible"): the lines of the rates
%   of MAIN_AxialBar & MAIN_SteadyConduction1D (README, open item 22 (d)),
%   with the rule of P3 for the rates that are not given; the order of the
%   Newton iteration of P4 is the same formula on the errors of one step &
%   of the next
%   run in GNU Octave 8.4.0 on 24/11/2026
%   run in MATLAB by Razin, reported on 25/11/2026 (runValidation: 28 of 28 test files passed in 178.1 s)
% Last updated: 24/11/2026 by Claude
% Double checked on: -

function Rate = PL_Rate(e,h,tol)
% UDF for the observed order between one entry of a sequence & the next,
%   Rate = log(e(k)/e(k+1))/log(h(k)/h(k+1))
% e: the errors (or residuals), the sequence along the 2nd dimension (a vector
% is one sequence, taken as a row); h: the element lengths, one per entry of
% that dimension, or a sequence of the size of e; tol: no rate (NaN) where
% e(k) < tol, the error of the coarser entry being already at the level of rounding
% Two uses: the rate of a mesh (h = L./Ne, the element length); the order p of
% an iteration, e(k+1) ~ C*e(k)^p (e = the errors from the 2nd step on, h =
% those from the 1st: then log(e(k)/e(k+1)) over log(e(k-1)/e(k)))

if isvector(e)
    e = e(:)';
    h = h(:)';
end
e1 = e(:,1:end-1,:);
e2 = e(:,2:end,:);
if isequal(size(h),size(e))
    h1 = h(:,1:end-1,:);
    h2 = h(:,2:end,:);
else
    h = h(:)';
    h1 = h(1:end-1);
    h2 = h(2:end);
end
Rate = log(e1./e2)./log(h1./h2);
Rate(e1 < tol) = NaN;

end
