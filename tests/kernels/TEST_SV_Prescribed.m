function TEST_SV_Prescribed()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'core','solver'));
K = [2 -2 0; -2 5 -3; 0 -3 3];
[U,R] = SV_Prescribed(K,[2; 5; 1],[3; 1],[2; -1]);
assert(norm(U-[-1; 1.8; 2],inf)<1e-12)
assert(norm(R-[-0.4; -7.6],inf)<1e-12)

[U,R] = SV_Prescribed(2,4,[],[]);
assert(U==2 && isempty(R))
[U,R] = SV_Prescribed(2,4,1,3);
assert(U==3 && R==2)

rejected = false;
try
    SV_Prescribed(2,4,[1 1],[3 3]);
catch err
    rejected = ~isempty(strfind(err.message,'unique'));
end
assert(rejected,'Duplicate constrained DOFs should be rejected.')
fprintf('TEST_SV_Prescribed passed.\n');
end
