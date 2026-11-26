% SPDX-License-Identifier: MIT
function K = AS_Global(K,ke,dofs)
% Add one element matrix to K using one-based global DOF indices.
for j = 1:numel(dofs)
    for i = 1:numel(dofs)
        K(dofs(i),dofs(j)) = K(dofs(i),dofs(j)) + ke(i,j);
    end
end
end
