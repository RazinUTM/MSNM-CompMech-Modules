% SPDX-License-Identifier: MIT
function fig = PL_AxialBar(xNodes,U,xElement,N,xExact,Uexact,Nexact,caption)
% Compare the nodal displacements and element forces with the closed form.
fig = figure('Name','P2 Axial Bar','Position',[100 100 800 600]);
subplot(2,1,1)
plot(xExact,Uexact,'k-',xNodes,U,'ro')
xlabel('x (m)'); ylabel('Displacement (m)')
legend('Closed form','FEM at nodes','Location','northwest')
title(caption)
subplot(2,1,2)
exactLine = plot(xExact,Nexact,'k-');
hold on
elementLines = plot(xElement,N,'r-');
hold off
xlabel('x (m)'); ylabel('Axial force (N)')
legend([exactLine; elementLines(1)],'Closed form','FEM by element','Location','northeast')
end
