function  Plot3D(truevalue, pf, marker)
pause(0.01);
plot3(truevalue(:,1), truevalue(:,2), truevalue(:,3), '.');
hold on;
plot3(pf(:,1), pf(:,2), pf(:,3), marker,'MarkerFace', 'r');
xlim([0 max(truevalue(:,1)) + 0.1]);
ylim([0 max(truevalue(:,2)) + 0.1]);
zlim([0 max(truevalue(:,3)) + 0.1]);
% xlim([0 max(pf(:,1)) + 0.1]);
% ylim([0 max(pf(:,2)) + 0.1]);
% zlim([0 max(pf(:,3)) + 0.1]);
xlabel('f_1', 'FontSize', 14);ylabel('f_2', 'FontSize', 14);zlabel('f_3', 'FontSize', 14);
view(135, 30);
hold off;
end