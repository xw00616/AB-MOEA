function  Plot3D1(truevalue, pf0,pf, marker)
pause(0.01);
plot3(truevalue(:,1), truevalue(:,2), truevalue(:,3), '.');
hold on;
pf1=pf0(1:5,:);
pf2=pf0(6:10,:);
plot3(pf1(:,1), pf1(:,2), pf1(:,3), marker,'MarkerFace', 'r');
hold on;
plot3(pf2(:,1), pf2(:,2), pf2(:,3), marker,'MarkerFace', 'r');
hole on;
plot3(pf(:,1), pf(:,2), pf(:,3), marker,'MarkerFace', 'r');
xlim([0 max(truevalue(:,1)) + 0.1]);
ylim([0 max(truevalue(:,2)) + 0.1]);
zlim([0 max(truevalue(:,3)) + 0.1]);
xlabel('f_1', 'FontSize', 14);ylabel('f_2', 'FontSize', 14);zlabel('f_3', 'FontSize', 14);
view(135, 30);
hold off;
end