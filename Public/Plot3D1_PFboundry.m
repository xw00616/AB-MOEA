function  Plot3D1_PFboundry(truevalue, pf0,pf, marker)
pause(0.01);
plot3(truevalue(:,1), truevalue(:,2), truevalue(:,3), '.');
hold on;
pf1=pf0(1:5,:);
pf2=pf0(6:10,:);
plot3(pf1(:,1), pf1(:,2), pf1(:,3), 'rs','MarkerFace', 'r');
hold on;
plot3(pf2(:,1), pf2(:,2), pf2(:,3), 'kd','MarkerFace', 'k');
hold on;
plot3(pf(1:end-10,1), pf(1:end-10,2), pf(1:end-10,3), 'bo','MarkerFace', 'b');
xm=max(max(pf(:,1)),max(pf0(:,1)));
ym=max(max(pf(:,2)),max(pf0(:,2)));
zm=max(max(pf(:,3)),max(pf0(:,3)));
xlim([0 xm + 0.1]);
ylim([0 ym + 0.1]);
zlim([0 zm + 0.1]);
xlabel('f_1', 'FontSize', 14);ylabel('f_2', 'FontSize', 14);zlabel('f_3', 'FontSize', 14);
view(135, 30);
hold off;
end