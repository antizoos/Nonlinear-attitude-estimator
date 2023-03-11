% this function is for plotting simulation results conviniently
function []=plot_fig(RR,RR_)
global dt tf N sigmaRR;
phii=zeros(1,N);
for i=1:N
    R_tilde=RR(:,:,i)*RR_(:,:,i)';
    phii(i)=acos((trace(R_tilde)-1)/2)*180/pi;% representation of attitude error
end
% figure('Name','attitude error','Interpreter','latex','fontsize',10);
figure
plot(0:dt:tf,phii,'g','LineWidth',1.4);
% xlabel('$$\int_0^x\!\int_y dF(u,v)$$','Interpreter','latex');
xlabel('time(s)','Interpreter','latex','fontsize',10);
ylabel('$ \vartheta (deg)  $','Interpreter','latex','fontsize',10);
legend('proposed observer','Interpreter','latex');
title('attitude error','Interpreter','latex','fontsize',10);
% % we need to test so we need to plot the sigmaR
% figure
% plot(0:dt:tf,sigmaRR,'b','LineWidth',1.4);
% title('innovation term','Interpreter','latex');
% end