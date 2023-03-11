% Nonlinear Attitude Estimation Using Intermittent and Multi-Rate Vector Measurements
%% initilization
clear all
clc
global ko sigmaR dt tf N sigmaRR;
R0=eye(3);% initial real attitude
v=[0.8;0.6;0];
R_0=eye(3)+sin(pi/2)*CrossMatrix(v)+(1-cos(pi/2))*CrossMatrix(v)^2;
% omega0=2*[0;sin(pi/3);1];% rad/s
f1=10; % measurement vector frequency Hz
f2=20; % Hz
f3=50; % Hz
N1=1000/f1;
N2=1000/f2;
N3=1000/f3;
dt=0.001;% integration step s
t0=0;% initial time
tf=14; % s
N=tf/dt+1;% N is the amount of state points
r1=[sqrt(2)/2;sqrt(2)/2;0]; % reference vector
r2=[sqrt(2)/2;-sqrt(2)/2;0];
r3=[0;0;-1];
% sigma=sqrt(0.08);% Gaussian white noise covariance matrix
sigma=0.08;
rho1=0.2;% innovation parameter
rho2=0.3;
rho3=0.5;
% rho1=0;% innovation parameter
% rho2=0;
% rho3=0;
% rho1=0.03;% innovation parameter
% rho2=0.045;
% rho3=0.075;
ko=15;
kr=0.1;
%% real state integration using 4th runge-kutta in Rotation Matrix
RR=zeros(3,3,N);
RR(:,:,1)=R0;
oomega=zeros(3,N);
tt=zeros(1,N);
tt(1)=t0;
for i=1:N-1
    oomega(:,i)=2*[sin(0.1*tt(i));sin(0.1*tt(i)+pi/3);cos(0.5*tt(i))];
    RR(:,:,i+1)=RR(:,:,i)*expm(CrossMatrix(oomega(:,i))*dt);
    tt(i+1)=tt(i)+dt;
end
%% Hybrid observer
RR_=zeros(3,3,N);
sigmaRR=zeros(1,N);
RR_(:,:,1)=R_0;
% RR_(:,:,1)=R0;
r1_=r1;
r2_=r2;
r3_=r3;
sigmaR=rho1*CrossMatrix(r1_)*r1+rho2*CrossMatrix(r2_)*r2+rho3*CrossMatrix(r3_)*r3;
sigmaRR(:,1)=norm(sigmaR);
for i=1:N-1
   % calculate innovation term
    sigmaR=rho1*CrossMatrix(r1_)*r1+rho2*CrossMatrix(r2_)*r2+rho3*CrossMatrix(r3_)*r3;
    sigmaRR(:,i+1)=norm(sigmaR);
    % generate vector measurements
    if mod(i-1,N1)==0
        b1=RR(:,:,i)' *r1+normrnd(zeros(3,1),sigma);% add Gaussian noise with standard deviation sigma
        b1=b1/sqrt(b1'*b1);
%         b1=RR(:,:,i)'*r1;
        
        r1_=r1_+kr*(RR_(:,:,i)*b1-r1_);
    else
        [r1_,~]=runge_kutta(@dr,r1_,tt(i),dt,dt);
    end
    if mod(i-1,N2)==0
        b2=RR(:,:,i)'*r2+normrnd(zeros(3,1),sigma);
        b2=b2/sqrt(b2'*b2);
%         b2=RR(:,:,i)'*r2;
        r2_=r2_+kr*(RR_(:,:,i)*b2-r2_);
    else
        [r2_,~]=runge_kutta(@dr,r2_,tt(i),dt,dt);
    end
    if mod(i-1,N3)==0
        b3=RR(:,:,i)'*r3+normrnd(zeros(3,1),sigma);
        b3=b3/sqrt(b3'*b3);
%         b3=RR(:,:,i)'*r3;
        r3_=r3_+kr*(RR_(:,:,i)*b3-r3_);
    else
        [r3_,~]=runge_kutta(@dr,r3_,tt(i),dt,dt);
    end
    
    
    
    
    RR_(:,:,i+1)=RR_(:,:,i)*expm(CrossMatrix((oomega(:,i)+ko*RR_(:,:,i)'*sigmaR)*dt));
end
%% plot figure
plot_fig(RR,RR_);
% phii=zeros(1,N);
% for i=1:N
%     R_tilde=RR(:,:,i)*RR_(:,:,i)';
%     phii(i)=acos((trace(R_tilde)-1)/2)*180/pi;% representation of attitude error
% end
% figure('Name','attitude error');
% plot(0:dt:tf,phii,'g','LineWidth',1.4);
% xlabel('time(s)','Interpreter','latex','fontsize',10);
% ylabel('$ \vartheta (deg)  $','Interpreter','latex','fontsize',10);
% legend('proposed observer')
% % we need to test so we need to plot the sigmaR
% figure('Name','innovation term');
% plot(0:dt:tf,sigmaRR,'b','LineWidth',1.4);