% 4th order runge kutta method
function [y,t]=runge_kutta(dif,y0,t0,dt,tf)% dif is the differential function y0 is the initial state vector dt is the integration step tf is the integration length
% y is the final state vector t is the time serial
N=floor(tf/dt)+1;
t=zeros(1,N);
y=y0;
t(1)=t0;
for i=1:N   
    yk=y;
    k1=dif(yk);
    k2=dif(yk+k1*dt/2);
    k3=dif(yk+k2*dt/2);
    k4=dif(yk+k3*dt);
    y=yk+dt*(k1/6+k2/3+k3/3+k4/6);
    y=y/sqrt(y'*y);
    t(i+1)=t(i)+dt;
end
