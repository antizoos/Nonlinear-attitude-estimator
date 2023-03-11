R0=eye(3);
R_0=eye(3)+sin(pi/2)*CrossMatrix(v)+(1-cos(pi/2))*CrossMatrix(v)^2;
rho1=0.2;% innovation parameter
rho2=0.3;
rho3=0.5;
rhoo=[rho1,rho2,rho3];
initial_val={R0,R_0,rhoo};
A=initial_val{1}