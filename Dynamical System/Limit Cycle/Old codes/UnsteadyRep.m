function Bifurcation()
clc
clear;

%########################################################################################################
%Bifurcation diagrams
% ########################################################################################################
hold on;



for U = 10:.1:15
    [t,x]=ode45(@(t,x)func(t,x,U),[0 600],[0 0 0.01 0 0 0]);
    ind = size(x(:,1));
    index = 3950*floor(ind/4000);
    m_p = max(x(index:ind,1));
    rms_p = m_p*((sqrt(2)/2));
    plot (U,rms_p,'b.');
end
title('Bifurcation diagram for pitch');
hold off;

end

function zdot = func(t,z,U)
W = {0.6 0.0325 -0.5 1.204 1.0662 3.836 0.0004438 0.0115 0.011 0.942 895.10 0.5 1 0.165 0.0455 0.335 0.3};
[span,b,a,rho,MW,MT,Ialpha,Calpha,Ch,Kalpha_0,Kh,Xalpha,C_0,C_1,C_2,C_3, C_4] = deal(W{:});
Kalpha_1 = 3.95;
Kalpha_2 = 107;
Salpha = MW*Xalpha*b;1
M_1 = [MT+pi*rho*b.^2 Salpha-a*pi*rho*b.^3 0; Salpha-a*pi*rho*b.^3 Ialpha+pi*(1/8+a^2)*rho*b.^4 0; 0 0 1];

D_1= [Ch+2*pi*rho*b*U*(C_0-C_1-C_3) (1+(C_0-C_1-C_3)*(1-2*a))*pi*rho*b.^2*U 2*pi*rho*U.^2*b*(C_1*C_2+C_3*C_4); -2*pi*(a+ 1/2)*rho*b.^2*(C_0-C_1-C_3)*U Calpha+(1/2-a)*(1-(C_0-C_1-C_3)*(1+2*a))*pi*rho*b.^3*U -2*pi*rho*b.^2*U.^2*(a+ 1/2)*(C_1*C_2+C_3*C_4); -1/b a-1/2 (C_2+C_4)*(U/b)];

K_1 = [Kh 2*pi*rho*b*U.^2*(C_0-C_1-C_3) 2*pi*rho*U.^3*C_2*C_4*(C_1+C_3); 0 Kalpha_0-2*pi*(0.5+a)*rho*(C_0-C_1-C_3)*b.^2*U.^2 -2*pi*rho*b*U.^3*(a+0.5)*C_2*C_4*(C_1+C_3); 0 -U/b C_2*C_4*(U.^2)/(b.^2)];

inv_M_1 = inv(M_1); 

D_1_star = -inv_M_1*D_1; K_1_star = -inv_M_1*K_1;

zdot(1) = z(2); zdot(3)= z(4); zdot(5)= z(6);
A = D_1_star*[z(2) z(4) z(6)]' + K_1_star*[z(1) z(3) z(5)]' - inv_M_1*[0 Kalpha_1*(z(3)^2) 0]' - inv_M_1* [0 Kalpha_2*(z(3)^3) 0]' ;
zdot(2)=A(1); zdot(4)= A(2); zdot(6)= A(3);
zdot = zdot';
end

