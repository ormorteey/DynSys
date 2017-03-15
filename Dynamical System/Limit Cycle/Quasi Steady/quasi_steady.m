function xdot = func(t,x)
clc; clear;
W = {0.6 0.0325 -0.5 1.204 1.0662 3.836 0.0004438 0.0115 0.011 0.942 895.10 0.5};
[span,b,a,rho,MW,MT,Ialpha,Calpha,Ch,Kalpha_0,Kh,Xalpha] = deal(W{:});
figure(5)
hold on
for U = 0:0.1:15
    Salpha = MW*Xalpha*b;
    M = [MT+pi*rho*b.^2 Salpha-a*pi*rho*b.^3; Salpha-a*pi*rho*b.^3 Ialpha+pi*(1/8+a^2)*rho*b.^4];
    D = [Ch+2*pi*rho*b*U 2*(1 - a)*pi*rho*b.^2*U; -2*pi*(a + 0.5)*rho*b.^2*U Calpha+a*(2*a-1)*pi*rho*b.^3*U];
    K = [Kh 2*pi*rho*b*U.^2; 0 Kalpha_0-2*pi*(0.5+a)*rho*b.^2*U.^2];
    inv_M = inv(M);
    D_star = -inv_M*D; K_star = -inv_M*K;
    B_U = [0 1 0 0; K_star(1,1) D_star(1,1) K_star(1,2) D_star(1,2); 0 0 0 1; K_star(2,1) D_star(2,1) K_star(2,2) D_star(2,2)]; 
    e = eig(B_U);
    [sorted, idx] = sort(real(e));
    e = e(idx);
    
    plot(U,real(e(3)),'b.');
    xlabel('U'); ylabel('Re {\lambda}_{3}');
    %plot(U,imag(e(1)),'b.');
    %xlabel('U'); ylabel('Im {\lambda}_{3}');  
end
hold off
% Kalpha_1 = 3.95;
% Kalpha_2 = 107;
% xdot(1) = x(2); xdot(3)= x(4);
% A = D_star*[x(2) x(4)]' + K_star*[x(1) x(3)]' - inv_M*[0 Kalpha_1*(x(3))^2]' - inv_M* [0 Kalpha_2*(x(3))^3]' ;
% xdot(2)=A(1); xdot(4)= A(2); 
% xdot = xdot';