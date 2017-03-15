function xdot = func(t,x)
clc;
W = {0.6 0.0325 -0.5 1.204 1.0662 3.836 0.0004438 0.0115 0.011 0.942 895.10 0.5};
[span,b,a,AD,MW,MT,Ialpha,Calpha,Ch,Kalpha_0,Kh,Xalpha] = deal(W{:});
hold on
'Constant array declaration';
 for U = 0:20   %U =10.9; %input ('Free Stream Velocity: ');
    %Kalpha_1 = 3.95;%input('Enter the quadratic torsional spring coeff: ');
    %Kalpha_2 = 107;%input('Enter the cubic torsional spring coeff: ');
    Salpha = MW*Xalpha*b;
    M = [ MT+(pi*AD*b^2) Salpha-a*pi*AD*b^3; Salpha-a*pi*AD*b^3 Ialpha+pi*((1/8)+a^2)*AD-b^4];
    D = [ Ch + 2*pi*AD*b*U 2*(1 - a)*pi*AD*b^2*U; -2*pi*(a + 0.5)*AD*b^2*U Calpha + a*(2*a)*pi*AD*b^3*U ];
    K = [ Kh 2*pi*AD*b*U^2; 0 Kalpha_0-2*pi*(a + 0.5)*AD*b^2*U^2 ];
    inv_M = inv(M);
    D_star = -inv_M*D; K_star = -inv_M*K;
    B_U = [0 1 0 0; digits (K_star(1,1)) D_star(1,1) K_star(1,2) D_star(1,2); 0 0 0 1; K_star(2,1) D_star(2,1) K_star(2,2) D_star(2,2)]; 
    e = eig(B_U);
    U
    real(e(1))
    plot (U,real(e(1)), 'b.');
 end
 
 hold off


