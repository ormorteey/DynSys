function Unsteadyload()
% Unsteady representation of the aerodynamic load

clc; %clf('reset');
W = {0.5 0.127 -0.5 1.225 1.716 3.53 0.684 80 30 0.6 1 0.165 0.0455 0.335 0.3};
[span,b,a,rho,MW,MT,Ralpha,Walpha, Wh, Xalpha,C_0,C_1,C_2,C_3, C_4] = deal(W{:});
Ialpha = MW*b^2*(Ralpha)^2; Calpha = 1.4069; Ch = 0.4579; Kh = (Wh^2)*MT; Kalpha_0 = (Walpha^2)*Ialpha; Salpha = MW*Xalpha*b;
figure(5)
hold on;

for U = 0:.1:45
    M_1 = [MT+pi*rho*b.^2 Salpha-a*pi*rho*b.^3 0; Salpha-a*pi*rho*b.^3 Ialpha+pi*(1/8+a^2)*rho*b.^4 0; 0 0 1];

    D_1= [Ch+2*pi*rho*b*U*(C_0-C_1-C_3) (1+(C_0-C_1-C_3)*(1-2*a))*pi*rho*b.^2*U 2*pi*rho*U.^2*b*(C_1*C_2+C_3*C_4); -2*pi*(a+ 1/2)*rho*b.^2*(C_0-C_1-C_3)*U Calpha+(1/2-a)*(1-(C_0-C_1-C_3)*(1+2*a))*pi*rho*b.^3*U -2*pi*rho*b.^2*U.^2*(a+ 1/2)*(C_1*C_2+C_3*C_4); -1/b a-1/2 (C_2+C_4)*(U/b)];
    
    K_1 = [Kh 2*pi*rho*b*U.^2*(C_0-C_1-C_3) 2*pi*rho*U.^3*C_2*C_4*(C_1+C_3); 0 Kalpha_0-2*pi*(0.5+a)*rho*(C_0-C_1-C_3)*b.^2*U.^2 -2*pi*rho*b*U.^3*(a+0.5)*C_2*C_4*(C_1+C_3); 0 -U/b C_2*C_4*(U.^2)/(b.^2)];

    inv_M_1 = inv(M_1); 
    D_1_star = -inv_M_1*D_1; K_1_star = -inv_M_1*K_1;
    B_U = [0 1 0 0 0 0; K_1_star(1,1) D_1_star(1,1) K_1_star(1,2) D_1_star(1,2) K_1_star(1,3) D_1_star(1,3); 0 0 0 1 0 0; K_1_star(2,1) D_1_star(2,1) K_1_star(2,2) D_1_star(2,2) K_1_star(2,3) D_1_star(2,3); 0 0 0 0 0 1; K_1_star(3,1) D_1_star(3,1) K_1_star(3,2) D_1_star(3,2) K_1_star(3,3) D_1_star(3,3)]; 
    e = eig(B_U);
% %###################################################################    
% %     Sorting the eigenvalues in order
%     [sorted, idx] = sort(real(e));
% %###################################################################       
%     e = e(idx);
    plot(U,real(e(5)),'b.');
    xlabel('U'); ylabel('Re {\lambda}_{5}');
    title('Unsteady representation of the aerodynamic load');
end

hold off;

