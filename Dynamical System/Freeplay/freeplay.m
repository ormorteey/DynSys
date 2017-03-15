function freeplay = func()
B_l = 2.12; B_U = 2.12; P = 0;
% 
 figure(8);
hold on;
% 
%tanh approximation
for e = 0:1:100
    for B = -4:.1:4
        f = 0.5*(1-tanh(e*(B+B_l)))*(B+B_l) + 0.5*(1+tanh(e*(B-B_U)))*(B-B_U)+P;
        plot(B,f, 'r.')
    end
end
d = 2.12;
% 
%freeplay aproximation
for B = -4:.1:4
    if abs(B)<= d
        f=0;
    else
        if B < d
            f=d+B;
        else
            f=B-d;
        end
    end
    plot(B,f,'b.')    
end

%Adapted Paidous constraints w/ k=1 n=1
% for B = -4:.1:4
%     f = (B - 0.5*(abs(B+d)-abs(B-d)));
%     plot(B,f, 'ro')
% end
% for B = -4:.1:4
%     f = B - 0.5*((B+d)-(B-d));
%     plot(B,f, 'r*')
% end

%Initial Paidous constraint
% for B = -4:.1:4
%     if abs(B)<= d
%         f=0;
%     elseif B>=d
%         f=B-d;
%     elseif B<=-d
%         f=B+d;
%     end
%     plot(B,f,'bo') 
% end

%trilinear Paidous constraint
% k=[1,270000,5600000,0,1000000000];
% d = [2.12,0.05,0.044,0,0.031];

% for i = 1
%     for B = -4:.1:4
%         f = (B - 0.5*(abs(B+d)-abs(B-d)))^2;
%         plot(B,f, 'ro')
%     end
% end
% hold off;
% end
