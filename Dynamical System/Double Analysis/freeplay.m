function f = freeplay(n,d,B)

Kalpla_0 = 0.942;
Kalpla_1 = 3.95;
Kalpla_2 = 107;

% freeplay
if n == 1
    if abs(B)<= d
        f=0;
    else
        if B < -d
            f=(d+B);
        else
            f=(B-d);
        end
    end   
end

% linear
if n == 2
    if abs(B)<= d
        f=0;
    else
        if B < -d
            f=(d+B)*(Kalpla_0);
        else
            f=Kalpla_0*(B-d);
        end
    end   
end



% % cubic nonlinear

if n == 3
    if abs(B)<= d
        f=0;
    else
        if B < -d
            f=(Kalpla_0+Kalpla_1*B+Kalpla_2*B^2)*(d+B);
        else
            f=(Kalpla_0+Kalpla_1*B+Kalpla_2*B^2)*(B-d);
        end
    end   
end

end
