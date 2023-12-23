function [roundx] = round_JY(x,n)
    roundx = round(x*10^n)/(10^n);
end