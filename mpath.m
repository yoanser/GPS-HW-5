function [MP1,CMC1] = mpath(firstPseudo, firstCarrier, f1, secondCarrier, f2)
%{
AUTHOR: Yoan Serafimog


INPUTS:
firstPseudo
firstCarrier
f1
secondCarrier
f2


OUTPUTS:
MP1 - Multipath observable
CMC1 - Code minus carrier
%}

c = 299792458; %m/s

lambda1 = c./f1;
lambda2 = c./f2;

MP1 = firstPseudo - ((f1^2 + f2^2)./(f1^2 - f2^2)).*firstCarrier.*lambda1  + ((2*f2^2)./(f1^2 - f2^2)).*secondCarrier.*lambda2;

CMC1 = firstPseudo - firstCarrier.*lambda1;

end 

