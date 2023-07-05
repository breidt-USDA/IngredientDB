function pKa = AdjustpKaMonoprotic(pKo, I, TempC)
%FUNCTION pKa = AdjustpKaMonoprotic(pKo, I, TempC)
%
%Returns the modified pKa of a monoprotic acid based on the
%ionic strength and a given temperature in degrees C. 
%See page 97 of Butler and Cogley ("Ionic equilibrium", 1998, 
%Wiley and Sons, Inc., NY). 
%
%Arguments:
%   pKo     The pKo in pure water at 25C or a vector of pKo's
%   I       Ionic strength of the solution 
%   TempC   The temperature in degrees C
%   Z       Vector of ionic charges
%   b       Value of adjustment parameter for Davies equation = 0.3
%-------------------------------------------------------------------------          
b = 0.3;         %0.3 adjustment factor from p97 of Butler and Cogley
Epsilon = 78.3808;      %dialectric constant of water
degK = TempC + 273;     %calc temp in K
A = 1.825 * 10^6 * (Epsilon * degK)^(-3/2); %from p45 of Butler and Cogley
%Note: A = 0.5112 at 25C
temp = (sqrt(I)/(1 + sqrt(I)))-(b*I);   %equation 9 p97 B & C
pKa = pKo - (2*A*temp);                 %equation 9 p97 B & C
end

