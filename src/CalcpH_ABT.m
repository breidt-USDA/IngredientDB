function pH = CalcpH_ABT(ABtable, IS, AdjC)
%FB Version 03062016
%PURPOSE 
%   Determine the pH for a given mixture of weak organic acid in solution
%NOTE
%   Function uses Newton's method to minimize equation 1 (see below)
%PARAMETERS
%   ABmatrix = struct array (Nx3) of:  conc (in moles/L), pKa values
%       Table arrays are essentially "struct arrays"?
%       a/b list values are the third col, with chars 'a' or 'b'. 
%       NOTE: only one acid or base (conc, pKa pair) is needed.
%       NOTE: if char is 'x' or other char the buffer is not calculated...
%   PercentNaCl = scalor, percent sodium chloride (g/100 ml)
%   AdjC = scalor, M/L of HCl (positive value) or M/L NaOH (mult. by -1)
%INTERNAL PARAMS
%   H = scalor, initial guess for the algorithm (pH 0..14)
%   tolerance = scalor, minimum value for approach to zero (eg. 10^-12)
%   maxevals = scalor, maximum number of function iterations (eg. 100)   
%NOTE: AdjC could also be the sum of salt concentrations for weak acids
%       or bases, essentially equal to cation or anion conc. i.e. 
%       [Na+Ace-] or [NH4+Cl-]. If [Na+] mult. conc. (M/L) by -1
%       if Cl- just use the concentration.
%EXAMPLE CALL 
%   pH = CalcpH_AB_Opt(ABtable,0.342,0) %for 2% (IS=0.342), no adjustment
% where 
% ABtable = (N x 3) of (conc, pKa, and a/b chars):
%     0.0230    2.6290      a
%     0.0039    4.0207      a
%     0.0241    4.8453      a
%     0.0055    6.9396      a
%     0.0133    9.6388      b
%     0.0503   12.0759      b
%RETURNS 
%   pHval = the predicted pH for the acid/base mixture
%   Note: pH = -1 if error
%Equations:
%     Concentrations:  
%     Ca = [HA] + [A-]
%     Cb = [B] + [BH+]  
%     C1 = [NaOH] i.e. Na+
%     C2 = [HCl]  i.e. Cl-    
%     Charge balance equation:
%     [Na+] +[BH+] + [H+] = [A-] + [OH-] + [Cl-]
%                                                                SALT or
%       ACID TERM            BASE TERM           OH       H    HCl or NaOH
%Equation 1:  f(H) = 0: 
%  0 = CaKa/(Ka + [H+]) - Cb[H+]/([H+] + Kb) + Kw/[H+] - [H+] + (C2 - C1)
%   OR:   0 = sum(acid) - sum(base) + Kw/[H+] - [H+] + adjC
%Equation 2:  f'(H) = 0, ACID TERM, BASE TERM = -CK/(H+K)^2
%   SO:   0 = sum(-CK/([H+]+K)^2) - Kw/[H+]^2 - 1
%   Using Newton's method with 1st derivative:
%      [H+]next = [H+] - f([H+])/f'([H+]), where the function result for
%      f([H+]) is minimized by repeated calls with [H+]next. 
%Reference: Butler and Cogely, 1998. Ionic Equilibrium, Solubility and pH 
%   Calculations. John Wiley and Sons, Inc. NY. 
%-----------------------------------------------------------------------

%initialize
H = 10^(-1*14);                             %startval (fixed at pH 14)
tolerance = 10^-15;                         %apx zero at 10^-15
maxevals = 1000;                            %Maximum of 1000 evaluations
%IS = NaClpercent/(5.84);                    %convert percent to IS
Kw = 1e-14;                                 %water equilib. constant
pKvals = ABtable{:,2};                      %get pK vals
conc = ABtable{:,1};                        %get conc in M
ab_chars = char(ABtable{:,3});              %get acid/base char list
ABmatrix = [conc pKvals];                   %make ab matrix (conc, pk)
ABmatrix(:,2) = adjpka(ABmatrix(:,2),IS);   %adjust pKa values for IS
[len,~] = size(ABmatrix);                   %number of conc,AB pairs
pH = -1;                                    %initialize pH to error code

%optimize
for itr = 0:maxevals                        %loop till maxevals (or tol)
   Fx = GetFx(len, ABmatrix, Kw, H, AdjC, ab_chars);  %eval current H: f(H)
   Fp = GetFp(len,ABmatrix,Kw,H);           %eval current H: f'(H)
   H = H - (Fx/Fp);                         %get 1/2 next [H+] value
   if ((abs(Fx)) < tolerance)               %see if tol is reached, H pos
      break;                                %if so breakout      
   end
end

%record pH only if optimum value found, otherwise, pH = -1
if (H >= 0) && (itr < maxevals)                             
   pH = -1 * log10(H);                      %record predicted pH   
end

end

function zval = GetFx(len, ConcpK, Kw, H, AdjC, ab_chars)
%Evaluate f(H), function to minimize (see Equation 1 above)
%Args:
%   len = scalor (number of buffers)
%   ConcpK = Nx2 (conc, pK matrix)
%   Kw val = scalor (10^-14)
%   H = scalor [H+] concentration
%   AdjC = scalor (adjustment, see above)
%Returns:
%   zval = scalor (zero when [H+] = pH of solution)
%-----------------------------------------------------------------------
Fxvec = zeros(1,len);                       %1xN for each conc, pK result
ConcK = ConcpK;                             %ConcK for pK to K
ConcK(:,2) = 10.^(-1*ConcpK(:,2));          %convert pKa to Ka
for i=1:len                                 %for each conc, pK pair...
    if ab_chars(i) == 'b'                      %check for base
        Fxvec(i) = -1*(ConcK(i,1)*H)/(H+ConcK(i,2)); %prooess base
    elseif ab_chars(i) == 'a'               %else if it is an acid
        Fxvec(i) = (ConcK(i,1)*ConcK(i,2))/(H+ConcK(i,2)); %process acid 
    end  
end
zval = sum(Fxvec) + (Kw/H) - H  + AdjC;     %sum equation 1 (see above)
end 

function zval = GetFp(len, ConcpK, Kw, H)
%Evaluate f'(H), first derivative of f(H)
%Args:
%   len = scalor (number of buffers)
%   ConcpK = Nx2 (conc, pK matrix)
%   Kw val = scalor (10^-14)
%   H = scalor [H+] concentration
%Returns:
%   rval = scalor (evaluation of f'(H), slope of tangent line)
%NOTE: acid or base have same derivative, so no need to check...
%-----------------------------------------------------------------------
Fpvec = zeros(1,len);                       %1xN for each conc, pK result
ConcK = ConcpK;                             %ConcK for pK to K
ConcK(:,2) = 10.^(-1*ConcpK(:,2));          %Convernt pKa to Ka
for j=1:len                                 %for each conc, pK pair
    Fpvec(j) = -1*((ConcK(j,1)*ConcK(j,2))/((ConcK(j,2) + H)^2)); %eval  
end                                    
zval = sum(Fpvec) - (Kw/(H^2)) - 1;         %sum equation 2 (see above)
end 

function pKa = adjpka(pKo, IS)
%Modify pKa of a monoprotic buffer due to IS or temp. 
%See page 97 of Butler and Cogley ("Ionic equilibrium", 1998, 
%Wiley and Sons, Inc., NY). 
%Args:
%   pKo = scalor (the pKo in pure water at 25C or a vector of pKo's)
%   IS = scakir (ionic strength of the solution)
%Returns:
%   pKa = scalor (adjusted pka value)
%NOTE: we found that temperature does not significantly affect pKa
%by printing IS, temperature, and predicted pKa graph, so temperature 
%was not included.
%-----------------------------------------------------------------------          
b = 0.3;                                    %from p97 of Butler and Cogley
A = 0.51;                                   %from p45 of B & C
SqrIS = sqrt(IS);                           %Square root of IS
temp = (SqrIS/(1 + SqrIS))-b*IS;            %equation 9 p97 B & C
pKa = pKo - (2*A*temp);                     %equation 9 p97 B & C
end


