function polynomialVector = polynomialBase_3States(x,y, z, n)
% This function produces all polynomials of x and y and z
% until degree n. 
% For instance for n = 2:
% polynomialVector = [1 x y z x^2 y^2 z^2 x*y x*z y*z ]
% Inputs:
% x and y and z as numbers (double)
% n: degree of the polynomial
polynomialVector = [];
counter = 0;
for i = 0 : n
    for j = 0:i %0: i
        for k = 0: i
            m = i - j - k;
            if m >= 0
                counter = counter + 1;
                %polynomialVector = [polynomialVector j k m];
                polynomialVector = [polynomialVector x^j*y^k*z^m];
            end
        end
    end
end

