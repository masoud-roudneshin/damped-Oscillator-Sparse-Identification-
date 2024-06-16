function polynomialVector = polynomialBase(x,y, n)
% This function produces all polynomials of x and y 
% until degree n. 
% For instance for n = 3:
% polynomialVector = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3]
% Inputs:
% x and y as numbers (double)
% n: degree of the polynomial
polynomialVector = [];
for i = 0 : n
    for j = 0: i
        k = i - j;
        polynomialVector = [polynomialVector x^j*y^k];
    end
end