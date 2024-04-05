function coef = CalcPolyCoeffs(i)
% Coefficients for order 6 and order 12 terms in Clayff potentials

    syms x y;  % Declare x and y as symbolic variables

    equation = ((x+y)/2)^i;

    % Expand the equation and retrieve the coefficients
    expanded_equation = expand(equation);

    % Get the coefficients
    coef = coeffs(expanded_equation, [x, y]);
    
    coef = double(coef);
end
