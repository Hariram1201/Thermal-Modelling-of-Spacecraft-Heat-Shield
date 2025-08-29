function x = tdm(a,b,c,d)
%  Function to DetermineTri-diagonal Matrix Solution to Determine the
%  Next Temperature Values for the Backward Differencing and
%  Crank-Nicholson Methods
    n = length(b);

    % Eliminate a terms
    for i = 2:n
        factor = a(i) / b(i-1);
        b(i) = b(i) - factor * c(i-1);
        d(i) = d(i) - factor * d(i-1);
    end

    x(n) = d(n) / b(n);

    % Loop backwards to find other x values by back-substitution
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end