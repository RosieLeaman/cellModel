% THIS FUNCTION TAKES IN TWO COLUMN VECTORS!!!

function quad = irregularSimpson(x,y)

    quad = 0;

    index = 1;
    
    % make the number of subintervals even if numel(x) is even
    
    if mod(numel(x),2) == 0
        
        % for the first four points, in the case where we have an odd
        % number of subintervals, we use the 3/8 rule which is solving a
        % cubic, not a quadratic
        
        V = vander(x(index:index+3));
        
        m = V\y(index:index+3);
        
        quad = quad + m(1)*(x(index+3)^4-x(index)^4)/4;
        quad = quad + m(2)*(x(index+3)^3-x(index)^3)/3;
        quad = quad + m(3)*(x(index+3)^2-x(index)^2)/2;
        quad = quad + m(4)*(x(index+3)-x(index));

        index = index + 3;
        
        % we now have an even number of subintervals remaining, so can
        % merrily go on our way as if we had an even number of subintervals
        
        % this probably makes our error term slightly worse (and does so in
        % tests goes from 10^-15 to 10^-11)
    end

    while index < numel(x)
        % until we hit the end we just composite simpson it

        % construct vandermonde
        V = vander(x(index:index+2));

        % solve for the coefficients of the quadratic
        m = V\y(index:index+2);

        quad = quad + m(1)*(x(index+2)^3-x(index)^3)/3;
        quad = quad + m(2)*(x(index+2)^2-x(index)^2)/2;
        quad = quad + m(3)*(x(index+2)-x(index));

        index = index + 2;
    end

end
