function quad = irregularSimpson(x,y)

    quad = 0;

    index = 1;

    while index < numel(x)-1
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

    % if there are an odd number of subintervals we have a missing subinterval

    if mod(size(x,1),2) == 0

        % for the last portion we use that our function is periodic and make an
        % extra copy

        newX = [x(index:end);x(1)];
        newY = [y(index:end);y(1)];

        % construct vandermonde
        V = vander(newX);

        % solve for the coefficients of the quadratic
        m = V\newY;

        % this time we use y(i+1) rather than y(i+2) as the final point of the
        % quadratic as this is our endpoint

        quad = quad + m(1)*(y(index+1)^3-y(index)^3);
        quad = quad + m(2)*(y(index+1)^2-y(index)^2);
        quad = quad + m(3)*(y(index+1)-y(index));


    end
