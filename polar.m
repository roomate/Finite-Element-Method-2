function R = polar(x)

r = sqrt(x(1)^2 + x(2)^2);

if (x(1) > 0)
    theta = atan(x(2)/x(1));
elseif (x(1) < 0)
    theta =  pi - atan(x(2)/abs(x(1)));
elseif (x(1) == 0 && x(2) > 0)
    theta = pi/2;
elseif( x(1) == 0 && x(2) < 0)
    theta = 3*pi/2;
end
R = [r , theta];


