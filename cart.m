function X = cart(R,theta)
    X(1) = R(1)*cos(theta) - sin(theta)*R(2);
    X(2) = R(1)*sin(theta) + cos(theta)*R(2); 