function dx = tamoxifen_model_ODEs(t,x,p,model_choice)


switch model_choice
    case 1
% tom+ cells

%p = [0; p];
p(1:3) = p(1:3).^2;

dx(1) = -p(2)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2);
dx(3) = p(3)*x(2) - p(4)*x(3);

% tom- cells

dx(4) = p(1) - p(2)*x(4);
dx(5) = p(2)*x(4) - p(3)*x(5);
dx(6) = p(3)*x(5) - p(4)*x(6);

    case 2
% tom+ cells

p(1:3) = p(1:3).^2;

dx(1) = -p(2)*x(1) + p(5)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2) + p(6)*x(2);
dx(3) = p(3)*x(2) + p(4)*x(3);

% tom- cells

dx(4) = p(1) - p(2)*x(4) + p(5)*x(4);
dx(5) = p(2)*x(4) - p(3)*x(5) + p(6)*x(5);
dx(6) = p(3)*x(5) + p(4)*x(6);

    case 3
        
       p = [p;0.236235871371677;0.218972744958695;0.0700200409797736];
       p(2:3) = p(2:3).^2;
       dx(1) = p(1) - p(2)*x(1);
       dx(2) = p(2)*x(1) - p(3)*x(2);
       dx(3) = p(3)*x(2) - p(4)*x(3);
       
    case 4
        
        p = p.^2;
dx(1) = p(1) - p(2)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2);
dx(3) = p(3)*x(2) - p(4)*x(3);

    case 5
        
        % tom+ cells

p(1:3) = p(1:3).^2;
%x = [x(1)+x(2) x(3) x(4)+x(5) x(6)];

dx(1) = -p(2)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2);

% tom- cells

dx(3) = p(1) - p(2)*x(3);
dx(4) = p(2)*x(3) - p(3)*x(4);

    case 6
        
        %x = [x(1)+x(2) x(3) x(4)+x(5) x(6)];
        
        % tom+ cells
        p(1:2) = p(1:2).^2;
        p(4) = p(4)^2;

dx(1) = -p(2)*x(1) + p(4)*x(1);
dx(2) = p(2)*x(1) + p(3)*x(2);

% tom- cells

dx(3) = p(1) - p(2)*x(3) + p(4)*x(3);
dx(4) = p(2)*x(3) + p(3)*x(4);

    case 7
        
        dx(1) = -p(2)*x(1) + p(4)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2);

% tom- cells

dx(3) = p(1) - p(2)*x(3) + p(4)*x(3);
dx(4) = p(2)*x(3) - p(3)*x(4);

    case 8
        
        dx(1) = -p(2)*x(1) + p(4)*x(1);
dx(2) = p(2)*x(1) - p(3)*x(2);

% tom- cells

dx(3) = p(1) - p(2)*x(3) + p(4)*x(3);
dx(4) = p(2)*x(3) - p(3)*x(4);
       

end

dx = dx';

end