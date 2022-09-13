function out = cr3bp_pseudoPotential(c3sys, xyz)
%CR3BP_PSEUDOPOTENTIAL Finds Lagrange points given system.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure containing at least:  c3sys.mu
%                                where mu is the characteristic mass ratio
%       2. xyz        [3x1]      Position Vector [x; y; z]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: structure 'out' containing fields:
%
%       1. .U         [1x1]      Pseudo-Potential Function
%       2. .jU        [3x1]      1st Partial of Pseudo-Potential Function
%       3. .jUjU      [6x6]      2nd Partial of Pseudo-Potential Function
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % System Properties
    out = struct;
    mu  = c3sys.mu;
    x1  = -mu;
    y1  = 0.0;
    z1  = 0.0;
    x2  = 1-mu;
    y2  = 0.0;
    z2  = 0.0;

    % State Properties
    x  = xyz(1);
    y  = xyz(2);
    z  = xyz(3);
    r1 = sqrt((x-x1)^2 + (y-y1)^2 + (z-z1)^2);
    r2 = sqrt((x-x2)^2 + (y-y2)^2 + (z-z2)^2);



    % Pseudo-Potential Function
    out.U = ((1/2)*(x^2 + y^2)) + ((1-mu)/r1) + (mu/r2);



    % Pseudo-Potential First Partial w.r.t x,y,z
    jUjx = x - (((1-mu)*(x+mu))/(r1^3)) - ((mu*(x-1+mu))/(r2^3));
    jUjy = y - (((1-mu)*y)/(r1^3))      - ((mu*y)/(r2^3));
    jUjz =   - (((1-mu)*z)/(r1^3))      - ((mu*z)/(r2^3));

    out.jU   = [jUjx; jUjy; jUjz];



    % Pseudo-Potential Second Partial Matrix w.r.t x,y,z
    Uxx = 1 - ((1-mu)/(r1^3)) - (mu/(r2^3)) + ...
        ((3*(1-mu)*((x+mu)^2))/(r1^5)) + ((3*(mu)*((x-1+mu)^2))/(r2^5));
    
    Uxy = ((3*(1-mu)*(x+mu)*y)/(r1^5)) + ((3*(mu)*(x-1+mu)*y)/(r2^5));
    Uyx = Uxy;
    
    Uxz = ((3*(1-mu)*(x+mu)*z)/(r1^5)) + ((3*(mu)*(x-1+mu)*z)/(r2^5));
    Uzx = Uxz;
    
    Uyy = 1 + ((3*(1-mu)*(y^2))/(r1^5)) - ((1-mu)/(r1^3)) - ((mu)/(r2^3)) ...
        + ((3*mu*(y^2))/(r2^5));
    
    Uyz = ((3*(1-mu)*y*z)/(r1^5)) + ((3*mu*y*z)/(r1^3));
    Uzy = Uyz;
    
    Uzz = ((3*(1-mu)*(z^2))/(r1^5)) - ((1-mu)/(r1^3)) + ((3*mu*(z^2))/(r2^5)) ...
        - (mu/(r2^3));
    
    out.jUjU = [Uxx Uxy Uxz;
            Uyx Uyy Uyz;
            Uzx Uzy Uzz];

    

    % Checking work with MATLAB Partials Computation
    %{
        clear x y z mu r1 r2;
        syms  x y z mu;
        
        
        r1 = sqrt((x-x1)^2 + (y-y1)^2 + (z-z1)^2);
        r2 = sqrt((x-x2)^2 + (y-y2)^2 + (z-z2)^2);
        
        jUjx = x - (((1-mu)*(x+mu))/(r1^3)) - ((mu*(x-1+mu))/(r2^3));
        jUjy = y - (((1-mu)*y)/(r1^3))      - ((mu*y)/(r2^3));
        jUjz =   - (((1-mu)*z)/(r1^3))      - ((mu*z)/(r2^3));
        
        Uxx = (diff(jUjx, x));
        Uxy = (diff(jUjx, y));    
        Uxz = (diff(jUjx, z));    
        Uyy = (diff(jUjy, y));
        Uyz = (diff(jUjy, z));
        Uzz = (diff(jUjz, z));
        
        Uxx = subs(Uxx,{x,y,z,mu},{xt,yt,zt,mut});
        Uxy = subs(Uxy,{x,y,z,mu},{xt,yt,zt,mut});
        Uyx = Uxy;
        Uxz = subs(Uxz,{x,y,z,mu},{xt,yt,zt,mut});
        Uzx = Uxz;
        
        Uyy = subs(Uyy,{x,y,z,mu},{xt,yt,zt,mut});
        Uyz = subs(Uyz,{x,y,z,mu},{xt,yt,zt,mut});
        Uzy = Uyz;
        
        Uzz = subs(Uzz,{x,y,z,mu},{xt,yt,zt,mut});
        
        
        U = double([Uxx Uxy Uxz;
             Uyx Uyy Uyz;
             Uzx Uzy Uzz]);
        
        U = round(U,8);
    %}

end