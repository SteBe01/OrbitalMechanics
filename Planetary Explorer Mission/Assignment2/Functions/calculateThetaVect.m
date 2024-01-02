function [theta_vect] = calculateThetaVect(mu, a, e, T_size)

% Calculates the correct theta_vect with time relations
% (correct theta_vect can be used to generate time related orbital points)
% The funcion also incorporates the functions Newton and vectorRectifier
%
% Usage:
% [theta_vect] = calculateThetaVect(mu, a, e, T_size)
%
% Input arguments:
% ----------------------------------------------------------------
% mu            [1x1]                   gravitational parameters        [Km^3/s^2]
% a             [1x1]                   semi-major axis                 [Km]
% e             [1x1]                   eccentricity                    [-]
% T_size        [1x1]                   lenght of theta_vect            [-]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% theta_vect    [T_sizex1]/[1xT_size]   theta vect                      [s]


    T_tot=2*pi*sqrt(a^3/mu);
    T_vect=linspace(0,T_tot,T_size);
    
    theta_vect=[];
    for k=1:length(T_vect)
        t=T_vect(k);
        tp=T_tot/2;
        
        M=sqrt(mu/(a^3))*(t-tp);
    
        f=@(x) -M+x-e*sin(x);
        df=@(x) 1-e.*cos(x);
    
        [xvect,~]=newton(2,100,1e-8,f,df);
        E=wrapTo2Pi(xvect(end));
    
        theta = 2*atan(tan(E/2)/sqrt((1-e)/(1+e)));
    
        theta_vect=[theta_vect; theta];
    end
    theta_vect=vectorRectifier(theta_vect);
end


% other functions

function [xvect,it]=newton(x0,nmax,toll,fun,dfun, mol)

    if (nargin == 5)
        mol = 1;
    end
    
    err = toll + 1;
    it = 0;
    xvect = x0;
    xv = x0;
    
    while (it< nmax && err> toll)
       dfx = dfun(xv);
       if dfx == 0
          error(' Arresto per azzeramento di dfun');
       else
          xn = xv - mol*fun(xv)/dfx;
          err = abs(xn-xv);
          xvect = [xvect; xn];
          it = it+1;
          xv = xn;
       end
    end
    
    % only critical output
    if (it > nmax)
        fprintf('Max iteration number k reached! : %d \n',it);
    end

end


function [rr_vect] = vectorRectifier(rr_vect)

    inversion=0;
    if size(rr_vect,1) < size(rr_vect,2)
        rr_vect = rr_vect';
        inversion = 1;
    end
    
    dim = size(rr_vect,1);
    
    if rem(dim,2) == 0
        for i = 1:dim
            if i<=dim/2
                dummy(i,:) = rr_vect(dim/2+i,:);
            else
                dummy(i,:) = rr_vect(i-dim/2,:);
            end
        end
        dummy=[dummy; dummy(1,:)];
    else
        for i = 1:dim
            if i <= floor(dim/2)
                dummy(i,:) = rr_vect(floor(dim/2)+i,:);
            else
                dummy(i,:) = rr_vect(i-floor(dim/2),:);
            end
        end
    end
    
    if inversion
        rr_vect = dummy';
    else
        rr_vect = dummy;
    end

end