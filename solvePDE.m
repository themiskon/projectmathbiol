function [t, C] = solvePDE(fn,f,fz,ff, v,vf,D,Df, tRange, xc, C0,lat,sensi)

% ----------------------------------------------
%  Grid parameters:
% ----------------------------------------------
n = length(xc);              % no. of grid cells

% Derived grid properties
dx = xc(2)-xc(1);     % Grid spacing; assume uniform spacing
xi = xc - dx/2;       % Points of interfaces between cells

% ----------------------------------------------
% Solve the system of ODEs:
% ----------------------------------------------
option=odeset('nonnegative',1:length(C0));
[t, C] = ode45(@derivative, tRange, C0,option,lat,sensi);

if nargout==0
    clf
    surf(t,xc,C')
    xlabel('t')
    ylabel('x')
    zlabel('C')
    shading interp
end

    function dCdt = derivative(t,C,lat,sensi)
        P=C(1:n);
        Z=C(n+1:2*n);
        F=C(2*n+1:+3*n);
        N=C(3*n+1:end);
        % Initialize
        dCdt = zeros(4*n,1);  % Rate of change of average concentration in each cells

        %
        % Transport: Loop over all interfaces
        %
        for i = 1:(n-1)
            %
            % Compute flux from cell i to cell i+1
            %
            %
            % Advective flux: Upwind method
            vi = v(t,xi(i));
            vif = vf(t,xi(i));
            if vi>0
                Jap = vi*C(i);
                Jaz = vi*C(n+i);
            else
                Jap = vi*C(i+1);
                Jaz = vi*C(n+i+1); % tou toum tss
            end
             if  vif>0
                Jaf = vif*C(2*n+i);
            else
                Jaf = vif*C(2*n+i+1);
            end

            Jdp = D(t,xi(i))*(C(i)-C(i+1))/dx;   % diffusive flux phytoplankton
            Jdz = D(t,xi(i))*(C(n+i)-C(n+i+1))/dx; %diffusive flux zooplankton
            Jdf = Df(t,xi(i))*(C(2*n+i)-C(2*n+1+i))/dx; %diffusive fish
            Jdn= D(t,xi(i))*(C(3*n+i)-C(3*n+1+i))/dx; %diffuse nutrients
            Jz=Jaz + Jdz;
            Jp = Jap  + Jdp;         % Total flux
            Jf=Jdf+Jaf; %fish spread
            %
            % Remove transport from cell i and add to cell i+1
            %
            %phytoplankton advection-diffusion
            dCdt(i) = dCdt(i) - Jp/dx;
            dCdt(i+1) = dCdt(i+1) + Jp/dx;
            %zooplankton advection-diffusion
            dCdt(n+i)=dCdt(n+i) - Jz/dx;
            dCdt(n+i+1)=dCdt(n+i+1) + Jz/dx;
            %fish diffusion
            dCdt(2*n+i+1)=dCdt(2*n+i+1)+Jf/dx;
            dCdt(2*n+i)=dCdt(2*n+i)-Jf/dx;
            %nutrient diffusion
            dCdt(3*n+1+i)=dCdt(3*n+1+i)+Jdn/dx;
            dCdt(3*n+i)=dCdt(3*n+i)-Jdn/dx;

        end
        dCdt(end)=0; %boundary condition for bottom nutrients (nutrient input)
        %
        % Reaction: Loop over all cells
        tt=t;
        ii=0;
        while tt>365
            tt=tt-365;
            ii=ii+1;
        end
        depth=dx*ones(1,n)';
        int=0.05.*cumsum(P).*dx;
        M=cumsum(depth).*(0.5+0.45*cos(2*pi*t/365))./dx;
        Iabs=(exp(-0.0375.*cumsum(depth)-int));
       
        if ii>=5
            if tt>=60 && tt<=61
                dCdt(2*n+1)=1*(sensi-1);
            end
            for i = 1:n
                g=Iabs(i)*exp(-0.025*M(i))*(0.5-0.6*sin(pi*lat/180)*cos(2*pi*t/365));
                dCdt(i,1) = dCdt(i) + f(t,xc(i),P(i),Z(i),F(i),N(i),g);
                dCdt(n+i,1)= dCdt(n+i)+fz(t,xc(i),Z(i),P(i),F(i));
                dCdt(2*n+i,1)=dCdt(2*n+i)+ff(t,xc(i),Z(i),F(i),P(i));
                dCdt(3*n+i,1)=dCdt(3*n+i)+min(0,fn(t,xc(i),P(i),N(i),g));
            end
            if tt>=270
                dCdt(2*n+1:3*n)=-C(2*n+1:+3*n);
            end
        else
            for i = 1:n
                g=exp(-0.025*M(i))*(0.5-0.6*sin(pi*lat/180)*cos(2*pi*t/365))*Iabs(i);
                dCdt(i,1) = dCdt(i) + f(t,xc(i),P(i),Z(i),F(i),N(i),g);
                dCdt(n+i,1)= dCdt(n+i)+fz(t,xc(i),Z(i),P(i),F(i));
                dCdt(2*n+i,1)=dCdt(2*n+i)+ff(t,xc(i),Z(i),F(i),P(i));
                dCdt(3*n+i,1)=dCdt(3*n+i)+fn(t,xc(i),P(i),N(i),g);
            end
        end
    end
end
