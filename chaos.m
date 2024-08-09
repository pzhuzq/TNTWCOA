function result = chaos(index, N, dim)

switch index
    case 1
        % tent %
        tent=1.5;  
        Tent=rand(N,dim);
        for i=1:N
            for j=2:dim
                if Tent(i,j-1)<tent
                    Tent(i,j)=Tent(i,j-1)/tent;
                elseif Tent(i,j-1)>=tent
                    Tent(i,j)=(1-Tent(i,j-1))/(1-tent);
                end
            end
        end
        result = Tent;
        
    case 2
        % Logistic %
        miu=3;  
        Logistic=rand(N,dim);
        for i=1:N
            for j=2:dim
                Logistic(i,j)=miu.* Logistic(i,j-1).*(1-Logistic(i,j-1));
            end
        end
        result = Logistic;
        
    case 3
        % Cubic%
        cubic=1.3;
        Cubic=rand(N,dim);
        for i=1:N
            for j=2:dim
                Cubic(i,j)=cubic.*Cubic(i,j-1).*(1-Cubic(i,j-1).^2);
            end
        end
        result = Cubic;
        
    case 4
        %chebyshev%
        chebyshev=4;
        Chebyshev=rand(N,dim);
        for i=1:N
            for j=2:dim
                Chebyshev(i,j)=cos(chebyshev.*acos(Chebyshev(i,j-1)));
            end
        end
        result = Chebyshev;
        
    case 5
        %Piecewise%
        p=1;
        Piecewise=rand(N,dim);
        for i=1:N
            for j=2:dim
                if Piecewise(i,j-1)>0&&Piecewise(i,j-1)<p
                    Piecewise(i,j)=Piecewise(i,j-1)/p;
                elseif Piecewise(i,j-1)>=p&&Piecewise(i,j-1)<0.5
                    Piecewise(i,j)=(Piecewise(i,j-1)-p)/(0.5-p);
                elseif Piecewise(i,j-1)>=0.5&&Piecewise(i,j-1)<1-p
                    Piecewise(i,j)=(1-p-Piecewise(i,j-1))/(0.5-p);
                elseif Piecewise(i,j-1)>=1-p&&Piecewise(i,j-1)<1
                    Piecewise(i,j)=(1-Piecewise(i,j-1))/p;
                end
            end
        end
        result = Piecewise;
        
        
    case 6
        %sinusoidal%
        sinusoidal=2;
        Sinusoidal=rand(N,dim);
        for i=1:N
            for j=2:dim
                Sinusoidal(i,j)=sinusoidal*Sinusoidal(i,j-1).^2*(sin(pi*Sinusoidal(i,j-1)));
            end
        end
        result = Sinusoidal;
        
    case 7
        %Sine%
        sine=2;
        Sine=rand(N,dim);
        for i=1:N
            for j=2:dim
                Sine(i,j)=(4/sine)*sin(pi*Sine(i,j-1));
            end
        end
        result = Sine;
        
        
    case 8
        %         ICMIC%
        icmic=2;
        ICMIC=rand(N,dim);
        for i=1:N
            for j=2:dim
                ICMIC(i,j)=sin(icmic/ICMIC(i,j-1));
            end
        end
        result = ICMIC;
        
        
    case 9
        % Circle%
        a = 0.3; b=0.3;
        Circle=rand(N,dim);
        for i=1:N
            for j=2:dim
                Circle(i,j)=mod(Circle(i,j-1)+a-b/(2*pi)*sin(2*pi*Circle(i,j-1)),1);
            end
        end
        result = Circle;
    case 10
        %Bernoulli%
        lammda = 0.4;
        Bernoulli=rand(N,dim);
        for i=1:N
            for j=2:dim
                if Bernoulli(i,j-1) <  1-lammda
                    Bernoulli(i,j)= Bernoulli(i,j-1)/(1-lammda);
                else
                    Bernoulli(i,j)= (Bernoulli(i,j-1)-1+lammda)/lammda;
                end
            end
        end
        result = Bernoulli;
        
        
        
    case 11    %NO
        result =rand(N,dim);
        
        
end
end


