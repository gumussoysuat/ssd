classdef balinfo
    % BALINFO  Constructs data for balanced realization and reduction.
    %
    %    Rc    Cholesky matrix of controllability gramian.
    %    Ro    Cholesky matrix of observability gramian.
    %    sv    Hankel singular values.
    %    TL    Left matrix for balanced realization or reduction.
    %    TR    Right matrix for balanced realization or reduction.
    %
    % Example:
    % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
    %
    % [~,info1] = balreal(sys)  % using balreal
    % [~,info2] = balred(sys,1) % using balred

    properties (SetAccess=private)
        Rc {balinfo.checkR}
        Ro {balinfo.checkR}
        sv {balinfo.checksv}
        TL {balinfo.checkT}
        TR {balinfo.checkT}
    end
    properties (SetAccess=private, Hidden)
        options {mustBeA(options,'ssdoptions')} = ssdoptions()
    end
    methods
        function obj = balinfo(Rc,Ro,sv,TL,TR,opt)
            obj.Rc = Rc;
            obj.Ro = Ro;
            obj.sv = sv;
            obj.TL = TL;
            obj.TR = TR;
            obj.options = opt;

            % validate consistency across arguments
            if ~isequal(size(Rc),size(Ro))
                error('Cholesky matrices must have same dimensions');
            end
            if ~isequal(size(TL),size(TR))
                error('Left and right matrices must have same dimensions');
            end
            if ~isequal(size(Rc),size(TL))
                error('Cholesky and Left/Right matrices must have same dimensions');
            end
            if size(Rc,2)~=length(sv)
                error('The number of singular values must be equal to system size');
            end
        end
        function plot(obj)
            % Plots state contributions. 
            %
            % Example:
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % [~,info] = balreal(sys);
            % plot(info)

            bar(1:size(obj.Rc,2),obj.sv)
            xlabel('Order (Number of states)');
            ylabel('State contribution');
            title('Energies of states');
        end
    end
    methods (Static = true)
        function checkR(R)
            if ~(isnumeric(R) && isreal(R) && ismatrix(R))
                error('Cholesky matrix must be a real upper triangle matrix.');
            end
        end
        function checkT(T)
            if ~(isnumeric(T) && isreal(T) && ismatrix(T))
                error('Left and right matrices must be a real square matrix.');
            end
        end
        function checksv(sv)
            if ~(isnumeric(sv) && isreal(sv) && all(sv>=0))
                error('Singular values must be non-negative real numbers.');
            end
        end        
    end
end