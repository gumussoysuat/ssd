classdef ssdoptions
    % SSDOPTIONS   Constructs options for ssd object.
    %
    %    AbsTol         Absolute error tolerance for integral. Default is 1e-3.
    %    RelTol         Relative error tolerance for integral. Default is 1e-3.
    %    Sparse         Enables sparse matrix computation. Default is false.
    %    FreqInt        Frequency intervals on which the gramian is defined.
    %                   Two-column matrix where each row defines an interval.
    %                   Default is [].
    %
    % Example:
    % opt = ssdoptions('AbsTol',1e-2,'Sparse',true,'FreqInt',[2 4;5 8])

    properties
        AbsTol (1,1) {ssdoptions.validateAbsTol};
        RelTol (1,1) {ssdoptions.validateRelTol};
        Sparse (1,1) {boolean};
        FreqInt (:,2) {ssdoptions.validateFrequencyIntervals};
    end

    methods
        function obj = ssdoptions(NameValueArgs)
            arguments
                NameValueArgs.AbsTol = 1e-3;
                NameValueArgs.RelTol = 1e-3;
                NameValueArgs.Sparse = false;
                NameValueArgs.FreqInt = [];
                NameValueArgs.MaxOrder = [];
            end
            obj.AbsTol = NameValueArgs.AbsTol;
            obj.RelTol = NameValueArgs.RelTol;
            obj.Sparse = NameValueArgs.Sparse;
            obj.FreqInt = NameValueArgs.FreqInt;
        end
        function obj = set.FreqInt(obj,val)
            y = val';
            v = y(:);
            [~,ix] = unique(v,'stable');
            ixdub = setdiff(1:numel(v),ix);
            obj.FreqInt = reshape(setdiff(v,v(ixdub)),2,[])';
        end
    end

    methods (Static = true)
        function validateAbsTol(tol)
            if ~(isfloat(tol) && isreal(tol) && tol >= 0)
                error('AbsTol must be a non-negative number');
            end
        end
        function validateRelTol(x)
            if ~(isfloat(x) && isreal(x) && x >= 0)
                error('RelTol must be a non-negative number');
            end
        end
        function validateFrequencyIntervals(x)
            if ~isempty(x)
                if ~(isnumeric(x) && isreal(x) && (size(x,2)==2) && all(all(x>=0)) && ismatrix(x))
                    error('Frequency intervals must be a two-column matrix of non-negative numbers');
                end
                y = x';
                v = y(:);
                if any(diff(y)<=0) && diff(v)<0
                    error('Boundaries of frequency intervals must be an non-decreasing sequence.');
                end
            end
        end
    end
end

