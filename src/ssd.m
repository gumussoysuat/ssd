classdef ssd < handle & matlab.mixin.Copyable
    % SSD   Constructs state-space object with delays in the form,
    %       x'(t) = sum A_i*x(t-hA_i) + sum B_i*u(t-hB_i)
    %       y (t) = sum C_i*x(t-hC_i) + sum D_i*u(t-hD_i)
    %
    %       SYS = SSD(A, DelaysA, B, DelaysB, C, DelaysC, D, DelaysD)
    %
    %       A : State matrices of state equation, cat(3,A0,A1,...)
    %       hA: Unique and increasing delays of A matrices, [0,hA_1,...]
    %       B : Input matrices of state equation, cat(3,B0,B1,...)
    %       hB: Unique and increasing delays of B matrices, [0,hB_1,...]
    %       C : State matrices of output equation, cat(3,C0,C1,...)
    %       hC: Unique and increasing delays of C matrices, [0,hC_1,...]
    %       D : State matrices of output equation, cat(3,D0,D1,...)
    %       hD: Unique and increasing delays of D matrices, [0,hD_1,...]
    %
    % Example:
    % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0)

    properties (SetAccess=private)
        A       {mustBeReal, mustBeFinite}
        hA(1,:) {mustBeNonnegative, mustBeFinite}
        B       {mustBeReal, mustBeFinite}
        hB(1,:) {mustBeNonnegative, mustBeFinite}
        C       {mustBeReal, mustBeFinite}
        hC(1,:) {mustBeNonnegative, mustBeFinite}
        D       {mustBeReal, mustBeFinite}
        hD(1,:) {mustBeNonnegative, mustBeFinite}
    end
    properties (SetAccess=public)
        name    % variable name
    end
    properties (Hidden, SetAccess=private)
        n       % state dimension
        nu      % input dimension
        ny      % output dimension
        mA      % number of A matrices
        mB      % number of B matrices
        mC      % number of C matrices
        mD      % number of D matrices
        isss = false    % ss exist, true or false
        ss      % delayss object
        Wc      % controllability gramian
        Wo      % observability gramian
        h2      % h2 norm
        sparse = false; % true or false
        spA     % sparse A
        spB     % sparse B
        spC     % sparse C
        spD     % sparse D
        id      % id for object
    end
    methods
        function obj = ssd(A,hA,B,hB,C,hC,D,hD)
            % Instantiates a delay system
            obj.A = A;
            obj.B = B;
            obj.C = C;

            if nargin<7
                obj.D = [];
            else
                obj.D = D;
            end

            obj.n = size(obj.A,1);
            obj.nu = size(obj.B,2);
            obj.ny = size(obj.C,1);

            obj = checkMatrices(obj);

            obj.hA = hA;
            obj.hB = hB;
            obj.hC = hC;

            if nargin<8
                obj.hD = [];
            else
                obj.hD = hD;
            end

            obj = checkDelays(obj);

            obj.mA = length(obj.hA);
            obj.mB = length(obj.hB);
            obj.mC = length(obj.hC);
            obj.mD = length(obj.hD);

            checkMatricesDelaysConsistency(obj)
        end
        % plotting
        function bode(obj,varargin)
            % Same as standard bode function.
            % see also bode.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % bode(sys)
            names = {};
            for ct=1:numel(varargin)+1
                names{ct} = inputname(ct);
            end
            plothelperfcn(obj,varargin,names,@bode)
        end
        function bodemag(obj,varargin)
            % Same as standard bodemag function.
            % see also bodemag.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % bodemag(sys)            
            names = {};
            for ct=1:numel(varargin)+1
                names{ct} = inputname(ct);
            end
            plothelperfcn(obj,varargin,names,@bodemag)
        end
        function step(obj,varargin)
            % Same as standard step function.
            % see also step.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % step(sys)            
            names = {};
            for ct=1:numel(varargin)+1
                names{ct} = inputname(ct);
            end
            plothelperfcn(obj,varargin,names,@step)
        end
        function sigma(obj,varargin)
            % Same as standard sigma function.
            % see also sigma.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % sigma(sys)
            names = {};
            for ct=1:numel(varargin)+1
                names{ct} = inputname(ct);
            end
            plothelperfcn(obj,varargin,names,@sigma)
        end
        function nyquist(obj,varargin)
            % Same as standard nyquist function.
            % see also nyquist.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % nyquist(sys)
            names = {};
            for ct=1:numel(varargin)+1
                names{ct} = inputname(ct);
            end
            plothelperfcn(obj,varargin,names,@nyquist)
        end        
        % rd2d functions
        function varargout = gram(obj,type,opt)
            % GRAM  Computes the controllability and observability gramians
            % of delay system.
            %
            %        Wc = GRAM(SYS,'c') and Wo = GRAM(SYS,'o') computes
            %        controllability and observability gramians of delay
            %        system SYS (see SSD).
            %            
            %        [Wc,Wo] = GRAM(SYS,'co') computes both.
            %
            %        Wc = GRAM(SYS,'c',OPT) and Wo = GRAM(SYS,'o',OPT)
            %        computes controllability and observability gramians
            %        with options (see SSDOPTIONS).
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % [Wc,Wo] = gram(sys,'co')

            arguments
                obj {mustBeA(obj,'ssd')}
                type char {mustBeMember(type,{'c','o','co','oc'})}
                opt {mustBeA(opt,'ssdoptions'), mustBeScalarOrEmpty(opt)} = ssdoptions()
            end
            
            [abstol,reltol,freqs] = ssd.extractOptions(opt);
            [A,hA,mA,B,hB,mB,C,hC,mC,~,~,~,n,nu,ny] = extractSystemInfo(obj,opt);

            if opt.Sparse
                cfun = @cgramsparsefun;
                valsumc = sparse(n);
                ofun = @ogramsparsefun;
                valsumo = sparse(n);
            else
                cfun = @cgramfun;
                valsumc = zeros(n);
                ofun = @ogramfun;
                valsumo = zeros(n);
            end

            if any(strfind(type,'c'))
                for ct1=1:size(freqs,1)
                    valsumc = valsumc + integral(cfun,freqs(ct1,1),freqs(ct1,2),'ArrayValued',true,'AbsTol',abstol,'RelTol',reltol);
                end
                Wc = valsumc/pi;
            end            
            if any(strfind(type,'o'))
                for ct1=1:size(freqs,1)
                    valsumo = valsumo + integral(ofun,freqs(ct1,1),freqs(ct1,2),'ArrayValued',true,'AbsTol',abstol,'RelTol',reltol); %#ok<*PROPLC>
                end
                Wo = valsumo/pi;
            end

            switch type
                case 'c'
                    varargout{1} = Wc;                    
                case 'o'
                    varargout{1} = Wo;                    
                case 'co'
                    varargout{1} = Wc;
                    varargout{2} = Wo;
                case 'oc'
                    varargout{1} = Wo;
                    varargout{2} = Wc;
            end
            function Y = cgramfun(w)
                sumA = 1j*w*eye(n);
                for ct=1:mA
                    sumA = sumA - A(:,:,ct)*exp(-hA(ct)*1j*w);
                end
                sumB = zeros(n,nu);
                for ct=1:mB
                    sumB = sumB + B(:,:,ct)*exp(-hB(ct)*1j*w);
                end
                AB = sumA\sumB;
                Y = real(AB*AB');
            end
            function Y = cgramsparsefun(w)
                sumA = 1j*w*speye(n);
                for ct=1:mA
                    sumA = sumA - A{ct}*exp(-hA(ct)*1j*w);
                end
                sumB = sparse(n,nu);
                for ct=1:mB
                    sumB = sumB + B{ct}*exp(-hB(ct)*1j*w);
                end
                AB = sumA\sumB;
                Y = real(AB*AB');
            end            
            function Y = ogramfun(w)
                sumA = 1j*w*eye(n);
                for ct=1:mA
                    sumA = sumA - A(:,:,ct)*exp(-hA(ct)*1j*w);
                end
                sumC = zeros(ny,n);
                for ct=1:mC
                    sumC = sumC + C(:,:,ct)*exp(-hC(ct)*1j*w);
                end
                CA = sumC/sumA;
                Y = real(CA'*CA);
            end
            function Y = ogramsparsefun(w)
                sumA = 1j*w*speye(n);
                for ct=1:mA
                    sumA = sumA - A{ct}*exp(-hA(ct)*1j*w);
                end
                sumC = sparse(ny,n);
                for ct=1:mC
                    sumC = sumC + C{ct}*exp(-hC(ct)*1j*w);
                end
                CA = sumC/sumA;
                Y = real(CA'*CA);
            end
        end
        function info = computeBalancedInfo(obj,opt)
            % COMPUTEBALANCEDINFO Compute the necessary data for balanced
            % realization and reduction.
            %        
            %        INFO = COMPUTEBALANCEDINFO(SYS) computes the balanced
            %        info of original delay system SYS (see SSD). INFO
            %        object contains necessary data to do realization and
            %        reduction without computing gramians (see BALINFO).
            %        You can see energy of the states by PLOT(INFO).
            %
            %        INFO = COMPUTEBALANCEDINFO(SYS,OPT) computes with
            %        options (see SSDOPTIONS).
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % info = computeBalancedInfo(sys)

            arguments
                obj {mustBeA(obj,'ssd')}
                opt {mustBeA(opt,'ssdoptions'), mustBeScalarOrEmpty(opt)} = ssdoptions()
            end

            % compute grammians
            [Wc,Wo] = gram(obj,'co',opt);
            Rc = computeChol(Wc);
            Ro = computeChol(Wo);
            [U,S,V] = svd(full(Rc*Ro'));
            TL=(S^(-1/2)*V'*Ro);
            TR=Rc'*U*S^(-1/2);
            
            info = balinfo(Rc,Ro,diag(S),TL,TR,opt);

            function R = computeChol(W)
                % modified from https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
                try
                    W = (W+W')/2;
                    [~,Sw,Vw] = svd(W);
                    W = (W+Vw*Sw*Vw')/2;
                    W = (W+W')/2;
                    R = chol(W);
                catch
                    % not positive definite due to round off
                    % if errors, increase offset to make it positive
                    % definite
                    offset=eps;
                    W = W + (offset+max(eps,abs(min(eig(W)))))*eye(size(W));
                    try
                        R = chol(W);
                    catch
                        error('Set offset value to make gramian positive definite above code.')
                    end
                end
            end
        end
        function [sysb,info] = balreal(obj,arg)
            % BALREAL Transform delay system into the balanced realization.
            %
            %        [SYSB,INFO] = BALREAL(SYS) computes the balanced
            %        realization SYSB of original delay system SYS (see
            %        SSD) and creates INFO for future computations
            %        capturing gramian information (see BALINFO).
            %
            %        [SYSB,INFO] = BALREAL(SYS,OPT) computes with
            %        options (see SSDOPTIONS).
            %            
            %        SYSB = BALREAL(SYS,INFO) computes the balanced
            %        realization SYSB of original delay system SYS  without
            %        computing gramians and using the information in
            %        INFO (see BALINFO).
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % sysb = balreal(sys)

            opt = ssdoptions();
            if nargin<2
                info = computeBalancedInfo(obj);
            elseif isa(arg,'ssdoptions')
                opt = arg;
                info = computeBalancedInfo(obj,opt);
            elseif isa(arg,'balinfo')
                info = arg;
            else
                error('Second argument must be ssdoptions or previously computed balinfo for this delay system.');
            end            

            TL = info.TL;
            TR = info.TR;

            sysb = applyTransformation(obj,TL,TR,info.options);
        end
        function [sysr,info] = balred(obj,nr,varargin)
            % BALRED Reduces delay system in balanced gramian sense.
            %
            %        [SYSR,INFO] = BALRED(SYS,NR) computes the reduced
            %        delay system RSYS with order NR of original delay
            %        system SYS (see SSD). INFO contains data to reduce
            %        system without recomputing gramians (see BALINFO).
            %        State contributions can seen by PLOT(INFO).
            %         
            %        [SYSR,INFO] = BALRED(SYS,NR,OPT) computes with options
            %        (see SSDOPTIONS).
            %
            %        SYSRs = BALRED(SYS,NR,INFO) computes the reduced
            %        delay system without recomputing gramians.
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % sysr = balred(sys,1)

            if ~(isequal(size(nr),[1 1]) && isnumeric(nr) && isreal(nr) && (rem(nr,1)==0) && (nr>=0))
                error('Reduction order must be a non-negative integer number.');
            end
            if (nr>obj.n)
                error('Reduction order must be smaller or equal to the original delay system order.')
            end

            if nargin>2 % parse opt or balinfo
                if isa(varargin{1},'balinfo')
                    info = varargin{1};
                elseif isa(varargin{1},'ssdoptions')
                    opt = varargin{1};
                    info = computeBalancedInfo(obj,opt);
                else
                    error('Third argument must be ssdoptions or balinfo.')
                end
            else
                info = computeBalancedInfo(obj);                
            end

            TL = info.TL;
            TR = info.TR;

            sysr = applyTransformation(obj,TL(1:nr,:),TR(:,1:nr),info.options);
        end
        function val = h2norm(obj,opt)
            % H2NORM Computes the H2 norm of exponentially stable delay system. H2 norm is
            % equivalent to L2 norm for this case.
            %
            %        H2 = H2NORM(SYS) computes the H2 norm of delay system
            %        SYS (see SSD).
            %
            %        H2 = H2NORM(SYS,OPT) computes the H2 norm of delay system
            %        with options (see SSDOPTIONS).
            %
            % Example: 
            % sys = ssd(cat(3,[-2 -1;-1.5 -.5],[0 .5;1 0]),[0 1],[1;-1],0,[2 .2],0,[],0);
            % h2norm(sys)
            arguments
                obj {mustBeA(obj,'ssd')}
                opt {mustBeA(opt,'ssdoptions'), mustBeScalarOrEmpty(opt)} = ssdoptions()
            end

            if norm(obj.D(:,:,1))>0 || (obj.mD>1)
                warning('H2 norm is infinite for non-zero D matrices.')
                val = Inf;
            else
                [abstol,reltol,freqs] = ssd.extractOptions(opt);
                [A,hA,mA,B,hB,mB,C,hC,mC,~,~,~,n,nu,ny] = extractSystemInfo(obj,opt);

                if opt.Sparse
                    fun = @l2normsparsefun;
                else
                    fun = @l2normfun;
                end

                valsum = 0;
                for ct=1:size(freqs,1)
                    valsum = valsum + integral(fun,freqs(ct,1),freqs(ct,2),'ArrayValued',true,'AbsTol',abstol,'RelTol',reltol);
                end
                val = sqrt(valsum/pi);
            end
          
            function Y = l2normfun(w)
                sumB = zeros(n,nu);
                for ct2=1:mB
                    sumB = sumB + B(:,:,ct2)*exp(-hB(ct2)*1j*w);
                end
                sumA = 1j*w*eye(n);                
                for ct2=1:mA
                    sumA = sumA - A(:,:,ct2)*exp(-hA(ct2)*1j*w);
                end
                sumC = zeros(ny,n);
                for ct2=1:mC
                    sumC = sumC + C(:,:,ct2)*exp(-hC(ct2)*1j*w);
                end
                CAB = (sumC/sumA)*sumB;
                Y = trace(real(CAB'*CAB));
            end
            function Y = l2normsparsefun(w)
                sumB = sparse(n,nu); %#ok<*CPROPLC> 
                for ct2=1:mB
                    sumB = sumB + B{ct2}*exp(-hB(ct2)*1j*w);
                end
                sumA = 1j*w*speye(n);                
                for ct2=1:mA
                    sumA = sumA - A{ct2}*exp(-hA(ct2)*1j*w);
                end
                sumC = sparse(ny,n);
                for ct2=1:mC
                    sumC = sumC + C{ct2}*exp(-hC(ct2)*1j*w);
                end
                CAB = (sumC/sumA)*sumB;
                Y = trace(real(CAB'*CAB));
            end              
        end
    end
    methods (Hidden)
        % system functions
        function obj = setss(obj)
            if ~obj.isss
                AllDelays = [obj.hA obj.hB obj.hC obj.hD];
                AllDelays = sort(unique(AllDelays(AllDelays~=0)));
                DelayT = struct('delay',[],'a',[],'b',[],'c',[],'d',[]); %#ok<*PROP>
                for ct=1:length(AllDelays)
                    curDelay = AllDelays(ct);
                    DelayT(ct).delay = curDelay;
                    ixa = find(obj.hA==curDelay);
                    if isempty(ixa)
                        DelayT(ct).a = [];
                    else
                        DelayT(ct).a = obj.A(:,:,ixa);
                    end
                    ixb = find(obj.hB==curDelay);
                    if isempty(ixb)
                        DelayT(ct).b = [];
                    else
                        DelayT(ct).b = obj.B(:,:,ixb);
                    end
                    ixc = find(obj.hC==curDelay);
                    if isempty(ixc)
                        DelayT(ct).c = [];
                    else
                        DelayT(ct).c = obj.C(:,:,ixc);
                    end
                    ixd = find(obj.hD==curDelay);
                    if isempty(ixd)
                        DelayT(ct).d = [];
                    else
                        DelayT(ct).d = obj.D(:,:,ixd);
                    end
                end
                obj.ss = delayss(obj.A(:,:,1),obj.B(:,:,1),obj.C(:,:,1),obj.D(:,:,1),DelayT);
                obj.ss.name = obj.name;
            end
        end
        function obj = constructSparse(obj)
            obj.sparse = true;
            A = obj.A; B = obj.B; C = obj.C; D = obj.D;
            for ct=1:obj.mA
                Asp{ct} = sparse(A(:,:,ct));
            end
            for ct=1:obj.mB
                Bsp{ct} = sparse(B(:,:,ct));
            end
            for ct=1:obj.mC
                Csp{ct} = sparse(C(:,:,ct));
            end
            for ct=1:obj.mD
                Dsp{ct} = sparse(D(:,:,ct));
            end
            obj.spA = Asp; obj.spB = Bsp; obj.spC = Csp; obj.spD = Dsp;
        end
        function [A,hA,mA,B,hB,mB,C,hC,mC,D,hD,mD,n,nu,ny] = extractSystemInfo(obj,opt)
            n = obj.n; nu = obj.nu; ny = obj.ny;
            hA = obj.hA; hB = obj.hB; hC = obj.hC; hD = obj.hD;
            mA = obj.mA; mB = obj.mB; mC = obj.mC; mD = obj.mD;

            if opt.Sparse
                if ~obj.sparse
                    constructSparse(obj);
                end
                A = obj.spA; B = obj.spB; C = obj.spC; D = obj.spD;
            else
                A = obj.A; B = obj.B; C = obj.C; D = obj.D;
            end
        end
        function sys = applyTransformation(obj,TL,TR,opt)
            opt.Sparse = false;
            [A,hA,mA,B,hB,mB,C,hC,mC,D,hD,~,n,~,~] = extractSystemInfo(obj,opt);

            if ~isequal(size(TL,2),n)
                error('Balinfo dimensions are inconsistent with delay system.');
            end

            nr = size(TL,1);

            Anew = zeros([nr nr mA]);
            for ct=1:mA
                Anew(:,:,ct) = TL*A(:,:,ct)*TR;
            end
            Bnew = zeros([nr obj.nu mB]);
            for ct=1:mB
                Bnew(:,:,ct) = TL*B(:,:,ct);
            end
            Cnew = zeros([obj.ny nr mC]);
            for ct=1:mC
                Cnew(:,:,ct) = C(:,:,ct)*TR;
            end

            sys = ssd(Anew,hA,Bnew,hB,Cnew,hC,D,hD);
        end
        % validation
        function obj = checkMatrices(obj)
            if size(obj.B,1)~=obj.n
                error('B matrices are not compatible with A matrices');
            end
            if size(obj.C,2)~=obj.n
                error('C matrices are not compatible with A matrices');
            end
            if isempty(obj.D)
                obj.D = zeros(obj.ny,obj.nu);
            elseif size(obj.D,2)~=obj.nu
                error('D matrices are not compatible with B matrices');
            elseif size(obj.D,1)~=obj.ny
                error('D matrices are not compatible with C matrices');
            end
        end
        function obj = checkDelays(obj)
            obj = setZeroIfEmpty(obj,'hA');
            obj = setZeroIfEmpty(obj,'hB');
            obj = setZeroIfEmpty(obj,'hC');
            obj = setZeroIfEmpty(obj,'hD');
            checkUnique(obj,'A');
            checkUnique(obj,'B');
            checkUnique(obj,'C');
            checkUnique(obj,'D');
        end
        function checkMatricesDelaysConsistency(obj)
            if size(obj.A,3)~=obj.mA
                error('The number of A matrices (%d) and DelaysA (%d) is not compatible',size(obj.A,3),obj.mA);
            end
            if size(obj.B,3)~=obj.mB
                error('The number of B matrices (%d) and DelaysB (%d) is not compatible',size(obj.B,3),obj.mB);
            end
            if size(obj.C,3)~=obj.mC
                error('The number of C matrices (%d) and DelaysC (%d) is not compatible',size(obj.C,3),obj.mC);
            end
            if size(obj.D,3)~=obj.mD
                error('The number of D matrices (%d) and DelaysD (%d) is not compatible',size(obj.D,3),obj.mD);
            end
        end
        function obj = setZeroIfEmpty(obj,field)
            if isempty(obj.(field))
                obj.(field) = 0;
            end
        end
        function checkUnique(obj,str)
            val = obj.(['h' str]);
            if ~isequal(unique(val),val)
                error('h%s must be non-repeating.',str)
            end
        end
        % plotting
        function checkss(obj)
            if ~obj.isss
                if (obj.mA==1) && (obj.mB==1) && (obj.mC==1) && (obj.mD==1)
                    obj.ss = ss(obj.A(:,:,1),obj.B(:,:,1),obj.C(:,:,1),obj.D(:,:,1)); %#ok<CPROP>
                    obj.ss.name = obj.name;
                else                    
                    setss(obj);
                end
            end
        end
        function plothelperfcn(obj,args,names,fun)
            checkss(obj)
            ix = find(cellfun(@(x) isa(x,'ssd'),args));
            if isempty(obj.name)
                obj.name = names{1};
            end
            for ct=1:length(ix)
                if isempty(args{ix(ct)}.name)
                    args{ix(ct)}.name = names{ix(ct)+1};
                end
            end
            if isempty(ix)
                fun(obj.ss,args{:})
            else
                vararg_cur = args(1:ix(1)-1);
                fun(obj,vararg_cur{:})
                hold on;
                vararg_later = args(ix(1):end);
                fun(vararg_later{:})
            end
        end
    end
    methods (Static,Hidden)
        function [atol, rtol, freqs] = extractOptions(opt)
                atol = opt.AbsTol;
                rtol = opt.RelTol;
                if isempty(opt.FreqInt)
                    freqs = [0 Inf];
                else
                    freqs = opt.FreqInt;
                end
        end
    end
end