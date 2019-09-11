function [u,v] = decompose(gI,MAX_SIZE,OVERLAP,lambda,tol,en_disp,Confi)
% structure-texture image decomposition, described in Wotao Yin, "A comparison
%           of three total variation based texture extraction models", 2007
% inputs:
%           gI:         gray image to decompose
%           MAX_SIZE:	maximum of patch size
%           OVERLAP:    overlap
%           lambda:     balance factor
%           tol:        tolerance of solutions
%           en_disp:    display option
%           Confi:      probability that an edge is caused by pattern on
%                       object
% outputs:
%           u:          structure component
%           v:          texture component
    
    if nargin == 1
        MAX_SIZE = 50;
        OVERLAP = 10;
        lambda = 0.5;
        tol = 1e-8;
        en_disp = 0;
        Confi = ones(size(gI));
    elseif nargin == 2
        OVERLAP = 10;
        lambda = 0.5;
        tol = 1e-8;
        en_disp = 0;
        Confi = ones(size(gI));
    elseif nargin == 3
        lambda = 0.5;
        tol = 1e-8;
        en_disp = 0;
        Confi = ones(size(gI));
    elseif nargin == 4
        tol = 1e-8;
        en_disp = 0;
        Confi = ones(size(gI));
    elseif nargin == 5
        en_disp = 0;
        Confi = ones(size(gI));
    elseif nargin == 6
         Confi = ones(size(gI));
  %      Confi = zeros(size(gI));
    end
    
    gI = double(gI);
    [H,W] = size(gI);
    addpath('C:\Program Files\Mosek\8\toolbox\r2014aom');
    
    index_y = 1:MAX_SIZE-OVERLAP:H;
    index_x = 1:MAX_SIZE-OVERLAP:W;
    u = zeros(H,W); cnt = zeros(H,W);
%     disp(['image size = ',num2str(H),'¡Á',num2str(W)]);
%     disp(['patch size = ',num2str(MAX_SIZE),'¡Á',num2str(MAX_SIZE),', overlap = ',num2str(OVERLAP)]);
    for y = index_y
        for x = index_x
            f = gI(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W));
            weight = Confi(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W));
            [h,w] = size(f);
            %disp(['position = (',num2str(y),',',num2str(x),'), size = ',num2str(h),'¡Á',num2str(w)]);
            [r, res] = mosekopt('symbcon echo(0)');
            % Specify the non-conic part of the problem.
            % 5*h*w+1 variables in total: t_{i,j} ,Dx_{i,j} ,Dy_{i,j}, u_{i,j} ,s_{i,j} ,s
            prob.c = [reshape(weight,h*w,1);zeros(4*h*w,1);lambda];
            prob.a = sparse(4*h*w+1,5*h*w+1);
            for j = 1:w
                for i = 1:h
                    index = (j-1)*h+i;
                    % the first h*w rows: -Dx_{i,j} + u_{i+1,j} - u_{i,j} = 0
                    prob.a(index,h*w+index) = -1;
                    if i ~= h
                        prob.a(index,3*h*w+index+1) = 1;
                        prob.a(index,3*h*w+index) = -1;
                    end
                    % the second h*w rows: -Dy_{i,j} + u_{i,j+1} - u_{i,j} = 0
                    prob.a(h*w+index,2*h*w+index) = -1;
                    if j ~= w
                        prob.a(h*w+index,3*h*w+index+h) = 1;
                        prob.a(h*w+index,3*h*w+index) = -1;
                    end
                    % the third h*w rows: u_{i,j} + s_{i,j} >= f_{i,j}
                    prob.a(2*h*w+index,3*h*w+index) = 1;
                    prob.a(2*h*w+index,4*h*w+index) = 1;
                    % the forth h*w rows: u_{i,j} - s_{i,j} <= f_{i,j}
                    prob.a(3*h*w+index,3*h*w+index) = 1;
                    prob.a(3*h*w+index,4*h*w+index) = -1;
                end
            end
            % the last row: sum(s_{i,j}) - s <= 0
            prob.a(4*h*w+1,:) = [zeros(1,4*h*w),ones(1,h*w),-1];
            prob.blc = [zeros(2*h*w,1);reshape(f,h*w,1);-inf*ones(h*w+1,1)];
            prob.buc = [zeros(2*h*w,1);inf*ones(h*w,1);reshape(f,h*w,1);0];
            % blx = -inf, and bux =inf
            % prob.blx = -inf*ones(1,5*h*w+1);
            % prob.bux = inf*ones(1,5*h*w+1);

            % Specify the cones.
            prob.cones.type   = repmat([res.symbcon.MSK_CT_QUAD],1,h*w);
            prob.cones.sub = ones(1,3*h*w);
            for j = 1:w
                for i = 1:h
                    index = (j-1)*h+i;
                    prob.cones.sub(3*index-2:3*index) = index + [0,h*w,2*h*w];
                end
            end
            prob.cones.subptr = 1:3:3*h*w;

            % sepcify parameters
            param = [];
            % Primal feasibility tolerance for the primal solution
            param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = tol;
            % Dual feasibility tolerance for the dual solution
            param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = tol;
            % Relative primal-dual gap tolerance.
            param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = tol;
            % maximum iterations
            param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 1000;

            % Optimize the problem. 
            if en_disp == 0
                [r,res]=mosekopt('minimize echo(0)',prob,param);
            elseif en_disp == 1
                [r,res]=mosekopt('minimize echo(1)',prob,param);
            elseif en_disp == 2
                [r,res]=mosekopt('minimize echo(2)',prob,param);
            else
                [r,res]=mosekopt('minimize echo(3)',prob,param);
            end
            
            
            % update the struct component
            old_u = u(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W));
            change_u = reshape(res.sol.itr.xx(3*h*w+1:4*h*w),h,w);
            old_cnt = cnt(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W));
            new_u = (old_u.*old_cnt+change_u)./(old_cnt+1);
            u(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W)) = new_u;
            cnt(y:min(y+MAX_SIZE-1,H),x:min(x+MAX_SIZE-1,W)) = old_cnt + 1;
        end
    end
    
    v = gI - u;

end