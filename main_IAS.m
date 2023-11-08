function main_IAS(inputFolder,outputFolder,categoryNbr)


SHOW = 0;
path(path,'MiscCodes/')
path(path,'../../../')

Nel = 32; %number of electrodes
z   = 1e-6*ones(Nel,1); %contact impedances
load([inputFolder '/ref.mat']) %reference data from water chamber

rmind       = 1:2*(categoryNbr - 1);

%undersampling of the data for the difficulty levels
vincl = true(31,76); %which measurements are used in the inversion - 31 voltage values for each of the 76 current injections

for ii=1:76
    if any(Injref(rmind,ii) ~= 0)
        vincl(:,ii) = false; %remove all voltage measurements for this current injection
    end
    vincl(rmind,ii) = false; %remove all voltage measurements collected using these electrodes
end
vincl = vincl(:);

load('Mesh_sparse.mat')



files=dir(fullfile(inputFolder,'*.mat')); %one of these is a reference data, the others are from targets with inclusions

for objectno = 1:1 %length(files)-1
    close all
    load([inputFolder '/data' num2str(objectno) '.mat'])
    load(['GroundTruths/true' num2str(objectno) '.mat']);

    % Defining the difference data
    U_coeffs = Uel - Uelref; % DeltaU
    
    solver = EITFEM(Mesh, Mesh2, Injref, Mpat, vincl, []);
    
    %our
    CondVertices = getfield(Mesh,'g');
    CondVertices = CondVertices';
    CondTopol    = getfield(Mesh,'H');
    nc           = size(CondVertices,2); % Number of nodes in the conductivity mesh
    sigma_0      = ones(nc,1); %linearization point
    
    % Finding interior vertices
    d2         = sum((CondVertices(1,:).^2 + CondVertices(2,:).^2),1);
    bound_d2   = (max(max(CondVertices)))^2;
    I_int_vert = find(d2 < bound_d2 * 0.99);
    
    % Construct a grid Laplacian with zero boundary conditions
    % Adjacency matrix
    rows = [];
    cols = [];
    vals = [];
    Edges  = GetEdges(CondTopol);
    for j = 1:nc
        % for the jth node, find all neighbors
        Ij = find(Edges(:,1) == j);
        rows = [rows;j*ones(length(Ij),1)];
        cols = [cols;Edges(Ij,2)];
        vals = [vals;ones(length(Ij),1)];
        Jj = find(Edges(:,2) == j);
        rows = [rows;j*ones(length(Jj),1)];
        cols = [cols;Edges(Jj,1)];
        vals = [vals;ones(length(Ij),1)];
    end
    Adj   = sparse(rows,cols,vals);
    dd    = sum(Adj,2);
    Lapl  = Adj - diag(dd);
    Lapl  = diag(1./dd)*Lapl;
    
    % Discarding the rows and columns of the grid laplacian corresponding to
    % boundary nodes
    Lapl_red = Lapl(I_int_vert,I_int_vert);%
    
    % Background conductivity. Find by fitting to data
    sigma_0 = 0.79*sigma_0;
    s_c     = zeros(nc,1); % Initial value of the perturbation
    sigma_c = sigma_0 + s_c;
    
    %get linear model using provided codes
%     b_0     = solver.SolveForward(sigma_c,z); % input sigma_c
    A       = solver.Jacobian(sigma_c,z);
    
    % Estimated noise level. Result of testing
    omega = 0.004; %0.006; %omega_0;
    
    % Setting up the likelihood model
    y = U_coeffs(vincl(:)); %!!!!!!
    % Whitening
    yw = (1/omega)*y;
    Aw = (1/omega)*A;
    % Assuming that the boundary nodal values are known, we remove those nodes
    Ared   = Aw(:,I_int_vert);
    
    % Phase I:
    % Solve the linearized problem using second order smoothness prior
    
    delta = 1000;
    sc_0 = (Ared'*Ared + delta*Lapl_red'*Lapl_red)\(Ared'*yw);
    
    sig_temp = zeros(nc,1);
    sig_temp(I_int_vert) = sc_0;
    sigma_c = sigma_0+sig_temp;
    
    %%%%
    sigma_pixgrid = interpolateRecoToPixGrid(sig_temp,Mesh);

    if SHOW == 1
        %plot the reconstruction
        sgplot = sigmaplotter(Mesh,5,'parula');
        sgplot.basic2Dplot(sigma_c,{['Prel. reconstruction ' num2str(objectno)]});
        %interpolate the reconstruction into a pixel image
        figure(1), imagesc(sigma_pixgrid), colorbar, axis image
        set(gcf,'Units','normalized','OuterPosition',[0.3 0.6 0.3 0.4])
    end
    %treshold the image histogram using Otsu's method
    [level,x] = Otsu2(sigma_pixgrid(:),256,7);
    reconstruction = zeros(size(sigma_pixgrid));
    ind1 = find(sigma_pixgrid < x(level(1)));
    ind2 = find(sigma_pixgrid >= x(level(1)) & sigma_pixgrid <= x(level(2)));
    ind3 = find(sigma_pixgrid > x(level(2)));
    
    %Check which class is the background (assumed to be the one with the
    %most pixels in it).
    [bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
    switch bgclass
        case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
            reconstruction([ind2]) = 2;
            reconstruction([ind3]) = 2;
        case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
            reconstruction([ind1]) = 1;
            reconstruction([ind3]) = 2;
        case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
            reconstruction([ind1]) = 1;
            reconstruction([ind2]) = 1;
    end
    
    if SHOW == 1
        %plot the segmented reconstruction in the pixel grid
        figure(2),
        drawnow
        imagesc(reconstruction)
        colormap gray
        colorbar, axis image
        title('segmented prel. reconstruction')
        set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
        % True target
        figure(100)
        drawnow
        imagesc(truth)
        colormap gray
        colorbar, axis image
        title('original target')
        set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
    end
    % For the IAS algorithm:
    % Finite difference matrix mapping nodal values to edge values
    L_full = GetIncrementMatrix(CondVertices,CondTopol,Edges);
    L      = L_full(:,I_int_vert);
    % Some of the edges are between boundary nodes. When leaving out boundary nodes,
    % a zero row in L results. These edges need to be removed.
    aux          = sum(abs(L'));
    I_free_edges = find(aux>0); % Index to free edges
    L            = L(I_free_edges,:);
    n_v          = length(I_int_vert); % Number of interior nodes, or degrees of freedom
    n_e          = size(L,1);          % Number of free edges
    
    % 2. IAS algorithm
    % ================
    
    % Parameters for phase 1 (gamma hypermodel)
    eta        = 0.001;
    maxiter_1  = 10;
    max_CGLS_1 = 200;
    stopcond_1 = 0;%5e-3;
    
    % Choose phase 2 model, and set iteration parameters. Observe: If r = 1,
    % phase 2 will be skipped
    r          = 1/2;  % Choose from 1,1/2,-1/2 and -1
    maxiter_2  = 10;
    max_CGLS_2 = 200;
    stopcond_2 = 0;%5e-3;
    
    % Compute parameters for phase 2
    if r ~= 1
        if r == 1/2
            m    = 1+3/(2*eta);
            beta = (6*m+1 + sqrt(48*m+1))/(2*m-2);
        elseif r == -1/2
            m    = 1+3/(2*eta);
            beta = (6+3*m + sqrt(m^2+80*m))/(2*m-2);
        elseif r == -1
            beta = 1 + 5*eta/3;
        end
        %         vartheta = theta_star*eta/((beta - 3/(2*r))^(1/r));
    end
    
    % sensitivities
    C           = .1;%0.1;%10; %manually set
    [Q_L,R_L]   = qr(L,0);
    diagRL = diag(R_L);
    R_L    = diag(diagRL+1e-2);
    AL_2   = Ared*(R_L\Q_L');
    theta_star  = C./(sum(AL_2.^2,1).^1.2);
    theta_star  = theta_star';
    if  r ~= 1
        vartheta    = theta_star*eta/((beta - 3/(2*r))^(1/r));
    end
    
    % Initializing
    iters = 1; % dummy
    z     = L*sc_0;
    theta = theta_star.*(eta/2*ones(n_e,1) + sqrt(eta^2/4*ones(n_e,1) + (z.^2)./(2*theta_star.*ones(n_e,1))));
    %theta = theta_star;
    w0      = zeros(n_e,1);
    count   = 0;
    d_th    = 1e2;
    % Running the IAS with r = 1
    while count < maxiter_1 && d_th > stopcond_1
        count = count + 1;
        % Updating z
        L_th = diag(1./sqrt(theta))*L;
        [Q,R] = qr(L_th,0); % Economy size QR
        AL = Ared*(R\Q');
        [w,n_CGLS] = CGLS_update_x0(yw,AL,w0,max_CGLS_1);
        z = diag(sqrt(theta))*w;
        % Updating theta
        theta_old = theta;
        theta     = theta_star.*(eta/2*ones(n_e,1) + sqrt(eta^2/4*ones(n_e,1) + (z.^2)./(2*theta_star.*ones(n_e,1))));
        d_th      = norm(theta - theta_old,'inf')/norm(theta_old,'inf');
        
        sigma_pixgrid = interpolateRecoToPixGrid(sig_temp,Mesh);
        sig_temp = zeros(nc,1);
        sig_temp(I_int_vert) = L\z;
        
        if SHOW == 1
            fprintf('\n LIN = %3d , r = %2.1f , COUNT = %4d , N_CGLS = %4d , D_TH = %6.5f',iters,1,count,n_CGLS,d_th)
            
            % Plot the current reconstruction
            sgplot = sigmaplotter(Mesh,5,'parula');
            sgplot.basic2Dplot(sigma_c,{['IAS r=1 reconstruction ' num2str(objectno)]});
            %interpolate the reconstruction into a pixel image
            figure(1), imagesc(sigma_pixgrid), colorbar, axis image
            set(gcf,'Units','normalized','OuterPosition',[0.3 0.6 0.3 0.4])
        end
        %treshold the image histogram using Otsu's method
        [level,x] = Otsu2(sigma_pixgrid(:),256,7);
        reconstruction = zeros(size(sigma_pixgrid));
        ind1 = find(sigma_pixgrid < x(level(1)));
        ind2 = find(sigma_pixgrid >= x(level(1)) & sigma_pixgrid <= x(level(2)));
        ind3 = find(sigma_pixgrid > x(level(2)));
        %Check which class is the background (assumed to be the one with the
        %most pixels in it).
        [bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
        switch bgclass
            case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
                reconstruction([ind2]) = 2;
                reconstruction([ind3]) = 2;
            case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
                reconstruction([ind1]) = 1;
                reconstruction([ind3]) = 2;
            case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
                reconstruction([ind1]) = 1;
                reconstruction([ind2]) = 1;
        end
        %plot the segmented reconstruction in the pixel grid
        if SHOW == 1
            figure(2+count),
            drawnow
            imagesc(reconstruction)
            colormap gray
            colorbar, axis image
            title('segmented prel. reconstruction')
            set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
        end
    end
    
    % Phase 2
    if r ~= 1
        count = 0;
        d_th = 1e2;
        while count < maxiter_2 && d_th > stopcond_2
            count = count + 1;
            % Updating z
            L_th  = diag(1./sqrt(theta))*L;
            [Q,R] = qr(L_th,0); % Economy size QR
            AL = Ared*(R\Q');
            [w,n_CGLS] = CGLS_update_x0(yw,AL,w0,max_CGLS_2);
            z = diag(sqrt(theta))*w;
            % Updating theta
            theta_old = theta;
            zeta   = z./sqrt(vartheta);
            zeta   = abs(zeta);
            xi     = GenGammaUpdate1D(zeta,r,beta);
            theta  = xi .* vartheta;
            d_th      = norm(theta - theta_old,'inf')/norm(theta_old,'inf');
            if SHOW == 1
                fprintf('\n LIN = %3d , r = %2.1f , COUNT = %4d , N_CGLS = %4d , D_TH = %6.5f',iters,r,count,n_CGLS,d_th)
            end
            % Plot the current reconstruction
            
            sig_temp = zeros(nc,1);
            sig_temp(I_int_vert) = L\z;
            sigma_pixgrid = interpolateRecoToPixGrid(sig_temp,Mesh);

            if SHOW == 1
                sgplot = sigmaplotter(Mesh,5,'parula');
                sgplot.basic2Dplot(sigma_c,{['IAS r=1 reconstruction ' num2str(objectno)]});
                %interpolate the reconstruction into a pixel image
                figure(1), imagesc(sigma_pixgrid), colorbar, axis image
                set(gcf,'Units','normalized','OuterPosition',[0.3 0.6 0.3 0.4])
            end
            %treshold the image histogram using Otsu's method
            [level,x] = Otsu2(sigma_pixgrid(:),256,7);
            reconstruction = zeros(size(sigma_pixgrid));
            ind1 = find(sigma_pixgrid < x(level(1)));
            ind2 = find(sigma_pixgrid >= x(level(1)) & sigma_pixgrid <= x(level(2)));
            ind3 = find(sigma_pixgrid > x(level(2)));
            %Check which class is the background (assumed to be the one with the
            %most pixels in it).
            [bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
            switch bgclass
                case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
                    reconstruction([ind2]) = 2;
                    reconstruction([ind3]) = 2;
                case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
                    reconstruction([ind1]) = 1;
                    reconstruction([ind3]) = 2;
                case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
                    reconstruction([ind1]) = 1;
                    reconstruction([ind2]) = 1;
            end
            %plot the segmented reconstruction in the pixel grid
            if SHOW == 1
                figure(2+count),
                drawnow
                imagesc(reconstruction)
                colormap gray
                colorbar, axis image
                title('segmented prel. reconstruction')
                set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
            end
        end
        
    end
    
%     
%     return
%     
%     
%     
%     % Solving for nodal values and plotting the result
%     
%     x               = L\z;
%     d_s             = zeros(nc,1);
%     d_s(I_int_vert) = x;
%     s_c             =  d_s;
%     sigma = sigma_0 + s_c;
%     OUT_SIG(:,iters) = sigma;
%     
%     nc_elem = size(CondTopol,1);
%     
%     
%     sigma_c = sigma;
%     
%     %end
%     count = count+1;
%     
%     figure(4+count)
%     trisurf(CondTopol,CondVertices(1,:),CondVertices(2,:),sigma)
%     axis('square')
%     shading flat
%     hold on
%     view(0,90)
%     colorbar
%     set(gca,'FontSize',15)
%     hold off
%     drawnow
%     
%     %plot the reconstruction
%     sgplot = sigmaplotter(Mesh,5,'parula');
%     sgplot.basic2Dplot(sigma,{['IAS reconstruction ' num2str(objectno)]});
%     
%     %interpolate the reconstruction into a pixel image
%     %sigma_pixgrid = interpolateRecoToPixGrid(sigma,Mesh);
%     sigma_pixgrid = interpolateRecoToPixGrid(s_c,Mesh);
%     figure(6), imagesc(sigma_pixgrid), colorbar, axis image
%     set(gcf,'Units','normalized','OuterPosition',[0.3 0.6 0.3 0.4])
%     
%     
%     %treshold the image histogram using Otsu's method
%     [level,x] = Otsu2(sigma_pixgrid(:),256,7);
%     
%     reconstruction = zeros(size(sigma_pixgrid));
%     
%     ind1 = find(sigma_pixgrid < x(level(1)));
%     ind2 = find(sigma_pixgrid >= x(level(1)) & sigma_pixgrid <= x(level(2)));
%     ind3 = find(sigma_pixgrid > x(level(2)));
%     
%     %Check which class is the background (assumed to be the one with the
%     %most pixels in it).
%     [bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
%     switch bgclass
%         case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
%             reconstruction([ind2]) = 2;
%             reconstruction([ind3]) = 2;
%         case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
%             reconstruction([ind1]) = 1;
%             reconstruction([ind3]) = 2;
%         case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
%             reconstruction([ind1]) = 1;
%             reconstruction([ind2]) = 1;
%     end
%     
%     %plot the segmented reconstruction in the pixel grid
%     figure(30+count),
%     drawnow
%     imagesc(reconstruction)
%     colormap gray
%     colorbar, axis image
%     title('segmented IAS reconstruction')
%     set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
%     
%     figure(100)
%     drawnow
%     imagesc(truth)
%     colormap gray
%     colorbar, axis image
%     title('original target')
%     set(gcf,'Units','normalized','OuterPosition',[0.6 0.2 0.3 0.4])
%     
%     
    
    
%     s = scoringFunction(truth, reconstruction);
%     
%     
%     disp(['Score from target ' num2str(objectno) ' = ' num2str(s)])
%     
    
    save([outputFolder '/' num2str(objectno) '.mat'],'reconstruction')
end
end








    
