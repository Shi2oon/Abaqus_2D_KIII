function [uX,uY,uZ,Flag,Tolerance,MaxIter] = removeoutliers(method,thresholdfactor,mask,uX,uY,varargin)
%REMOVEOUTLIERS Removes displacement outliers
% [uX,uY,uZ]=REMOVEOUTLIERS(method,thresholdfactor,mask,uX,uY) Removes
%   displacement outliers using a finite element analysis (FEA) framework.
%   REMOVEOUTLIERS applies the displacements contained in FIELD as boundary
%   conditions to compute the nodal forces. Under the assumption that
%   internal forces should be zero, any non-zero nodal forces are used to
%   identify outliers. REMOVEOUTLIERS uses an iterative process to identify
%   and remove outlier displacements to avoid convergence to a
%   linear-elastic solution.
%
%   REMOVEOUTLIERS('method','thresholdfactor') specifies the outlier
%   identification method. 'Method' can be set to 'median' or 'grubbs'.
%   'thresholdfactor' specifies the outlier detection thresholds; for
%   'grubbs' this a scalar between 0 and 1. For 'median' this is a
%   nonnegative scalar. See MATLAB’s 'isoutlier' documentation for more
%   information.
%
%   'mask' defines a mask of ones and 'nan' values that has the same size
%   as uX. Data points set to nan are ignored and assumed not to contribute
%   to the continuity of the solution. 'uX', 'uY' (and 'uZ') are the
%   respective displacement data. 'mask', 'uX', 'uY' (and 'uZ') mussed be
%   in a 'meshgrid' format
%
%   REMOVEOUTLIERS(...,'MaxIter',n,'Tol',tol) specifies the exit criterion.
%   'MaxIter' specifies the maximum number of iterations. 'Tol' specifies
%   the tolerance between each iteration. Default values are n=10, tol=0.1.
%
%   [uX,uY,uZ,Flag,Tolerance,Iter]=REMOVEOUTLIERS provides exit conditions.
%   Flag, 1. convergence criteria was met, 2. a zero median was obtained,
%   3. maximum number of iterations reached, and 4. all data has been
%   replaced. Tolerance is the final tolerance value, and ‘Iter’ the total
%   number of iterations.

%	History: 
%   15 May 2020, Thorsten Becker. Created. 
%   18 June 2020, Thorsten Becker. Included 'isoulier' functionality. 
%   7 June 2021, Thorsten Becker. Simplified code.


%% Check inputs
% DataSize
DataSize=size(mask);
% Defaults
defaultTol=0.1;
defaultMaxIter=10;
defaultuZ=zeros(DataSize);
% Input parser
p=inputParser;
addRequired(p,'method',@ischar)
addRequired(p,'ThresholdFactor',@(x) all([isnumeric(x),numel(x)==1]))
addRequired(p,'mask',@(x) all([isnumeric(x),size(x)==DataSize]))
addRequired(p,'uX',@(x) all([isnumeric(x),size(x)==DataSize]))
addRequired(p,'uY',@(x) all([isnumeric(x),size(x)==DataSize]))
addOptional(p,'uZ',defaultuZ,@(x) all([isnumeric(x),size(x)==DataSize]))
addParameter(p,'Tol',defaultTol,@(x) all([isnumeric(x),numel(x)==1]));
addParameter(p,'MaxIter',defaultMaxIter,@(x) all([isnumeric(x),numel(x)==1]));
parse(p,method,thresholdfactor,mask,uX,uY,varargin{:});
% Data
mask=p.Results.mask;
uX=p.Results.uX;
uY=p.Results.uY;
uZ=p.Results.uZ;
% Variables
Method=lower(p.Results.method);
ThresholdFactor=p.Results.ThresholdFactor;
Tolerance=p.Results.Tol;
MaxIter=p.Results.MaxIter;
% Check data type
if numel(DataSize)==2 && range(uZ(:))==0
    type=1;     % 2D-DIC data  
elseif numel(DataSize)==2 && range(uZ(:))>0
    type=2;     % Stereo-DIC data
elseif numel(DataSize)==3
    type=3;     % DVC data
end

%% Outlier filter
% Pad array
if type==1      
    u.x=padarray(uX,[1,1],nan);
    u.y=padarray(uY,[1,1],nan);
    f.x=padarray(nan(DataSize),[1,1],0);
    f.y=padarray(nan(DataSize),[1,1],0);
    m=padarray(mask,[1,1],'replicate','both');
    edges=padarray(zeros(DataSize),[1,1],1);
elseif type==2  
    u.x=padarray(uX,[1,1,1],nan);
    u.y=padarray(uY,[1,1,1],nan);
    u.z=padarray(uZ,[1,1,1],nan);
    f.x=padarray(nan(DataSize),[1,1,1],0);
    f.y=padarray(nan(DataSize),[1,1,1],0);
    f.z=padarray(nan(DataSize),[1,1,1],0);
    m=padarray(mask,[1,1,1],'replicate','both');
    edges=padarray(zeros(DataSize),[1,1,1],1);
elseif type==3 
    u.x=padarray(uX,[1,1,1],nan);
    u.y=padarray(uY,[1,1,1],nan);
    u.z=padarray(uZ,[1,1,1],nan);
    f.x=padarray(nan(DataSize),[1,1,1],0);
    f.y=padarray(nan(DataSize),[1,1,1],0);
    f.z=padarray(nan(DataSize),[1,1,1],0);
    m=padarray(mask,[1,1,1],'replicate','both');
    edges=padarray(zeros(DataSize),[1,1,1],1);
end
DataSize=size(u.x);
% Edge and masked indices.
indedge=find(edges==1);
indmask=find(isnan(m));
% Index identifiers.
indnofilter=unique([indmask;indedge]);
inddofilter=setdiff(1:prod(DataSize),indnofilter)';
matfiltered=false(DataSize);
numfilt=zeros(1,MaxIter+1);
maxnumfilt=numel(inddofilter);
% Initialize
flag=0;
kk=0;
fbias=zeros(1,MaxIter+1);
fvar=zeros(1,MaxIter+1);
ftol=zeros(1,MaxIter+1);
% Show output
fprintf(' Iter | Points replaced | Bias (N) | Var. (N) | Tol.   \n');
% Iteratively remove outliers.
while flag==0
    % Reset filters
    indfiltered=[];
    filtered=false(DataSize);
    % Count
    kk=kk+1;
    % Remove edge & masked regions.
    u.x(indnofilter)=nan;
    u.y(indnofilter)=nan;
    if type>1, u.z(indnofilter)=nan; end
    f.x(indnofilter)=0;
    f.y(indnofilter)=0;
    if type>1, f.z(indnofilter)=0; end
    % Find outliers.
    switch Method
        % All data.
        case 'median'
            fx=isoutlier(f.x(inddofilter),'median','ThresholdFactor',ThresholdFactor);
            fy=isoutlier(f.y(inddofilter),'median','ThresholdFactor',ThresholdFactor);
            if type>1, fz=isoutlier(f.z(inddofilter),'median','ThresholdFactor',ThresholdFactor); end
            if type==1, filtered(inddofilter)=fy|fx; else, filtered(inddofilter)=fy|fx|fz; end
        case 'grubbs'
            fx=isoutlier(f.x(inddofilter),'grubbs','ThresholdFactor',ThresholdFactor);
            fy=isoutlier(f.y(inddofilter),'grubbs','ThresholdFactor',ThresholdFactor);
            if type>1, fz=isoutlier(f.z(inddofilter),'grubbs','ThresholdFactor',ThresholdFactor); end
            if type==1, filtered(inddofilter)=fy|fx; else, filtered(inddofilter)=fy|fx|fz; end
        otherwise
            error('Method not recognised.')
    end
    % Outlier indicies.
    indfiltered=setdiff(find(filtered),indnofilter);
    matfiltered(indfiltered)=true;
    % Remove outliers
    u.x(indfiltered)=nan;
    u.y(indfiltered)=nan;
    if type>1, u.z(indfiltered)=nan; end
    f.x(:)=nan; f.x(indfiltered)=0; f.x(indnofilter)=0;
    f.y(:)=nan; f.y(indfiltered)=0; f.y(indnofilter)=0;
    if type>1, f.z(:)=nan; f.z(indfiltered)=0; f.z(indnofilter)=0; end
    % FEA
    [u,f]=gridFEA(u,f,1,m);
    f.x=round(f.x,10); f.y=round(f.y,10); if type>1, f.z=round(f.z,10); end
    % Force vector length.
    if type>1
        fv=(f.x(inddofilter).^2+f.y(inddofilter).^2+f.z(inddofilter).^2).^(1/2);
    else
        fv=(f.x(inddofilter).^2+f.y(inddofilter).^2).^(1/2);
    end
    % Stats.
    numfilt(kk)=numel(indfiltered);
    fbias(kk)=nanmedian(fv(fv~=0));
    fvar(kk)=nanstd(fv(fv~=0));
    if kk==1
        ftol(kk)=nan;
        ntol(kk)=nan;
    else
        ftol(kk)=abs((fbias(kk-1)-fbias(kk))/fbias(1));
        ntol(kk)=numfilt(kk)-numfilt(kk-1);
    end
    % Show output.
    fprintf(' %4.0f | %7.0f/%7.0f | %7.2e | %7.2e | %7.4f \n',kk-1,numfilt(kk),maxnumfilt,fbias(kk),fvar(kk),ftol(kk));
    % Convergence criteria.
    if kk~=1
        if (ftol(kk)<=Tolerance || ntol(kk)==0)     % C2
            flag=1;
            if nargout<=3, disp(' Function tolerance reached.'); end
        elseif round(fbias(kk),10)==0   % C1
            flag=2;
            if nargout<=3, disp(' Bias equals zero.'); end
        elseif kk==MaxIter  % C3
            flag=3;
            if nargout<=3, disp(' Maxmimum iterations reached.'); end
        elseif numfilt(kk)==maxnumfilt
            flag=4;
            if nargout<=3, disp(' Total number of points filtered equals the total number of analysis points.'); end
        end
    end
end
% Output.
if type==1
    uX=u.x(2:end-1,2:end-1).*mask;
    uY=u.y(2:end-1,2:end-1).*mask;
    uZ=[];
else
    uX=u.x(2:end-1,2:end-1,2:end-1).*mask;
    uY=u.y(2:end-1,2:end-1,2:end-1).*mask;
    uZ=u.z(2:end-1,2:end-1,2:end-1).*mask; 
end
MaxIter=kk;
Tolerance=ftol(kk);
Flag=flag;
end

%% Support functions
function varargout = gridFEA(varargin)
%GRIDFEA Grid based Finite Element Analysis (FEA) tool.
%   GRIDFEA is an efficient grid based FEA tool with nodes equally spaced
%   using MESHGRID. GRIDFEA accepts 2D planar (plane stress), surface
%   (membrane - feature not coded yet) and 3D (volume) data to folve the
%   equation [K]{u}={f}.
%
%   GRIDFEA requires structured boundary conditions for displacement (u.x,
%   u.y or u.x, u.y, u.z), force (f.x, f.y or f.x, f.y, f.z), a mask (size
%   u or f), and gridspacing (scalar as distance between data points). u,f,
%   and mask are defined using either values or NaNs (i.e. not given).
%
%   [u,f]=GRIDFEA(u,f,GridSpacing,Mask) computes displacements and forces
%   for boundary conditions defined in u and f with a grid spacing of
%   GridSpacing and using a mask. GRIDFEA assumes a default Young's modulus
%   of 10e3 and Poisson's ratio of 0.25.
%
%   [u,f,e,s]=GRIDFEA(u,f,GridSpacing,mask,'E',E,'v',v,...
%       'ComputeDisplacementGradients','On','ComputeStrain','On',...
%       'ComputeStress','On') defines Young's modulus and Poisson's ratio
%       and outputs the symmetric strain and stress tensor respectively.
%
%   Plane stress example:
%       % Define domain.
%       GridSpacing=1;
%       DataSize=[20,40];
%       [p.x,p.y]=meshgrid(0:GridSpacing:DataSize(2)-1,0:GridSpacing:DataSize(1)-1);
%       mask=ones(DataSize); mask(1:9,21)=nan;
%       % Boundary conditions, displacement.
%       u.x=nan(DataSize); u.x(:,1)=0;
%       u.y=nan(DataSize); u.y(:,1)=0;
%       % Boundary conditions, force.
%       f.x=nan(DataSize); f.x(2:end-1,end)=100; f.x([1,end],end)=100/2;
%       f.y=nan(DataSize);
%       % girdFEA
%       [u,f,~,s]=gridFEA(u,f,GridSpacing,mask,'ComputeStress','On');
%       % plot
%       mask(mask==0)=NaN;
%       figure; hold on;
%       surf(p.x+u.x,p.y+u.y,zeros(DataSize).*mask,s.xx)
%       quiver(p.x+u.x,p.y+u.y,f.x,f.y,'k','LineWidth',2,'AutoScaleFactor',1)
%       plot(reshape(p.x+u.x,1,[]),reshape(p.y+u.y,1,[]),'.k')
%       % Set figure properties.
%       view(2); axis equal; axis tight; box on; grid on;
%       title('Stress in xx')
%       xlabel('position [mm]'); ylabel('position [mm]')
%       hcb=colorbar; title(hcb,'[MPa]')
%
%   Volume example:
%       GridSpacing=1;
%       DataSize=[20,40,5];
%       [p.x,p.y,p.z]=meshgrid(0:GridSpacing:DataSize(2)-1,0:GridSpacing:DataSize(1)-1,0:GridSpacing:DataSize(3)-1);
%       mask=ones(DataSize); mask(1:9,21,:)=nan;
%       % Boundary conditions, displacement.
%       u.x=nan(DataSize); u.x(:,1,:)=0;
%       u.y=nan(DataSize); u.y(:,1,:)=0;
%       u.z=nan(DataSize); u.z(:,1,:)=0;
%       % Boundary conditions, force.
%       f.x=nan(DataSize); f.x(2:end-1,end,:)=10; f.x([1,end],end,:)= 10/2;
%       f.y=nan(DataSize);
%       f.z=nan(DataSize);
%       % girdFEA
%       [u,f,~,s]=gridFEA(u,f,GridSpacing,mask,'ComputeStress','On');
%       % plot
%       mask(mask==0)=NaN;
%       figure; hold on;
%       slice(p.x,p.y,p.z,s.xx,20,10,2)
%       quiver3(p.x,p.y,p.z,f.x,f.y,f.z,'k','LineWidth',2,'AutoScaleFactor',1)
%       plot3(reshape(p.x,1,[]),reshape(p.y,1,[]),reshape(p.z,1,[]),'.k')
%       % Set figure properties.
%       view(3); axis equal; axis tight; box on; grid on;
%       title('Stress in xx')
%       xlabel('position [mm]'); ylabel('position [mm]')
%       hcb=colorbar; title(hcb,'[MPa]')
%
%   Version:
%       Created by Thorsten Becker, June 2018

%% Function preamble.
% Check format of input data.
if numel(varargin)<4;  error('Not enough input arguments.'); end
if ~isstruct(varargin{1}); error('u input structure not correct'); end
if ~all(isfield(varargin{1},{'x','y'})); error('u input structure not correct'); end
if ~isstruct(varargin{2}); error('f input structure not correct'); end
if ~all(isfield(varargin{2},{'x','y'})); error('f input structure not correct'); end
if ~isnumeric(varargin{3}); error('GridSpacing input not correct'); end
if ~numel(varargin{3})==1; error('GridSpacing input not correct'); end
if ~isnumeric(varargin{4}); error('Mask input not correct'); end
% Input data dimensions.
if ~isfield(varargin{1},'z')&&~isfield(varargin{2},'z')&&...
        ndims(varargin{1}.x)&&ndims(varargin{1}.y)&&ndims(varargin{2}.x)&&ndims(varargin{2}.y)==2
    % Check if gridded .
    if ~isMeshGrid(varargin{1}.x,varargin{1}.y); error('u not in gridded format.'); end
    if ~isMeshGrid(varargin{2}.x,varargin{2}.y); error('f not in gridded format.'); end
    % 2D data.
    UX=varargin{1}.x;
    UY=varargin{1}.y;
    FX=varargin{2}.x;
    FY=varargin{2}.y;
    MASK=varargin{4};
    % Datasize.
    DataSize=[size(UX),1];
    DataType='PlaneStress';
    GridSpacing=varargin{3};
elseif isfield(varargin{1},'z')&&isfield(varargin{2},'z')&&...
        ndims(varargin{1}.x)&&ndims(varargin{1}.y)&&ndims(varargin{2}.x)&&ndims(varargin{2}.y)==2
    % Check if gridded .
    if ~isMeshGrid(varargin{1}.x,varargin{1}.y); error('u not in gridded format.'); end
    if ~isMeshGrid(varargin{1}.x,varargin{1}.z); error('u not in gridded format.'); end
    if ~isMeshGrid(varargin{2}.x,varargin{2}.y); error('f not in gridded format.'); end
    if ~isMeshGrid(varargin{2}.x,varargin{2}.z); error('f not in gridded format.'); end
    % Membrane data.
    UX=varargin{1}.x;
    UY=varargin{1}.y;
    UZ=varargin{1}.z;
    FX=varargin{2}.x;
    FY=varargin{2}.y;
    FZ=varargin{2}.z;
    MASK=varargin{4};
    % Datasize.
    DataSize=[size(UX),1];
    DataType='Membrane';
    GridSpacing=varargin{3};
elseif isfield(varargin{1},'z')&&isfield(varargin{1},'z')&&...
        ndims(varargin{1}.x)&&ndims(varargin{1}.y)&&ndims(varargin{2}.x)&&ndims(varargin{2}.y)==3
    % Check dimensions
    if ~isequal(size(varargin{1}.x),size(varargin{1}.y),size(varargin{1}.z),...
            size(varargin{2}.x),size(varargin{2}.y),size(varargin{2}.z),...
            size(varargin{4}))
        error('Input dimensions not equal.')
    end
    % Volume data.
    UX=varargin{1}.x;
    UY=varargin{1}.y;
    UZ=varargin{1}.z;
    FX=varargin{2}.x;
    FY=varargin{2}.y;
    FZ=varargin{2}.z;
    MASK=varargin{4};
    % Datasize.
    DataSize=size(UX);
    DataType='Volume';
    GridSpacing=varargin{3};
else
    error('Input format not correct.')
end
% Function defaults.
default.E=300e3;
default.v=0.3;
default.compstress='Off';
default.compstrain='Off';
default.Seem=ones(DataSize);
p=inputParser;
addOptional(p,'Seem',default.Seem,@(x) isnumeric(x));
addParameter(p,'E',default.E,@isnumeric);
addParameter(p,'v',default.v,@isnumeric);
addParameter(p,'ComputeStress',default.compstress,@ischar);
addParameter(p,'ComputeStrain',default.compstrain,@ischar);
% Parse input arguments.
parse(p,varargin{5:end});
% Material properties.
E=p.Results.E;
v=p.Results.v;
% Ouput options.
SEEM=p.Results.Seem;
compstress=p.Results.ComputeStress;
compstrain=p.Results.ComputeStrain;

%% FEA discretisation.
switch DataType
    case 'PlaneStress'
        [K,u,f,dofids,elids,elmask]=FEA2D(DataSize,UX,UY,FX,FY,GridSpacing,MASK,SEEM,E,v);
    case 'Membrane'
        error('Membrane feature not coded. On the to-do list.')
    case 'Volume'
        [K,u,f,dofids,elids,elmask]=FEA3D(DataSize,UX,UY,UZ,FX,FY,FZ,GridSpacing,MASK,SEEM,E,v);
end

%% Boundary conditions.
% Check that BCs are defined correctly.
if any(ismember(find(~isnan(u)),find(~isnan(f)))); error('Multiple boundary conditions defined'); end
% Fixed (displacment) BCs.
fixddofs=find(~isnan(u));
% Free (force) BCs.
freedofs=setdiff(1:numel(u),fixddofs)';
% Set non-defined freedofs to zero. Change to apply gravity, etc.
f(isnan(f))=0;

%% Solve FEA.
% Solve displacements.
if numel(dofids)>1e5
    % For large problem use pcg operator.
    [u(freedofs),flag]=pcg(K(freedofs,freedofs),(f(freedofs)-K(freedofs,fixddofs)*u(fixddofs)),...
        1e-10,1e4,diag(diag(K(freedofs,freedofs))));
    if flag~=0, warning('pcg solver did not converge.'); end
else
    % For small problem use backslash operator.
    u(freedofs)=K(freedofs,freedofs)\(f(freedofs)-K(freedofs,fixddofs)*u(fixddofs));
end
% Solve forces.
f(fixddofs)=K(fixddofs,freedofs)*u(freedofs)+K(fixddofs,fixddofs)*u(fixddofs);

%% Gridded discretisation.
switch DataType
    case 'PlaneStress'
        U=reshape(u(dofids),2,[]);
        varargout{1}.x=reshape(U(1,:),DataSize);
        varargout{1}.y=reshape(U(2,:),DataSize);
        F=reshape(f(dofids),2,[]);
        varargout{2}.x=reshape(F(1,:),DataSize);
        varargout{2}.y=reshape(F(2,:),DataSize);
    case 'Membrane'
        error('Membrane feature not coded. On the to-do list.')
    case 'Volume'
        U=reshape(u(dofids),3,[]);
        varargout{1}.x=(reshape(U(1,:),DataSize));
        varargout{1}.y=(reshape(U(2,:),DataSize));
        varargout{1}.z=(reshape(U(3,:),DataSize));
        F=reshape(f(dofids),3,[]);
        varargout{2}.x=(reshape(F(1,:),DataSize));
        varargout{2}.y=(reshape(F(2,:),DataSize));
        varargout{2}.z=(reshape(F(3,:),DataSize));
end

%% Calculate displacement gradients.
% Displacement derivatives.
UXX=zeros(DataSize); UXY=zeros(DataSize); UXZ=zeros(DataSize);
UYX=zeros(DataSize); UYY=zeros(DataSize); UYZ=zeros(DataSize);
UZX=zeros(DataSize); UZY=zeros(DataSize); UZZ=zeros(DataSize);
% Number of elements neighbouring each node.
elnN=zeros(DataSize);
switch DataType
    case 'PlaneStress'
        e1=[-1,1,1,-1]';
        e2=[-1,-1,1,1]';
        Bx=[e2-1,-e2+1,e2+1,-e2-1]/(2*GridSpacing);
        By=[e1-1,-e1-1,e1+1,-e1+1]/(2*GridSpacing);
        % Loop through elements.
        for elnum=1:size(elids,1)
           if elmask(elnum)==1
                % Displacement at nodes [x1,y1,x2,y2,x3,y3,x4,y4].
                elu=u(reshape([dofids(elids(elnum,:)*2-1)';dofids(elids(elnum,:)*2)'],[],1));
                % Displacment derivatives at nodes (row=[uxx,uyy,uxy,uyx],col=nodes).
                eldu=[Bx*elu(1:2:end-1),By*elu(1:2:end-1),...
                      Bx*elu(2:2:end)  ,By*elu(2:2:end)  ]';
                % Track number of elements neighbouring each node.
                elnN(elids(elnum,:))=elnN(elids(elnum,:))+1;
                % Gridded displacement derivatives.
                UXX(elids(elnum,:))=UXX(elids(elnum,:))+eldu(1,:);
                UXY(elids(elnum,:))=UXY(elids(elnum,:))+eldu(2,:);
                UYX(elids(elnum,:))=UYX(elids(elnum,:))+eldu(3,:);
                UYY(elids(elnum,:))=UYY(elids(elnum,:))+eldu(4,:);
           end
        end
        % Average nodal data by number of elements neighbouring each node.
        UXX=UXX./elnN; UXY=UXY./elnN;
        UYX=UYX./elnN; UYY=UYY./elnN; 
        % Exclude Mask and Seem nodes.
        varargout{1}.xx=UXX.*MASK; 
        varargout{1}.xy=UXY.*MASK;
        varargout{1}.yx=UYX.*MASK; 
        varargout{1}.yy=UYY.*MASK;
        % Displacment derivative in ZZ.
        UZZ=v/(v-1)*(UXX+UYY);
        varargout{1}.zz=UZZ;
    case 'Membrane'
        error('Membrane feature not coded. On the to-do list.')
    case 'Volume'
        e1=[-1, 1, 1,-1,-1, 1, 1,-1]';
        e2=[-1,-1, 1, 1,-1,-1, 1, 1]';
        e3=[-1,-1,-1,-1, 1, 1, 1, 1]';
        Bx= 1/(4*GridSpacing)*[-(e2-1).*(e3-1),(e2-1).*(e3-1),-(e2+1).*(e3-1),(e2+1).*(e3-1),...
            (e2-1).*(e3+1),-(e2-1).*(e3+1),(e2+1).*(e3+1),-(e2+1).*(e3+1)];
        By= 1/(4*GridSpacing)*[(e1-1).*(e3-1),-(e1+1).*(e3-1),(e1+1).*(e3-1),-(e1-1).*(e3-1),...
            -(e1-1).*(e3+1),(e1+1).*(e3+1),-(e1+1).*(e3+1),(e1-1).*(e3+1)];
        Bz= 1/(4*GridSpacing)*[-(e1-1).*(e2-1),(e1+1).*(e2-1),-(e1+1).*(e2+1),(e1-1).*(e2+1),...
            (e1-1).*(e2-1),-(e1+1).*(e2-1),(e1+1).*(e2+1),-(e1-1).*(e2+1)];
        for elnum=1:size(elids,1)
            if elmask(elnum)==1
                % Displacement at nodes [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4].
                elu=u(reshape([dofids(elids(elnum,:)*3-2)';
                               dofids(elids(elnum,:)*3-1)';
                               dofids(elids(elnum,:)*3  )'],[],1));
                % Displacment derivatives at nodes (row=[uxx,uyy,uxy,uyx],col=nodes).
                eldu=[Bx*elu(1:3:end-2),By*elu(1:3:end-2),Bz*elu(1:3:end-2),...
                    Bx*elu(2:3:end-1),By*elu(2:3:end-1),Bz*elu(2:3:end-1),...
                    Bx*elu(3:3:end)  ,By*elu(3:3:end)  ,Bz*elu(3:3:end)  ]';
                % Track number of elements neighbouring each node.
                elnN(elids(elnum,:))=elnN(elids(elnum,:))+1;
                % Gridded displacement derivatives.
                UXX(elids(elnum,:))=UXX(elids(elnum,:))+eldu(1,:);
                UXY(elids(elnum,:))=UXY(elids(elnum,:))+eldu(2,:);
                UXZ(elids(elnum,:))=UXZ(elids(elnum,:))+eldu(3,:);
                UYX(elids(elnum,:))=UYX(elids(elnum,:))+eldu(4,:);
                UYY(elids(elnum,:))=UYY(elids(elnum,:))+eldu(5,:);
                UYZ(elids(elnum,:))=UYZ(elids(elnum,:))+eldu(6,:);
                UZX(elids(elnum,:))=UZX(elids(elnum,:))+eldu(7,:);
                UZY(elids(elnum,:))=UZY(elids(elnum,:))+eldu(8,:);
                UZZ(elids(elnum,:))=UZZ(elids(elnum,:))+eldu(9,:);
           end
        end
        % Average nodal data by number of elements neighbouring each node.
        UXX=UXX./elnN; UXY=UXY./elnN; UXZ=UXZ./elnN;
        UYX=UYX./elnN; UYY=UYY./elnN; UYZ=UYZ./elnN;
        UZX=UZX./elnN; UZY=UZY./elnN; UZZ=UZZ./elnN;
        % Exclude Mask and Seem nodes.
        varargout{1}.xx=UXX.*MASK; 
        varargout{1}.xy=UXY.*MASK; 
        varargout{1}.xz=UXZ.*MASK;
        varargout{1}.yx=UYX.*MASK;
        varargout{1}.yy=UYY.*MASK;
        varargout{1}.yz=UYZ.*MASK;
        varargout{1}.zx=UZX.*MASK; 
        varargout{1}.zy=UZY.*MASK;
        varargout{1}.zz=UZZ.*MASK;
end

%% Calculate strain and stress tensor
% Symmetric strain tensor.
if compstress=="On" || compstrain=="On" || nargout>=3
    EXX=UXX;
    EYY=UYY;
    EXY=(UXY+UYX)/2;
    if DataType =="PlaneStress"
        EXZ=zeros(DataSize);
        EYZ=zeros(DataSize);
        EZZ=v/(v-1)*(EXX+EYY);
    elseif DataType =="Membrane"
        error('Membrane feature not coded. On the to-do list.')
    else
        EZZ=UZZ;
        EXZ=(UXZ+UZX)/2;
        EYZ=(UYZ+UZY)/2;
    end
    % Store data.
    if compstrain=="On" || nargout>=3
        varargout{3}.xx=EXX;
        varargout{3}.yy=EYY;
        varargout{3}.zz=EZZ;
        varargout{3}.xy=EXY;
        varargout{3}.xz=EXZ;
        varargout{3}.yz=EYZ;
    end
end
% Stress tensor.
if compstress=="On" || nargout==4
    SXX=E/(1-v^2)*(EXX+v*(EYY+EZZ));
    SYY=E/(1-v^2)*(EYY+v*(EXX+EZZ));
    SZZ=E/(1-v^2)*(EZZ+v*(EXX+EYY));
    SXY=E/(1-v^2)*(1-v)*EXY;
    SXZ=E/(1-v^2)*(1-v)*EXZ;
    SYZ=E/(1-v^2)*(1-v)*EYZ;
    % Store data.
    varargout{4}.xx=SXX;
    varargout{4}.yy=SYY;
    varargout{4}.zz=SZZ;
    varargout{4}.xy=SXY;
    varargout{4}.xz=SXZ;
    varargout{4}.yz=SYZ;
end

end

function [K,u,f,dofids,elids,elmask]=FEA2D(DataSize,ux,uy,fx,fy,gridspacing,mask,seem,E,v)
% FEA2D Computes stiffness matrix (K), displacment fector (u), force vector
% (f) and degree of freedom ids (dofids) under plane stress assumption.

% Number of elements.
nelx=DataSize(2)-1;
nely=DataSize(1)-1;
nel=nelx*nely;
% Element node ids.
[elrow,elcol]=meshgrid(0:nelx-1,0:nely-1);
elid4=elrow*(nely+1)+(nely+1-(elcol));
elid3=elid4+(nely+1);
elid2=elid4+nely;
elid1=elid4-1;
elids=[elid1(:),elid2(:),elid3(:),elid4(:)];
% Dof ids.
xdof=reshape(flipud(reshape((1:prod(DataSize))*2-1,DataSize)),1,[]);
ydof=reshape(flipud(reshape((1:prod(DataSize))*2,DataSize)),1,[]);
dofids=reshape([xdof;ydof],[],1);
% Mask & Seem elements.
nodmask=find(isnan(mask(:)));
nodseem=find(isnan(seem(:)));
elmask=zeros(nel,4);
elseem=zeros(nel,4);
for elnum=1:prod(nel)
    elmask(elnum,:)=ismember(elids(elnum,:),nodmask);
    elseem(elnum,:)=ismember(elids(elnum,:),nodseem);
end
elmask=double(~any(elmask,2));
elseem=double(~elseem(:,4));
elmask(elmask==0|elseem==0)=1e-12;
% Element stiffness matrix.
A11=[12, 3,-6,-3; 3,12, 3, 0;-6, 3,12,-3;-3, 0,-3,12];
A12=[-6,-3, 0, 3;-3,-6,-3,-6; 0,-3,-6, 3; 3,-6, 3,-6];
B11=[-4, 3,-2, 9; 3,-4,-9, 4;-2,-9,-4,-3; 9, 4,-3,-4];
B12=[ 2,-3, 4,-9;-3, 2, 9,-2; 4, 9, 2, 3;-9,-2, 3, 2];
KE=1/(1-v^2)/24*([A11,A12;A12',A11]+v*[B11,B12;B12',B11])*gridspacing;
% Global stiffness matrix.
nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat=repmat(edofVec,1,8)+repmat([0,1,2*nely+[2,3,0,1],-2,-1],nelx*nely,1);
iK=reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK=reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK=reshape(KE(:)*E*elmask',64*nelx*nely,1);
K=sparse(iK,jK,sK);
K=(K+K')/2;
% Displacement & force dofs.
u(dofids)=reshape(reshape([reshape(ux,1,[]);reshape(uy,1,[])],[],1),1,[]); u=u';
f(dofids)=reshape(reshape([reshape(fx,1,[]);reshape(fy,1,[])],[],1),1,[]); f=f';
end

function [K,u,f,dofids,elids,elmask]=FEA3D(DataSize,ux,uy,uz,fx,fy,fz,gridspacing,mask,seem,E,v)
% FEA3D Computes stiffness matrix (K), displacment fector (u), force vector
% (f) and degree of freedom ids (dofids) for 3D data. Adabted from "An
% efficient 3D topology optimization code written in Matlab", Kai Liu,
% Andrés Tovar, 25 June 2014 see https://top3dapp.com/

% Number of elements.
nelx=DataSize(2)-1;
nely=DataSize(1)-1;
nelz=DataSize(3)-1;
nel=nelx*nely*nelz;
% Element node ids.
elzid=(nelx+1)*(nely+1);
[elii,eljj,elkk]=meshgrid(0:nelx-1,0:nely-1,0:nelz-1);
elid1=(elkk)*(nelx+1)*(nely+1)+(elii)*(nely+1)+(nely+1-(eljj));
elid2=elid1+(nely+1);
elid3=elid1+nely;
elid4=elid1-1;
elid5=elid1+elzid;
elid6=elid2+elzid;
elid7=elid3+elzid;
elid8=elid4+elzid;
elids=[elid1(:),elid2(:),elid3(:),elid4(:),elid5(:),elid6(:),elid7(:),elid8(:)];
% Dof ids.
xdof=reshape(flipud(reshape((1:prod(DataSize))*3-2,DataSize)),1,[]);
ydof=reshape(flipud(reshape((1:prod(DataSize))*3-1,DataSize)),1,[]);
zdof=reshape(flipud(reshape((1:prod(DataSize))*3,DataSize)),1,[]);
dofids=reshape([xdof;ydof;zdof],[],1);
% Mask & Seem elements.
nodmask=find(isnan(mask(:)));
nodseem=find(isnan(seem(:)));
elmask=zeros(nel,8);
elseem=zeros(nel,8);
for elnum=1:prod(nel)
    elmask(elnum,:)=ismember(elids(elnum,:),nodmask);
    elseem(elnum,:)=ismember(elids(elnum,1),nodseem);
end
elmask=double(~any(elmask,2));
elseem=double(~any(elseem,2));
elmask(elmask==0|elseem==0)=1e-12;
% Element stiffness matrix
KE=lk_H8(v)*gridspacing;
nodegrd=reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids=reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz=0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids=repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec=3*nodeids(:)+1;
edofMat=repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]],nel,1);
% Global stiffness matrix assembly.
iK=reshape(kron(edofMat,ones(24,1))',24*24*nel,1);
jK=reshape(kron(edofMat,ones(1,24))',24*24*nel,1);
sK=reshape(KE(:)*E*elmask(:)',24*24*nel,1);
K=sparse(iK,jK,sK);
K=(K+K')/2;
% Displacement & force dofs.
u(dofids)=reshape(reshape([reshape(ux,1,[]);reshape(uy,1,[]);reshape(uz,1,[])],[],1),1,[]); u=u';
f(dofids)=reshape(reshape([reshape(fx,1,[]);reshape(fy,1,[]);reshape(fz,1,[])],[],1),1,[]); f=f';
end

function KE = lk_H8(v)
A=[32 6 -8   6 -6 4 3 -6 -10   3 -3 -3 -4 -8;
    -48 0  0 -24 24 0 0  0  12 -12  0 12 12 12];
k=1/144*A'*[1; v];
K1=[k(1)  k(2)  k(2)  k(3)  k(5)  k(5);
    k(2)  k(1)  k(2)  k(4)  k(6)  k(7);
    k(2)  k(2)  k(1)  k(4)  k(7)  k(6);
    k(3)  k(4)  k(4)  k(1)  k(8)  k(8);
    k(5)  k(6)  k(7)  k(8)  k(1)  k(2);
    k(5)  k(7)  k(6)  k(8)  k(2)  k(1)];
K2=[k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3=[k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4=[k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5=[k(1)  k(2)  k(8)  k(3)  k(5)  k(4);
    k(2)  k(1)  k(8)  k(4)  k(6)  k(11);
    k(8)  k(8)  k(1)  k(5)  k(11) k(6);
    k(3)  k(4)  k(5)  k(1)  k(8)  k(2);
    k(5)  k(6)  k(11) k(8)  k(1)  k(8);
    k(4)  k(11) k(6)  k(2)  k(8)  k(1)];
K6=[k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE=1/((v+1)*(1-2*v))*...
    [K1  K2  K3  K4;
    K2' K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end

%function isMeshGrid
function isMG=isMeshGrid(X,Y)
% Adapted from Matlab.
if ~ismatrix(X)||isempty(X)||~isequal(size(X),size(Y))
    isMG=false;
elseif (~isnumeric(X)&&~islogical(X))||(~isnumeric(Y)&&~islogical(Y))
    isMG=false;
else
    isMG=true;
end
end