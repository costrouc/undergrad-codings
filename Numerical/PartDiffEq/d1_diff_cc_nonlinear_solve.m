%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is the approximation of 1-d convection diffusion problems
%   -d/dx( a(x) du/dx)  + b(x) u du/dx + c(x) u = f(x) , xleft < x < xright
%            u(xleft) = uleft,   u(xright) = uright
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global xpts nnds
global Global_r  Global_s  Global_u
global rad_bas_type  str_bas_type  vel_bas_type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Define the problem and the solution method


%
rad_bas_type = 'd1_CtsQuad' ;  %% used for temperature

xleft = 0.0 ; %% left end of interval
xright = 2.0 ; %% right end of the interval

uleft = 0.0 ; %% temperature at the left end
uright = 8.0 ; %% temperature at the right end

epsilon = 10^-6 ; %% maximum error in calculation allowed
maxiterations = 22;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Define the partition of the domain (interval)

% We assume that the subroutine will pass back xpts a 
% vector of dimension nnds x 1 containing the nodal values.
%
nnds = 20 ;
xpts = linspace(xleft, xright, nnds)' ;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Set up and initialize the Global Solution vectors.
%

if strcmp(rad_bas_type, 'd1_CtsLin') == 1 
   Global_r = zeros(nnds, 1) ;
   prev_sol = zeros(nnds, 1) + 1 ;
   % Set up initial guess for the temperature.
   % prev_sol is fudged to enter loop but records the prev sol
     
elseif  strcmp(rad_bas_type, 'd1_CtsQuad') == 1 
   Global_r = zeros(2*nnds - 1, 1) ;
   prev_sol = zeros(2*nnds - 1, 1) + 1; 
   % Set up initial guess for the temperature.
   % prev_sol is fudged to enter loop but records the prev sol    
elseif  strcmp(rad_bas_type, 'd1_CtsCub') == 1 
   Global_r = zeros(3*nnds - 2, 1) ;
   prev_sol = zeros(3*nnds - 2, 1) + 1;
   % Set up initial guess for the radial function.
   % prev_sol is fudged to enter loop but records the prev sol

end

   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Now to solve the problem !  

% Set up the approximating linear system.
% Since it is non-linear we itteratre until errror epsilon is below 
% tolerance
%% Guess an allocation for the coefficient matrix and RHS
rdim = size(Global_r,1) ;
Acoeff = spalloc(rdim, rdim, 5*rdim) ;
RHSvec = zeros(rdim,1) ;

iter = 0;
while norm(Global_r - prev_sol) > epsilon && iter < maxiterations

for isub = 1 : nnds-1 % loop over the subintervals
      
      Alocal = [ ] ;
      RHSloc = [ ] ;
      % Identify the global unknown coefficients
      if strcmp(rad_bas_type, 'd1_CtsLin') == 1
         GlTrg_r = [isub ; isub + 1] ; 
         
      elseif strcmp(rad_bas_type, 'd1_CtsQuad') == 1
         GlTrg_r = [2*isub - 1; 2*isub; 2*isub + 1] ;
          
      elseif strcmp(rad_bas_type, 'd1_CtsCub') == 1
         GlTrg_r = [3*isub - 2; 3*isub - 1; 3*isub; 3*isub + 1] ;
           
      end
      
      %%%%%%%%%%%%%%%%
      % Evaluate the integrals.
      %% First the local LHS matrix.
      
      %% The diffusion term
      Alocal =  ...
        d1_ip_Dten0_Dten0(isub, 'd1_quad_5','Afun', rad_bas_type, rad_bas_type) ;
    
      %% The convective term
      Alocal =  Alocal + ...
        d1_ip_ten0_ten0_Dten0(isub, 'd1_quad_5','Bfun', rad_bas_type, rad_bas_type)  ;
   
      %% The source term
      Alocal =  Alocal + ...
        d1_ip_ten0_ten0(isub, 'd1_quad_5','Cfun', rad_bas_type, rad_bas_type)  ;
       
      %% We compute the right hand side "vector" by two calls.
         

      RHSloc =  d1_ip_ten0(isub, 'd1_quad_5','RHSfun', rad_bas_type)  ;
      
      %% Distribute the matrix and the RHS vector to the system.
         
      Acoeff(GlTrg_r , GlTrg_r) = Acoeff(GlTrg_r , GlTrg_r) + Alocal ;
      RHSvec(GlTrg_r) = RHSvec(GlTrg_r) + RHSloc ;
      
end
   
%% Impose any boundary conditions.
%% For r(xleft) = uleft and r(right) = uright we do the following
ndim = size(Acoeff,1) ;
RHSvec(1) = uleft ;
Acoeff(1,:) = 0 ;
Acoeff(1,1) = 1 ;
RHSvec(ndim) = uright ;
Acoeff(ndim,:) = 0 ;
Acoeff(ndim,ndim) = 1 ;

%% Solve the linear system
%% Scale the system   
[Acoeff , RHSvec] = ScaleM(Acoeff , RHSvec) ;
   
prev_sol = Global_r;
Global_r = Acoeff \ RHSvec ;
   
iter = iter + 1;
end

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Postprocess the approximation.
disp(sprintf('Resulting Error: %f\n',norm(Global_r - prev_sol)));
disp(sprintf('Number of Iterations: %d\n',iter));

% Identify number of points in the basis
if strcmp(rad_bas_type, 'd1_CtsLin') == 1
   basis_pts = linspace(xleft, xright, nnds); 
elseif strcmp(rad_bas_type, 'd1_CtsQuad') == 1
   basis_pts = linspace(xleft, xright, 2 * nnds - 1);  ;       
elseif strcmp(rad_bas_type, 'd1_CtsCub') == 1
   basis_pts=linspace(xleft, xright, 3 * nnds - 2) ;    
end

plot(basis_pts,Global_r)

% Call a subroutine to graph the various quantities.


      
