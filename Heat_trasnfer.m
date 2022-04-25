%%
times = time_new;
total_time = max(times);

thickness = 7.5; %mm
steps = 20;

dx = thickness/steps;
dt = total_time/steps;
nx = (thickness/dx)+1;
nt = steps+1;
T = [802.547307533049,965.661881333826,1202.44641386309,1501.21104397401,1752.69850759634,1859.38619473850,1957.95903713806,1859.71122077243,1716.06689614075,1425.14947190476,1249.33751792153,893.458623715963,612.250131554634,415.560082030761,292.681768973805,225.141898014106];

%%
thermal_diffusivity = 0.013; % slide 4 lecture 19 for TiC

lambda = thermal_diffusivity.*dt./(dx.^2);

b = lambda;
c = b;
a = 2*(1+lambda);

%% https://matlabgeeks.weebly.com/uploads/8/0/4/8/8048228/crank-example_with_matlab_code-v3__doc_.pdf

%Initial Conditions
Uo(1)=T_traj(1); Uo(2:nx-1)=T_traj(1); Uo(nx)=T_traj(1);

%boundary conditions
Un(1)= 
Un(nx)= 50;

% Store results for future use
UUU(1,:)=Uo;
% Loop over time
for k=2:nt

 for ii=1:nx-2
 if ii==1
 d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2)+c*Un(1);
 elseif ii==nx-2
 d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2)+b*Un(nx);
 else
 d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2);
 end
 end % d is a row vector
 % Transform a, b, c scalar constants in column vectors:
 bb=b*ones(nx-3,1);
 cc=bb;
 aa=a*ones(nx-2,1);
 % Use column vectors to construct diagonal matrices
 AA=diag(aa)+ diag(-bb,1)+ diag(-cc,-1); %AA is one triadiagonal Matrix
 % Find the solution for interior nodes i=2,3,4,5
 UU=AA\d'; % UU is temp at interior nodes only
 % Build the whole solution by including BCs
 Un=[Un(1),UU',Un(nx)]; % row vector
 % Store results for future use
 UUU(k,:)=Un;
 % to start over
 Uo=Un;
end
UUU % Output
toc


