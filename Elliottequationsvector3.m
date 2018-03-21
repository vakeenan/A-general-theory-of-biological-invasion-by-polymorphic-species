% MATLAB code for solving Elliot equations in 1-D for any number of morphs 
% using a Crank-Nicholson scheme. Composed by Vincent Keenan. All of the 
% calculated values are stored in an array "a". The numerical allocation of 
% the array is 1=n, 2=R, 3=D, 4=C, 5=MU, 6=r, 7=CNmatrix1, 8=CNmatrix2, 
% 9=kinetics, 10=x, 11=t, 12=x and t values for approximating the speed.

%% Initialise mesh
%L is the length of space which acts as the environment. T is the time
%over which the simulation models.
tic                   %Begins timer to determine run time.
L=28000;
Lmax=L;
T=7200;
a=cell(13,1);
%temporal and spatial step sizes.
dx=0.05;
dt=0.05;
%creating vectors for time and space to be used within the calculations.
%Other parameters are included for ease later in the calculation.
a{10}=(0:dx:L);      %Intialise space discretisation.
a{11}=(0:dt:T);      %Intitialise time discretisation.

%% Initialise parameters 
%parameters of the equation which can be selected by the user. Ri are
%reproductive parameters. Di are dispersal parameters. MUi are mutation
%parameters. Theta is set at 1/2 for use within the Crank-Nicholson scheme 
%for solving the PDE. 
theta=1/2;

D=[8.5, 8, 4.5, 3, 0.5];         %Vector containing diffusion constant.
R=[0.05, 0.3, 0.45, 0.8, 0.85];  %Vector containing growth rates.
a{2}=R;
a{3}=D;
MU=0.005*ones(length(a{2}));     %Matrix containing mutation rates.
C=ones(length(a{2}));            %Matrix containing competition coeff.
a{4}=C;
a{5}=MU;

%Vector that contains r(i) values that are important for stability of
%finitite difference methods - analagous to the CFD condition. For the
%Crank-Nicholson method this is not an issue as it is always stable for all
%values of r(i). Not to be confused with per-capita growth rate. See
%literature for numerical methods if unsure.
a{6}=zeros(length(a{2}),1);
for ii=1:length(a{2})
    a{6}(ii)=(a{3}(ii)*dt)/(dx^2);
end



%% Initialise Initial conditions
%Creating the initial vectors for each morph and concatanating them for
%ease within the calculation.
a{1}=zeros(length(a{10}),length(a{2})); %Solution vector.I'C's overwritten.

%Loop that fills the matrix with the initial conditions for the specified 
%spatial location and pop. density values.

    for jj=1:(L/dx/100)          %spatial entries to be filled.
    a{1}(jj,:)=1/length(a{2});   %pop density values
    end


%% Main solution loop for solving the equation.

a=nmorphsol3(a,dt,theta,dx,Lmax);
    
%% Data Extraction

%For calculating the speed of the travelling wave using "v=d/t".
space=diff(a{12}(1:end-1,1));
time=diff(a{12}(1:end-1,2));
speed=space./time;

%Here we store all the appropriate information in b{} required for
%reporting in the manuscript. Alternatively for all values the user could
%choose to save a{} instead.
a{13}=mean(speed);
b=cell(5,1);
b{1}=a{1};
b{2}=a{12};
b{3}=a{13};
toc;
soltime=toc;
b{4}=soltime;
b{5}=a{10};
save(['Enter your directory here/final_at_x _',num2str(a{12}(end,1)),...
'_t_',num2str(a{11}(end)),'.mat'], 'b')


