%Function for solving the solution of the reaction diffusion equations. The
%function uses the Crank-Nicholson scheme. The kinetic terms should be 
%adjusted for the users need. The input is a cellular array, the input
%for the a{1} must be an nxn matrix.
function X=nmorphsol3(a,dt,theta,dx,Lmax)

%allocates the Crank-Nicolson matrices used for the calculation. 
a{7}=CNmatrix1(a,theta);
a{8}=CNmatrix2(a,theta);
%Numerical solver which obtains the solution.
[L,U,P] = ilu(a{8},struct('type','ilutp','droptol',1e-6));
%Generate cells and vectors for use in the function.
b=cell(3,1);

k=1;        %Creates a count for saving the solution at time intervals.
a{12}=zeros(101,2);

%hvec is used later for a solution shift in space.
hvec=length(a{10}(floor(length(a{10})/2):end));

%% Main solution loop

for ii=1:length(a{11})
    %Compresses the kinetics calculations into an array
    a{9}=kinetics(a,dt);
    %First step of the Crank-Nicolson scheme
    x=a{7}*a{1}(:)+a{9};
    %Second step of the Crank-Nicolson scheme. Using an LU decomposition
    %for faster calculation.
    q=U\(L\(P*x));
    
    %This part of the code has been included in an attempt to remove the
    %problem of small numbers within the solution propagting changes in
    %decay rate through the solution due to round-off error.
    q(q < realmin)=0;
    %Reshapes the solution into a matrix corresponding to each strain.
    a{1}=reshape(q,length(a{1}),length(a{2}));
    
    %Finds the midpoint of the wave.
    xrow=xwavefront(a);
    %IF the midpoint is 9/10*L then the solution is shifted by a half. This
    %keeps the boundary the same size, but shifts space for a continuing
    %solution. This is a step to reduce computation time.
%     if xrow(end) >= 9*(length(a{10})/10)
%        %This part shifts the solution.
%        nzeros=zeros(size(a{1}));
%        nzeros(1:length(hvec:end),:)=a{1}(hvec:end,:); %%%%%%HEre is where you stopped
%        a{1}=nzeros;
%        %This part shifts the space.
%        a{10}=a{10}(hvec):dx:(a{10}(hvec)+hvec*dx);
%     end
    
    %If loop for saving equally spaced time points throughout the
    %simulation. 
    if rem(ii,floor(length(a{11})/100))==0
       %Here we store the wavefront at it's time and position within the
       %simulations in a{12}.
       xrow=xwavefront(a);
       a{12}(k,1)=a{10}(xrow);
       a{12}(k,2)=a{11}(ii);
       %Seperate array for saving the essential information required for
       %assessing speeds, i.e. the numerical solution, b{1}, at the time
       %and position, b{2}.
       b{1}=a{1};
       b{2}=a{12};
       b{3}=a{10};
       %Solution and plot from the solution is saved.
       save(['//home/vakeenan/paperfigs2/5morph/sim',num2str(k),'at_x_',num2str(a{10}(xrow)),...
       '_t_',num2str(a{11}(ii)),'.mat'], 'b')

       k=k+1;
    end
    
%IF the midpoint is 9/10*L then the solution is shifted by a half. This
%keeps the boundary the same size, but shifts space for a continuing
%solution. This is a step to reduce computation time.

%     if xrow(end) >= 9*(length(a{10})/10)
%        %This part shifts the solution.
%        nzeros=zeros(size(a{1}));
%        nzeros(1:length(hvec:end),:)=a{1}(hvec:end,:); %%%%%%HEre is where you stopped
%        a{1}=nzeros;
%        %This part shifts the space.
%        a{10}=a{10}(hvec):dx:(a{10}(hvec)+Lmax);
%     end

% There are two shifts, both commented out. This is because the user may
% want the shift done before or after saving and is up to them.
end

X=a;


end
    