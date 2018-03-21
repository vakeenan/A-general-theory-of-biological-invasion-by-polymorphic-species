%This is a function for calculating the kinetics to be inserted into the
%main calculation loop. I take each individual term within the kinetics and
%calculate these seperately, and concatenate them at the end. This is the 
%function for nearest-neighbour mutation.

function B=kinetics(a,dt)

%Initialise the arrays to be used within the loop.
b=cell(3,1);
d=cell(length(a{2}));
%This cell is for storing the mutation matrix with 0 on the diagonals for
%calculating the final terms within the kinetics easier.
mu=(a{5}(1,2))/2;

%Creates empty vectors for use in the kinetic terms.
for i=1:3
    b{i}=zeros(length(a{2}),1);
end

%This loop calculates the kinetics for each individual species.
for ii=1:length(a{2})
    % This term is the competition coefficient part.
    b{1}=-1*sum(bsxfun(@times,a{4}(ii,:),a{1}),2);
    %These terms are the incoming and leaving mutations for each species.
    if ii==1                          %case for 1 due to indexing
        b{2}=a{1}(:,2)*mu;
        b{3}=-1*a{1}(:,1)*mu;
    elseif ii==length(a{2})           %case for N due to indexing
       b{2}=a{1}(:,end-1)*mu;
       b{3}=-1*a{1}(:,end)*mu;
    else
    b{2}=(a{1}(:,ii-1)+a{1}(:,ii+1))*mu;
    b{3}=-1*a{1}(:,ii)*2*mu;
    end
    %This term is the full kinetic equation using the above.
    d{ii}=(a{2}(ii).*a{1}(:,ii).*(1+b{1})+b{2}+b{3})*dt;
end
%Concatenating the matrix into the form used for solving numerically.
B=cat(1,d{:});

end



