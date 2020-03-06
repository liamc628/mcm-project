% Simulating the 2-D Diffusion equation by the Finite Difference
...Method 
% Numerical scheme used is a first order upwind in time and a second
...order central difference in space (Implicit and Explicit)

%%
%Specifying parameters
Lx=2;
Ly=2;
nx=50;                           %Number of steps in space(x)
ny=50;                           %Number of steps in space(y)       
nt=0;                           %Number of time steps 
dt=0.01;                         %Width of each time step
dx=Lx/(nx-1);                     %Width of space step(x)
dy=Lx/(ny-1);                     %Width of space step(y)
x=0:dx:Lx;                        %Range of x(0,2) and specifying the grid points
y=0:dy:Lx;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
vis=0.01;                         %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)

Pmax=2;
sigmax=0.4;
sigmay=0.4;

%%
%Initial Conditions
for i=1:nx
    for j=1:ny
        %{
        if ((0.75<=y(j))&&(y(j)<=1.25)&&(0<=x(i))&&(x(i)<=0.5))
            u(i,j)=4;
        else
            u(i,j)=0;
        end
        %}
        u(i,j)=(Pmax^2/(sigmax*sigmay*2*pi))*exp(-(x(i)-Lx/2)^2/(2*sigmax^2)-(y(j))^2/(2*sigmay^2));
    end
end

%%
%B.C vector
bc=zeros(nx-2,ny-2);
%bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
%bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
bc(:,1)=-UnS/dy; bc(:,nx-2)=UnN/dy;  %Neumann B.Cs
%B.Cs at the corners:
bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;

%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
%Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
%Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;

%%
%Calculating the field variable for each time step
i=2:nx-1;
j=2:ny-1;
for it=0:nt
    
    
    %UnW=-0.0000013*it;
    %UnE=-UnW;
    
    UnN=0.00013*abs((power(it,0.5)-10));
    UnS=UnN;
    
    bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
    bc(:,1)=-UnS/dy; bc(:,nx-2)=UnN/dy;  %Neumann B.Cs
%}
    un=u;
    h=surf(x,y,u','EdgeColor','none');       %plotting the field variable
    shading interp
    axis ([0 2 0 2 0 5])
    title({['2-D Diffusion with '],['time (\itt) = ',num2str((it*dt))]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('fish population density (\rho) \rightarrow')
    drawnow; 
    refreshdata(h)
    %Uncomment as necessary
    
    %Implicit method:
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
    %Boundary conditions
    %Dirichlet:
    %{
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %}
    %Neumann:
    
    u(1,:)=u(2,:)-UnW*dx;
    u(nx,:)=u(nx-1,:)+UnE*dx;
    u(:,1)=u(:,2)-UnS*dy;
    u(:,ny)=u(:,ny-1)+UnN*dy;
    
    %Explicit method:
    %{
    u(i,j)=un(i,j)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(dy*dy));
    %Boundary conditions
    %Dirichlet:
    %{
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %}
    %Neumann:
    u(1,:)=u(2,:)-UnW*dx;
    u(nx,:)=u(nx-1,:)+UnE*dx;
    u(:,1)=u(:,2)-UnS*dy;
    u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
end