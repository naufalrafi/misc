%%NAUFAL RAFI ANTARES======================================================
%FOR PROFESSOR RYAN ORSZULIK, ESSE 4350 ASSIGNMENT 2
%FEBRUARY 15, 2020
%==========================================================================
clear all

%USER INPUT PROMPTS SECTION------------------------------------------------
mp = input('For MANUAL PROMPTED input, please enter 0\nFor PRESET input, please enter 1\n\nPlease enter 0 or 1: ');

if mp == 1
%PRESET INPUTS, ((EDIT BEFORE RUNNING THE CODE))---------------------------
    %number of nodes
    nnum = 5;
    
    %node coordinates [x1 y1;...;xn yn]
    nodes = [0 0;
        1.8 3.6;
        3.6 0;
        5.4 3.6;
        7.2 0];
    x = nodes (:,1);
    y = nodes(:,2);
    
    %number of elements
    enum = 7;
    
    %element lenghts (1..enum) [m]
    L = [3.6 3.6 3.6 3.6 3.6 3.6 3.6];
    
    %angles elements make c-clkwise from x-axis [degree]
    T = [60 0 -60 60 0 -60 0];
    
    %Cross section area of elements 1..enum [m^2)]
    A = [0.0030 0.0030 0.0030 0.0030 0.0030 0.0030 0.0030];
    
    %Presence of boundary conditions (1 = node axis fixed)
    BC = [1 1 0 0 0 0 0 0 0 1];
    
    %Young's modulus of elements 1..enum [Pa]
    E = [2.1e11 2.1e11 2.1e11 2.1e11 2.1e11 2.1e11 2.1e11 2.1e11];
    
    %elemet nodal connections
    elements = [1 2; %element 1
        1 3; %element 2
        2 3;  %element 3
        3 4; %element 4
        3 5; %element 5
        4 5; %element 6
        2 4]; %element 7
    
    %Force applied at nodes (directional x1,y1,...,xn,y)[N]
    F = [0 -280000 0 0 0 -210000 0 -80000 0 -280000];
%--------------------------------------------------------------------------

elseif mp == 0    
%MANUAL PROMPTED INPUT-----------------------------------------------------
nnum = input('How many nodes are there in the truss? ');

while mod(nnum,1) ~= 0
    fprintf('Please input an integer\n')
    nnum = input('Try again: ');
end

for i = 1:1:nnum
    j = (2*(i-1))+1;
    k = 2*i;

    promptx = ['Please insert node '  num2str(i)  ' x-coordinate: '];
    x(i) = input(promptx);
    prompty = ['Please insert node '  num2str(i)  ' y-coordinate: '];
    y(i) = input(prompty);
    nodes(i,:)= [x(i) y(i)];
    
    fp = input('Is there applied force on this node (1 for yes, 0 for no)? ');
    if fp == 0
        F(j) = 0;
        F(k) = 0;
    elseif fp == 1
        F(j) = input('What is the applied force in the global x-direction? (N) ');
        F(k) = input('What is the applied force in the global y-direction? (N) ');
    else
        fprintf('Neither 0 or 1 was entered, assuming no applied force')
        F(j) = 0;
        F(k) = 0;
    end
    
    bcp = input('Is this node fixed (1 for yes, 0 for no)? ');
    if bcp == 0
        BC(j) = 0;
        BC(k) = 0;
    elseif bcp == 1
        BC(j) = input('is it fixed in the global x-direction (1 for yes, 0 for no)? ');
        BC(k) = input('is it fixed in the global y-direction (1 for yes, 0 for no)? ');
        if BC(j)==0 && BC(k)==0
            fprintf('WARNING you have actually entered no boundary conditions\n')
        elseif xor(BC(j),BC(k)) == 1
            fprintf(['Pin-roller added to node ' num2str(i) '\n'])
        elseif BC(j)==1 && BC(k)==1
            fprintf(['Node ' num2str(i) 'pinned\n'])
        end
    else 
        fprintf('Neither 0 or 1 was entered, assuming a free moving node')
        BC(j) = 0;
        BC(k) = 0;
    end
end

enum = input('How many elements are there in the truss? ');

while mod(enum,1) ~= 0
    fprintf('Please input an integer\n')
    nnum = input('Try again: ');
end

for i = 1:1:enum
    promptx = ['Which two nodes does element '  num2str(i)  ' connect? \nFirst node number: '];
    en1(i) = input(promptx);
    prompty = ['Second node number: '];
    en2(i) = input(prompty);
    elements(i,:) = [en1(i) en2(i)];
    promptE = ['What is the Young''s modulus value for element '  num2str(i)  '? (Pa) '];
    E(i) = input(promptE);
    promptA = ['What is the cross section area of element '  num2str(i)  '? (m^2) '];
    A(i) = input(promptA);
    promptL = ['What is the length of element '  num2str(i)  '? (m) '];
    L(i) = input(promptL);
    promptT = ['What is the angle that element '  num2str(i)  ' makes counter-clockwise from the positive global x-axis? (degrees) '];
    T(i) = input(promptT);
end
else
    fprintf('[ERROR] BETTER LUCK NEXT TIME AT READING INSTRUCTIONS!\n\n')
    return
end
%--------------------------------------------------------------------------
%END OF INPUT SECTION------------------------------------------------------
%--------------------------------------------------------------------------

%ASSEMBLING THE GLOBAL STIFFNESS MATRIX------------------------------------
Ke = zeros(4,4,enum); %Setting up element stiffness matrix size
Kg = zeros(2*nnum,2*nnum); %Setting up global stiffness matrix size

for i = 1:1:enum
    Ke(:,:,i) = elstiffmat(E(i),A(i),L(i),T(i));
    
    j = 2*(elements(i,1)-1)+1; %index of u of first node of the element
    k = j + 1; %index of v of first node of the element
    l = 2*(elements(i,2)-1)+1; %index of u of second node of the element
    m = l + 1; %index of v of second node of the element
    
    Kg(j,j) = Kg(j,j) + Ke(1,1,i);
    Kg(k,k) = Kg(k,k) + Ke(2,2,i);
    Kg(l,l) = Kg(l,l) + Ke(3,3,i);
    Kg(m,m) = Kg(m,m) + Ke(4,4,i);
    Kg(j,k) = Kg(j,k) + Ke(1,2,i);
    Kg(k,j) = Kg(j,k); %symmetry
    Kg(j,l) = Kg(j,l) + Ke(1,3,i);
    Kg(l,j) = Kg(j,l); %symmetry
    Kg(j,m) = Kg(j,m) + Ke(1,4,i);
    Kg(m,j) = Kg(j,m); %symmetry
    Kg(k,l) = Kg(k,l) + Ke(2,3,i);
    Kg(l,k) = Kg(k,l); %symmetry
    Kg(k,m) = Kg(k,m) + Ke(2,4,i);
    Kg(m,k) = Kg(k,m); %symmetry
    Kg(l,m) = Kg(l,m) + Ke(3,4,i);
    Kg(m,l) = Kg(l,m); %symmetry
end
%--------------------------------------------------------------------------
%REDUCING THE EQUATIONS USING BCs------------------------------------------
BCi = find(BC); %indices fixed
ui = find(~BC); %indices free
Kgr = Kg;

for i = 1:1:length(BCi)
    Kgr(BCi(i)-(i-1),:) = [];
    Kgr(:,BCi(i)-(i-1)) = [];
end

for i = 1:1:length(ui)
    Ff(i) = F(ui(i));
end
%--------------------------------------------------------------------------
%SOLVING REDUCED EQUATIONS-------------------------------------------------
u = linsolve(Kgr,transpose(Ff));

U = zeros(2*nnum,1);
for i = 1:1:length(ui)
    U(ui(i)) = u(i);
end
%--------------------------------------------------------------------------
%SOLVING REACTION FORCES---------------------------------------------------
Kgrx = Kg;
for i = 1:1:length(BCi)
    Kgrx(:,BCi(i)-(i-1)) = [];
    Frx(i) = F(BCi(i));
end

for i = 1:1:length(ui)
    Kgrx(ui(i)-(i-1),:) = [];
end

R = Kgrx*u - Frx;
%--------------------------------------------------------------------------
%COMPUTING ELEMENTAL STRESS------------------------------------------------
for i = 1:1:enum
    RotG2E = [cosd(T(i)) sind(T(i)) 0 0;
        -sind(T(i)) cosd(T(i)) 0 0;
        0 0 cosd(T(i)) sind(T(i));
        0 0 -sind(T(i)) cosd(T(i))];
    
    j = 2*(elements(i,1)-1)+1; %index of u of first node of the element
    k = j + 1; %index of v of first node of the element
    l = 2*(elements(i,2)-1)+1; %index of u of second node of the element
    m = l + 1; %index of v of second node of the element
    
    Ue = RotG2E*[U(j);U(k);U(l);U(m)];
    s(i) = (210e9)*(Ue(3) - Ue(1))/3.6;
end
%--------------------------------------------------------------------------
%PLOTS---------------------------------------------------------------------
nodeplot = plot(nodes(:,1),nodes(:,2),'ko','LineWidth',3);
hold on; grid on; grid minor

for i = 1:1:nnum
    j = (2*(i-1))+1;
    k = 2*i;
    Up(i,1) = 100*U(j);
    Up(i,2) = 100*U(k);
end
    

exy = zeros(2,2,enum);
for i = 1:1:enum
    j = elements(i,1); %index of x and y of first node of the element
    k = elements(i,2); %index of x and y of second node of the element

    exy(:,:,i) = [x(j) y(j); x(k) y(k)];
    before = plot(exy(:,1,i),exy(:,2,i),'r','LineWidth',2); 
    
    exyd(:,:,i) = [x(j)+Up(j,1) y(j)+Up(j,2); x(k)+Up(k,1) y(k)+Up(k,2)];
    after = plot(exyd(:,1,i),exyd(:,2,i),'b','LineWidth',2);
end

%-----------------------MAKING THE PLOT PRETTY-----------------------------
plot([0 0;(max(x)+2) 0],[0 0; 0 (max(y)+2)],'k','LineWidth',0.5)

arrowx = 0.1;
arrowy = 0.05;
patch(...
    [(max(x)+2)-arrowx -arrowy; (max(x)+2)-arrowx +arrowy; (max(x)+2) 0.0], ...
    [arrowy (max(y)+2)-arrowx; -arrowy (max(y)+2)-arrowx; 0 (max(y)+2)], 'k', 'clipping', 'off')

text(max(x)+2+arrowy, 0, 'X', 'fontsize', 14)
text(-arrowy, max(y)+2.1+arrowx, 'Y', 'fontsize', 14)

legend([nodeplot before after],'Original node locations','Bar elements before deformation','Bar elements after deformation')
title('Naufal Rafi Antares - Truss Deformation Visualizer')
set(gcf,'position',[200,100,1500,850])
%--------------------------------------------------------------------------


%FUNCTION GENERATING ELEMENT STIFFNESS MATRIX------------------------------
%MATLAB VERSION R2018a REQUIRES FUNCTIONS DECLARED AT THE END OF FILE------
function K = elstiffmat(E, A, L, T)
l = cosd(T);
m = sind(T);

k = [l^2 l*m -(l^2) -l*m;
    l*m m^2 -l*m -(m^2);
    -(l^2) -l*m l^2 l*m;
    -l*m -(m^2) l*m m^2];

K = (E.*A./L)*k;    
end
%--------------------------------------------------------------------------
