%% Calibration
clear
clc

% Loading the calibration image
I = imread("IMG_2506.jpg");
imshow(I)

% Inuput the 2D and 3D (M in mm) points
% m = ginput(11)';
% 
% M = 20 * [1, 0, 2
%       0, 1, 1
%       8, 0, 1
%       7, 0, 8
%       0, 8, 2
%       0, 7, 8
%       0, 1, 7
%       0, 4, 4
%       3, 0, 6
%       5, 0, 3
%       0, 6, 5]';

% Save matrices

% writematrix(m, 'm_s.txt');
% writematrix(M, 'm_b.txt');

m = readmatrix('m_s.txt');
M = readmatrix('m_b.txt');

P = dlt(m, M);

% Checking if everything is ok
size(I)/2;
[K, ~, ~] = ud_krt(P);

% Recontruction of the 2D points
m1 = htx(P,M);

% Plotting the reconstructed points
figure
imshow(I)
hold on
    Checking(m, m1)

%% Epipolar Lines

S = imread("bottle_1.jpg");
% figure
imshow(S)

% m_S1 = ginput(15)';
% writematrix(m_S1, "m_s1.txt")

S2 = imread("bottle_2.jpg");
% figure
imshow(S2)

% m_S2 = ginput(15)';
% writematrix(m_S2, "m_s2.txt")

m_S1 = readmatrix("m_s1.txt");
m_S2 = readmatrix("m_s2.txt");

F = eight_pts(m_S2, m_S1)

% Homogeneous Coordinates
m_S1_o = [m_S1; ones(1, size(m_S1, 2))];
m_S2_o = [m_S2; ones(1, size(m_S2, 2))];

for i=1:size(m_S1_o,2)
    m_S2_o(:,i)'*F*m_S1_o(:,i)
end


% Epipolar Lines and plotting them
line_1 = F * m_S1_o;
line_2 = F' * m_S2_o;

figure
imshow(S2)
hold on
    Epipolar_Line(line_1, S2);
    plot(m_S2(1,:), m_S2(2,:), '+g','MarkerSize',15);  

figure
imshow(S)
hold on
    Epipolar_Line(line_2, S);
    plot(m_S1(1,:), m_S1(2,:), '+g','MarkerSize',15);

%% Factorisation of E and Triangulation

E = essential_lin(m_S2, m_S1, K, K);

[R_out,t_out] = relative_lin(m_S2, m_S1, K, K);
[R_out,t_out] = relative_nonlin(R_out, t_out, m_S2, m_S1, K, K);

P_m{1} = K * [eye(3), zeros(3,1)];
P_m{2} = K * [R_out, t_out];

for i = 1:size(m_S1, 2)
    M_S(:, i) = ud_triang(P_m, {m_S1(:, i), m_S2(:, i)});
end

% Projection on the images

m_S1_comp = ud_htx(P_m{1}, M_S);
m_S2_comp = ud_htx(P_m{2}, M_S);

% tiledlayout(2,4)

figure
imshow(S)
size(S)
hold on
    Checking(m_S1, m_S1_comp)

figure
imshow(S2)
hold on
    Checking(m_S2, m_S2_comp)

%% Aux Functions 

function Checking(m, m1)
     plot(m(1,:), m(2,:), '+g','MarkerSize',15);
     plot(m1(1,:), m1(2,:),   'or','MarkerSize',15);

 for ii = 1:size(m, 2)
                text(m(1,ii)+4, m(2,ii)+4,num2str(ii),'Color','r', 'FontSize', 25)
 end

end

function Epipolar_Line(epline, S)
    
    for i = 1:size(epline, 2)
    
        A = epline(1,i);
        B = epline(2,i);
        C = epline(3,i);
    
        X1 = 0;
        X2 = size(S, 1);
    
        x = linspace(X1, X2);
    
        y = -(A.* x./B) - (C/B);
        line(x, y)

    end
end