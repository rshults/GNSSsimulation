function GNSS_sim1(iteration, dt)
    dt = 1;
    iteration = 1000;
    file_name = 'all_data.txt';
    receiver = 2;

    k=1;
    while k < 3
        k = menu ('GNSS Deformation Network Simulation','Import file data','Iterations','Calculate');
        fprintf('Name file data = %s\n',file_name);
        fprintf('Number of iterations = %g\n',iteration);
        fprintf('Number of receivers = %g\n',iteration);
        fprintf('Calculate\n\n');
        if k == 1
            file_name = input([sprintf('Current Value file name = %s',file_name), '\nNew Value file name = ']);
        elseif k == 2
            iteration = input([sprintf('Current Value iteration = %g',iteration), '\nNew Value iteration =']);                       
        %elseif k == 3
            %receiver = input([sprintf('Current Value receiver = %g',receiver), 'New Value receiver =']); 
        end
    end
     
    fid = fopen(file_name, 'r');
    % �������� ������������ ��������
    if fid == -1
        error('File is not opened');
        return
    end

    cell_arr = textscan(fid, '%f %f %f %s','Delimiter',';','CommentStyle',{'/*','*/'});
    receiver = cell_arr{1,1}; % ���������� ������� ����
    nbaselines = cell_arr{1,2}; % ���������� ������� �����
    n_xr_lines = cell_arr{1,3}; % ���������� ��������� � ������� �����
    xr_name = cell_arr{1,4}; % ��� �������� ������
    
    % ������ ������ ��������� ������� � �������: id ������ - ��� ������ - x - y - z - b - l;
    rec_data = cell(receiver,7);
    for i = 1:1:receiver
        rec_data(i,1) = num2cell(0);
        rec_data(i,2:7) = textscan(fid, '%s %n %n %n %n %n',1,'Delimiter',';','CommentStyle',{'/*','*/'});
        
        rec_data(i,2) = cellstr(rec_data{i,2});
    end

    %������� ����� ��� ���������
    %������: id � ������� - id � ��������� - � ������� - � ��������� - Dx - Dy - Dz
    base_lines = cell(nbaselines,7);
    for i = 1:1:nbaselines
        base_lines(i,1:2) = [num2cell(0), num2cell(0)];
        base_lines(i,3:7) = textscan(fid, '%s %s %n %n %n',1,'Delimiter',';','CommentStyle',{'/*','*/'});
        
        base_lines(i,3:4) = [cellstr(base_lines{i,3}), cellstr(base_lines{i,4})];
    end
    fclose(fid);
    
    % ��������� �������� ����� ��� ����� � ��������, ��������� ������ �� ������� � �������� �� ���������������� �������
    index_xr = find(strcmp(rec_data(:,2), xr_name));
    ind = 1;
    rec_data{index_xr,1} = ind;
    for i = 1:1:receiver 
        if rec_data{i,1} == 0
            ind = ind + 1;
            rec_data{i,1} = ind;
        end
    end
    
    % ����������� ������ ��������� ������� �� ����������� ��������, ����� ������� ����� � ��� ���������� ������
    rec_data = sortrows(rec_data);
    
    % ��������� ������������������ ������� � ���������� 
    seq_receiver = rec_data(:,2)';

    % ��������� ������� ����� ������� � ���������
    for i = 1:1:nbaselines
        base_lines{i,1} = find(strcmp(seq_receiver, base_lines(i,3)));
        base_lines{i,2} = find(strcmp(seq_receiver, base_lines(i,4)));
    end
    
    % ����������� ������ ��������� ������� ����� �� ����������� ��������,
    % ����� ��������� �� ������� ����� ���� �������
    base_lines = sortrows(base_lines);
    base = conv_cell_to_matrix([base_lines(:,5); base_lines(:,6); base_lines(:,7)]); % ������� ����� ��� ���������

    A = zeros(nbaselines*3,receiver*3); % ������� ������������� ��������� ��������
    for i = 1:1:nbaselines
        ind = find(strcmp(seq_receiver, base_lines(i,3)));
        j = ind+2*(ind-1);
        A(i,j) = -1; % �� ������ ������� x
        A(i+nbaselines,j+1) = -1; % �� ������ ������� y
        A(i+2*nbaselines,j+2) = -1; % �� ������ ������� z

        ind = find(strcmp(seq_receiver, base_lines(i,4)));
        j = ind+2*(ind-1);
        A(i,j) = 1; % �� ������ ��������� x
        A(i+nbaselines,j+1) = 1; % �� ������ ��������� y
        A(i+2*nbaselines,j+2) = 1; % �� ������ ��������� z
    end
    
    % ��������� xyz ���������� ������� ��������� �� ������� ������ rec_data
    coord = conv_cell_to_matrix(rec_data(:,3:5)); 

    coordmean = mean(coord); % ������ ������� �������� �� �������� coord
    x = coord - coordmean; % �������������� ���������� ����� ��� ���������

    coord_rad = conv_cell_to_matrix(rec_data(:,6:7)); % ������������� ���������� ������� � �������� b � l

    Cov = zeros(iteration,1); % ������ ������ �������� ����������� ������� ����� �� ���� ���������
    
    mux = zeros(1,iteration); % ������ �������� ����������� ��������� �� ������������� �� x
    muy = zeros(1,iteration); % ������ �������� ����������� ��������� �� ������������� �� y
    muz = zeros(1,iteration); % ������ �������� ����������� ��������� �� ������������� �� z
    
    DevX = zeros(iteration,receiver);
    DevY = zeros(iteration,receiver);
    DevZ = zeros(iteration,receiver);

    for t = 1:dt:iteration
        sx = zeros(nbaselines,1); sy = zeros(nbaselines,1); sz = zeros(nbaselines,1);
        lx = abs(base(1:nbaselines,1));
        ly = abs(base((nbaselines+1):nbaselines*2,1));
        lz = abs(base((nbaselines*2+1):nbaselines*3,1));

        a = [0.003; 0.003; 0.0035];
        b = 0.0000001;
        m = 0;
        for n = 1:1:nbaselines
            sx(n,1) = (a(1,1) + b*lx(n,1))^2; % ��������� ��������� ������� ����� x
            sy(n,1) = (a(2,1) + b*ly(n,1))^2; % ��������� ��������� ������� ����� y
            sz(n,1) = (a(3,1) + b*lz(n,1))^2; % ��������� ��������� ������� ����� z
        end
        s = [sx;sy;sz];
        r = m+sqrt(s).*randn(nbaselines*3,1); % r - ��������������� ������ ��������� ������� �����

        lm = base + r; % ���������� ������� ����� (��������� + ������)

        xxo = x(1,1)+lm(1:n_xr_lines,1); % ������������ ���������� x ������������ �����
        yyo = x(1,2)+lm((nbaselines+1):(nbaselines+n_xr_lines),1); % ������������ ���������� y ������������ �����
        zzo = x(1,3)+lm((nbaselines*2+1):(nbaselines*2+n_xr_lines),1); % ������������ ���������� z ������������ �����
 
        xo = zeros(receiver*3,1); % ������ ������� ������������ �������� ���������� (���������)
        xo(1:3,1) = x(1,1:3)';
        xo(4:3:receiver*3,1) = xxo;
        xo(5:3:receiver*3,1) = yyo;
        xo(6:3:receiver*3,1) = zzo;
        
        deltab = A*xo - lm; % ������ ������� ��������� ������ ������� �����

        wx = zeros(nbaselines,1); wy = zeros(nbaselines,1); wz = zeros(nbaselines,1);
        lx = abs(lm(1:nbaselines,1));
        ly = abs(lm((nbaselines+1):nbaselines*2,1));
        lz = abs(lm((nbaselines*2+1):nbaselines*3,1));
        
        for n = 1:1:nbaselines
            wx(n,1) = 0.01/(a(1,1) + b*lx(n,1)); % ���� ������� ����� x
            wy(n,1) = 0.01/(a(2,1) + b*ly(n,1)); % ���� ������� ����� y
            wz(n,1) = 0.01/(a(3,1) + b*lz(n,1)); % ���� ������� ����� z
        end

        W = [wx; wy; wz]; %������-������� ����� ���������
        P = diag(W); %������� ����� ��������� (� ��������� ����������� �� ����� x1x2x3...y1y2y3...z1z2z3...)
        
        r = zeros(nbaselines*3,1);
        %����������� ������-������� ����� ��������� �� ����� x1y1z1, x2y2z2
        r(1:3:nbaselines*3,1) = wx(1:nbaselines,1);
        r(2:3:nbaselines*3,1) = wy(1:nbaselines,1);
        r(3:3:nbaselines*3,1) = wz(1:nbaselines,1);
        R = diag(r); % ������������� ������� ����� ��������� (� ��������� �������� �� ����� x1y1z1, x2y2z2 � �.�. 

        Rnew = zeros(nbaselines*3); % ����� ��������� ������� ����� ������������� �� ����� x1y1z1, x2y2z2 � �.�.
        for k = 1:1:nbaselines
            % ��������� ������ ������� ����� �������
            ind = find(strcmp(seq_receiver, base_lines(k,3)));
           
            % ��������� ���������� (b,l) ��� ����� �������� �� �� ��������� ���������� ������� ��������
            rot = calculate_rotation(coord_rad(ind,1),coord_rad(ind,2));
            C = rot';
            
            i = 3*k-2; 
            j = k+2*k;
            Rnew(i:j,i:j) = C*R(i:j,i:j)*C';
        end

        Pnew = zeros(nbaselines*3);
        for n = 1:1:nbaselines
          i = n+2*(n-1);

          % ������ ����
          Pnew(n,n) = Rnew(i,i);
          Pnew(n,nbaselines+n) = Rnew(i+1,i);
          Pnew(n,nbaselines*2+n) = Rnew(i+2,i);

          % ������ ����
          Pnew(nbaselines+n,n) = Rnew(i,i+1);
          Pnew(nbaselines+n,nbaselines+n) = Rnew(i+1,i+1);
          Pnew(nbaselines+n,nbaselines*2+n) = Rnew(i+2,i+1);

          % ������ ����
          Pnew(nbaselines*2+n,n) = Rnew(i,i+2);
          Pnew(nbaselines*2+n,nbaselines+n) = Rnew(i+1,i+2);
          Pnew(nbaselines*2+n,nbaselines*2+n) = Rnew(i+2,i+2);
        end
        
        Pnew = (receiver/2)*Pnew;

        supp = eye(3); % ��������������� ��������� ������� ��� ���������� ������� D
        D = zeros(7,receiver*3);
        for n = 1:1:receiver
            i = n + 2*(n-1);
            D(1:3,i:i+2) = supp; 
            D(4:6,i:i+2) = format_coord_D(xo(i:i+2));
            D(7,i:i+2) = xo(i:i+2)';
        end

        N = A'*Pnew*A;
        QQ = N + D'*D;
        Q11 = QQ\N/QQ;
        hat = -Q11*A'*Pnew*deltab; %������ �������� � ������������ ��������� ����������
        xhat = xo - hat; %������ ����������
        
        dx = x(:,1) - xhat(1:3:receiver*3,1); %������ ������� ���������� ������ ��������� � �������� ���������
        dy = x(:,2) - xhat(2:3:receiver*3,1); %������ ������� ���������� ������ ��������� � �������� ���������
        dz = x(:,3) - xhat(3:3:receiver*3,1); %������ ������� ���������� ������ ��������� � �������� ���������
        
        mux(1,t) = sqrt((dx'*dx)/receiver); %������ �������� ����������� ��������� �� �������������
        muy(1,t) = sqrt((dy'*dy)/receiver); %������ �������� ����������� ��������� �� �������������
        muz(1,t) = sqrt((dz'*dz)/receiver); %������ �������� ����������� ��������� �� �������������
        
        v = A*xhat - lm; %������ ������� ��������
        
        Cov(t,1) = sqrt((v'*P*v)/60); %������ �������� ����������� ������� �����
        
        DevX(t,:) = dx(:,1)';
        DevY(t,:) = dy(:,1)';
        DevZ(t,:) = dz(:,1)';
    end
    
    vecx = zeros(receiver,1);
    vecy = zeros(receiver,1);
    vecz = zeros(receiver,1);
    for i = 1:1:receiver
        vecx(i,1) = sqrt((DevX(:,i)'*DevX(:,i))/iteration);
        fprintf('mx%d = %.4f\n',i,vecx(i,1));
    
        vecy(i,1) = sqrt((DevY(:,i)'*DevY(:,i))/iteration);
        fprintf('my%d = %.4f\n',i,vecy(i,1));
        
        vecz(i,1) = sqrt((DevZ(:,i)'*DevZ(:,i))/iteration);
        fprintf('mz%d = %.4f\n',i,vecz(i,1));
        
        fprintf('\n');
    end

    close all;
    t = 1:dt:iteration;
    
    figure;
    plot3(coord(:,1),coord(:,2),coord(:,3),'g^'), title('Network points')

    figure;
    quiver(coord(:,1),coord(:,2),vecx,vecy, 'r')
    %pause

    figure;
    plot(t,Cov);
    grid;
    xlabel('Attempt');
    ylabel('Position Error (meter)');
    title('Figure 1 - Position Error (meter)')

    figure;
    plot(t,mux);
    grid;
    xlabel('Attempt');
    ylabel('Position Error X (meter)');
    title('Figure 2 - Position Error X (meter)')

    figure;
    plot(t,muy);
    grid;
    xlabel('Attempt');
    ylabel('Position Error Y (meter)');
    title('Figure 3 - Position Error Y (meter)')

    figure;
    plot(t,muz);
    grid;
    xlabel('Attempt');
    ylabel('Position Error Z (meter)');
    title('Figure 4 - Position Error Z (meter)')

    figure;
    plot(t,DevX(:,1));
    grid;
    xlabel('Attempt');
    ylabel('Deviation X1 (meter)');
    title('Figure 5 - Deviation (meter)')

    figure;
    plot(t,DevY(:,1));
    grid;
    xlabel('Attempt');
    ylabel('Deviation Y1 (meter)');
    title('Figure 6 - Deviation (meter)')

    figure;
    plot(t,DevZ(:,1));
    grid;
    xlabel('Attempt');
    ylabel('Deviation Z1 (meter)');
    title('Figure 7 - Deviation (meter)')
end

% ������������ ������� �������� (����������� 3�3) �� ��������� �������� ��������� b � l
function rotation = calculate_rotation(b, l)
    rotation = [-sin(l)          cos(l)         0;
             -sin(b)*cos(l)     -sin(b)*sin(l)  cos(b);
              cos(b)*cos(l)     cos(b)*sin(l)   sin(b)];
end

% ����������� coord ([x y z]) � ������� 3�3 ���������� ���� ��� ���������� ������� D
function format = format_coord_D(coord)
    format = [ 0,             -coord(3,1),    coord(2,1);
                coord(3,1),     0,              -coord(1,1);
              	-coord(2,1),    coord(1,1),     0           ];
end

% ����������� ������ ����� � �������
function matrix = conv_cell_to_matrix(cells)
  	n = size(cells, 1);
    m = size(cells, 2);
    matrix = zeros(n,m);
    
    for i = 1:1:n
        for j = 1:1:m
            matrix(i,j) = cells{i,j};
        end
    end
end