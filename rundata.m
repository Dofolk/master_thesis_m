%run all data
%fft version and without empty checking
%less commment

clear
close all
% record output messages
diary on

% parameters
dirnum = '主機2020.01';
start_date = '2020_1_1';
end_date   = '2020_1_1';
Sample_Rate = 400;
Sample_Period = 300;
Amp_max = 0.01;
r_normal = 0.001;
f0 = (1/Sample_Period:1/Sample_Period:(Sample_Rate/2))';

% program start
dir1 = [dirnum '/'];

parentfolder = [pwd filesep dir1];
if(~exist([parentfolder 'Problems'],'dir'))
    mkdir(parentfolder,'Problems');
end


for dd = datenum(start_date):datenum(end_date)
    
    close all
    pause(3)
    
    ds = datestr(dd,29);
    dL = find(ds == '-');
    date_string = [ds(1:(dL(1)-1)) '_' num2str(str2num(ds((dL(1)+1):(dL(2)-1)))) '_' num2str(str2num(ds((dL(2)+1):end)))];

    % list all files with channel X
    
    Fi = dir([dir1 sprintf('*Ch_X*%s_*.txt',date_string)]);
    FX2=[]; FY2=[]; FZ2=[];
    XA = []; YA = []; ZA = []; tA = []; 
    Q3X=[]; Q3Y=[]; Q3Z = [];
    Q1X=[]; Q1Y=[]; Q1Z = [];
    NA =[]; FA = {}; tt = [];
    EX = []; EY = []; EZ = [];
    Amp_meanX = []; Amp_varX = [];
    Amp_meanY = []; Amp_varY = [];
    Amp_meanZ = []; Amp_varZ = [];
    N_original = [];
    if ~isempty(Fi)
        wbar = waitbar(0, 'wbar', 'name', 'Read files');
        for k=1:length(Fi)
            Xfile = Fi(k).name;
            waitbar(k/length(Fi),wbar,sprintf('%s: %d / %d', datestr(date_string,29), k, length(Fi)),'title','Reading files');
            [num2str(k) ':' Xfile];
            n = strfind(Xfile, 'Ch_X');
            Yfile = [Xfile(1:n+2) 'Y' Xfile((n+4):end)];
            Zfile = [Xfile(1:n+2) 'Z' Xfile((n+4):end)];

            B = textsplit(Xfile,'_');
            BS = textsplit(B{18},'.');
            if B{15}(1) == '下'
                c = 12;
            else
                c = 0;
            end
            t1 = datenum([B{12} '-' B{13} '-' B{14} ' ' num2str(mod(str2num(B{16}),12)+c) ':' B{17} ':' BS{1}]);
            [t_x,X] = read_txt([dir1 Xfile]); t_x = double(t_x);
            
            if exist([dir1 Yfile],'file') > 0
                [t_y,Y] = read_txt([dir1 Yfile]); t_y = double(t_y);
            else
                t_y = t_x; Y = zeros(size(X));
            end
            if exist([dir1 Zfile],'file') > 0
                [t_z,Z] = read_txt([dir1 Zfile]); t_z = double(t_z);
            else
                t_z = t_x; Z = zeros(size(X));
            end
    
    % index error catch===================================================
            if ((min(t_x) > Sample_Period*Sample_Rate) && (max(t_x) < 1e+9) )
                tm = min(t_x);
                t_x = t_x - tm;
                fprintf('%s has time index problem. \n',Xfile);
            end
            if ((min(t_y) > Sample_Period*Sample_Rate) && (max(t_y) < 1e+9) )
                tm = min(t_y);
                t_y = t_y - tm;
                fprintf('%s has time index problem. \n',Yfile);
            end
            if ((min(t_z) > Sample_Period*Sample_Rate) && (max(t_z) < 1e+9) )
                tm = min(t_z);
                t_z = t_z - tm;
                fprintf('%s has time index problem. \n',Zfile);
            end
    %=====================================================================
            t_x(t_x>Sample_Period*Sample_Rate+10)=[];
            t_y(t_y>Sample_Period*Sample_Rate+10)=[];
            t_z(t_z>Sample_Period*Sample_Rate+10)=[];
            tmin = min([length(t_x),length(t_y),length(t_z)]);
            t_x = t_x(1:tmin);t_y = t_y(1:tmin);t_z = t_z(1:tmin);
            X = X(1:length(t_x));Y = Y(1:length(t_y));Z = Z(1:length(t_z));
            
            N_original = [N_original;tmin];
            Amp_meanX = [Amp_meanX; mean(X)]; Amp_varX = [Amp_varX; var(X)];
            Amp_meanY = [Amp_meanY; mean(Y)]; Amp_varY = [Amp_varY; var(Y)];
            Amp_meanZ = [Amp_meanZ; mean(Z)]; Amp_varZ = [Amp_varZ; var(Z)];
            
            tA = [tA; double(t_x)/86400/Sample_Rate+t1];
            XA = [XA; X]; Q3X = [Q3X; quantile(X,1-r_normal)]; Q1X = [Q1X; quantile(X,r_normal)];
            YA = [YA; Y]; Q3Y = [Q3Y; quantile(Y,1-r_normal)]; Q1Y = [Q1Y; quantile(Y,r_normal)];
            ZA = [ZA; Z]; Q3Z = [Q3Z; quantile(Z,1-r_normal)]; Q1Z = [Q1Z; quantile(Z,r_normal)];
            
            N = min([length(t_x),length(t_y),length(t_z)]); NA = [NA; N];
            f = 1:1:(Sample_Rate/2);
            
            FX = zeros(size(X));
            FX = nufft(X,(t_x./Sample_Rate),f) / (length(t_x)/2);
            FA(k).X = abs(FX);
            FX2 = [FX2, abs(FX)];
            EX = [EX;sqrt(sum(abs(FX).^2))];
            
            FY = zeros(size(Y));
            FY = nufft(Y,(t_y./Sample_Rate),f) / (length(t_y)/2);
            FA(k).Y = abs(FY);
            FY2 = [FY2, abs(FY)];
            EY = [EY;sqrt(sum(abs(FY).^2))];
            
            FZ = zeros(size(Z));
            FZ = nufft(Z,(t_z./Sample_Rate),f) / (length(t_z)/2);
            FZ2 = [FZ2, abs(FZ)];
            EZ = [EZ;sqrt(sum(abs(FZ).^2))];
            
            tt(k) = t1;
        end
        close(wbar);
        tt(tt == 0) = [];
        [tt,id] = sort(tt);
        NA  = NA(id); FA = FA(id);
        FX2 = FX2(:,id); Q3X = Q3X(id); Q1X = Q1X(id);
        FY2 = FY2(:,id); Q3Y = Q3Y(id); Q1Y = Q1Y(id);
        FZ2 = FZ2(:,id); Q3Z = Q3Z(id); Q1Z = Q1Z(id);
        EX = EX(id); EY = EY(id); EZ = EZ(id);
        [tA,id] = sort(tA);
        XA = XA(id); YA = YA(id); ZA = ZA(id);
    else
        fprintf('No file in %s with date %s\n',dirnum, date_string);
    end
    
    % figure 1-3, time-frequency for X, Y, Z
    if ~isempty(tt)
        %{
        figure(1);
        imagesc(tt,f,FX2);
        %imagesc(tt,f0,FX2);
        datetick('x','HH:MM');
        caxis([0,Amp_max]);
        colorbar;
        set(gca,'YDir','normal');
        xlabel('time');
        ylabel('frequency');
        title(sprintf('X-axis Date: %s',datestr(tt(1),29)));
        print([dirnum '-' date_string '_X.jpg'],'-djpeg')
        
        figure(2);
        imagesc(tt,f,FY2);
        %imagesc(tt,f0,FY2);
        datetick('x','HH:MM');
        caxis([0,Amp_max]);
        colorbar;
        set(gca,'YDir','normal');
        xlabel('time');
        ylabel('frequency');
        title(sprintf('Y-axis Date: %s',datestr(tt(1),29)));
        print([dirnum '-' date_string '_Y.jpg'],'-djpeg')
        
        figure(3);
        imagesc(tt,f,FZ2);
        %imagesc(tt,f0,FZ2);
        datetick('x','HH:MM');
        caxis([0,Amp_max]);
        colorbar;
        set(gca,'YDir','normal');
        xlabel('time');
        ylabel('frequency');
        title(sprintf('Z-axis Date: %s',datestr(tt(1),29)));
        print([dirnum '-' date_string '_Z.jpg'],'-djpeg')
        
        % Figure 4: time-amplitude for X, Y, Z
        figure(4);
        subplot(3,1,1);
        plot(tA,XA,'b.');
        xlabel('time');
        ylabel('X');
        datetick('x','HH:MM');
        ax = axis;
        Q3 = quantile(XA,1-r_normal);
        Q1 = quantile(XA,r_normal);
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), Q1-1.5*(Q3-Q1), Q3+1.5*(Q3-Q1)]);
        title(sprintf('Date: %s',datestr(tt(1),29)));
        
        subplot(3,1,2);
        plot(tA,YA,'g.');
        xlabel('time');
        ylabel('Y');
        datetick('x','HH:MM');
        ax = axis;
        mY = mean(YA);
        Q3 = quantile(YA,1-r_normal);
        Q1 = quantile(YA,r_normal);
        if Q1==Q3
            Q3 = 0.5;
            Q1 = -0.5;
        end
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), Q1-1.5*(Q3-Q1), Q3+1.5*(Q3-Q1)]);
        
        subplot(3,1,3);
        plot(tA,ZA,'r.');
        xlabel('time');
        ylabel('Z');
        datetick('x','HH:MM');
        ax = axis;
        mZ = mean(Z);
        Q3 = quantile(ZA,1-r_normal);
        Q1 = quantile(ZA,r_normal);
        if Q1==Q3
            Q3 = 0.5;
            Q1 = -0.5;
        end
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), Q1-1.5*(Q3-Q1), Q3+1.5*(Q3-Q1)]);
        print([dirnum '-' date_string '_T.jpg'],'-djpeg')
        
        %
        figure(5);
        plot(tt,N_original/(Sample_Period*Sample_Rate),'.');
        xlabel('time');
        ylabel('Effective Data Rate');
        datetick('x','HH:MM');
        ax = axis;
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), 0, 1.2]);     
        print([dirnum '-' date_string '_N.jpg'],'-djpeg')
        %}
        intfreqX = FX2;
        intfreqY = FY2;
        intfreqZ = FZ2;
        %save all data 
        %eval(sprintf('save %s-%s-%s.mat tt f0 FX2 FY2 FZ2 tA XA YA ZA',dirnum,date_string,'basicinfo'));
        eval(sprintf('save %s-%s-%s.mat tt EX EY EZ Amp_meanX Amp_meanY Amp_meanZ Amp_varX Amp_varY Amp_varZ intfreqX intfreqY intfreqZ N_original',dirnum,date_string,'features'));
    end
end
diary off
