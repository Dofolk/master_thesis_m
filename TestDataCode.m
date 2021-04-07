clear
close all
% record output messages
diary on

% parameters
dirnum = 'autotest';
start_date = '2020_2_1';
end_date   = '2020_2_1';
Sample_Rate = 400;
Sample_Period = 300;
Amp_max = 0.01;
r_normal = 0.001;
f0 = (1/Sample_Period:1/Sample_Period:(Sample_Rate/2))';
% MX = 1.9100e-04;
% MY = 1.4870e-04;
% MZ = 1.5211e-04;

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
%     UE = []; LE = [];
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
            if B{15}(1) == '¤U'
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
            
%             if isempty(X)
%                 fprintf('%s has data problem.',Xfile);
%                 eval(sprintf('movefile %s Problem_data',Xfile(1:75)));
%                 continue;
%             end
%             if isempty(Y)
%                 fprintf('%s has data problem.',Yfile);
%                 eval(sprintf('movefile %s Problem_data',Yfile(1:75)));
%                 continue
%             end
%             if isempty(Z)
%                 fprintf('%s has data problem.',Zfile);
%                 eval(sprintf('movefile %s Problem_data',Zfile(1:75)));
%                 continue
%             end

            epy_control = 0;
            if(isempty(t_x))
                fprintf('%s has no data. \n',[B{13} '/' B{14} ' ' num2str(mod(str2num(B{16}),12)+c) ':' B{17} ':' BS{1} '-X']);
                epy_control = 1;
%                 movefile(Xfile,'Problem_data');
            end
            if(isempty(t_y))
                fprintf('%s has no data. \n',[B{13} '/' B{14} ' ' num2str(mod(str2num(B{16}),12)+c) ':' B{17} ':' BS{1} '-Y']);
                epy_control = 1;
%                 movefile(Yfile,'Problem_data');
            end
            if(isempty(t_z))
                fprintf('%s has no data. \n',[B{13} '/' B{14} ' ' num2str(mod(str2num(B{16}),12)+c) ':' B{17} ':' BS{1} '-Z']);
                epy_control = 1;
%                 movefile(Zfile,'Problem_data');
            end
            
            if epy_control
                continue;
            end


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
            
            t = interp1(t_x,t_x,(min(t_x):max(t_x))');         
            X = interp1(t_x,X,(min(t_x):max(t_x))');
            Y = interp1(t_y,Y,(min(t_y):max(t_y))');
            Z = interp1(t_z,Z,(min(t_z):max(t_z))');
            
            tA = [tA; double(t(:))/86400/Sample_Rate+t1];
            XA = [XA; X]; Q3X = [Q3X; quantile(X,1-r_normal)]; Q1X = [Q1X; quantile(X,r_normal)];
            YA = [YA; Y]; Q3Y = [Q3Y; quantile(Y,1-r_normal)]; Q1Y = [Q1Y; quantile(Y,r_normal)];
            ZA = [ZA; Z]; Q3Z = [Q3Z; quantile(Z,1-r_normal)]; Q1Z = [Q1Z; quantile(Z,r_normal)];
            
%             [u,l] = envelope(Z,5,'peak');
%             UE = [UE;u];
%             LE = [LE;l];        

            T = double(t)/Sample_Rate/86400+t1;
            N = length(t); NA = [NA; N];
            TT = N/Sample_Rate;
            f  = 1/TT:1/TT:(Sample_Rate/2); FA(k).t = t1; FA(k).f = f;
            j  = find(f0>1/TT,1,'first');
            
            FX  = zeros(size(f0));          
            F = fft(X)/(N/2);       FA(k).X = abs(F(2:(1+length(f))));
            FX(j:end) = interp1(f,abs(F(2:(1+length(f)))),f0(j:end),'nearest');
            FX(isnan(FX)) = 0;
            FX2 = [FX2, FX];
            EX = [EX;sqrt(sum(FX.^2))];
            
            FY = zeros(size(f0));
            F = fft(Y)/(N/2);       FA(k).Y = abs(F(2:(1+length(f))));
            FY(j:end) = interp1(f,abs(F(2:(1+length(f)))),f0(j:end),'nearest');
            FY(isnan(FY)) = 0;
            FY2 = [FY2, FY];
            EY = [EY;sqrt(sum(FY.^2))];
            
            FZ = zeros(size(f0));
            F = fft(Z)/(N/2);       FA(k).Z = abs(F(2:(1+length(f))));
            FZ(j:end) = interp1(f,abs(F(2:(1+length(f)))),f0(j:end),'nearest');
            FZ(isnan(FZ)) = 0;
            FZ2 = [FZ2, FZ];
            EZ = [EZ;sqrt(sum(FZ.^2))];
            
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
    
%     FX2(FX2>MX) = FX2(FX2>MX)-MX;
%     FX2(FX2<=MX) = 0;
%     FY2(FY2>MY) = FY2(FY2>MY)-MY;
%     FY2(FY2<=MY) = 0;
%     FZ2(FZ2>MZ) = FZ2(FZ2>MZ)-MZ;
%     FZ2(FZ2<=MZ) = 0;
    
    % make all figures,
    % figure 1-3, time-frequency for X, Y, Z
    if ~isempty(tt)
        %{
        figure(1);
        imagesc(tt,f0,FX2);
        datetick('x','HH:MM');
        caxis([0,Amp_max]);
        colorbar;
        set(gca,'YDir','normal');
        xlabel('time');
        ylabel('frequency');
        title(sprintf('X-axis Date: %s',datestr(tt(1),29)));
        print([dirnum '-' date_string '_X.jpg'],'-djpeg')
        
        figure(2);
        imagesc(tt,f0,FY2);
        datetick('x','HH:MM');
        caxis([0,Amp_max]);
        colorbar;
        set(gca,'YDir','normal');
        xlabel('time');
        ylabel('frequency');
        title(sprintf('Y-axis Date: %s',datestr(tt(1),29)));
        print([dirnum '-' date_string '_Y.jpg'],'-djpeg')
        
        figure(3);
        imagesc(tt,f0,FZ2);
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
        plot(tA,XA,'b.'); %,tt,Q3X,'k-',tt,Q1X,'k--');
        xlabel('time');
        ylabel('X');
        datetick('x','HH:MM');
        ax = axis;
        Q3 = quantile(XA,1-r_normal);
        Q1 = quantile(XA,r_normal);
        %axis([ax(1), ceil(ax(1)+1/86400), Q1-1.5*(Q3-Q1), Q3+1.5*(Q3-Q1)]);
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), Q1-1.5*(Q3-Q1), Q3+1.5*(Q3-Q1)]);
        %xlim([ax(1), ceil(ax(1)+1/86400)]);
        %ylim([-2.5 2.5]);
        title(sprintf('Date: %s',datestr(tt(1),29)));
        
        subplot(3,1,2);
        plot(tA,YA,'g.'); %,tt,Q3Y,'k-',tt,Q1Y,'k--');
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
        %xlim([ax(1), min([ax(2),ceil(ax(1)+1/86400)])]);
        %ylim([-2.5 2.5]);
        
        subplot(3,1,3);
        plot(tA,ZA,'r.'); %,tt,Q3Z,'k-',tt,Q1Z,'k--');
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
        %xlim([ax(1), min([ax(2),ceil(ax(1)+1/86400)])]);
        %ylim([-2.5 2.5]);
        
        print([dirnum '-' date_string '_T.jpg'],'-djpeg')
        
        %
        figure(5);
        plot(tt,N_original/(Sample_Period*Sample_Rate),'.');
        xlabel('time');
        ylabel('Effective Data Rate');
        datetick('x','HH:MM');
        ax = axis;
        axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), 0, 1.2]); 
        %axis([ax(1), min([ax(2),ceil(ax(1)+1/86400)]), 0, 1]);    
        print([dirnum '-' date_string '_N.jpg'],'-djpeg')
        %}
        intfreqX = [];intfreqY = [];intfreqZ = [];
        for i = 1:199
            intfreqX = [intfreqX; mean(FX2((i*300-10:i*300+10),:))];
            intfreqY = [intfreqY; mean(FY2((i*300-10:i*300+10),:))];
            intfreqZ = [intfreqZ; mean(FZ2((i*300-10:i*300+10),:))];
        end
        intfreqX = [intfreqX; mean(FX2((end-10:end),:))];
        intfreqY = [intfreqY; mean(FY2((end-10:end),:))];
        intfreqZ = [intfreqZ; mean(FZ2((end-10:end),:))];
        
        intfreqX2 = [];intfreqY2 = [];intfreqZ2 = [];
        intfreqX2 = [intfreqX2; FX2(300:300:60000,:)];
        intfreqY2 = [intfreqY2; FY2(300:300:60000,:)];
        intfreqZ2 = [intfreqZ2; FZ2(300:300:60000,:)];
        
        % save all data 
        %eval(sprintf('save %s-%s-%s.mat tt f0 FX2 FY2 FZ2 tA XA YA ZA',dirnum,date_string,'basicinfo'));
        eval(sprintf('save %s-%s-%s.mat EX EY EZ Amp_meanX Amp_meanY Amp_meanZ Amp_varX Amp_varY Amp_varZ intfreqX intfreqY intfreqZ N_original',dirnum,date_string,'features'));
    end
end
diary off
