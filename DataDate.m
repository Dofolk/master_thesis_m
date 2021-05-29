T = [];
id = 1;
for j = 1:5
    s = num2str(j);
    foldertmp = ['D:' filesep '震動資料' filesep '主機2020.0' s]; 
    folder = [foldertmp filesep];
    Fi = dir([folder sprintf('*Ch_X*.txt')]);
    for i = 1:length(Fi)
        B = textsplit(Fi(i).name,'_');
        BS = textsplit(B{18},'.');
        if B{15}(1) == '下'
            c = 12;
        else
            c = 0;
        end
        t1 = datenum([B{12} '-' B{13} '-' B{14} ' ' num2str(mod(str2num(B{16}),12)+c) ':' B{17} ':' BS{1}]);
        T(id) = t1;
        id = id + 1;
    end
end

function str_output = textsplit(str_input,delimiter)
loc = find(str_input==delimiter);
s = 1;
k = 1;
for i=1:length(loc)
    str_output{k} = str_input(s:loc(i)-1);
    k = k + 1;
    s = loc(i)+1;
end
str_output{k} = str_input(s:end);
end