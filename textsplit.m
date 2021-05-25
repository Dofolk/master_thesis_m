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