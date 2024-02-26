
function arr = str_to_0s_and_1s(str1)
    % Explanation about ASCII code: https://www.ascii-code.com/ -- uses 7
    % bits to store information on charcaters
    arr_utf_16 = uint16(str1);
    char_arr = dec2bin(arr_utf_16);
    size_char_arr = size(char_arr);
    str_arr = '';
    for i = (1:size_char_arr(1))
        str_arr = append(str_arr, char_arr(i,:));
    end
    arr = [];
    for i = (1:length(str_arr))
        if str_arr(i) == '1'
            arr(end+1) = 1;
        elseif str_arr(i) == '0'
            arr(end+1) = 0;
        else
            fprintf("str_arr should be only 1s and 0s")
        end
    end
end