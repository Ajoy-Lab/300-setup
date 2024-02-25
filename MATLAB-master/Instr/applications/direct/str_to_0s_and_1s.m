
function arr = str_to_0s_and_1s(str1)
    % Explanation about ASCII code: https://www.ascii-code.com/ -- uses 7
    % bits to store information on charcaters
    arr_utf_16 = uint16(str1);
    char_arr = dec2bin(arr_utf_16);
    size_char_arr = size(char_arr);
    length_char_arr = size_char_arr(1) * size_char_arr(2);
    arr = [];
    for i = (1:length_char_arr)
        arr(end+1) = str2double(char_arr(i));
    end
end