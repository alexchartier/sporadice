function [lat, lon] = maidenhead_to_ll(code)

Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
kv = 0:25;


cr = cellstr(upper(code)')';


%% Longitude
Field = kv(Alphabet == cr{1}) * 20;
Square = str2double(cr{3}) * 2;

if length(code) == 6
    SubSquareLow = kv(Alphabet == cr{5}) * (2/24);
    SubSquareHigh = SubSquareLow + (2/24);
elseif length(code) == 4
    SubSquareLow = 0;
    SubSquareHigh = 0;
end
StartLon = Field + Square + SubSquareLow - 180;
EndLon = Field + Square + SubSquareHigh - 180;

lon = mean([StartLon, EndLon]);


%% Latitude

Field = kv(Alphabet == cr{2}) * 10;
Square = str2double(cr{4});
if length(code) == 6
    SubSquareLow = kv(Alphabet == cr{6}) * (1/24);
    SubSquareHigh = SubSquareLow + (1/24);
elseif length(code) == 4
    SubSquareLow = 0;
    SubSquareHigh = 0;
end

StartLat = Field + Square + SubSquareLow - 90;
EndLat = Field + Square + SubSquareHigh - 90;
lat = mean([StartLat, EndLat]);


% def main(strMaidenHead = __MH__):
