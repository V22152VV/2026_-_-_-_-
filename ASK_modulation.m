% パラメータ設定
bits_per_word = 9; 
bit_rate = 50e3;
fs = 5e6;
samples_per_bit = fs / bit_rate;

num_words = 130;   % セット数(ここを毎回チェック)
wave_filenames = cell(num_words, 1);  % ← これを忘れずに

carrier_freq = 200e3;   % キャリア周波数
mod_index = 0.5;        % 変調度（0〜1） (ここを毎回チェック)

amp0 = 1 - mod_index;   % ビット0の振幅
amp1 = 1;               % ビット1の振幅（基準）

all_bitstreams = zeros(num_words, bits_per_word);  % ← 追加
for w = 1:num_words
    bit_data = randi([0 1], 1, bits_per_word);

    modulated_wave = [];
    all_bitstreams(w, :) = bit_data;

    for i = 1:bits_per_word
        t = (0:samples_per_bit-1)/fs;
        amp = amp1 * bit_data(i) + amp0 * (1 - bit_data(i));
        wave = amp * sin(2 * pi * carrier_freq * t);
        modulated_wave = [modulated_wave, wave];
    end

    % トリガパルス生成
    trigger_freq = 300e3;
    trigger_duration = 20e-6;
    trigger_amplitude = 0.8;
    trigger_len = round(trigger_duration * fs);
    t_trigger = (0:trigger_len-1) / fs;
    trigger = trigger_amplitude * sin(2 * pi * trigger_freq * t_trigger);

    full_wave = [trigger, modulated_wave];

    disp('生成されたビット列 (各行が1ワード):');
    disp(bit_data);

    % 書き出し処理（そのまま）
    wave_duration = bits_per_word * (1 / bit_rate);
    fs_awes = 2e7;
    awes_points = round(fs_awes * (wave_duration + trigger_duration));
    max_points = 16384; %4096 %16384
    if awes_points > max_points
        fprintf('最大サンプル数を超えているのでデコードできません\n');
    end

    resampled_wave = resample(full_wave, fs_awes, fs);
    resampled_wave = resampled_wave - mean(resampled_wave);
    resampled_wave = resampled_wave / max(abs(resampled_wave));
    resampled_wave = resampled_wave * 6500; %500 %6500

    filename = sprintf('ask_modidx%.2f_bit9_%02d_awe.csv', mod_index, w);
    wave_filenames{w} = filename;  % ←ここで名前を保存しておく
    fid = fopen(filename, 'w');
    fprintf(fid, 'Start:,0,\n');
    fprintf(fid, 'Length:,%d,\n', awes_points);
    fprintf(fid, 'Sample Rate:,%.1f,\n', fs_awes);

    for i = 1:awes_points
        fprintf(fid, '%.8f,\n', resampled_wave(i));
    end
    fclose(fid);
    disp(['保存完了: ', filename]);
end

filename = 'bit_data.xlsx';

if isfile(filename)
    [~, ~, existing_data] = xlsread(filename);
else
    existing_data = {};
end

% 追記データ準備（ビット列文字列）
data_to_save = cell(num_words, 1);
for i = 1:num_words
    bit_str = char(join(string(all_bitstreams(i, :)), ''));
    data_to_save{i} = bit_str;
end

% ヘッダー
headers = {'BitString', 'Filename'};

% ここで「波形ファイル名」を2列目にセット（ループ内で保存したファイル名を使う）
file_column = wave_filenames;

% file_columnが横ベクトルであれば、縦ベクトルに転置
if isrow(file_column)
    file_column = file_column';
end

% これで列方向の連結が可能になる
append_data = [data_to_save, file_column];

% 空行を1行空ける
if ~isempty(existing_data)
    existing_data(end+1, :) = {''};
end

% ヘッダーの有無判定
has_header = ~isempty(existing_data) && ischar(existing_data{1,1}) && strcmp(existing_data{1,1}, headers{1});

if has_header
    all_data = [existing_data; append_data];
else
    all_data = [headers; existing_data; append_data];
end

% Excelに書き込み
writecell(all_data, filename);

fprintf('Excelファイル "%s" にビット列＋対応波形ファイル名を1行空けて追記しました。\n', filename);