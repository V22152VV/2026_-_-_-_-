% パラメータ設定
bits_per_word = 7;             % 情報ビット数（7ビット）
num_words = 180;                % セット数
bit_rate = 50e3;               % ビットレート（50 kbps）
fs = 5e6;                      % 波形生成時のサンプリング周波数（5 MHz）
samples_per_bit = fs / bit_rate;

carrier_freq = 200e3;          % BPSKの搬送波周波数

for w = 1:num_words
    % --- 情報ビットの生成 ---
    info_bits = randi([0 1], 1, bits_per_word);

    % --- プリアンブルの挿入（固定: 0,1） ---
    preamble = [0 1];
    bit_data = [preamble, info_bits];  % 合計9ビット

    all_bitstreams(w, :) = bit_data;

% BPSK変調波形生成（cosベース）
modulated_wave = [];
for i = 1:length(bit_data)
    t = (0:samples_per_bit-1)/fs;
    phase = pi * bit_data(i);  % 0→0°, 1→180°
    wave = sin(2 * pi * carrier_freq * t + phase);
    modulated_wave = [modulated_wave, wave];
end

    % トリガパルス生成（正弦波）
    trigger_freq = 300e3;
    trigger_duration = 20e-6;
    trigger_amplitude = 0.8;
    trigger_len = round(trigger_duration * fs);
    t_trigger = (0:trigger_len-1) / fs;
    trigger = trigger_amplitude * sin(2 * pi * trigger_freq * t_trigger);

    % トリガと変調波形を結合
    full_wave = [trigger, modulated_wave];

    % 表示（ビット列）
    disp(['セット ', num2str(w), ' のビット列 (プリアンプル + 情報):']);
    disp(bit_data);

    % 波形の時間長さ (秒)
    wave_duration = length(bit_data) / bit_rate;

    % AWES用のサンプリング周波数とサンプル数
    fs_awes = 2e7;
    awes_points = round(fs_awes * (wave_duration + trigger_duration));
    max_points = 16384;

    if awes_points > max_points
        fprintf('最大サンプル数を超えているのでデコードできません\n');
    end

    % リサンプリング
    resampled_wave = resample(full_wave, fs_awes, fs);

    % スケーリング（-500〜+500）
    resampled_wave = resampled_wave - mean(resampled_wave);
    resampled_wave = resampled_wave / max(abs(resampled_wave));
    resampled_wave = resampled_wave * 6500;

    % ファイル保存
    filename = sprintf('bpsk_bit9_%02d_awe.csv', w);
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

% file_columnが縦ベクトルか確認し、違うなら転置する
if size(file_column, 1) == 1 && size(file_column, 2) == num_words
    file_column = file_column';  % 転置して縦ベクトルにする
end

% data_to_saveは縦ベクトルなので問題なし

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