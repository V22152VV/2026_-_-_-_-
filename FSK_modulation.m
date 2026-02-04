% パラメータ設定
bits_per_word = 9;           % 各セットのビット数 %9 %11
bit_rate = 50e3;             % ビットレート（50 kbps）
fs = 5e6;                    % 波形生成時のサンプリング周波数（5 MHz）
samples_per_bit = fs / bit_rate;

% 周波数設定(ここを毎回チェック)
f0 = 200e3;                  % '0'の周波数
f1 = 225e3;                  % '1'の周波数
delta_f_khz = round((f1 - f0) / 1e3);  % 差分をkHzにして四捨五入

%セット数(ここを毎回チェック)
num_words = 150;     

wave_filenames = cell(num_words, 1); 

% 変調波形生成
modulated_wave = [];
for w = 1:num_words
    bit_data = randi([0 1], 1, bits_per_word);
    modulated_wave = [];
   
    all_bitstreams(w, :) = bit_data;

    for i = 1:bits_per_word
        t = (0:samples_per_bit-1)/fs;
        freq = f0 * (bit_data(i)==0) + f1 * (bit_data(i)==1);

        % 位相の連続性を考慮しない波形
        phase = 2 * pi * freq * t;
        wave = sin(phase);

        modulated_wave = [modulated_wave, wave];
    end

    % トリガパルス生成（正弦波バージョン）
    trigger_freq = 350e3;           % トリガ正弦波の周波数（例：2MHz）
    trigger_duration = 20e-6;     % トリガパルスの継続時間（例：20μs）
    trigger_amplitude = 0.8;      % トリガパルスの振幅
    trigger_len = round(trigger_duration * fs);  % トリガパルスの長さ（サンプル数）

    t_trigger = (0:trigger_len-1) / fs;  % 時間ベクトル
    trigger = trigger_amplitude * sin(2 * pi * trigger_freq * t_trigger);  % 正弦波生成

    % トリガを波形に追加
    full_wave = [trigger, modulated_wave];

    % 生成したビットデータを表示
    disp('生成されたビット列 (各行が1ワード):');
    disp(bit_data);

    % 波形の時間長さ (秒)
    wave_duration = bits_per_word * (1 / bit_rate);  % 5ビット分の波形の時間 (秒)

    % 保存時のサンプリング周波数（AWES用）
    fs_awes = 2e7;         % 保存用サンプリング周波数（20 MHz）

    % サンプル数を計算
    awes_points = round(fs_awes * (wave_duration + trigger_duration));  % 波形の総サンプル数を計算

    % 最大サンプル数の制限 (例えば4096サンプル)
    max_points = 16384; %4096 %16384

    if awes_points > max_points
        fprintf('最大サンプル数を超えているのでデコードできません');
    end

    % FFTベースのリサンプリング
    N = length(full_wave);  % 元のサンプル数
    freq_original = fft(full_wave);  % FFTで周波数領域に変換

    % 新しいサンプリング周波数に合わせた新しいサンプル数
    N_new = round(N * fs_awes / fs);  % 新しいサンプル数

    % 周波数領域でサンプリングを変更するための準備
    freq_new = zeros(1, N_new);
    half_old = floor(N / 2);
    half_new = floor(N_new / 2);

    % 元の周波数成分を新しい周波数成分にコピー
    freq_new(1:half_old+1) = freq_original(1:half_old+1);
    freq_new(end-half_old:end) = freq_original(end-half_old:end);

    % 新しい周波数領域で逆FFTを実行
    resampled_wave = ifft(freq_new, 'symmetric');  % 逆FFTで時間領域に戻す

    % 時間長さを調整
    t_new = (0:1/fs_awes:(N_new-1)/fs_awes);  % 新しい時間軸

    % スケーリング（-500〜+500）
    resampled_wave = resampled_wave - mean(resampled_wave);  % DC除去
    resampled_wave = resampled_wave / max(abs(resampled_wave));  % 振幅を正規化
    resampled_wave = resampled_wave * 6500; %500 %6500

    %ファイルの名前設定(ここを毎回チェック)
    filename = sprintf('fsk_m%dk_bit9_%02d_awe.csv', delta_f_khz, w);
    wave_filenames{w} = filename;  % ←ここで名前を保存しておく

    % 書き出し
    fid = fopen(filename, 'w');
    fprintf(fid, 'Start:,0,\n');
    fprintf(fid, 'Length:,%d,\n', N_new);  % 保存用のサンプル数
    fprintf(fid, 'Sample Rate:,%.1f,\n', fs_awes);  % 保存用のサンプリング周波数を記載

    for i = 1:N_new
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