% パラメータ設定
bits_per_word = 6;            % 情報ビット数（偶数にしてください、例：8ビット）
num_words = 200;               % セット数
bit_rate = 50e3;              % ビットレート（50 kbps）
fs = 5e6;                    % 波形生成時のサンプリング周波数（5 MHz）
samples_per_bit = fs / bit_rate;
samples_per_symbol = samples_per_bit * 2;

carrier_freq = 200e3;         % QPSKの搬送波周波数

% --- パラメータ ---
bits_per_word = 6;  % 必ず偶数にする

if mod(bits_per_word, 2) ~= 0
    error('bits_per_word は偶数である必要があります。');
end

% --- ビット列保存用配列（各行が1ワード分） ---
all_bitstreams = zeros(num_words, bits_per_word + 2);  % プリアンブル2ビットを含む
append_data = strings(num_words, 2);

for w = 1:num_words
    % --- 情報ビットの生成（偶数長でランダム） ---
    info_bits = randi([0 1], 1, bits_per_word);
    preamble = [0 0];

    % --- プリアンブル + 情報ビット結合 ---
    bit_data = [preamble, info_bits];
    all_bitstreams(w, :) = bit_data;

    % --- 2ビットずつペアにしてQPSKシンボルに変換 ---
    bit_pairs = reshape(bit_data, 2, []).';  % 各行が [bit1 bit2]
    symbols = bit_pairs(:,1)*2 + bit_pairs(:,2);  % 0~3の整数

    bit1 = floor(symbols / 2);     % 上位ビット（MSB）
    bit2 = mod(symbols, 2);        % 下位ビット（LSB）

    % --- QPSK変調波形生成 ---
    modulated_wave = [];
    for i = 1:length(symbols)
        t = (0:samples_per_symbol-1) / fs;
        phase = (pi/2) * symbols(i);  % 0→0°, 1→90°, 2→180°, 3→270°
        wave = sin(2 * pi * carrier_freq * t + phase);
        modulated_wave = [modulated_wave, wave];
    end
    
    % --- トリガパルス生成（正弦波） ---
    trigger_freq = 300e3;
    trigger_duration = 20e-6;
    trigger_amplitude = 0.8;
    trigger_len = round(trigger_duration * fs);
    t_trigger = (0:trigger_len-1) / fs;
    trigger = trigger_amplitude * sin(2 * pi * trigger_freq * t_trigger);
    
    % --- トリガと変調波形結合 ---
    full_wave = [trigger, modulated_wave];

    % --- 表示 ---
    disp(['セット ', num2str(w), ' のビット列（2ビット単位でQPSKシンボルに変換）:']);
    disp(bit_data);  
    disp('対応するQPSKシンボル列（0〜3）:');
    disp(symbols(:)');  % 縦ベクトルでも横ベクトルに変換
    
    % --- 波形時間長さ（秒） ---
    wave_duration = length(symbols) * 2 / bit_rate;  % 1シンボル=2ビット
    
    % --- AWES用設定 ---
    fs_awes = 2e7;
    awes_points = round(fs_awes * (wave_duration + trigger_duration));
    max_points = 4096;
    
    if awes_points > max_points
        fprintf('最大サンプル数 %d を超えているのでデコードできません\n', max_points);
    end
    
    % --- リサンプリング ---
    resampled_wave = resample(full_wave, fs_awes, fs);
    
    % --- スケーリング (-500〜+500) ---
    resampled_wave = resampled_wave - mean(resampled_wave);
    resampled_wave = resampled_wave / max(abs(resampled_wave));
    resampled_wave = resampled_wave * 6500;
    
    % --- ファイル保存 ---
    filename = sprintf('qpsk_8bit_%02d_awe.csv', w);
    excel_filenames{w} = filename;  % ←ここで名前を保存しておく
    fid = fopen(filename, 'w');
    fprintf(fid, 'Start:,0,\n');
    fprintf(fid, 'Length:,%d,\n', awes_points);
    fprintf(fid, 'Sample Rate:,%.1f,\n', fs_awes);
    for i = 1:awes_points
        fprintf(fid, '%.8f,\n', resampled_wave(i));
    end
    fclose(fid);
    disp(['保存完了: ', filename]);

    bit_str = num2str(bit_data);
    bit_char = char(join(string(bit_data), ""));
    append_data(w,1) = string(bit_char);
    append_data(w,2) = string(filename);
end

excel_filename = "bit_data.xlsx";
% まず既存データを読み込む（または初期化）
if isfile(excel_filename)
    existing_data = readcell(excel_filename);
    existing_data = string(existing_data);  % string配列に変換
else
    existing_data = strings(0, 2);  % 新規ファイル用の空データ
end

% 空行を1行空ける（既存データがある場合）
if ~isempty(existing_data)
    existing_data(end+1, :) = {"", ""};  % ← 空セル2列分に修正（エラー対策）
end

% ヘッダー
headers = ["BitString", "Filename"];

% 既存データの読み込み（string配列に変換）
if isfile(excel_filename)
    existing_data = readcell(excel_filename);
    existing_data = string(existing_data);  % string配列に変換
else
    existing_data = strings(0, 2);
end

% 空行を1行空ける（既存データがある場合）
if ~isempty(existing_data)
    existing_data = [existing_data; strings(1, 2)];
end

% append_data を string 配列に変換
append_data_str = string(append_data);

% ヘッダー有無の判定
has_header = ~isempty(existing_data) && existing_data(1,1) == headers(1);

if has_header
    all_data = [existing_data; append_data_str];
else
    all_data = [headers; existing_data; append_data_str];
end

% all_data は string 配列なので cell 配列に変換
all_data_cell = cellstr(all_data);

% 書き込み
writecell(all_data_cell, excel_filename);

fprintf('Excelファイル "%s" にビット列とファイル名を追記しました。\n', excel_filename);