% パラメータ設定
bits_per_word = 4;           % 情報ビット数（必ず4の倍数に）
num_words = 300;
bit_rate = 50e3;
fs = 5e6;
samples_per_bit = fs / bit_rate;
samples_per_symbol = samples_per_bit * 4;  % 1シンボル=4ビット

carrier_freq = 200e3;

% --- 16QAM用レベルマッピング
level_map = [-3, -1, +1, +3];  % 00, 01, 10, 11

% ビット列保存用配列
all_bitstreams = zeros(num_words, bits_per_word + 4);  % プリアンブル4ビット
append_data = strings(num_words, 2);

for w = 1:num_words
    % 情報ビット生成（4の倍数）
    info_bits = randi([0 1], 1, bits_per_word);
    preamble = [0 0 0 0];
    
    bit_data = [preamble, info_bits];
    all_bitstreams(w, :) = bit_data;

    % --- 4ビットずつ → 16QAMシンボル
    bit_quads = reshape(bit_data, 4, []).';  % 各行が [b1 b2 b3 b4]

    % 各シンボルについて、I = f(b1,b2), Q = f(b3,b4)
    I_bits = bit_quads(:, 1:2);
    Q_bits = bit_quads(:, 3:4);

    I_idx = I_bits(:,1)*2 + I_bits(:,2) + 1;
    Q_idx = Q_bits(:,1)*2 + Q_bits(:,2) + 1;

    I_vals = level_map(I_idx);
    Q_vals = level_map(Q_idx);

    % --- 変調波形生成
    modulated_wave = [];
    for i = 1:length(I_vals)
        t = (0:samples_per_symbol-1) / fs;
        I_comp = I_vals(i) * cos(2*pi*carrier_freq*t);
        Q_comp = Q_vals(i) * sin(2*pi*carrier_freq*t);
        wave = I_comp + Q_comp;
        modulated_wave = [modulated_wave, wave];
    end

    % --- トリガパルス生成
    trigger_freq = 300e3;
    trigger_duration = 20e-6;
    trigger_amplitude = 3.0;
    trigger_len = round(trigger_duration * fs);
    t_trigger = (0:trigger_len-1) / fs;
    trigger = trigger_amplitude * sin(2 * pi * trigger_freq * t_trigger);

    % --- トリガと変調波形結合
    full_wave = [trigger, modulated_wave];

    % --- 表示
    disp(['セット ', num2str(w), ' のビット列（4ビット単位で16QAMシンボルに変換）:']);
    disp(bit_data);  

    disp('対応するIレベル:');
    disp(I_vals');
    disp('対応するQレベル:');
    disp(Q_vals');

    % 波形長（秒）
    wave_duration = length(I_vals) * 4 / bit_rate;

    % --- AWES用設定
    fs_awes = 2e7;
    awes_points = round(fs_awes * (wave_duration + trigger_duration));
    max_points = 4096;

    if awes_points > max_points
        fprintf('最大サンプル数 %d を超えているのでデコードできません\n', max_points);
    end

    % リサンプリング
    resampled_wave = resample(full_wave, fs_awes, fs);

    % スケーリング (-500 ～ +500)
    resampled_wave = resampled_wave - mean(resampled_wave);
    resampled_wave = resampled_wave / max(abs(resampled_wave));
    resampled_wave = resampled_wave * 6500;

    % ファイル保存
    filename = sprintf('16qam_4bit_%02d_awe.csv', w);
    excel_filenames{w} = filename;
    fid = fopen(filename, 'w');
    fprintf(fid, 'Start:,0,\n');
    fprintf(fid, 'Length:,%d,\n', awes_points);
    fprintf(fid, 'Sample Rate:,%.1f,\n', fs_awes);
    for i = 1:awes_points
        fprintf(fid, '%.8f,\n', resampled_wave(i));
    end
    fclose(fid);
    disp(['保存完了: ', filename]);

    % ビット列保存
    bit_char = char(join(string(bit_data), ""));
    append_data(w,1) = string(bit_char);
    append_data(w,2) = string(filename);
end

% --- Excelファイル追記処理（元コードと同様）
excel_filename = "bit_data.xlsx";
if isfile(excel_filename)
    existing_data = readcell(excel_filename);
    existing_data = string(existing_data);
else
    existing_data = strings(0, 2);
end

if ~isempty(existing_data)
    existing_data = [existing_data; strings(1, 2)];
end

headers = ["BitString", "Filename"];
append_data_str = string(append_data);
has_header = ~isempty(existing_data) && existing_data(1,1) == headers(1);

if has_header
    all_data = [existing_data; append_data_str];
else
    all_data = [headers; existing_data; append_data_str];
end

all_data_cell = cellstr(all_data);
writecell(all_data_cell, excel_filename);
fprintf('Excelファイル "%s" にビット列とファイル名を追記しました。\n', excel_filename);