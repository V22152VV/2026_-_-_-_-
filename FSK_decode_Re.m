% === ファイル一覧取得 ===
files = dir('fsk_ost2140_m25_*.csv');
num_files = length(files);      % 読み込んだファイルの数を取得
disp(num_files); % ファイルの数を表示
decode_results = {};  

for k = 1:length(files)
    filename = files(k).name;
    fprintf('\n==== %s の処理 ====\n', filename);

    % === ファイル読み込み ===
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    opts.DataLines = 2;
    data = readmatrix(filename, opts);

    time = data(:,1);
    ch1 = data(:,2);
    ch2 = data(:,3);    

    % === サンプリング周波数推定 ===
    fs = 1 / mean(diff(time));
    fprintf('fs = %.3e\n', fs);
    f0 = 200e3;  % '0'の周波数
    f1 = 225e3;  % '1'の周波数
    num_bits = 9;
    bit_rate = 50e3;
    samples_per_bit = round(fs / bit_rate);

    % === CH2変調開始点の検出 ===
    start_ch2 = detect_fsk_start(ch2, fs, f0, f1, samples_per_bit);
    if isnan(start_ch2)
        warning("CH2: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    % トリガ終了点の最小値をデコード開始点に設定
    decode_start_idx = start_ch2;
    decode_end_idx = decode_start_idx + samples_per_bit * num_bits - 1;

    % === 結果表示 ===
    %{
    fprintf('変調開始点 (CH1): %d\n', start_ch1);
    %}
    fprintf('変調開始点 (CH2): %d\n', start_ch2);
    fprintf('変調終了点: %d\n', decode_end_idx);

    % === 対象波形を切り出し ===
    wave1 = ch1(decode_start_idx : decode_end_idx);
    wave2 = ch2(decode_start_idx : decode_end_idx);

    figure;
    subplot(2,1,1)
    plot(wave1);
    title(filename);
    xlabel('Time [ms]');
    ylabel('Amplitude');
    grid on;

    subplot(2,1,2)
    plot(wave2);
    title(filename);
    xlabel('Time [ms]');
    ylabel('Amplitude');
    grid on;

    % === FSK復調 ===
    bitstream1 = decode_fsk(wave1,f0, f1, fs, num_bits, samples_per_bit);
    bitstream2 = decode_fsk(wave2,f0, f1, fs, num_bits, samples_per_bit);
    
    % === 通常ビット列とBER ===
    fprintf('\n--- CH1 ビット列 ---\n');
    disp(bitstream1);
    fprintf('\n--- CH2 ビット列 ---\n');
    disp(bitstream2);
    bit_errors = sum(bitstream1 ~= bitstream2);
    ber = bit_errors / num_bits;
    fprintf('\nBER (非回転): %.2f%% (%d / %d)\n', ber * 100, bit_errors, num_bits);

    % デコード結果を行ごとに追加

    % ファイル名だけを取得（パスを除去）
    [~, filename, ext] = fileparts(files(k).name);  % 'files(i).name'でファイル名を取得
    csv_name = [filename, ext];  % 拡張子を付けたファイル名

    % デコード結果をセル配列に格納
    decode_results{k, 1} = datetime('now');           % デコード日時
    decode_results{k, 2} = csv_name;                  % ファイル名（例: 'test_18_01.csv'）
    decode_results{k, 3} = sprintf('%d', bitstream1); % 出力ビット列（例: '0110101'）
    decode_results{k, 4} = sprintf('%d', bitstream2); % 入力ビット列（例: '0110111'）
    decode_results{k, 5} = ber;                        % ビット誤り率
    decode_results{k, 6} = '';             % 備考欄

    % 平均ビット誤り率を追加()
    average_ber = mean([decode_results{:, 5}]);
    decode_results{k + 1, 1} = datetime('now');   % デコード日時
    decode_results{k + 1, 2} = '平均';            % ラベル
    decode_results{k + 1, 3} = '';             
    decode_results{k + 1, 4} = '';                
    decode_results{k + 1, 5} = average_ber;       % 平均BER
    decode_results{k + 1, 6} = '';     % 備考欄
  
end

if isempty(mfilename)
    % スクリプトとして実行されていて mfilename が使えない場合
    [~, script_name] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    script_name = mfilename;
end

% 最初に - または _ が出る位置を探す
token = regexp(script_name, '^[^_-]+', 'match');

if ~isempty(token)
    sheet_name = token{1};
else
    sheet_name = script_name;  % 記号がなければ全体を使う
end

excel_filename = 'decode.xlsx';

% 結果を書き込む
write_to_excel(excel_filename, sheet_name, decode_results);
disp('Excelに書き込みました');

% === 結果をエクセルに書き込む関数 ===
function write_to_excel(filename, sheet_name, data)
    % Excelファイルが存在するか確認
    if isfile(filename)
        % 既存のファイルがある場合は、シートが存在するか確認
        try
            % シートが存在するかチェック
            opts = detectImportOptions(filename);
            opts.Sheet = sheet_name;
            % 読み込み処理を試みる
            readtable(filename, opts);
            % シートが存在する場合、最初の空行にデータを書き込む
            start_row = height(readtable(filename, opts)) +1;
        catch
            % シートが存在しない場合、新しいシートを作成
            start_row = 2;  % ヘッダーを除いて次の行から書き込む
            % ヘッダーを追加
            header = {'デコード日時','デコードファイル' '出力ビット', '入力ビット', 'ビット誤り率', '備考欄'};
            writecell(header, filename, 'Sheet', sheet_name, 'Range', 'A1');
        end
    else
        % ファイルが存在しない場合、新規作成
        start_row = 2;
        header = {'デコード日時','デコードファイル' '出力ビット', '入力ビット', 'ビット誤り率', '備考欄'};
        writecell(header, filename, 'Sheet', sheet_name, 'Range', 'A1');
    end
    
    % データの書き込み
    writecell(data, filename, 'Sheet', sheet_name, 'Range', sprintf('A%d', start_row));
end

function bits = decode_fsk(wave, f0, f1, fs, num_bits, samples_per_bit)
    bits = zeros(1, num_bits);
    
    for k = 1:num_bits
        idx_start = (k-1)*samples_per_bit + 1;
        idx_end = k*samples_per_bit;
        segment = wave(idx_start:idx_end);

        % 平均を引いて直流成分除去
        segment = segment - mean(segment);

        % ビット区間の長さで窓関数作成＆掛ける
        w = hanning(samples_per_bit);
        windowed_segment = segment .* w;

        % Goertzel計算（非整数k対応版）
        power0 = goertzel_nonint_k(windowed_segment, fs, f0);
        power1 = goertzel_nonint_k(windowed_segment, fs, f1);

        % 判定（power0 > power1なら0、そうでなければ1）
        if power0 > power1
            bits(k) = 0;
        else
            bits(k) = 1;
        end
    end
end

function power = goertzel_nonint_k(bit_wave, fs, f)
    N = length(bit_wave);
    w = 2 * pi * f / fs;  % ここでkは使わず周波数そのものからw計算
    coeff = 2 * cos(w);

    s_prev = 0;
    s_prev2 = 0;

    for n = 1:N
        s = bit_wave(n) + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end

    power = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
end

function power = goertzel_power(bit_wave, fs, f)
    N = length(bit_wave);
    k = round(f / fs * N);
    w = 2 * pi * k / N;
    coeff = 2 * cos(w);

    s_prev = 0;
    s_prev2 = 0;

    for n = 1:N
        s = bit_wave(n) + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end

    power = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
end

function start_idx = detect_fsk_start(signal, fs, f0, f1, samples_per_bit)
    trigger_freq = 350e3;
    trigger_threshold = 0.6;  % 正規化した振幅で判定
    threshold_fsk = 0.4;
    step = round(samples_per_bit / 4);

    found_trigger = false;
    trigger_ended = false;
    low_count = 0;
    low_threshold = 2;

    max_p_trigger = 0;
    start_idx = NaN;

    i = 1;
    while i <= length(signal) - samples_per_bit
        segment = signal(i : i + samples_per_bit - 1);
        segment = segment - mean(segment);
        segment = segment / max(abs(segment));

        p_trigger = goertzel_power(segment, fs, trigger_freq);
        p0 = goertzel_power(segment, fs, f0);
        p1 = goertzel_power(segment, fs, f1);

        if p_trigger > max_p_trigger
            max_p_trigger = p_trigger;
        end

        norm_p_trigger = p_trigger / max(max_p_trigger, 1e-10);

        if ~found_trigger
            if norm_p_trigger > trigger_threshold && p_trigger > p0 && p_trigger > p1
                found_trigger = true;
            end
        elseif found_trigger && ~trigger_ended
            if norm_p_trigger < 0.3
                low_count = low_count + 1;
                if low_count >= low_threshold
                    trigger_ended = true;
                end
            else
                low_count = 0;
            end
        elseif trigger_ended
            if p0 > threshold_fsk || p1 > threshold_fsk
                % ★ 補正：少し前にずらす（最大でもsamples_per_bit/2だけ）
                corrected_idx = i + floor(samples_per_bit / 12);
                %corrected_idx = i;
                start_idx = max(1, corrected_idx);  % マイナス防止
                break;
                
            end
        end

        i = i + step;
    end
end