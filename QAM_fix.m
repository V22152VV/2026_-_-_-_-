% === ファイル一覧取得 ===
files = dir('ost_qam_last_*.csv');
num_files = length(files);
disp(num_files);
decode_results = {};  
row_idx = 1;

for k = 1:length(files)
    filename = files(k).name;
    fprintf('\n==== %s の処理 ====\n', filename);

    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    opts.DataLines = 2;
    data = readmatrix(filename, opts);

    time = data(:,1);
    ch1 = data(:,2);
    ch2 = data(:,3);

    fs = 1 / mean(diff(time));
    fprintf('fs = %.3e\n', fs);
    carrier_freq = 200e3;
    num_symbols = 2;  % 1シンボル4ビットで2シンボル = 8ビット
    num_bits = num_symbols * 4;
    bit_rate = 50e3;
    samples_per_bit = round(fs / bit_rate);

    start_ch2 = detect_bpsk_start(ch2, fs, carrier_freq, samples_per_bit);
    if isnan(start_ch2)
        warning("CH2: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    decode_start_idx =  start_ch2;
    decode_end_idx = decode_start_idx + samples_per_bit * num_bits - 1;

    fprintf('変調開始点 (CH2): %d\n', start_ch2);
    fprintf('変調終了点: %d\n', decode_end_idx);

    wave1 = ch1(decode_start_idx : decode_end_idx);
    wave2 = ch2(decode_start_idx : decode_end_idx);

    bits_per_symbol = 1;
    num_bits = 8;
    fc = 200e3;
    bw = 20e3;
    wave1_true = freq_domain_smooth_qam(wave1, samples_per_bit, num_bits, fs, fc, bw, bits_per_symbol);

    % --- デコード（16QAM） ---

    % --- まだ±1,±3などにはしない、あくまでIとQの複素数にするだけ ---%
    symbol_pre1 = demod_16qam_symbols(wave1_true, fs, carrier_freq, num_symbols, samples_per_bit);
    symbol_2 = demod_16qam_symbols(wave2, fs, carrier_freq, num_symbols, samples_per_bit);

    % 新コード
    num_subsegments = 4;
    %symbol_pre1 = demod_16qam_symbols_averaged(wave1_true, fs, carrier_freq, num_symbols, samples_per_bit, num_subsegments);
    %symbol_2 = demod_16qam_symbols_averaged(wave2, fs, carrier_freq, num_symbols, samples_per_bit, num_subsegments);

    % --- 先頭シンボル（プリアンブル用）を取り出す --- %
    preamble_symbol_ch1 = symbol_pre1(1);
    preamble_symbol_ch2 = symbol_2(1);
    
    % --- ch1の振幅補正・位相補正 --- %
    symbol_1 = correct_phase_amplitude(preamble_symbol_ch1, preamble_symbol_ch2, symbol_pre1);

    % --- 16QAMのデコード、ビットで出てくる ---%
    bitstream1 = cluster16QAM_threshold(symbol_1);
    bitstream2 = cluster16QAM_threshold(symbol_2);  % CH2はそのまま

    % === 通常ビット列とBER ===
    fprintf('\n--- CH1 ビット列 ---\n');
    disp(bitstream1);
    fprintf('\n--- CH2 ビット列 ---\n');
    disp(bitstream2);
    bit_errors = sum(bitstream1 ~= bitstream2);
    ber = bit_errors / num_bits;
    fprintf('\nBER: %.2f%% (%d / %d)\n', ber * 100, bit_errors, num_bits);

    % デコード結果を行ごとに追加

    % ファイル名だけを取得（パスを除去）
    [~, filename, ext] = fileparts(files(k).name);  % 'files(i).name'でファイル名を取得
    csv_name = [filename, ext];  % 拡張子を付けたファイル名

    if  isnan(start_ch2)
        % デコードできなかった場合、BER = NaN としてスキップ（平均には含めない）
        decode_results{row_idx, 1} = datetime('now');
        decode_results{row_idx, 2} = csv_name;
        decode_results{row_idx, 3} = '---';
        decode_results{row_idx, 4} = '---';
        decode_results{row_idx, 5} = NaN;
        decode_results{row_idx, 6} = 'デコード失敗';
        row_idx = row_idx + 1;
        continue;  % 次のファイルへ
    end

    % デコード成功時の結果格納
    decode_results{row_idx, 1} = datetime('now');
    decode_results{row_idx, 2} = csv_name;
    decode_results{row_idx, 3} = sprintf('%d', bitstream1);
    decode_results{row_idx, 4} = sprintf('%d', bitstream2);
    decode_results{row_idx, 5} = ber;
    decode_results{row_idx, 6} = '';  % 備考空白
    row_idx = row_idx + 1;
end

valid_bers = [decode_results{:, 5}];
valid_bers = valid_bers(~isnan(valid_bers));
average_ber = mean(valid_bers);

% 平均BER追加
decode_results{row_idx, 1} = datetime('now');
decode_results{row_idx, 2} = '平均';
decode_results{row_idx, 3} = '';
decode_results{row_idx, 4} = '';
decode_results{row_idx, 5} = average_ber;
decode_results{row_idx, 6} = '';

% Excel書き出し
write_to_excel('decode.xlsx', 'decode16QAM', decode_results);
fprintf('\nBER: %.2f%%\n', average_ber * 100);
disp('Excelに書き込みました');

% --- 16QAMのIQ複素数(マッピングはまだ) --- %
function rx_symbols = demod_16qam_symbols(signal, fs, carrier_freq, num_symbols, samples_per_bit)
    samples_per_symbol = samples_per_bit * 4;
    rx_symbols = zeros(1, num_symbols);

    for i = 1:num_symbols
        idx_start = (i-1)*samples_per_symbol + 1;
        idx_end   = idx_start + samples_per_symbol - 1;
        if idx_end > length(signal)
            break;
        end

        x = signal(idx_start:idx_end);
        t = (0:(samples_per_symbol - 1)) / fs;

        ref_I = cos(2*pi*carrier_freq*t);
        ref_Q = -sin(2*pi*carrier_freq*t);

        I = sum(x .* ref_I(:));
        Q = sum(x .* ref_Q(:));

        rx_symbols(i) = I + 1j * Q;
    end
end

% --- ch1の位相補正と振幅補正 --- % 
function sym_all_corrected = correct_phase_amplitude(sym1, sym2, sym_all)
    % sym1, sym2 : 複素数シンボル列 (同じ長さ)
    % bitstream1 : CH1 のビット列（補正対象の全データ）
    % I2, Q2     : CH2 のプリアンブル I, Q
    % 出力       : 位相補正＋振幅補正済みの CH1 シンボル列

    rotations = [0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi, 7*pi/6, 4*pi/3, 3*pi/2, 5*pi/3, 11*pi/6];
    % rotations = (0:35) * (pi/18);  % 0, 10°, 20°, ..., 350°

    best_match = 0;
    min_error = inf;

    for km = 1:length(rotations)
        rotated_sym1 = sym1 * exp(-1j * rotations(km));  % ch1を回転させる
        % 誤差（MSE風）
        error = sum(abs(rotated_sym1 - sym2).^2);
        if error < min_error
            min_error = error;
            best_match = rotations(km);
        end
    end

    % --- CH1全体にこの回転を適用して補正する ---
    sym_all_corrected = sym_all * exp(-1j * best_match);

end

% --- 16QAMクラスタリング --- %
function bits = cluster16QAM_threshold(sym_corrected)
    % --- I/Q を取り出す ---
    I_vals = real(sym_corrected);
    Q_vals = imag(sym_corrected);

    % --- 最大振幅に基づく閾値決定 ---
    max_I = max(abs(I_vals));
    max_Q = max(abs(Q_vals));

    % ±1 / ±3 の境界を比率で設定（0.6: ±3、0.2: ±1）
    th_I = [max_I*0.2, max_I*0.6];
    th_Q = [max_Q*0.2, max_Q*0.6];

    % --- I のクラスタリング ---
    I_cluster = zeros(size(I_vals));
    for k = 1:length(I_vals)
        I = I_vals(k);
        if I > th_I(2)
            I_cluster(k) = 3;
        elseif I > th_I(1)
            I_cluster(k) = 1;
        elseif I > -th_I(1)
            I_cluster(k) = -1;
        else
            I_cluster(k) = -3;
        end
    end

    % --- Q のクラスタリング ---
    Q_cluster = zeros(size(Q_vals));
    for k = 1:length(Q_vals)
        Q = Q_vals(k);
        if Q > th_Q(2)
            Q_cluster(k) = 3;
        elseif Q > th_Q(1)
            Q_cluster(k) = 1;
        elseif Q > -th_Q(1)
            Q_cluster(k) = -1;
        else
            Q_cluster(k) = -3;
        end
    end

    % --- I/Q をビット列に変換（Gray符号なし） ---
    % I_cluster, Q_cluster ∈ {-3,-1,1,3} → 2ビットずつ
    N = length(I_cluster);
    bits = zeros(N, 4);
    for k = 1:N
        bits(k,:) = levels_to_bits(I_cluster(k), Q_cluster(k));
    end

    % --- フラットな1行ベクトルに変換 ---
    bits = reshape(bits.', [], 1).';
end

function b = levels_to_bits(I_level, Q_level)
    % -3 -> 0, -1 -> 1, 1 -> 2, 3 -> 3
    decI = (I_level + 3)/2;
    decQ = (Q_level + 3)/2;

    % 2ビットに変換
    bI = [bitget(decI,2), bitget(decI,1)];
    bQ = [bitget(decQ,2), bitget(decQ,1)];

    b = [bI, bQ];
end

function write_to_excel(filename, sheet_name, data)
    if isfile(filename)
        try
            opts = detectImportOptions(filename);
            opts.Sheet = sheet_name;
            readtable(filename, opts);
            start_row = height(readtable(filename, opts)) + 1;
        catch
            start_row = 2;
            header = {'デコード日時','デコードファイル','出力ビット', '入力ビット', 'ビット誤り率', '備考欄'};
            writecell(header, filename, 'Sheet', sheet_name, 'Range', 'A1');
        end
    else
        start_row = 2;
        header = {'デコード日時','デコードファイル','出力ビット', '入力ビット', 'ビット誤り率', '備考欄'};
        writecell(header, filename, 'Sheet', sheet_name, 'Range', 'A1');
    end
    writecell(data, filename, 'Sheet', sheet_name, 'Range', sprintf('A%d', start_row));
end

% === Goertzel関数 ===
function power = goertzel_exact(x, fs, f_target)
    N = length(x);
    k = f_target * N / fs;
    w = 2 * pi * k / N;
    coeff = 2 * cos(w);
    s_prev = 0; s_prev2 = 0;
    for n = 1:N
        s = x(n) + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end
    power = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
end

function [start_idx, trigger_range] = detect_bpsk_start(signal, fs, carrier_freq, samples_per_bit)
    trigger_freq = 300e3;
    threshold_fsk = 0.3;
    step = round(samples_per_bit/4);

    found_trigger = false;
    trigger_ended = false;
    i = 1;

    p_trigger_values = [];
    time_values = [];
    trigger_start = NaN;
    trigger_end = NaN;

    while i <= length(signal) - samples_per_bit
        segment = signal(i : i + samples_per_bit - 1);
        segment = segment - mean(segment);
        segment = segment / max(abs(segment));

        p_trigger = goertzel_exact(segment, fs, trigger_freq);
        p0 = goertzel_exact(segment, fs, carrier_freq);

        if isempty(p_trigger_values)
            p_trigger_values = [p_trigger_values, p_trigger];
        else
            p_trigger_values = [p_trigger_values, p_trigger / max(p_trigger_values)];
        end
        time_values = [time_values, i/fs];

        if ~found_trigger
            if p_trigger > p0 * 10
                found_trigger = true;
                trigger_start_idx = i;  % トリガ開始位置記録
            end
        elseif found_trigger && ~trigger_ended
            if p_trigger < p0 / 10
                trigger_ended = true;
                trigger_end_idx = i;  % トリガ終了位置記録
            end
        elseif trigger_ended
            if p0 > threshold_fsk || p1 > threshold_fsk
                % ★ 補正：少し前にずらす（最大でもsamples_per_bit/2だけ）
                corrected_idx = i - floor(samples_per_bit / 10);
                %corrected_idx = i;
                start_idx = max(1, corrected_idx);  % マイナス防止
                break;
            end
        end

        i = i + step;
    end

    if ~found_trigger
        start_idx = NaN;
        trigger_range = [NaN NaN];
    else
        trigger_range = [trigger_start, trigger_end];
    end
end


function waveform_out = freq_domain_smooth_qam(waveform, samples_per_bit, num_bits, fs, fc, bw, bits_per_symbol)
% QPSKや16QAMなど、複数ビット/シンボルに対応した周波数領域平滑化
% 改良版：reshapeエラー対策、入力チェック、複素対応、波形長の自動調整（長い場合のみ）

% 入力ベクトルを列ベクトルに（行ベクトルが来ても対応）
if isrow(waveform)
    waveform = waveform(:);
end

% 基本パラメータ
samples_per_symbol = samples_per_bit * bits_per_symbol;

% num_symbolsは整数であるべき
if mod(num_bits, bits_per_symbol) ~= 0
    error('num_bits が bits_per_symbol の倍数ではありません。num_bits=%d, bits_per_symbol=%d', num_bits, bits_per_symbol);
end
num_symbols = num_bits / bits_per_symbol;

% 必要なサンプル数
total_needed = samples_per_symbol * num_symbols;

% 長さチェック
len_wave = length(waveform);
if len_wave < total_needed
    error('waveform が短すぎます: length(waveform)=%d, 必要=%d', len_wave, total_needed);
elseif len_wave > total_needed
    warning('waveform が期待より長いので先頭 %d サンプルに切り詰めます (length=%d -> %d)。', total_needed, len_wave, total_needed);
    waveform = waveform(1:total_needed);
end

% reshape（ここでエラーが出るはずがない）
segments = reshape(waveform, samples_per_symbol, num_symbols);

% 周波数軸（シンボル長に合わせる）
f = (0:samples_per_symbol-1)' * fs / samples_per_symbol;

% 各シンボルのFFT（複素用で確保）
fft_segments = complex(zeros(size(segments)));
for k = 1:num_symbols
    fft_segments(:,k) = fft(segments(:,k));
end

% 信号帯域マスク（正負周波数両方を明示的に保護）
mask_pos = (f >= (fc - bw)) & (f <= (fc + bw));
mask_neg = (f >= (fs - (fc + bw))) & (f <= (fs - (fc - bw)));
mask = mask_pos | mask_neg;

% 帯域外を平均化（ビット/シンボル平均）
fft_smoothed = fft_segments;
for i = 1:samples_per_symbol
    if ~mask(i)
        % 全シンボルで平均した値で置き換え
        fft_smoothed(i,:) = mean(fft_segments(i,:));
    end
end

% IFFTで時間波形に戻す（列ごとに実行される）
segments_out = ifft(fft_smoothed, [], 1);  % 明示的に dim=1 指定
% 'symmetric' オプションは実数化のために使うこともあるが、複素信号を扱う場合は不要
% segments_out = ifft(fft_smoothed, 'symmetric');

% 連結して出力
waveform_out = reshape(segments_out, [], 1);

end

function [optimized_start, best_snr] = optimize_symbol_timing(signal, fs, carrier_freq, initial_start, samples_per_symbol, num_symbols)
    % === シンボルタイミング最適化（AWGN対策）===
    % 原理: 複数のタイミング候補でSNRを計算し、最良点を選択
    % 
    % 入力:
    %   signal              : 受信信号
    %   fs                  : サンプリング周波数
    %   carrier_freq        : 搬送波周波数
    %   initial_start       : 初期開始点（detect_bpsk_startの結果）
    %   samples_per_symbol  : 1シンボルのサンプル数
    %   num_symbols         : シンボル数
    % 出力:
    %   optimized_start     : 最適化された開始点
    %   best_snr            : その時のSNR推定値
    
    % === パラメータ ===
    search_range = round(samples_per_symbol * 0.25);  % ±25%探索
    search_step = 2;  % 2サンプルずつ探索
    
    best_snr = -inf;
    optimized_start = initial_start;
    
    fprintf('\n=== シンボルタイミング最適化 ===\n');
    
    % === 探索ループ ===
    for offset = -search_range:search_step:search_range
        test_start = initial_start + offset;
        
        % 範囲チェック
        test_end = test_start + samples_per_symbol * num_symbols - 1;
        if test_start < 1 || test_end > length(signal)
            continue;
        end
        
        % このタイミングで復調
        snr_estimate = calculate_demod_snr(signal, fs, carrier_freq, test_start, samples_per_symbol, num_symbols);
        
        % 最良点を更新
        if snr_estimate > best_snr
            best_snr = snr_estimate;
            optimized_start = test_start;
        end
    end
    
    % === 結果表示 ===
    timing_shift = optimized_start - initial_start;
    fprintf('初期開始点: %d\n', initial_start);
    fprintf('最適開始点: %d (シフト: %+d サンプル)\n', optimized_start, timing_shift);
    fprintf('推定SNR: %.2f dB\n', best_snr);
    fprintf('シフト量: %.2f シンボル\n', timing_shift / samples_per_symbol);
end

function snr_db = calculate_demod_snr(signal, fs, carrier_freq, start_idx, samples_per_symbol, num_symbols)
    % === 復調後のSNR推定 ===
    % 原理: コヒーレント積分後の信号電力 vs ノイズ電力
    
    signal_powers = zeros(1, num_symbols);
    
    % 各シンボルの信号電力を計算
    for i = 1:num_symbols
        idx_start = start_idx + (i-1)*samples_per_symbol;
        idx_end = idx_start + samples_per_symbol - 1;
        
        if idx_end > length(signal)
            snr_db = -inf;
            return;
        end
        
        x = signal(idx_start:idx_end);
        t = (0:(samples_per_symbol-1)) / fs;
        
        % コヒーレント積分
        ref_I = cos(2*pi*carrier_freq*t);
        ref_Q = -sin(2*pi*carrier_freq*t);
        
        I = sum(x .* ref_I(:));
        Q = sum(x .* ref_Q(:));
        
        signal_powers(i) = I^2 + Q^2;
    end
    
    % === SNR推定 ===
    % 仮定: 信号はある程度一定、ばらつきはノイズ由来
    mean_power = mean(signal_powers);
    noise_power = var(signal_powers);  % 分散をノイズと見なす
    
    if noise_power < 1e-10
        snr_db = inf;
    else
        snr_db = 10 * log10(mean_power / noise_power);
    end
end