% === ファイル一覧取得 ===
files = dir('qpsk_10m_*.csv');
num_files = length(files);      % 読み込んだファイルの数を取得
disp(num_files); % ファイルの数を表示
decode_results = {};  
row_idx = 1;          % 行番号を管理

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
    carrier_freq = 200e3;
    num_symbols = 4;
    num_bits = num_symbols * 2 ;
    bit_rate = 50e3;
    samples_per_bit = round(fs / bit_rate);

    % === CH2変調開始点の検出 ===
    start_ch2 = detect_bpsk_start(ch2, fs, carrier_freq, samples_per_bit);
    if isnan(start_ch2)
        warning("CH2: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    % トリガ終了点の最小値をデコード開始点に設定
    decode_start_idx =  start_ch2;
    decode_end_idx = decode_start_idx + samples_per_bit * num_bits - 1;

    % === 結果表示 ===
    %{
    fprintf('変調開始点 (CH1): %d\n', start_ch1);
    %}
    fprintf('変調開始点 (CH2): %d\n', start_ch2);
    fprintf('変調終了点: %d\n', decode_end_idx);

    % CH2はそのまま切り出し
    wave1 = ch1(decode_start_idx : decode_end_idx);
    wave2 = ch2(decode_start_idx : decode_end_idx);

    bits_per_symbol = 1;
    num_bits = 8;
    fc = 200e3;
    bw = 1e3;
    wave1_true = freq_domain_smooth_qam(wave1, samples_per_bit, num_bits, fs, fc, bw, bits_per_symbol);

    % QPSK復調
    bitstream1 = decode_qpsk(wave1, fs, carrier_freq, num_symbols, samples_per_bit);
    bitstream2 = decode_qpsk(wave2, fs, carrier_freq, num_symbols, samples_per_bit);

    % --- 受信プリアンブル抽出 ---
    preamble1 = bitstream1(1:2);
    preamble2 = bitstream2(1:2);

    % --- 各プリアンブルをシンボル番号（0〜3）に変換 ---
    sym1 = preamble1(1)*2 + preamble1(2);  % CH1のシンボル
    sym2 = preamble2(1)*2 + preamble2(2);  % CH2のシンボル

    % --- CH1をCH2の位相に合わせる補正量（90度単位） ---
    phase_offset = mod(sym1 - sym2, 4);

    fprintf("CH1の位相を -%d×90° 回転してCH2に揃えます\n", phase_offset);

    % === 位相補正処理（CH1を補正）===
    bitstream1_corrected = zeros(size(bitstream1));

    for i = 1:2:length(bitstream1)
        % CH1の2ビット → シンボル番号
        sym = bitstream1(i)*2 + bitstream1(i+1);
    
        % 位相補正（CH2側に揃えるため逆回転）
        sym_corr = mod(sym - phase_offset, 4);
    
        % シンボル番号 → ビット列に戻す
        bitstream1_corrected(i)   = floor(sym_corr / 2);
        bitstream1_corrected(i+1) = mod(sym_corr, 2);
    end

    % === ビット列 === %
    info_bits1 = bitstream1_corrected;
    info_bits2 = bitstream2;
    fprintf('\n--- CH1 ビット列 ---\n');
    disp(info_bits1);
    fprintf('\n--- CH2 ビット列 ---\n');
    disp(info_bits2);
    % --- ビット列 → シンボル（2ビットずつ → 0〜3） ---
    % 長さが奇数だった場合に備えて切り捨て
    len = floor(min(length(info_bits1), length(info_bits2)) / 2) * 2;
    bit_pairs1 = reshape(info_bits1(1:len), 2, []).';  % 各行が [b1 b2]
    bit_pairs2 = reshape(info_bits2(1:len), 2, []).';

    symbols1 = bit_pairs1(:,1)*2 + bit_pairs1(:,2);  % シンボル列（0〜3）
    symbols2 = bit_pairs2(:,1)*2 + bit_pairs2(:,2);

    fprintf('\n--- CH1 シンボル列 ---\n');
    disp(symbols1.');
    fprintf('\n--- CH2 シンボル列 ---\n');
    disp(symbols2.');


    % --- BER計算（情報ビット部分のみで） ---
    bit_errors = sum(info_bits1 ~= info_bits2);
    ber = bit_errors / length(info_bits2);
    fprintf('\nBER : %.2f%% (%d / %d)\n', ber * 100, bit_errors,length(info_bits2));
  
    shift = 1;

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
    decode_results{row_idx, 3} = sprintf('%d', info_bits1);
    decode_results{row_idx, 4} = sprintf('%d', info_bits2);
    decode_results{row_idx, 5} = ber;
    decode_results{row_idx, 6} = '';  % 備考空白
    row_idx = row_idx + 1;
  
end

% 全処理後に一括で平均BER追加（NaN除外）
valid_bers = [decode_results{:, 5}];
valid_bers = valid_bers(~isnan(valid_bers));  % NaNを除く
average_ber = mean(valid_bers);

% 平均行を最後に追加
decode_results{row_idx, 1} = datetime('now');
decode_results{row_idx, 2} = '平均';
decode_results{row_idx, 3} = '';
decode_results{row_idx, 4} = '';
decode_results{row_idx, 5} = average_ber;
decode_results{row_idx, 6} = '';
    

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
fprintf('BER = %.2f%%\n',average_ber*100);
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

function corrected_ch1 = align_phase_by_preamble(ch1, ch2, fs, carrier_freq, samples_per_bit)

    % === プリアンブル部分だけ取り出し（1シンボル＝2ビット） ===
    N = samples_per_bit * 2;
    x1 = ch1(1:N);
    x2 = ch2(1:N);
    t = (0:N-1)' / fs;

    % === 同期検波して複素包絡を取得 ===
    z1 = x1 .* exp(-1j * 2 * pi * carrier_freq * t);
    z2 = x2 .* exp(-1j * 2 * pi * carrier_freq * t);

    % === 平均複素ベクトルの位相差を計算（よりロバストな方法） ===
    delta_phi = angle(mean(z2) * conj(mean(z1)));  % z2/z1 の角度
    fprintf("補正前のQPSK位相差: %.2f度\n", rad2deg(delta_phi));

    % === CH1全体に補正を適用 ===
    t_full = (0:length(ch1)-1)' / fs;
    z = ch1 .* exp(-1j * 2 * pi * carrier_freq * t_full);     % 複素包絡
    z_corr = z * exp(1j * delta_phi);                         % 位相補正
    corrected_ch1 = real(z_corr .* exp(1j * 2 * pi * carrier_freq * t_full)); % 実波形へ戻す
end

function bit_array = decode_qpsk(signal, fs, carrier_freq, num_symbols,samples_per_bit)
    bit_array = zeros(1, num_symbols*2);
    samples_per_symbol = samples_per_bit * 2;

    I_vals = zeros(1, num_symbols);
    Q_vals = zeros(1, num_symbols);

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

        x = x(:);       % 80×1
        ref_I = ref_I(:); % 80×1 に変更
        ref_Q = ref_Q(:); % 80×1 に変更

        I = sum(x .* ref_I);
        Q = sum(x .* ref_Q);

        I_vals(i) = I;
        Q_vals(i) = Q;

        if I > 0 && Q > 0
            bits = [0 0];
        elseif I < 0 && Q > 0
            bits = [0 1];
        elseif I < 0 && Q < 0
            bits = [1 0];
        else
            bits = [1 1];
        end
        bit_array(2*i-1:2*i) = bits;
    end
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
                corrected_idx = i - floor(samples_per_bit / 8);
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
