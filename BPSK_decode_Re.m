% === ファイル一覧取得 ===
files = dir('gr_bp_3_*.csv');
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
    num_bits = 9;
    bit_rate = 50e3;
    samples_per_bit = round(fs / bit_rate);

    % === CH1変調開始点の検出 ===
    
    start_ch1 = detect_bpsk_start(ch1, fs, carrier_freq, samples_per_bit);
    if isnan(start_ch1)
        warning("CH1: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    % === CH2変調開始点の検出 ===
    start_ch2 = detect_bpsk_start(ch2, fs, carrier_freq, samples_per_bit);
    if isnan(start_ch2)
        warning("CH2: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    % トリガ終了点の最小値をデコード開始点に設定
    decode_start_idx = start_ch2;
    decode_end_idx = decode_start_idx + samples_per_bit * num_bits - 1;

    % === 結果表示 ===
  
    fprintf('変調開始点 (CH1): %d\n', start_ch1);
    fprintf('変調開始点 (CH2): %d\n', start_ch2);
    fprintf('変調終了点: %d\n', decode_end_idx);
    
    % 波形切り出し
    wave1 = ch1(decode_start_idx : decode_end_idx);
    wave2 = ch2(decode_start_idx : decode_end_idx);

    % === BPSK復調 ===
    bitstream1 = decode_bpsk(wave1, fs,  carrier_freq, num_bits, samples_per_bit);
    bitstream2 = decode_bpsk(wave2, fs,  carrier_freq, num_bits, samples_per_bit);

    % --- プリアンブル（最初の2ビット）で反転を判定 ---
    preamble1 = bitstream1(1:2);
    preamble2 = bitstream2(1:2);

    % 一致しない場合（反転してると仮定）
    if sum(preamble1 ~= preamble2) >= 1
        fprintf("プリアンブルのビットが反転していたため、bitstream1 を反転します\n");
        bitstream1 = 1 - bitstream1;  % ビットを反転
    end

    % --- BER計算（情報ビット部分のみで） ---
    info_bits1 = bitstream1;
    info_bits2 = bitstream2;
    fprintf('\n--- CH1 ビット列 ---\n');
    disp(info_bits1);
    fprintf('\n--- CH2 ビット列 ---\n');
    disp(info_bits2);
    bit_errors = sum(info_bits1 ~= info_bits2);
    ber = bit_errors / length(info_bits2);
    fprintf('\nBER : %.2f%% (%d / %d)\n', ber * 100, bit_errors,length(info_bits2));

    shift = 1;

     num_bits);

    % === 全体波形と範囲表示 ===
    
    fig2 = figure;
    subplot(2,1,1);
    plot(time, ch1, 'b'); hold on;
    area(time(decode_start_idx:decode_end_idx), ch1(decode_start_idx:decode_end_idx), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    title('CH1 全体波形とデコード範囲',filename);
    xlabel('Time [s]'); ylabel('Amplitude');

    subplot(2,1,2);
    plot(time, ch2, 'b'); hold on;
    area(time(decode_start_idx:decode_end_idx), ch2(decode_start_idx:decode_end_idx), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    title('CH2 全体波形とデコード範囲',filename);
    xlabel('Time [s]'); ylabel('Amplitude');

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

% ✅ 全処理後に一括で平均BER追加（NaN除外）
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

% 結果を書き込む
write_to_excel(excel_filename, sheet_name, decode_results);
fprintf('\nBER : %.2f%%\n', average_ber * 100);
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
        header = {'デコード日時','デコードファイル' ,'出力ビット', '入力ビット', 'ビット誤り率', '備考欄'};
        writecell(header, filename, 'Sheet', sheet_name, 'Range', 'A1');
    end
    
    % データの書き込み
    writecell(data, filename, 'Sheet', sheet_name, 'Range', sprintf('A%d', start_row));
end

function corrected_ch1 = align_phase_by_preamble(ch1, ch2, fs, carrier_freq, samples_per_bit)
    % 2ビットのプリアンブル用信号区間を取得
    N = 2 * samples_per_bit;

    if length(ch1) < N || length(ch2) < N
        error("波形が短すぎてプリアンブルが検出できません");
    end

    % プリアンブル区間を抽出
    x1 = ch1(1:N);
    x2 = ch2(1:N);

    % 時間ベクトル
    t = (0:(N-1)) / fs;

    % 複素信号への変換（同期検波に相当）
    ref_cos = cos(2 * pi * carrier_freq * t);
    ref_sin = sin(2 * pi * carrier_freq * t);

    I1 = sum(x1 .* ref_cos);
    Q1 = sum(x1 .* ref_sin);
    angle1 = atan2(Q1, I1);

    I2 = sum(x2 .* ref_cos);
    Q2 = sum(x2 .* ref_sin);
    angle2 = atan2(Q2, I2);

    delta_phi = wrapToPi(angle2 - angle1);  % [-π, π]

    %fprintf("補正前の位相差: %.2f度\n", rad2deg(delta_phi));

    % プリアンブルが反転している場合（約±180°ずれてる）
    if abs(delta_phi) > pi/2
        delta_phi = wrapToPi(delta_phi + pi);
        fprintf("プリアンブルが反転していたため、180度回転を追加\n");
    end

    % === ch1全体に複素回転を適用 ===
    t_full = (0:length(ch1)-1)' / fs;
    z = ch1 .* exp(-1j * 2 * pi * carrier_freq * t_full);  % 同期検波
    z_corr = z .* exp(1j * delta_phi);  % 位相補正
    corrected_ch1 = real(z_corr .* exp(1j * 2 * pi * carrier_freq * t_full));  % 元の周波数に復調
end

% === BPSKデコード関数===
function bit_array = decode_bpsk(signal, fs,  carrier_freq, num_bits, samples_per_bit)
    bit_array = zeros(1, num_bits);  % 出力配列の初期化

    for i = 1:num_bits
        idx_start = round((i-1)*samples_per_bit) + 1;
        idx_end   = round(i*samples_per_bit);

        if idx_end > length(signal)
            warning("信号が短いため、%dビット目は処理されません", i);
            break;
        end

        x = signal(idx_start:idx_end);
        t = (0:(samples_per_bit - 1)) / fs;  % ← 長さを明示的にsamples_per_bitに固定
        ref = sin(2 * pi * carrier_freq * t);   % 搬送波

        % 念のため両方のサイズを一致させる（列ベクトルでも可）
        x = x(:);         % 縦ベクトルに変換
        ref = ref(:);     % 縦ベクトルに変換

        corr = sum(x .* ref);        % 相関計算（スカラー）
        bit_array(i) = corr < 1;  % ← 反転判定    % BPSK復調（1:反転, 0:正位相）

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

function start_idx = detect_bpsk_start(signal, fs, carrier_freq,samples_per_bit)
    trigger_freq = 300e3;
    threshold_fsk = 0.3;
    step = round(samples_per_bit/8);

    found_trigger = false;
    trigger_ended = false;
    i = 1;

    % トリガパワーを格納するための配列
    p_trigger_values = [];
    time_values = [];  % 時間軸（またはサンプルインデックス）

    while i <= length(signal) - samples_per_bit
        segment = signal(i : i + samples_per_bit - 1);
        segment = segment - mean(segment);  % DC成分を取り除く
        segment = segment / max(abs(segment));  % 正規化

        % トリガ周波数と変調周波数のパワーを計算
        p_trigger = goertzel_exact(segment, fs, trigger_freq);
        p0 = goertzel_exact(segment, fs, carrier_freq);


        % 初回はそのまま格納、2回目以降は最大値で正規化
        if isempty(p_trigger_values)
            p_trigger_values = [p_trigger_values, p_trigger];  % 初回はそのまま
        else
            p_trigger_values = [p_trigger_values, p_trigger / max(p_trigger_values)];  % 2回目以降は正規化
        end

        if ~found_trigger
            % トリガ周波数が最初に現れたらフラグを立てる
            if p_trigger > p0 * 3
                found_trigger = true;
            end
        elseif found_trigger && ~trigger_ended
            % トリガが終わったら次の段階へ
            if p_trigger < 0.3 * max(p_trigger_values)  % 比率で判定
                trigger_ended = true;
            end
        elseif trigger_ended
            if p0 > threshold_fsk
                % ★ 補正：少し前にずらす（最大でもsamples_per_bit/2だけ）
                corrected_idx = i + floor(samples_per_bit / 12);
                % corrected_idx = i;
                start_idx = max(1, corrected_idx);  % マイナス防止
                break;
                
            end
        end

        i = i + step;
    end

    % トリガ開始点が検出できなかった場合
    if ~found_trigger
        start_idx = NaN;
    end
end

function waveform_out = freq_domain_smooth(waveform, samples_per_bit, num_bits, fs, fc, bw)
% waveform: 入力波形（列ベクトル）
% samples_per_bit: 1ビットあたりサンプル数
% num_bits: ビット数
% fs: サンプリング周波数 (Hz)
% fc: 信号周波数 (Hz)
% bw: 信号帯域幅 (Hz, ±)

total_samples = samples_per_bit * num_bits;
segments = reshape(waveform(1:total_samples), samples_per_bit, num_bits);
fft_segments = zeros(size(segments));

% 周波数軸
f = (0:samples_per_bit-1)' * fs / samples_per_bit;

% 各ビットのFFT
for k = 1:num_bits
    fft_segments(:,k) = fft(segments(:,k));
end

% 信号帯域マスク
mask = (f >= (fc-bw)) & (f <= (fc+bw));
mask = mask | (f >= (fs-(fc+bw)) & f <= (fs-(fc-bw))); % 対称周波数も保護

% 信号帯域以外をビット平均に置き換え
fft_smoothed = fft_segments;
for i = 1:samples_per_bit
    if ~mask(i)
        fft_smoothed(i,:) = mean(fft_segments(i,:));
    end
end

% IFFT で時間領域に戻す
segments_out = ifft(fft_smoothed, 'symmetric');

% 結合
waveform_out = reshape(segments_out, [], 1);
end
