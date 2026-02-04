% === ファイル一覧取得 ===
files = dir('ask_08_10m_*.csv');
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
    carrier_freq = 200e3;
    num_bits = 9;
    bit_rate = 50e3;
    samples_per_bit = round(fs / bit_rate);

    % ====スペクトル減算====
    [ch1, ch2, do_ss, noise_psd_ch1, alpha]  = apply_spectral_subtraction(ch1, ch2, fs);
    
    target_freq = 200e3;      % 注目周波数 200 kHz
    [amp200, mean_other, ratio] = analyze_fft(ch1, fs, target_freq);

    % === CH2変調開始点の検出 ===
    start_ch2 = detect_fsk_start(ch2, fs, carrier_freq, samples_per_bit);
    if isnan(start_ch2)
        warning("CH2: 変調開始点が検出できませんでした: %s", filename);
        continue;
    end

    % トリガ終了点の最小値をデコード開始点に設定
    decode_start_idx = start_ch2;
    decode_end_idx = decode_start_idx + samples_per_bit * num_bits - 1;

    % === 結果表示 ===
    fprintf('変調開始点 (CH2): %d\n', start_ch2);
    fprintf('変調終了点: %d\n', decode_end_idx);

    % === 対象波形を切り出し ===
    wave1 = ch1(decode_start_idx : decode_end_idx);
    wave2 = ch2(decode_start_idx : decode_end_idx);

    %バンドパスフィルタ(チャンネル1,2別・毎回ここをチェック)
    
    N = length(wave1);
    % === 周波数軸 ===
    f_axis = (-N/2:N/2-1) * fs / N;

    % === チャンネル1（例: 中心150kHz, ±50kHz） ===
    
    target_freq1 = 200e3;
    band_width1 =  200e3;
    mask1 = (abs(f_axis) >= (target_freq1 - band_width1)) & ...
            (abs(f_axis) <= (target_freq1 + band_width1));

    % === FFTしてフィルタリング ===
    fft_wave1 = fft(wave1);
    fft_wave1_shifted = fftshift(fft_wave1);
    filtered_fft1 = fft_wave1_shifted .* mask1.';
    filtered_fft1 = ifftshift(filtered_fft1);
    filtered_wave1 = real(ifft(filtered_fft1));

    % === ASK復調 ===
    bitstream1 = decode_ask_kmeans(wave1, samples_per_bit, num_bits);
    bitstream2 = decode_ask_kmeans(wave2, samples_per_bit, num_bits);

    % === 通常ビット列とBER ===
    fprintf('\n--- CH1 ビット列 ---\n');
    disp(bitstream1);
    fprintf('\n--- CH2 ビット列 ---\n');
    disp(bitstream2);
    bit_errors = sum(bitstream1 ~= bitstream2);
    ber = bit_errors / num_bits;
    fprintf('\nBER (非回転): %.2f%% (%d / %d)\n', ber * 100, bit_errors, num_bits);

    shift = 1;

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

% === ASKデコード関数（k-meansしきい値自動推定、安定版）===
function bits = decode_ask_kmeans(signal, samples_per_bit, num_bits)
    bits = zeros(1, num_bits);
    energies = zeros(1, num_bits);

    % 各ビット区間のエネルギーを計算
    for i = 1:num_bits
        segment = signal((i-1)*samples_per_bit + 1 : i*samples_per_bit);
        segment = segment - mean(segment);  % DCオフセット除去
        energies(i) = sum(segment .^ 2);     % エネルギー計算
    end

    % --- k-meansクラスタリング（2クラスタ）でしきい値を推定 ---
    rng(0);  % 乱数シード固定 → 結果を毎回同じにする
    opts = statset('MaxIter',1000);  % 収束改善のオプション
    [labels, centers] = kmeans(energies', 2, 'Replicates',10, 'Options', opts);

    % エネルギーの小さいクラスタを0、大きい方を1とする
    [~, low_cluster] = min(centers);
    bits = (labels ~= low_cluster)';  % 0 or 1 のビット列に変換（行ベクトル）
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

function start_idx = detect_fsk_start(signal, fs, carrier_freq, samples_per_bit)
    trigger_freq = 300e3;
    threshold_Ask = 0.3;
    step = round(samples_per_bit/4);

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
        time_values = [time_values, i/fs];  % 時間軸に変換（サンプル数から時間に変換）

        if ~found_trigger
            % トリガ周波数が最初に現れたらフラグを立てる
            if p_trigger > p0 * 10
                found_trigger = true;
            end
        elseif found_trigger && ~trigger_ended
            % トリガが終わったら次の段階へ
            if p_trigger < max(p_trigger_values)
                trigger_ended = true;
            end
        elseif trigger_ended
            % トリガが終わったあと、f0 または f1 が現れたら開始とみなす
            if p0 > threshold_Ask 
                % ★ 補正：少し前にずらす（最大でもsamples_per_bit/2だけ）
                corrected_idx = i + floor(samples_per_bit / 12);
                %corrected_idx = i;
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

% ==== スペクトル減算 ====
function [wave1_out, wave2_out, do_spectral_subtraction, noise_psd_ch1, alpha] = apply_spectral_subtraction(wave1, wave2, fs)
    % ノイズファイルを読み込み、CH1にスペクトル減算を行う
    % 入力:
    %   wave1, wave2 : 元波形
    %   fs           : 波形のサンプリング周波数
    % 出力:
    %   wave1_out, wave2_out : ノイズ除去後の波形
    %   do_spectral_subtraction : 処理が可能かどうか
    %   noise_psd_ch1 : 推定ノイズPSD
    %   alpha : 減算係数

    % === ノイズPSD推定モードを選択 ===
    % "average" : 平均値で計算
    % "max"     : 各周波数の最大値を使用
    % "random"  : ランダムに1つのノイズファイルを使用
    method = "random";   % ←ここを "average" や "random" に変更して比較

    % === ノイズファイル読み込み ===
    noise_file = dir('noise_l10_*.csv');
    alpha = 1.0; % デフォルト係数

    if isempty(noise_file)
        fprintf('ノイズファイルが見つかりません。スペクトル減算はスキップします。\n');
        wave1_out = wave1;
        wave2_out = wave2;
        do_spectral_subtraction = false;
        noise_psd_ch1 = [];
        return;
    end

    fprintf('合計 %d 個のノイズファイルを使用します（モード: %s）\n', numel(noise_file), method);

    % === PSDの計算 ===
    switch method
        case "average"
            noise_psd_sum = [];
            for k = 1:numel(noise_file)
                fname = noise_file(k).name;
                data = readmatrix(fname);
                ch1_noi = data(:,2);

                [S_noi, ~, ~] = stft(ch1_noi, fs, ...
                    'Window', hamming(128,'periodic'), ...
                    'OverlapLength', 64, ...
                    'FFTLength', 256);

                psd_this = mean(abs(S_noi).^2, 2);
                if isempty(noise_psd_sum)
                    noise_psd_sum = psd_this;
                else
                    noise_psd_sum = noise_psd_sum + psd_this;
                end
            end
            noise_psd_ch1 = noise_psd_sum / numel(noise_file);
            fprintf('平均ノイズPSDを計算しました。\n');

        case "max"
            noise_psd_max = [];
            for k = 1:numel(noise_file)
                fname = noise_file(k).name;
                data = readmatrix(fname);
                ch1_noi = data(:,2);

                [S_noi, ~, ~] = stft(ch1_noi, fs, ...
                    'Window', hamming(128,'periodic'), ...
                    'OverlapLength', 64, ...
                    'FFTLength', 256);

                psd_this = mean(abs(S_noi).^2, 2);
                if isempty(noise_psd_max)
                    noise_psd_max = psd_this;
                else
                    noise_psd_max = max(noise_psd_max, psd_this);
                end
            end
            noise_psd_ch1 = noise_psd_max;
            fprintf('全ノイズファイルの最大PSDを計算しました。\n');

        case "random"
            idx = randi(numel(noise_file));
            fname = noise_file(idx).name;
            data = readmatrix(fname);
            ch1_noi = data(:,2);

            [S_noi, ~, ~] = stft(ch1_noi, fs, ...
                'Window', hamming(128,'periodic'), ...
                'OverlapLength', 64, ...
                'FFTLength', 256);

            noise_psd_ch1 = mean(abs(S_noi).^2, 2);
            fprintf('ノイズファイル "%s" をランダムに選択して使用しました。\n', fname);

        otherwise
            error('method の設定が不正です。 "average", "max", "random" のいずれかを指定してください。');
    end

    % === CH1のみスペクトル減算 ===
    wave1_out = spectral_subtract_ch1(wave1, fs, noise_psd_ch1, alpha);
    wave2_out = wave2; % CH2はそのまま
    do_spectral_subtraction = true;
end

% ==== CH1スペクトル減算（短い信号対応版） ====
function wave1_denoised = spectral_subtract_ch1(wave1, fs, noise_psd_ch1, alpha)
    win_len = min(128, length(wave1));
    overlap = round(win_len/2);
    nfft = max(256, 2^nextpow2(win_len));
    win = hamming(win_len,'periodic');

    [S1, ~, ~] = stft(wave1, fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
    mag1 = abs(S1);
    phase1 = angle(S1);

    % ノイズPSDの補間（STFTサイズに合わせる）
    if length(noise_psd_ch1) ~= size(S1,1)
        f_noi = linspace(0, fs/2, length(noise_psd_ch1)).';
        f_sig = linspace(0, fs/2, size(S1,1)).';
        noise_psd_ch1_interp = interp1(f_noi, noise_psd_ch1, f_sig, 'linear', 'extrap');
    else
        noise_psd_ch1_interp = noise_psd_ch1;
    end

    % スペクトル減算
    mag1_denoised = sqrt(max(mag1.^2 - alpha * noise_psd_ch1_interp, 0));
    S1_denoised = mag1_denoised .* exp(1j * phase1);

    x1_denoised = istft(S1_denoised, fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);

    % 元の信号長に合わせる
    x1_denoised = x1_denoised(1:min(length(x1_denoised), length(wave1)));
    if length(x1_denoised) < length(wave1)
        x1_denoised = [x1_denoised; zeros(length(wave1)-length(x1_denoised),1)];
    end

    wave1_denoised = real(x1_denoised);
end

function [amp_target, mean_other, ratio] = analyze_fft(signal, fs, target_freq)
% === FFT解析関数 ===
% 指定した信号に対してFFTを実行し、
% 注目周波数（例: 200kHz）成分の振幅と、
% それ以外の平均振幅、および比率を計算・表示する。
%
% 使用例:
%   [amp200, mean_other, ratio] = analyze_fft(wave1, 2e6, 200e3);
%
% 引数:
%   signal       : 時間波形データ（列ベクトルでもOK）
%   fs           : サンプリング周波数 [Hz]
%   target_freq  : 注目する周波数 [Hz]
%
% 戻り値:
%   amp_target   : 注目周波数の振幅 [V]
%   mean_other   : その他の平均振幅 [V]
%   ratio        : amp_target / mean_other の比率

    % 信号を列ベクトルに変換
    signal = signal(:);
    N = length(signal);

    % === FFT計算 ===
    Y = fft(signal);
    f = (0:N-1)*(fs/N);     % 周波数軸

    % === 振幅スペクトル（片側）===
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f1 = f(1:N/2+1);

    % === 注目周波数の成分 ===
    [~, idx_target] = min(abs(f1 - target_freq));
    amp_target = P1(idx_target);

    % === 注目周波数以外の平均値 ===
    P1_exclude = P1;
    P1_exclude(idx_target) = [];
    mean_other = mean(P1_exclude);

    % === 比率計算 ===
    ratio = amp_target / mean_other;

    % === 結果表示 ===
    fprintf('--- FFT解析結果 ---\n');
    fprintf('サンプリング周波数: %.3f MHz\n', fs/1e6);
    fprintf('注目周波数: %.3f kHz\n', target_freq/1e3);
    fprintf('%.3f kHz成分振幅: %.6f V\n', target_freq/1e3, amp_target);
    fprintf('%.3f kHz以外の平均振幅: %.6f V\n', target_freq/1e3, mean_other);
    fprintf('比率 (%.0fkHz / 平均): %.2f 倍\n\n', target_freq/1e3, ratio);
end
