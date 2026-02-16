import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from 'recharts';
import { Calculator, Activity, AlertCircle, TrendingUp } from 'lucide-react';

/**
 * TAZ/PIPC投与設計支援アプリケーション
 * 1コンパートメントモデルとモンテカルロシミュレーションによるPTA算出
 */

// ===============================
// 定数定義
// ===============================
const PIPC_PHARMACOKINETICS = {
  Vd: 13.4,           // 分布容積 (L)
  baseCL: 8.74,       // 基本クリアランス (L/hr)
  ccrCoefficient: 0.0472,  // CCr係数
  ccrReference: 82.6, // 参照CCr値
  IIV_CV: 0.189,      // 個人間誤差（変動係数）18.9%
  freeFormFraction: 0.8  // 遊離形濃度係数（タンパク結合率20%）
};

const DEFAULT_MIC_VALUES = [1, 2, 4, 8, 16, 32, 64, 128];
const MONTE_CARLO_ITERATIONS = 10000;

// ===============================
// ユーティリティ関数
// ===============================

/**
 * 理想体重（IBW）の計算
 * @param {number} height - 身長 (cm)
 * @param {string} gender - 性別 ('male' or 'female')
 * @returns {number} IBW (kg)
 */
function calculateIBW(height, gender) {
  if (gender === 'male') {
    return 50 + 0.9 * (height - 152.4);
  } else {
    return 45.5 + 0.9 * (height - 152.4);
  }
}

/**
 * BMIの計算
 * @param {number} weight - 体重 (kg)
 * @param {number} height - 身長 (cm)
 * @returns {number} BMI
 */
function calculateBMI(weight, height) {
  const heightM = height / 100;
  return weight / (heightM * heightM);
}

/**
 * クレアチニンクリアランス（CCr）の計算
 * BMI > 30 または 実体重 > IBW * 1.2 の場合はSalazar-Corcoran式を使用
 * それ以外はCockcroft-Gault式を使用
 * 
 * @param {number} age - 年齢
 * @param {number} weight - 体重 (kg)
 * @param {number} height - 身長 (cm)
 * @param {number} scr - 血清クレアチニン (mg/dL)
 * @param {string} gender - 性別 ('male' or 'female')
 * @returns {object} { ccr, method }
 */
function calculateCCr(age, weight, height, scr, gender) {
  if (scr <= 0) {
    throw new Error('血清クレアチニン値が不正です');
  }

  const ibw = calculateIBW(height, gender);
  const bmi = calculateBMI(weight, height);
  
  let ccr, method;

  // Salazar-Corcoran式またはCockcroft-Gault式の選択
  if (bmi > 30 || weight > ibw * 1.2) {
    // Salazar-Corcoran式
    method = 'Salazar-Corcoran';
    if (gender === 'male') {
      ccr = ((137 - age) * ((0.285 * weight) + (12.1 * (height / 100) ** 2))) / (51 * scr);
    } else {
      ccr = ((146 - age) * ((0.287 * weight) + (9.74 * (height / 100) ** 2))) / (60 * scr);
    }
  } else {
    // Cockcroft-Gault式
    method = 'Cockcroft-Gault';
    ccr = ((140 - age) * weight) / (72 * scr);
    if (gender === 'female') {
      ccr *= 0.85;
    }
  }

  return { ccr: Math.max(0, ccr), method };
}

/**
 * 個別のクリアランス（CL）の計算
 * @param {number} ccr - クレアチニンクリアランス (mL/min)
 * @returns {number} CL (L/hr)
 */
function calculateCL(ccr) {
  return PIPC_PHARMACOKINETICS.baseCL + 
         PIPC_PHARMACOKINETICS.ccrCoefficient * (ccr - PIPC_PHARMACOKINETICS.ccrReference);
}

/**
 * 対数正規分布に従う乱数生成（Box-Muller法）
 * @param {number} mean - 平均値
 * @param {number} cv - 変動係数（CV）
 * @returns {number} 対数正規分布に従う乱数
 */
function generateLogNormal(mean, cv) {
  // Box-Muller法で標準正規分布の乱数を生成
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
  
  // 対数正規分布のパラメータ計算
  const sigma = Math.sqrt(Math.log(1 + cv * cv));
  const mu = Math.log(mean) - 0.5 * sigma * sigma;
  
  return Math.exp(mu + sigma * z);
}

/**
 * 1コンパートメントモデルによる血中濃度計算（定常状態）
 * 定常状態における濃度推移を2サイクル分計算
 * 
 * @param {number} dose - 1回投与量 (g)
 * @param {number} infusionTime - 点滴時間 (hr)
 * @param {number} interval - 投与間隔 (hr)
 * @param {number} cl - クリアランス (L/hr)
 * @param {number} vd - 分布容積 (L)
 * @returns {array} 時間と濃度の配列
 */
function calculateConcentrationProfile(dose, infusionTime, interval, cl, vd) {
  const doseG = dose; // グラム単位
  const doseMg = doseG * 1000; // mg単位に変換
  const k = cl / vd; // 消失速度定数 (1/hr)
  const infusionRate = doseMg / infusionTime; // mg/hr
  
  // 定常状態のCmin（トラフ値）とCmax（ピーク値）を先に計算
  // Cmin = (Dose/Tinf) * (1/CL) * (1 - exp(-k*Tinf)) * exp(-k*(τ-Tinf)) / (1 - exp(-k*τ))
  const cMin = (infusionRate / cl) * 
               (1 - Math.exp(-k * infusionTime)) * 
               Math.exp(-k * (interval - infusionTime)) / 
               (1 - Math.exp(-k * interval));
  
  // Cmax = Cmin * exp(k*(τ-Tinf)) + (Dose/Tinf) * (1/CL) * (1 - exp(-k*Tinf))
  const cMax = cMin * Math.exp(k * (interval - infusionTime)) + 
               (infusionRate / cl) * (1 - Math.exp(-k * infusionTime));
  
  const points = [];
  const timeStep = 0.1; // 0.1時間刻み
  
  // 2サイクル分の濃度推移を計算（定常状態を視覚化）
  for (let cycle = 0; cycle < 2; cycle++) {
    const cycleOffset = cycle * interval;
    
    // 点滴中の濃度上昇（0 ≤ t ≤ Tinf）
    // C(t) = Cmin * exp(k*t) + (Dose/Tinf) * (1/CL) * (1 - exp(-k*t))
    for (let t = 0; t <= infusionTime; t += timeStep) {
      const concentration = cMin * Math.exp(k * t) + 
                           (infusionRate / cl) * (1 - Math.exp(-k * t));
      points.push({ 
        time: cycleOffset + t, 
        concentration 
      });
    }
    
    // 点滴終了後の濃度減衰（Tinf < t ≤ τ）
    // C(t) = Cmax * exp(-k*(t-Tinf))
    for (let t = infusionTime + timeStep; t <= interval; t += timeStep) {
      const concentration = cMax * Math.exp(-k * (t - infusionTime));
      points.push({ 
        time: cycleOffset + t, 
        concentration 
      });
    }
  }
  
  return points;
}

/**
 * %fT > MICの計算（定常状態の1サイクル）
 * @param {array} concentrationProfile - 濃度プロファイル（2サイクル分含む）
 * @param {number} mic - MIC値 (μg/mL)
 * @param {number} interval - 投与間隔 (hr)
 * @returns {number} %fT > MIC
 */
function calculatePercentTimeMIC(concentrationProfile, mic, interval) {
  const freeFormFraction = PIPC_PHARMACOKINETICS.freeFormFraction;
  let timeAboveMIC = 0;
  
  // 定常状態の1サイクル分（0 ≤ t ≤ interval）のみを使用
  const oneCycleData = concentrationProfile.filter(point => point.time <= interval);
  
  for (let i = 0; i < oneCycleData.length - 1; i++) {
    const currentConc = oneCycleData[i].concentration * freeFormFraction;
    const nextConc = oneCycleData[i + 1].concentration * freeFormFraction;
    const dt = oneCycleData[i + 1].time - oneCycleData[i].time;
    
    if (currentConc >= mic && nextConc >= mic) {
      // 両方がMIC以上の場合、全時間を加算
      timeAboveMIC += dt;
    } else if (currentConc >= mic || nextConc >= mic) {
      // どちらか一方がMIC以上の場合、線形補間で正確な交差時間を計算
      const ratio = currentConc >= mic ? 
                   (currentConc - mic) / (currentConc - nextConc) :
                   (mic - currentConc) / (nextConc - currentConc);
      timeAboveMIC += dt * ratio;
    }
  }
  
  return (timeAboveMIC / interval) * 100;
}

/**
 * モンテカルロシミュレーションによるPTA計算
 * @param {number} dose - 1回投与量 (g)
 * @param {number} infusionTime - 点滴時間 (hr)
 * @param {number} interval - 投与間隔 (hr)
 * @param {number} meanCL - 平均クリアランス (L/hr)
 * @param {number} targetPercent - 目標%fT > MIC
 * @returns {array} MIC値ごとのPTA
 */
function runMonteCarloSimulation(dose, infusionTime, interval, meanCL, targetPercent) {
  const results = {};
  
  DEFAULT_MIC_VALUES.forEach(mic => {
    results[mic] = 0;
  });
  
  // 10,000人の仮想患者でシミュレーション
  for (let i = 0; i < MONTE_CARLO_ITERATIONS; i++) {
    // 対数正規分布でCLを生成（個人間誤差を考慮）
    const individualCL = generateLogNormal(meanCL, PIPC_PHARMACOKINETICS.IIV_CV);
    
    // 濃度プロファイル計算
    const profile = calculateConcentrationProfile(
      dose, 
      infusionTime, 
      interval, 
      individualCL, 
      PIPC_PHARMACOKINETICS.Vd
    );
    
    // 各MIC値での%fT > MICを計算
    DEFAULT_MIC_VALUES.forEach(mic => {
      const percentTime = calculatePercentTimeMIC(profile, mic, interval);
      if (percentTime >= targetPercent) {
        results[mic]++;
      }
    });
  }
  
  // PTA（%）に変換
  const ptaResults = DEFAULT_MIC_VALUES.map(mic => ({
    mic,
    pta: (results[mic] / MONTE_CARLO_ITERATIONS) * 100
  }));
  
  return ptaResults;
}

// ===============================
// メインコンポーネント
// ===============================
export default function TazPipcDosingApp() {
  // 患者情報の状態管理
  const [patientData, setPatientData] = useState({
    age: 65,
    gender: 'male',
    height: 170,
    weight: 70,
    scr: 1.0
  });

  // 投与設計の状態管理
  const [dosingRegimen, setDosingRegimen] = useState({
    dose: 4.5,          // TAZ/PIPC 4.5g
    infusionTime: 0.5,  // 30分点滴
    interval: 6         // 6時間ごと
  });

  // 目標設定
  const [targetFtMic, setTargetFtMic] = useState(50);

  // 計算結果の状態管理
  const [results, setResults] = useState(null);
  const [isCalculating, setIsCalculating] = useState(false);
  const [error, setError] = useState(null);

  // 入力変更ハンドラ
  const handlePatientChange = (field, value) => {
    setPatientData(prev => ({ ...prev, [field]: value }));
  };

  const handleDosingChange = (field, value) => {
    setDosingRegimen(prev => ({ ...prev, [field]: value }));
  };

  // 計算実行
  const handleCalculate = () => {
    setIsCalculating(true);
    setError(null);

    // 非同期処理でUIをブロックしない
    setTimeout(() => {
      try {
        // バリデーション
        if (patientData.scr <= 0) {
          throw new Error('血清クレアチニン値は0より大きい値を入力してください');
        }
        if (patientData.age <= 0 || patientData.height <= 0 || patientData.weight <= 0) {
          throw new Error('年齢、身長、体重は0より大きい値を入力してください');
        }

        // 腎機能計算
        const { ccr, method } = calculateCCr(
          patientData.age,
          patientData.weight,
          patientData.height,
          patientData.scr,
          patientData.gender
        );

        // クリアランス計算
        const cl = calculateCL(ccr);

        // 濃度プロファイル計算（平均的な患者）
        const concentrationProfile = calculateConcentrationProfile(
          dosingRegimen.dose,
          dosingRegimen.infusionTime,
          dosingRegimen.interval,
          cl,
          PIPC_PHARMACOKINETICS.Vd
        );

        // 定常状態のCminとCmax（最初のサイクルから取得）
        const cMin = concentrationProfile[0].concentration;
        const cMax = Math.max(...concentrationProfile
          .filter(p => p.time <= dosingRegimen.interval)
          .map(p => p.concentration));

        // モンテカルロシミュレーション
        const ptaResults = runMonteCarloSimulation(
          dosingRegimen.dose,
          dosingRegimen.infusionTime,
          dosingRegimen.interval,
          cl,
          targetFtMic
        );

        setResults({
          ccr,
          ccrMethod: method,
          cl,
          cMin,
          cMax,
          concentrationProfile,
          ptaResults
        });
      } catch (err) {
        setError(err.message);
      } finally {
        setIsCalculating(false);
      }
    }, 100);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-blue-50 to-slate-100 p-6">
      <div className="max-w-7xl mx-auto">
        {/* ヘッダー */}
        <div className="mb-8 text-center">
          <div className="flex items-center justify-center gap-3 mb-2">
            <Activity className="w-10 h-10 text-blue-600" />
            <h1 className="text-4xl font-bold bg-gradient-to-r from-blue-600 to-indigo-600 bg-clip-text text-transparent">
              TAZ/PIPC 投与設計支援システム
            </h1>
          </div>
          <p className="text-slate-600 text-lg">
            1コンパートメントモデル × モンテカルロシミュレーション
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* 左側：入力セクション */}
          <div className="space-y-6">
            {/* 患者情報カード */}
            <div className="bg-white rounded-2xl shadow-lg p-6 border border-slate-200">
              <h2 className="text-2xl font-bold text-slate-800 mb-4 flex items-center gap-2">
                <Calculator className="w-6 h-6 text-blue-600" />
                患者情報
              </h2>
              
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-sm font-semibold text-slate-700 mb-2">
                      年齢 (歳)
                    </label>
                    <input
                      type="number"
                      value={patientData.age}
                      onChange={(e) => handlePatientChange('age', parseFloat(e.target.value))}
                      className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                    />
                  </div>
                  
                  <div>
                    <label className="block text-sm font-semibold text-slate-700 mb-2">
                      性別
                    </label>
                    <select
                      value={patientData.gender}
                      onChange={(e) => handlePatientChange('gender', e.target.value)}
                      className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                    >
                      <option value="male">男性</option>
                      <option value="female">女性</option>
                    </select>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-sm font-semibold text-slate-700 mb-2">
                      身長 (cm)
                    </label>
                    <input
                      type="number"
                      value={patientData.height}
                      onChange={(e) => handlePatientChange('height', parseFloat(e.target.value))}
                      className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                    />
                  </div>
                  
                  <div>
                    <label className="block text-sm font-semibold text-slate-700 mb-2">
                      体重 (kg)
                    </label>
                    <input
                      type="number"
                      value={patientData.weight}
                      onChange={(e) => handlePatientChange('weight', parseFloat(e.target.value))}
                      className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                    />
                  </div>
                </div>

                <div>
                  <label className="block text-sm font-semibold text-slate-700 mb-2">
                    血清クレアチニン (mg/dL)
                  </label>
                  <input
                    type="number"
                    step="0.01"
                    value={patientData.scr}
                    onChange={(e) => handlePatientChange('scr', parseFloat(e.target.value))}
                    className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                  />
                </div>
              </div>
            </div>

            {/* 投与設計カード */}
            <div className="bg-white rounded-2xl shadow-lg p-6 border border-slate-200">
              <h2 className="text-2xl font-bold text-slate-800 mb-4 flex items-center gap-2">
                <TrendingUp className="w-6 h-6 text-indigo-600" />
                投与設計
              </h2>
              
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-semibold text-slate-700 mb-2">
                    1回投与量 (g)
                  </label>
                  <input
                    type="number"
                    step="0.5"
                    value={dosingRegimen.dose}
                    onChange={(e) => handleDosingChange('dose', parseFloat(e.target.value))}
                    className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                  />
                </div>

                <div>
                  <label className="block text-sm font-semibold text-slate-700 mb-2">
                    点滴時間 (時間)
                  </label>
                  <input
                    type="number"
                    step="0.25"
                    value={dosingRegimen.infusionTime}
                    onChange={(e) => handleDosingChange('infusionTime', parseFloat(e.target.value))}
                    className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                  />
                </div>

                <div>
                  <label className="block text-sm font-semibold text-slate-700 mb-2">
                    投与間隔 (時間)
                  </label>
                  <input
                    type="number"
                    step="1"
                    value={dosingRegimen.interval}
                    onChange={(e) => handleDosingChange('interval', parseFloat(e.target.value))}
                    className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                  />
                </div>

                <div>
                  <label className="block text-sm font-semibold text-slate-700 mb-2">
                    目標 %fT {'>'} MIC (%)
                  </label>
                  <input
                    type="number"
                    step="5"
                    value={targetFtMic}
                    onChange={(e) => setTargetFtMic(parseFloat(e.target.value))}
                    className="w-full px-4 py-2 border-2 border-slate-300 rounded-lg focus:border-blue-500 focus:outline-none transition-colors"
                  />
                </div>
              </div>
            </div>

            {/* 計算ボタン */}
            <button
              onClick={handleCalculate}
              disabled={isCalculating}
              className="w-full bg-gradient-to-r from-blue-600 to-indigo-600 text-white font-bold py-4 px-6 rounded-xl shadow-lg hover:shadow-xl transform hover:-translate-y-1 transition-all duration-200 disabled:opacity-50 disabled:cursor-not-allowed disabled:transform-none"
            >
              {isCalculating ? (
                <span className="flex items-center justify-center gap-2">
                  <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-white"></div>
                  計算中... (10,000シミュレーション実行中)
                </span>
              ) : (
                'シミュレーション実行'
              )}
            </button>

            {/* エラー表示 */}
            {error && (
              <div className="bg-red-50 border-2 border-red-200 rounded-xl p-4 flex items-start gap-3">
                <AlertCircle className="w-5 h-5 text-red-600 flex-shrink-0 mt-0.5" />
                <div className="text-red-800">{error}</div>
              </div>
            )}
          </div>

          {/* 右側：結果セクション */}
          <div className="space-y-6">
            {results ? (
              <>
                {/* 腎機能情報カード */}
                <div className="bg-white rounded-2xl shadow-lg p-6 border border-slate-200">
                  <h3 className="text-xl font-bold text-slate-800 mb-4">腎機能評価・定常状態濃度</h3>
                  <div className="grid grid-cols-2 gap-4">
                    <div className="bg-blue-50 rounded-lg p-4">
                      <div className="text-sm text-slate-600 mb-1">計算式</div>
                      <div className="text-lg font-bold text-blue-700">{results.ccrMethod}</div>
                    </div>
                    <div className="bg-indigo-50 rounded-lg p-4">
                      <div className="text-sm text-slate-600 mb-1">CCr</div>
                      <div className="text-lg font-bold text-indigo-700">{results.ccr.toFixed(1)} mL/min</div>
                    </div>
                    <div className="bg-purple-50 rounded-lg p-4">
                      <div className="text-sm text-slate-600 mb-1">クリアランス (CL)</div>
                      <div className="text-lg font-bold text-purple-700">{results.cl.toFixed(2)} L/hr</div>
                    </div>
                    <div className="bg-emerald-50 rounded-lg p-4">
                      <div className="text-sm text-slate-600 mb-1">Cmin (トラフ)</div>
                      <div className="text-lg font-bold text-emerald-700">{results.cMin.toFixed(2)} μg/mL</div>
                    </div>
                    <div className="bg-rose-50 rounded-lg p-4 col-span-2">
                      <div className="text-sm text-slate-600 mb-1">Cmax (ピーク)</div>
                      <div className="text-lg font-bold text-rose-700">{results.cMax.toFixed(2)} μg/mL</div>
                    </div>
                  </div>
                  <div className="mt-3 text-xs text-slate-600 bg-blue-50 rounded-lg p-2">
                    ※ 定常状態における総濃度。遊離形濃度 = 総濃度 × {PIPC_PHARMACOKINETICS.freeFormFraction}
                  </div>
                </div>

                {/* 血中濃度グラフ */}
                <div className="bg-white rounded-2xl shadow-lg p-6 border border-slate-200">
                  <h3 className="text-xl font-bold text-slate-800 mb-4">血中濃度推移（定常状態・2サイクル表示）</h3>
                  <ResponsiveContainer width="100%" height={300}>
                    <LineChart data={results.concentrationProfile}>
                      <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                      <XAxis 
                        dataKey="time" 
                        label={{ value: '時間 (hr)', position: 'insideBottom', offset: -5 }}
                        stroke="#64748b"
                      />
                      <YAxis 
                        label={{ value: '濃度 (μg/mL)', angle: -90, position: 'insideLeft' }}
                        stroke="#64748b"
                      />
                      <Tooltip 
                        contentStyle={{ backgroundColor: 'white', border: '2px solid #e2e8f0', borderRadius: '8px' }}
                        formatter={(value) => [value.toFixed(2) + ' μg/mL', '総濃度']}
                        labelFormatter={(label) => `時間: ${label.toFixed(2)} hr`}
                      />
                      <Legend />
                      <Line 
                        type="monotone" 
                        dataKey="concentration" 
                        stroke="#3b82f6" 
                        strokeWidth={3}
                        name="PIPC総濃度"
                        dot={false}
                      />
                      <ReferenceLine y={16} stroke="#ef4444" strokeDasharray="5 5" label={{ value: 'MIC 16', fill: '#ef4444', fontSize: 12 }} />
                      <ReferenceLine y={8} stroke="#f97316" strokeDasharray="5 5" label={{ value: 'MIC 8', fill: '#f97316', fontSize: 12 }} />
                      <ReferenceLine y={4} stroke="#eab308" strokeDasharray="5 5" label={{ value: 'MIC 4', fill: '#eab308', fontSize: 12 }} />
                    </LineChart>
                  </ResponsiveContainer>
                  <div className="mt-3 text-sm text-slate-600 bg-slate-50 rounded-lg p-3">
                    <div className="mb-1"><strong>注意:</strong> グラフは総濃度を示しています。遊離形濃度は総濃度 × {PIPC_PHARMACOKINETICS.freeFormFraction} です。</div>
                    <div className="text-blue-700"><strong>定常状態:</strong> グラフはCmin（トラフ値）から始まり、2サイクル繰り返すことで定常状態を視覚的に確認できます。</div>
                  </div>
                </div>

                {/* PTAテーブル */}
                <div className="bg-white rounded-2xl shadow-lg p-6 border border-slate-200">
                  <h3 className="text-xl font-bold text-slate-800 mb-4">
                    PTA (Probability of Target Attainment)
                  </h3>
                  <div className="text-sm text-slate-600 mb-4 bg-blue-50 rounded-lg p-3">
                    目標: {targetFtMic}% fT {'>'} MIC を達成する確率（10,000シミュレーション）
                  </div>
                  <div className="overflow-x-auto">
                    <table className="w-full">
                      <thead>
                        <tr className="bg-slate-100">
                          <th className="px-4 py-3 text-left text-sm font-bold text-slate-700 rounded-tl-lg">
                            MIC (μg/mL)
                          </th>
                          <th className="px-4 py-3 text-right text-sm font-bold text-slate-700 rounded-tr-lg">
                            PTA (%)
                          </th>
                        </tr>
                      </thead>
                      <tbody>
                        {results.ptaResults.map((result, index) => (
                          <tr 
                            key={result.mic} 
                            className={`border-b border-slate-200 ${
                              result.pta >= 90 ? 'bg-green-50' : 
                              result.pta >= 70 ? 'bg-yellow-50' : 
                              'bg-red-50'
                            }`}
                          >
                            <td className="px-4 py-3 text-sm font-semibold text-slate-700">
                              {result.mic}
                            </td>
                            <td className="px-4 py-3 text-right">
                              <span className={`text-sm font-bold ${
                                result.pta >= 90 ? 'text-green-700' : 
                                result.pta >= 70 ? 'text-yellow-700' : 
                                'text-red-700'
                              }`}>
                                {result.pta.toFixed(1)}%
                              </span>
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                  <div className="mt-4 text-xs text-slate-600 space-y-1">
                    <div className="flex items-center gap-2">
                      <div className="w-4 h-4 bg-green-50 border border-green-200 rounded"></div>
                      <span>PTA ≥ 90%: 高い達成確率</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <div className="w-4 h-4 bg-yellow-50 border border-yellow-200 rounded"></div>
                      <span>PTA 70-90%: 中等度の達成確率</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <div className="w-4 h-4 bg-red-50 border border-red-200 rounded"></div>
                      <span>PTA {'<'} 70%: 低い達成確率</span>
                    </div>
                  </div>
                </div>
              </>
            ) : (
              <div className="bg-white rounded-2xl shadow-lg p-12 border border-slate-200 text-center">
                <Activity className="w-16 h-16 text-slate-300 mx-auto mb-4" />
                <p className="text-slate-500 text-lg">
                  患者情報と投与設計を入力して<br />シミュレーションを実行してください
                </p>
              </div>
            )}
          </div>
        </div>

        {/* フッター */}
        <div className="mt-8 text-center text-sm text-slate-600 bg-white rounded-xl p-6 shadow-md border border-slate-200">
          <div className="font-semibold mb-2">薬物動態パラメータ（PIPC）</div>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-xs">
            <div>
              <div className="text-slate-500">分布容積 (Vd)</div>
              <div className="font-bold">{PIPC_PHARMACOKINETICS.Vd} L</div>
            </div>
            <div>
              <div className="text-slate-500">基本CL</div>
              <div className="font-bold">{PIPC_PHARMACOKINETICS.baseCL} L/hr</div>
            </div>
            <div>
              <div className="text-slate-500">個人間誤差 (IIV)</div>
              <div className="font-bold">{(PIPC_PHARMACOKINETICS.IIV_CV * 100).toFixed(1)}%</div>
            </div>
            <div>
              <div className="text-slate-500">遊離形分画</div>
              <div className="font-bold">{(PIPC_PHARMACOKINETICS.freeFormFraction * 100).toFixed(0)}%</div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}