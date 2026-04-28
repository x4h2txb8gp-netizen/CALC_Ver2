"use strict";

/* ══════════════════════════════════════════════════════════
   УТИЛИТЫ
══════════════════════════════════════════════════════════ */

var $ = function(id) { return document.getElementById(id); };

var TO_MBAR = {
    mbar:  function(v) { return v; },
    kPa:   function(v) { return v * 10; },
    mmHg:  function(v) { return v * 1.333224; },
    atm:   function(v) { return v * 1013.25; },
    kgcm2: function(v) { return v * 980.665; }
};

/* обратная конверсия: мбар → выбранная единица */
var FROM_MBAR = {
    mbar:  function(v) { return v; },
    kPa:   function(v) { return v / 10; },
    mmHg:  function(v) { return v / 1.333224; },
    atm:   function(v) { return v / 1013.25; },
    kgcm2: function(v) { return v / 980.665; }
};

/* русские названия единиц давления */
var UNIT_LABELS = {
    mbar:  "мбар",
    kPa:   "кПа",
    mmHg:  "мм рт.ст.",
    atm:   "атм",
    kgcm2: "кг/см²"
};

function getPmbar() {
    return TO_MBAR[$("pru").value](parseFloat($("prv").value));
}

function getDPmbar() {
    return TO_MBAR[$("pru").value](Math.abs(parseFloat($("prc").value) || 0));
}

/* ══════════════════════════════════════════════════════════
   ТОЧНОСТЬ (степпер)
══════════════════════════════════════════════════════════ */

var PRECISION = 3;
var PREC_MIN  = 0;
var PREC_MAX  = 8;

function updatePrecDisplay() {
    $("prec-display").textContent = PRECISION;
    $("dec-prec").disabled = (PRECISION <= PREC_MIN);
    $("inc-prec").disabled = (PRECISION >= PREC_MAX);
}

/* ══════════════════════════════════════════════════════════
   ДАВЛЕНИЕ НАСЫЩЕННОГО ПАРА
══════════════════════════════════════════════════════════ */

var GW = [
    -2.8365744e3,
    -6.028076559e3,
     1.954263612e1,
    -2.737830188e-2,
     1.6261698e-5,
     7.0229056e-10,
    -1.8680009e-13,
     2.7150305
];

function esWater(Tc) {
    var Tk = Tc + 273.15;
    var lnE =
        GW[0] / (Tk * Tk) +
        GW[1] / Tk +
        GW[2] +
        GW[3] * Tk +
        GW[4] * Tk * Tk +
        GW[5] * Tk * Tk * Tk +
        GW[6] * Tk * Tk * Tk * Tk +
        GW[7] * Math.log(Tk);
    return Math.exp(lnE) / 100;
}

function desWater(Tc) {
    var Tk = Tc + 273.15;
    var dlnE =
        GW[0] * (-2) / (Tk * Tk * Tk) +
        GW[1] * (-1) / (Tk * Tk) +
        GW[3] +
        GW[4] * 2 * Tk +
        GW[5] * 3 * Tk * Tk +
        GW[6] * 4 * Tk * Tk * Tk +
        GW[7] / Tk;
    return esWater(Tc) * dlnE;
}

function esBuckWater(Tc) {
    return 6.1121 * Math.exp((18.678 - Tc / 234.5) * Tc / (257.14 + Tc));
}

function desBuckWater(Tc) {
    var a = 18.678, b = 234.5, c = 257.14;
    var num = (a - 2 * Tc / b) * (c + Tc) - (a * Tc - Tc * Tc / b);
    var den = (c + Tc) * (c + Tc);
    return esBuckWater(Tc) * num / den;
}

function esW(Tc)  { return Tc >= 0 ? esWater(Tc)  : esBuckWater(Tc); }
function desW(Tc) { return Tc >= 0 ? desWater(Tc) : desBuckWater(Tc); }

var T0GG = 273.16;

function esIce(Tc) {
    var Tk = Tc + 273.15;
    var ratio = T0GG / Tk;
    var log10ei =
        -9.09718  * (ratio - 1) -
         3.56654  * Math.log10(ratio) +
         0.876793 * (1 - Tk / T0GG) +
         Math.log10(6.1071);
    return Math.pow(10, log10ei);
}

function desIce(Tc) {
    var Tk = Tc + 273.15;
    var d1 =  9.09718  * T0GG / (Tk * Tk);
    var d2 =  3.56654  / (Tk * Math.LN10);
    var d3 = -0.876793 / T0GG;
    return esIce(Tc) * Math.LN10 * (d1 + d2 + d3);
}

function esF(Tc)  { return esW(Tc); }
function desF(Tc) { return desW(Tc); }





/* ══════════════════════════════════════════════════════════
   ФАКТОР УСИЛЕНИЯ (ENHANCEMENT FACTOR) f(T, P)
   ─────────────────────────────────────────────────────────
   Hardy (1998), ITS-90. Учитывает неидеальность влажного
   воздуха: реальное равновесное давление пара в воздухе
   немного выше, чем над чистой водой:
       e_w*(T,P) = f(T,P) · es(T)

   ln(f) = α(t)·(1 - es/P) + β(t)·(P/es - 1)
     α(t) = A0 + A1·t + A2·t² + A3·t³
     β(t) = exp(B0 + B1·t + B2·t² + B3·t³)
   t — температура в °C, es и P — в одинаковых единицах.

   При нормальных условиях f ≈ 1.004 (поправка ~0.4 %).
══════════════════════════════════════════════════════════ */

/* Коэффициенты Hardy 1998 для воды (t в °C) */
var FA = [3.53624e-4, 2.9328363e-5, 2.6168979e-7, 8.5813609e-9];
var FB = [-1.07588e1, 6.3268134e-2, -2.5368934e-4, 6.3405286e-7];

/* Коэффициенты Hardy 1998 для льда */
var FA_ICE = [3.64449e-4, 2.9367585e-5, 4.8874766e-7, 4.3669918e-9];
var FB_ICE = [-1.0728876e1, 7.6215115e-2, -1.7540413e-4, 2.1689756e-6];

/* Универсальный расчёт f */
function _enhFactor(Tc, Pmbar, esmbar, A, B) {
    if (!isFinite(Pmbar) || Pmbar <= 0 || esmbar <= 0) return 1;
    var t = Tc;
    var alpha = A[0] + t*(A[1] + t*(A[2] + t*A[3]));
    var betaExp = B[0] + t*(B[1] + t*(B[2] + t*B[3]));
    var beta = Math.exp(betaExp);
    var lnf = alpha * (1 - esmbar / Pmbar) + beta * (Pmbar / esmbar - 1);
    return Math.exp(lnf);
}

/* f над водой */
function fWater(Tc, Pmbar) {
    return _enhFactor(Tc, Pmbar, esW(Tc), FA, FB);
}

/* f надо льдом */
function fIce(Tc, Pmbar) {
    return _enhFactor(Tc, Pmbar, esIce(Tc), FA_ICE, FB_ICE);
}

/* Эффективное давление насыщения в воздухе */
function esEffW(Tc, Pmbar)   { return fWater(Tc, Pmbar) * esW(Tc); }
function esEffIce(Tc, Pmbar) { return fIce(Tc, Pmbar)   * esIce(Tc); }

/* Производная esEff по T (для Ньютона в invertES).
   Численная — точности достаточно, df/dT мало. */
function desEffW(Tc, Pmbar) {
    var h = 1e-4;
    return (esEffW(Tc + h, Pmbar) - esEffW(Tc - h, Pmbar)) / (2 * h);
}
function desEffIce(Tc, Pmbar) {
    var h = 1e-4;
    return (esEffIce(Tc + h, Pmbar) - esEffIce(Tc - h, Pmbar)) / (2 * h);
}



/* ══════════════════════════════════════════════════════════
   ОБРАТНАЯ ЗАДАЧА: e → T
══════════════════════════════════════════════════════════ */

function invertES(eTarget, esFn, desFn, lo, hi) {
    if (!isFinite(eTarget) || eTarget <= 0) return NaN;
    var eLo = esFn(lo);
    var eHi = esFn(hi);
    if (eTarget < eLo || eTarget > eHi) return NaN;

    var i, mid;
    for (i = 0; i < 200; i++) {
        mid = (lo + hi) / 2;
        if (esFn(mid) < eTarget) lo = mid; else hi = mid;
        if (hi - lo < 1e-6) break;
    }

    var T = (lo + hi) / 2;
    var f, df, dt;
    for (i = 0; i < 50; i++) {
        f  = esFn(T) - eTarget;
        df = desFn(T);
        if (Math.abs(df) < 1e-30) break;
        dt = f / df;
        T -= dt;
        if (Math.abs(dt) < 1e-9) break;
    }
    return T;
}

function TdFromE(e, Pmbar) {
    return invertES(
        e,
        function(T) { return esEffW(T, Pmbar); },
        function(T) { return desEffW(T, Pmbar); },
        -80, 200
    );
}

function TfFromE(e, Pmbar) {
    var eMax = esEffIce(-0.01, Pmbar);
    if (e > eMax) return NaN;
    return invertES(
        e,
        function(T) { return esEffIce(T, Pmbar); },
        function(T) { return desEffIce(T, Pmbar); },
        -100, -0.01
    );
}

/* ══════════════════════════════════════════════════════════
   ВСПОМОГАТЕЛЬНЫЕ ФОРМУЛЫ
══════════════════════════════════════════════════════════ */

function absHum(e, Tc) {
    return 216.679 * e / (Tc + 273.15);
}

function mixRatio(e, P) {
    if (P <= e) return NaN;
    return 621.9907 * e / (P - e);
}

function enthalpy(Tc, Wgkg) {
    var w = Wgkg / 1000;
    return 1.006 * Tc + w * (2500.9 + 1.86 * Tc);
}

function wetBulb(Tc, e, P) {
    var A = 0.000662;
    var Tw = TdFromE(e, P);
    if (!isFinite(Tw)) Tw = Tc - 5;
    if (Tw > Tc) Tw = Tc;

    var i, F, dF, dt;
    for (i = 0; i < 80; i++) {
        F  = esEffW(Tw, P) - A * P * (Tc - Tw) - e;
        dF = desEffW(Tw, P) + A * P;
        if (Math.abs(dF) < 1e-30) break;
        dt = F / dF;
        Tw -= dt;
        if (Math.abs(dt) < 1e-9) break;
    }
    return Tw;
}

function primaryToE(key, val, T, P) {
    var A;
    switch (key) {
        case "RH": return esEffW(T, P) * val / 100;
        case "Td": return esEffW(val, P);
        case "Tf": return esEffIce(val, P);
        case "Av": return val * (T + 273.15) / 216.679;
        case "W":  return val * P / (621.9907 + val);
        case "e":  return val;
        case "Tw":
            A = 0.000662;
            return esEffW(val, P) - A * P * (T - val);
        default:   return NaN;
    }
}

/* ══════════════════════════════════════════════════════════
   ПЕРВИЧНЫЙ ПАРАМЕТР → e
══════════════════════════════════════════════════════════ */

function primaryToE(key, val, T, P) {
    var A;
    switch (key) {
        case "RH": return esEffW(T, P) * val / 100;
        case "Td": return esEffW(val, P);
        case "Tf": return esEffIce(val, P);
        case "Av": return val * (T + 273.15) / 216.679;
        case "W":  return val * P / (621.9907 + val);
        case "e":  return val;
        case "Tw":
            A = 0.000662;
            return esEffW(val, P) - A * P * (T - val);
        default:   return NaN;
    }
}

/* ══════════════════════════════════════════════════════════
   ПОЛНЫЙ РАСЧЁТ
══════════════════════════════════════════════════════════ */

function calcAll(key, val, T, P) {
    var e = primaryToE(key, val, T, P);
    if (!isFinite(e) || e < 0) return null;

    var esVal   = esF(T);                  /* чистое es(T) — для вывода */
    var esEff   = esEffW(T, P);            /* эффективное — для RH */
    var RH      = e / esEff * 100;
    var Td      = TdFromE(e, P);
    var Tf      = TfFromE(e, P);
    var Av      = absHum(e, T);
    var W       = mixRatio(e, P);
    var Tw      = wetBulb(T, e, P);
    var H       = isFinite(W) ? enthalpy(T, W) : NaN;

    return {
        RH: RH, Td: Td, Tf: Tf,
        Av: Av, W: W, es: esVal,
        e: e, Tw: Tw, H: H
    };
}

/* ══════════════════════════════════════════════════════════
   ПОГРЕШНОСТИ
   Центральные конечные разности 2-го порядка
   σf = sqrt( Σ (∂f/∂xi * σxi)² )
══════════════════════════════════════════════════════════ */

var RKEYS = ["RH", "Td", "Tf", "Av", "W", "es", "e", "Tw", "H"];




['pu', 'tu', 'prc'].forEach(function(id) {
    var input = document.getElementById(id);
    if (!input) return;
    
    input.setAttribute('maxlength', '3');
    
    input.addEventListener('input', function() {
        var val = this.value.replace(/[^\d.-]/g, '');
        val = val.replace(/(?!^)-/g, '');
        val = val.replace(/^(.+?)\.(.*)\./, '$1.$2');
        if (val.length > 3) val = val.slice(0, 3);
        this.value = val;
    });
});




function calcUnc(key, val, dVal, T, dT, P, dP) {
    var unc = {};
    var k, j, s, sumSq, h, rp, rm, deriv;

    var sources = [
        {
            dx: dVal, ref: val,
            fp: function(hh) { return calcAll(key, val + hh, T, P); },
            fm: function(hh) { return calcAll(key, val - hh, T, P); }
        },
        {
            dx: dT, ref: T,
            fp: function(hh) { return calcAll(key, val, T + hh, P); },
            fm: function(hh) { return calcAll(key, val, T - hh, P); }
        },
        {
            dx: dP, ref: P,
            fp: function(hh) { return calcAll(key, val, T, P + hh); },
            fm: function(hh) { return calcAll(key, val, T, P - hh); }
        }
    ];

    for (j = 0; j < RKEYS.length; j++) {
        k = RKEYS[j];
        sumSq = 0;

        for (s = 0; s < sources.length; s++) {
            if (sources[s].dx <= 0) continue;
            h  = Math.max(Math.abs(sources[s].ref) * 1e-5, 1e-7);
            rp = sources[s].fp(h);
            rm = sources[s].fm(h);
            if (rp && rm && isFinite(rp[k]) && isFinite(rm[k])) {
                deriv = (rp[k] - rm[k]) / (2 * h);
                sumSq += deriv * deriv * sources[s].dx * sources[s].dx;
            }
        }

        unc[k] = Math.sqrt(sumSq);
    }
    return unc;
}

/* ══════════════════════════════════════════════════════════
   ВАЛИДАЦИЯ
══════════════════════════════════════════════════════════ */

function validate() {
    var errors = [];
    var ids = ["pv", "tv", "prv"];
    var i;

    for (i = 0; i < ids.length; i++) {
        $(ids[i]).classList.remove("ferr");
    }

    var T    = parseFloat($("tv").value);
    var pVal = parseFloat($("prv").value);
    var pU   = $("pru").value;
    var Pkpa = TO_MBAR[pU](pVal) / 10;
    var pv   = parseFloat($("pv").value);

    if (!isFinite(T)) {
        errors.push("Введите температуру воздуха.");
        $("tv").classList.add("ferr");
    } else if (T < -100 || T > 200) {
        errors.push("Температура воздуха должна быть от −100 до +200 °С.");
        $("tv").classList.add("ferr");
    }

    if (!isFinite(pVal)) {
        errors.push("Введите атмосферное давление.");
        $("prv").classList.add("ferr");
    } else if (!isFinite(Pkpa) || Pkpa < 10 || Pkpa > 150) {
        errors.push("Давление должно быть от 10 до 150 кПа (в пересчёте).");
        $("prv").classList.add("ferr");
    }

    if (!isFinite(pv)) {
        errors.push("Введите значение первичного параметра.");
        $("pv").classList.add("ferr");
    }

    return errors;
}

/* ══════════════════════════════════════════════════════════
   ОТОБРАЖЕНИЕ РЕЗУЛЬТАТОВ
══════════════════════════════════════════════════════════ */

function setVal(id, v, digits) {
    var el = $(id);
    if (!el) return;
    if (!isFinite(v)) {
        el.textContent = "—";
        el.classList.add("na");
    } else {
        el.textContent = v.toFixed(digits).replace(".", ",");
        el.classList.remove("na");
    }
}

function setUnc(id, u, digits) {
    var el = $(id);
    if (!el) return;
    if (isFinite(u) && u > 0) {
        el.textContent = "\u00B1\u00A0" + u.toFixed(digits).replace(".", ",");
        el.style.display = "";
    } else {
        el.textContent = "";
        el.style.display = "none";
    }
}

function clearResults() {
    var j, k;
    for (j = 0; j < RKEYS.length; j++) {
        k = RKEYS[j];
        setVal("o_" + k, NaN, PRECISION);
        setUnc("u_" + k, 0, PRECISION);
    }
}

/* ══════════════════════════════════════════════════════════
   ГЛАВНАЯ ФУНКЦИЯ ПЕРЕСЧЁТА
══════════════════════════════════════════════════════════ */

function render() {
    var msg  = $("msg");
    var errs = validate();

    if (errs.length > 0) {
        msg.className = "msgbox err";
        msg.innerHTML = errs.join("<br>");
        clearResults();
        return;
    }

    var key  = $("pp").value;
    var val  = parseFloat($("pv").value);
    var T    = parseFloat($("tv").value);
    var P    = getPmbar();
    var dVal = Math.abs(parseFloat($("pu").value) || 0);
    var dT   = Math.abs(parseFloat($("tu").value) || 0);
    var dP   = getDPmbar();

    var esAtT = esEffW(T, P);   /* эффективное для проверки физичности */

    if (key === "RH" && val > 100 && esAtT >= P) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "При заданной температуре (" + T.toFixed(1) + " °С) и давлении (" +
            (P / 10).toFixed(2) + " кПа) относительная влажность не может " +
            "достигать " + val.toFixed(1) + " %.";
        clearResults();
        return;
    }

    var res = calcAll(key, val, T, P);

    if (!res) {
        msg.className = "msgbox err";
        msg.innerHTML = "Невозможно рассчитать параметры. Проверьте значения.";
        clearResults();
        return;
    }

    if (res.e > P) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "Парциальное давление пара (" + res.e.toFixed(2) +
            " мбар) превышает общее давление (" + P.toFixed(2) + " мбар).";
        clearResults();
        return;
    }

    if (res.RH > 100.05) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "Расчётная RH = " + res.RH.toFixed(2) +
            " % > 100 %. Воздух пересыщен, результаты могут быть некорректны.";
        clearResults();
        return;
    }

    msg.className = "msgbox";
    msg.innerHTML = "";

    var j, k;
    var pu     = $("pru").value;
    var convFn = FROM_MBAR[pu];
    var uLabel = UNIT_LABELS[pu];

    for (j = 0; j < RKEYS.length; j++) {
        k = RKEYS[j];
        if (k === "es" || k === "e") {
            setVal("o_" + k, convFn(res[k]), PRECISION);
        } else {
            setVal("o_" + k, res[k], PRECISION);
        }
    }

   
    if ($("unit_es")) $("unit_es").textContent = uLabel;
    if ($("unit_e"))  $("unit_e").textContent  = uLabel;

    /* ── Погрешности ── */
    var hasUnc = (dVal > 0 || dT > 0 || dP > 0);
    if (hasUnc) {
        var unc = calcUnc(key, val, dVal, T, dT, P, dP);
        for (j = 0; j < RKEYS.length; j++) {
            k = RKEYS[j];
            if (k === "es" || k === "e") {
                setUnc("u_" + k, convFn(unc[k]), PRECISION);
            } else {
                setUnc("u_" + k, unc[k], PRECISION);
            }
        }
    } else {
        for (j = 0; j < RKEYS.length; j++) {
            setUnc("u_" + RKEYS[j], 0, PRECISION);
        }
    }

    /* ── Подсветка входного параметра ── */
    var rows = document.querySelectorAll(".out-row[data-k]");
    for (j = 0; j < rows.length; j++) {
        if (rows[j].dataset.k === key) {
            rows[j].classList.add("is-input");
        } else {
            rows[j].classList.remove("is-input");
        }
    }

    /* ── График ── */
    updatePsyChart(T, res);

}

/* ══════════════════════════════════════════════════════════
   ВЫДЕЛЕНИЕ СТРОКИ (клик)
══════════════════════════════════════════════════════════ */

function initHighlight() {
    var section = $("out-section");
    if (!section) return;

    var rows = section.querySelectorAll(".out-row[data-k]");
    var i;

    for (i = 0; i < rows.length; i++) {
        rows[i].addEventListener("click", function() {
            if (this.classList.contains("hl")) {
                this.classList.remove("hl");
                return;
            }
            var all = section.querySelectorAll(".out-row.hl");
            for (var j = 0; j < all.length; j++) {
                all[j].classList.remove("hl");
            }
            this.classList.add("hl");
        });
    }
}

/* ══════════════════════════════════════════════════════════
   ПРИВЯЗКА СОБЫТИЙ
══════════════════════════════════════════════════════════ */

function initEvents() {
    var inputIds  = ["pv", "pu", "tv", "tu", "prv", "prc"];
    var selectIds = ["pp", "pru"];
    var i;

    for (i = 0; i < inputIds.length; i++) {
        (function(id) {
            var el = $(id);
            if (!el) return;
            el.addEventListener("input", function() { render(); });
        })(inputIds[i]);
    }

    for (i = 0; i < selectIds.length; i++) {
        (function(id) {
            var el = $(id);
            if (!el) return;
            el.addEventListener("change", function() { render(); });
        })(selectIds[i]);
    }

    // Единицы первичного параметра
    var PP_UNITS = {
        RH: "%", Td: "°С", Tf: "°С",
        Av: "г/м³", W: "г/кг", e: "мбар", Tw: "°С"
    };

    $("pp").addEventListener("change", function() {
        var u = $("pp-unit");
        if (u) u.textContent = PP_UNITS[this.value] || "";
    });

        /* Перерисовка графика при ресайзе */
    var resizeTimer;
    window.addEventListener('resize', function() {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(function() {
            if (psyChart) {
                psyChart.destroy();
                psyChart = null;
                render();
            }
        }, 250);
    });
}

/* ══════════════════════════════════════════════════════════
   СТЕППЕР ТОЧНОСТИ
══════════════════════════════════════════════════════════ */

function initPrecision() {
    updatePrecDisplay();

    $("dec-prec").addEventListener("click", function() {
        if (PRECISION > PREC_MIN) {
            PRECISION--;
            updatePrecDisplay();
            render();
        }
    });

    $("inc-prec").addEventListener("click", function() {
        if (PRECISION < PREC_MAX) {
            PRECISION++;
            updatePrecDisplay();
            render();
        }
    });
}


/* ═══════════════════════════════════════
   Психрометрический график es(T)
   Адаптивный: desktop / tablet / mobile
   ═══════════════════════════════════════ */

var psyChart = null;

function getChartConfig() {
    var w = window.innerWidth;

    if (w <= 380) {
        return {
            aspectRatio: 1.0,
            pointRadius: 6,
            pointHover: 8,
            lineWidth: 2,
            dashWidth: 1,
            axisFont: 9,
            titleFont: 10,
            tickFont: 8,
            legendFont: 9,
            legendPad: 8,
            stepX: 50,
            showAxisTitleY: false,
            showAxisTitleX: true,
            tooltipTitle: 11,
            tooltipBody: 10,
            tooltipPad: 6
        };
    }
    if (w <= 700) {
        return {
            aspectRatio: 1.1,
            pointRadius: 7,
            pointHover: 9,
            lineWidth: 2,
            dashWidth: 1.2,
            axisFont: 10,
            titleFont: 11,
            tickFont: 9,
            legendFont: 10,
            legendPad: 10,
            stepX: 40,
            showAxisTitleY: true,
            showAxisTitleX: true,
            tooltipTitle: 12,
            tooltipBody: 11,
            tooltipPad: 8
        };
    }
    if (w <= 1024) {
        return {
            aspectRatio: 1.5,
            pointRadius: 8,
            pointHover: 10,
            lineWidth: 2.5,
            dashWidth: 1.5,
            axisFont: 11,
            titleFont: 12,
            tickFont: 10,
            legendFont: 11,
            legendPad: 12,
            stepX: 30,
            showAxisTitleY: true,
            showAxisTitleX: true,
            tooltipTitle: 12,
            tooltipBody: 11,
            tooltipPad: 8
        };
    }
    return {
        aspectRatio: 1.8,
        pointRadius: 9,
        pointHover: 12,
        lineWidth: 2.5,
        dashWidth: 1.5,
        axisFont: 11,
        titleFont: 13,
        tickFont: 11,
        legendFont: 12,
        legendPad: 14,
        stepX: 20,
        showAxisTitleY: true,
        showAxisTitleX: true,
        tooltipTitle: 13,
        tooltipBody: 12,
        tooltipPad: 10
    };
}

function updatePsyChart(T, res) {
    if (!res) return;

    var eKpa  = res.e  / 10;
    var esKpa = res.es / 10;
    var TdVal = isFinite(res.Td) ? res.Td : null;
    var cfg   = getChartConfig();

    /* ── Кривая насыщения ── */
    var curve = [];
    var t;
    for (t = -100; t <= 200; t += 1) {
        curve.push({ x: t, y: esF(t) / 10 });
    }

    /* ── Датасеты ── */
    var datasets = [
        {
            label: 'Кривая насыщения es(T)',
            data: curve,
            showLine: true,
            borderColor: '#3b82f6',
            backgroundColor: 'rgba(59,130,246,0.08)',
            borderWidth: cfg.lineWidth,
            pointRadius: 0,
            pointHitRadius: 0,
            fill: true,
            tension: 0.4,
            order: 4
        }
    ];

    if (TdVal !== null) {
        datasets.push({
            label: '_hline',
            data: [
                { x: TdVal, y: eKpa },
                { x: T,     y: eKpa }
            ],
            showLine: true,
            borderColor: '#94a3b8',
            borderWidth: cfg.dashWidth,
            borderDash: [6, 4],
            pointRadius: 0,
            pointHitRadius: 0,
            fill: false,
            order: 3
        });
    }

    if (Math.abs(esKpa - eKpa) > 0.001) {
        datasets.push({
            label: '_vline',
            data: [
                { x: T, y: eKpa  },
                { x: T, y: esKpa }
            ],
            showLine: true,
            borderColor: '#94a3b8',
            borderWidth: cfg.dashWidth,
            borderDash: [4, 3],
            pointRadius: 0,
            pointHitRadius: 0,
            fill: false,
            order: 3
        });
    }

    datasets.push({
        label: 'Текущее (' + T.toFixed(1) + ' °C; ' + eKpa.toFixed(3) + ' кПа)',
        data: [{ x: T, y: eKpa }],
        pointRadius: cfg.pointRadius,
        pointHoverRadius: cfg.pointHover,
        pointBackgroundColor: '#ef4444',
        pointBorderColor: '#fff',
        pointBorderWidth: 2,
        showLine: false,
        order: 1
    });

    if (TdVal !== null) {
        datasets.push({
            label: 'Точка росы (' + TdVal.toFixed(1) + ' °C; ' + eKpa.toFixed(3) + ' кПа)',
            data: [{ x: TdVal, y: eKpa }],
            pointRadius: cfg.pointRadius,
            pointHoverRadius: cfg.pointHover,
            pointBackgroundColor: '#22c55e',
            pointBorderColor: '#fff',
            pointBorderWidth: 2,
            pointStyle: 'triangle',
            showLine: false,
            order: 1
        });
    }

    var yMax = Math.ceil(esKpa * 1.4);
    if (yMax < 1) yMax = 1;



// ★ Вычисляем адаптивные границы оси X ★
var userMin = T;  // самая холодная точка — текущая температура
if (TdVal !== null) userMin = Math.min(userMin, TdVal);  // точка росы может быть холоднее

var xMin = Math.floor(userMin - 30);  // левая граница: на 30 меньше самой холодной точки
var xMax = Math.ceil(userMin + 30);   // правая граница: на 30 больше самой тёплой (симметрично)

// Не даём графику сжаться слишком сильно
if (xMax - xMin < 60) {
    var center = (xMin + xMax) / 2;
    xMin = center - 30;
    xMax = center + 30;
}





    var chartOptions = {
        responsive: true,
        maintainAspectRatio: true,
        aspectRatio: cfg.aspectRatio,
        animation: { duration: 300 },
        interaction: {
            mode: 'nearest',
            intersect: true
        },
        plugins: {
            legend: {
                position: 'bottom',
                labels: {
                    usePointStyle: true,
                    padding: cfg.legendPad,
                    font: { size: cfg.legendFont },
                    boxWidth: cfg.legendFont,
                    filter: function(item) {
                        return item.text.charAt(0) !== '_';
                    }
                }
            },
            tooltip: {
                backgroundColor: 'rgba(15,23,42,0.92)',
                titleFont: { size: cfg.tooltipTitle, weight: '600' },
                bodyFont: { size: cfg.tooltipBody },
                padding: cfg.tooltipPad,
                cornerRadius: 8,
                callbacks: {
                    label: function(context) {
                        var p = context.raw;
                        return context.dataset.label + ': ' +
                               p.x.toFixed(1) + ' °C, ' +
                               p.y.toFixed(4) + ' кПа';
                    }
                }
            }
        },
        scales: {
            x: {
                title: {
                    display: cfg.showAxisTitleX,
                    text: 'Температура, °C',
                    font: { size: cfg.titleFont, weight: '600' },
                    color: '#374151',
                    padding: { top: 4, bottom: 0 }
                },
                min: xMin,
                max: xMax,
                ticks: {
                    stepSize: cfg.stepX,
                    font: { size: cfg.tickFont },
                    color: '#6b7280',
                    maxRotation: 0,
                    autoSkip: true,
                    autoSkipPadding: 8
                },
                grid: { color: 'rgba(0,0,0,0.05)' }
            },
            y: {
                title: {
                    display: cfg.showAxisTitleY,
                    text: 'Давление, кПа',
                    font: { size: cfg.titleFont, weight: '600' },
                    color: '#374151',
                    padding: { top: 0, bottom: 4 }
                },
                min: 0,
                max: yMax,
                ticks: {
                    font: { size: cfg.tickFont },
                    color: '#6b7280',
                    callback: function(v) { return v.toFixed(1); },
                    maxTicksLimit: 8
                },
                grid: { color: 'rgba(0,0,0,0.05)' }
            }
        }
    };

    /* ── Обновить или создать ── */
    if (psyChart) {
        psyChart.data.datasets = datasets;
        psyChart.options = chartOptions;
        psyChart.update();
        return;
    }

    var ctx = document.getElementById('psyChart');
    if (!ctx) return;

    psyChart = new Chart(ctx, {
        type: 'scatter',
        data: { datasets: datasets },
        options: chartOptions
    });
}


/* ══════════════════════════════════════════════════════════
   УНИВЕРСАЛЬНЫЕ КАСТОМНЫЕ ДРОПДАУНЫ (.cdd)
══════════════════════════════════════════════════════════ */

function initCustomDropdowns() {
    var dropdowns = document.querySelectorAll(".cdd");
    for (var d = 0; d < dropdowns.length; d++) {
        initSingleCdd(dropdowns[d]);
    }

    document.addEventListener("click", function() {
        var all = document.querySelectorAll(".cdd.open");
        for (var j = 0; j < all.length; j++) all[j].classList.remove("open");
    });
}

function initSingleCdd(dd) {
    var trigger  = dd.querySelector(".cdd-trigger");
    var menu     = dd.querySelector(".cdd-menu");
    var hidden   = dd.querySelector('input[type="hidden"]');
    var items    = menu.querySelectorAll("li");
    var focusIdx = -1;

    /* Тип триггера: двухколоночный (code+desc) или простой (label) */
    var tCode  = trigger.querySelector(".cdd-code");
    var tDesc  = trigger.querySelector(".cdd-desc");
    var tLabel = trigger.querySelector(".cdd-label");

    /* ── Открыть / закрыть ── */
    trigger.addEventListener("click", function(e) {
        e.stopPropagation();
        var wasOpen = dd.classList.contains("open");
        closeOtherCdd(dd);
        if (!wasOpen) {
            dd.classList.add("open");
            focusIdx = selectedIdx();
            scrollTo(focusIdx);
        } else {
            dd.classList.remove("open");
        }
    });

    /* ── Выбор по клику ── */
    for (var i = 0; i < items.length; i++) {
        (function(li) {
            li.addEventListener("click", function(e) {
                e.stopPropagation();
                pick(li);
            });
        })(items[i]);
    }

    /* ── Клавиатура ── */
    trigger.addEventListener("keydown", function(e) {
        if (e.key === "Enter" || e.key === " " || e.key === "ArrowDown") {
            e.preventDefault();
            if (!dd.classList.contains("open")) {
                dd.classList.add("open");
                focusIdx = selectedIdx();
                scrollTo(focusIdx);
            }
        }
    });

    dd.addEventListener("keydown", function(e) {
        if (!dd.classList.contains("open")) return;
        if (e.key === "ArrowDown") {
            e.preventDefault();
            focusIdx = Math.min(focusIdx + 1, items.length - 1);
            highlight(focusIdx); scrollTo(focusIdx);
        } else if (e.key === "ArrowUp") {
            e.preventDefault();
            focusIdx = Math.max(focusIdx - 1, 0);
            highlight(focusIdx); scrollTo(focusIdx);
        } else if (e.key === "Enter") {
            e.preventDefault();
            if (focusIdx >= 0) pick(items[focusIdx]);
        } else if (e.key === "Escape" || e.key === "Tab") {
            dd.classList.remove("open");
        }
    });

    /* ── Выбор элемента ── */
    function pick(li) {
        for (var j = 0; j < items.length; j++) items[j].classList.remove("selected");
        li.classList.add("selected");

        hidden.value = li.dataset.value;

        /* Обновить текст триггера */
        if (tCode && tDesc) {
            var c = li.querySelector(".cdd-code");
            var d = li.querySelector(".cdd-desc");
            tCode.textContent = c ? c.textContent : "";
            tDesc.textContent = d ? d.textContent : "";
        } else if (tLabel) {
            tLabel.textContent = li.textContent;
        }

        dd.classList.remove("open");
        trigger.focus();

        /* Dispatch change на hidden input */
        var evt;
        if (typeof Event === "function") {
            evt = new Event("change", {bubbles: true});
        } else {
            evt = document.createEvent("Event");
            evt.initEvent("change", true, true);
        }
        hidden.dispatchEvent(evt);
    }

    /* ── Утилиты ── */
    function closeOtherCdd(except) {
        var all = document.querySelectorAll(".cdd.open");
        for (var j = 0; j < all.length; j++) {
            if (all[j] !== except) all[j].classList.remove("open");
        }
    }

    function selectedIdx() {
        for (var j = 0; j < items.length; j++) {
            if (items[j].classList.contains("selected")) return j;
        }
        return 0;
    }

    function highlight(idx) {
        for (var j = 0; j < items.length; j++) {
            items[j].style.outline       = (j === idx) ? "2px solid #3b82f6" : "";
            items[j].style.outlineOffset = (j === idx) ? "-2px" : "";
        }
    }

    function scrollTo(idx) {
        if (idx >= 0 && items[idx]) {
            items[idx].scrollIntoView({block: "nearest"});
        }
    }
}


/* ══════════════════════════════════════════════════════════
   ЗАПУСК
══════════════════════════════════════════════════════════ */

document.addEventListener("DOMContentLoaded", function() {
    initCustomDropdowns();  
    initEvents();
    initHighlight();
    initPrecision();
    render();
});
