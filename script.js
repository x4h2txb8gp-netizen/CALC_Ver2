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

   T >= 0 °С  →  Hardy (1998) ITS-90
   T <  0 °С  →  Goff-Gratch (1946) WMO
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

function esF(Tc)  { return Tc >= 0 ? esWater(Tc) : esIce(Tc); }
function desF(Tc) { return Tc >= 0 ? desWater(Tc) : desIce(Tc); }

/* ══════════════════════════════════════════════════════════
   ОБРАТНАЯ ЗАДАЧА: e → T  (бисекция + Ньютон)
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

function TdFromE(e) {
    return invertES(e, esWater, desWater, -80, 200);
}

function TfFromE(e) {
    var eMax = esIce(-0.01);
    if (e > eMax) return NaN;
    return invertES(e, esIce, desIce, -100, -0.01);
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
    var Tw = TdFromE(e);
    if (!isFinite(Tw)) Tw = Tc - 5;
    if (Tw > Tc) Tw = Tc;

    var i, F, dF, dt;
    for (i = 0; i < 80; i++) {
        F  = esF(Tw) - A * P * (Tc - Tw) - e;
        dF = desF(Tw) + A * P;
        if (Math.abs(dF) < 1e-30) break;
        dt = F / dF;
        Tw -= dt;
        if (Math.abs(dt) < 1e-9) break;
    }
    return Tw;
}

/* ══════════════════════════════════════════════════════════
   ПЕРВИЧНЫЙ ПАРАМЕТР → e (парциальное давление)
══════════════════════════════════════════════════════════ */

function primaryToE(key, val, T, P) {
    var A;
    switch (key) {
        case "RH": return esF(T) * val / 100;
        case "Td": return esWater(val);
        case "Tf": return esIce(val);
        case "Av": return val * (T + 273.15) / 216.679;
        case "W":  return val * P / (621.9907 + val);
        case "e":  return val;
        case "Tw":
            A = 0.000662;
            return esF(val) - A * P * (T - val);
        default:   return NaN;
    }
}

/* ══════════════════════════════════════════════════════════
   ПОЛНЫЙ РАСЧЁТ ВСЕХ ПАРАМЕТРОВ
══════════════════════════════════════════════════════════ */

function calcAll(key, val, T, P) {
    var e = primaryToE(key, val, T, P);
    if (!isFinite(e) || e < 0) return null;

    var esVal = esF(T);
    var RH    = e / esVal * 100;
    var Td    = TdFromE(e);
    var Tf    = TfFromE(e);
    var Av    = absHum(e, T);
    var W     = mixRatio(e, P);
    var Tw    = wetBulb(T, e, P);
    var H     = isFinite(W) ? enthalpy(T, W) : NaN;

    return {
        RH: RH,
        Td: Td,
        Tf: Tf,
        Av: Av,
        W:  W,
        es: esVal,
        e:  e,
        Tw: Tw,
        H:  H
    };
}

/* ══════════════════════════════════════════════════════════
   ПОГРЕШНОСТИ
   Центральные конечные разности 2-го порядка
   σf = √( Σ (∂f/∂xi · σxi)² )
══════════════════════════════════════════════════════════ */

var RKEYS = ["RH", "Td", "Tf", "Av", "W", "es", "e", "Tw", "H"];

function calcUnc(key, val, dVal, T, dT, P, dP) {
    var unc = {};
    var k, j, s, sumSq, h, rp, rm, deriv;

    var sources = [
        {
            dx: dVal,
            ref: val,
            fp: function(hh) { return calcAll(key, val + hh, T, P); },
            fm: function(hh) { return calcAll(key, val - hh, T, P); }
        },
        {
            dx: dT,
            ref: T,
            fp: function(hh) { return calcAll(key, val, T + hh, P); },
            fm: function(hh) { return calcAll(key, val, T - hh, P); }
        },
        {
            dx: dP,
            ref: P,
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
    if (!isFinite(v)) {
        el.textContent = "—";
        el.className   = "rv na";
    } else {
        el.textContent = v.toFixed(digits).replace(".", ",");
        el.className   = "rv";
    }
}

function setUnc(id, u, digits) {
    var el = $(id);
    if (isFinite(u) && u > 0) {
        el.textContent = "\u00B1 " + u.toFixed(digits).replace(".", ",");
    } else {
        el.textContent = "";
    }
}

function clearResults() {
    var j, k;
    for (j = 0; j < RKEYS.length; j++) {
        k = RKEYS[j];
        setVal("o_" + k, NaN, PRECISION);
        setUnc("u_" + k, NaN, PRECISION);
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

    /* Проверка es > P */
    var esAtT = esF(T);

    if (key === "RH" && val > 100 && esAtT >= P) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "При заданной температуре (" + T.toFixed(1) + " °С) и давлении (" +
            (P / 10).toFixed(2) + " кПа) относительная влажность не может " +
            "достигать " + val.toFixed(1) + " %, так как давление насыщенного пара (" +
            esAtT.toFixed(2) + " мбар) превышает общее давление (" +
            P.toFixed(2) + " мбар).";
        clearResults();
        return;
    }

    /* Расчёт */
    var res = calcAll(key, val, T, P);

    if (!res) {
        msg.className = "msgbox err";
        msg.innerHTML = "Невозможно рассчитать параметры при заданных входных данных. Проверьте значения.";
        clearResults();
        return;
    }

    /* Проверка e > P */
    if (res.e > P) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "При заданной температуре и давлении парциальное давление пара (" +
            res.e.toFixed(2) + " мбар) превышает общее давление (" +
            P.toFixed(2) + " мбар). Расчёт невозможен.";
        clearResults();
        return;
    }

    /* Проверка RH > 100 */
    if (res.RH > 100.05) {
        msg.className = "msgbox warn";
        msg.innerHTML =
            "Расчётная относительная влажность (" + res.RH.toFixed(2) +
            " %) превышает 100 %. При заданных параметрах воздух пересыщен — " +
            "результаты могут быть некорректны.";
        clearResults();
        return;
    }

    /* Всё ок */
    msg.className = "msgbox";
    msg.innerHTML = "";

    var j, k, d;
    for (j = 0; j < RKEYS.length; j++) {
        k = RKEYS[j];
        d = PRECISION;
        setVal("o_" + k, res[k], d);
    }

    /* Погрешности */
    var hasUnc = (dVal > 0 || dT > 0 || dP > 0);
    if (hasUnc) {
        var unc = calcUnc(key, val, dVal, T, dT, P, dP);
        for (j = 0; j < RKEYS.length; j++) {
            k = RKEYS[j];
            d = PRECISION;
            setUnc("u_" + k, unc[k], d);
        }
    } else {
        for (j = 0; j < RKEYS.length; j++) {
            k = RKEYS[j];
            setUnc("u_" + k, 0, PRECISION);
        }
    }
}

/* ══════════════════════════════════════════════════════════
   ВЫДЕЛЕНИЕ СТРОКИ РЕЗУЛЬТАТА
══════════════════════════════════════════════════════════ */

function initHighlight() {
    var rt = $("rt");
    if (!rt) return;

    var rows = rt.querySelectorAll("tr[data-k]");
    var i;

    for (i = 0; i < rows.length; i++) {
        rows[i].addEventListener("click", function() {
            if (this.classList.contains("hl")) {
                this.classList.remove("hl");
                return;
            }
            var all = rt.querySelectorAll("tr.hl");
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
}

/* ══════════════════════════════════════════════════════════
   СТЕППЕР ТОЧНОСТИ — СОБЫТИЯ
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

    // Единицы первичного параметра
    var PP_UNITS = {
        RH: "%",
        Td: "°С",
        Tf: "°С",
        Av: "г/м³",
        W:  "г/кг",
        e:  "мбар",
        Tw: "°С"
    };

    $("pp").addEventListener("change", function() {
        var u = $("pp-unit");
        if (u) u.textContent = PP_UNITS[this.value] || "";
    });

/* ══════════════════════════════════════════════════════════
   ЗАПУСК
══════════════════════════════════════════════════════════ */

document.addEventListener("DOMContentLoaded", function() {
    initEvents();
    initHighlight();
    initPrecision();
    render();
});