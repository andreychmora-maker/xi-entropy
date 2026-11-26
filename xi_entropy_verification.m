/* ============================================================================
 * xi_entropy_verification.m
 * Verification & Monte Carlo Simulation Suite for:
 * "Asymptotics of xi-Entropy and Optimal Parameter Selection for
 *  Reed-Solomon Masked Cryptosystems" (Andrey L. Chmora, 2026)
 *
 * Compatible with Magma V2.29 and Magma Online Calculator (< 10 s runtime).
 * ============================================================================ */

t_start := Cputime();
SetSeed(20261126);

RR  := RealField(50);
LN2 := Log(RR!2);

/* --- Elementary Helpers & Safe Fixed-Precision Formatting --- */
function Log2(x)
    return Log(RR!x) / LN2;
end function;

function Fmt(x, d)
    if d eq 0 then
        return IntegerToString(Round(RR!x));
    end if;
    scale := 10^d;
    n_int := Round(RR!x * scale);
    sgn := "";
    if n_int lt 0 then
        sgn := "-";
        n_int := -n_int;
    end if;
    ip := n_int div scale;
    fp := n_int mod scale;
    fs := IntegerToString(fp);
    while #fs lt d do
        fs := "0" cat fs;
    end while;
    return sgn cat IntegerToString(ip) cat "." cat fs;
end function;

function h2(x)
    xr := RR!x;
    if xr le 0 or xr ge 1 then
        return RR!0;
    end if;
    return -xr * Log2(xr) - (RR!1 - xr) * Log2(RR!1 - xr);
end function;

function Ihyp(sigma, tau)
    sr := RR!sigma;
    tr := RR!tau;
    if sr le 0 or sr ge tr or (tr - sr) ge (RR!1 - tr) then
        return RR!0;
    end if;
    return h2(tr) - tr * h2(sr / tr) - (RR!1 - tr) * h2((tr - sr) / (RR!1 - tr));
end function;

function EcollBM(R, beta)
    delta := RR!1 - R;
    tau   := beta * delta / 2;
    theta := delta / 2;
    sigma := 2*tau - theta;
    if sigma le tau*tau or sigma ge tau then
        return RR!0;
    end if;
    return Ihyp(sigma, tau);
end function;

function EcollGS(R, beta)
    delta := RR!1 - R;
    tau   := beta * delta / 2;
    theta := RR!1 - Sqrt(R);
    sigma := 2*tau - theta;
    if sigma le tau*tau or sigma ge tau then
        return RR!0;
    end if;
    return Ihyp(sigma, tau);
end function;

/* --- Log-Factorial Table for Exact Finite-Length Sums (n <= 4096) --- */
MAX_N := 4096;
LogFact := [RR!0 : i in [1..MAX_N+1]]; // LogFact[k+1] = ln(k!)
for k in [1..MAX_N] do
    LogFact[k+1] := LogFact[k] + Log(RR!k);
end for;

function LogBinom(n, k)
    if k lt 0 or k gt n then
        return RR!-10^30;
    end if;
    return LogFact[n+1] - LogFact[k+1] - LogFact[n-k+1];
end function;

function StaticHxi(n, q, w_xi)
    return LogBinom(n, w_xi) / LN2 + (RR!w_xi) * Log2(RR!q - 1);
end function;

function OverThreshGrind(n, q, W, t)
    dt := W - t;
    log2_prob := (LogBinom(W, dt) - LogBinom(n, dt)) / LN2 - (RR!dt) * Log2(RR!q - 1);
    return -log2_prob;
end function;

/* Exact double sum (Eq. 14): returns H_pair (bits) and E[C | collision] */
function ExactPairCollision(n, q, w_xi, t)
    log_denom := LogBinom(n, w_xi);
    p       := RR!1 / (RR!q - 1);
    log_p   := Log(p);
    log_1mp := Log(RR!1 - p);

    total := RR!0;
    exp_c := RR!0;
    s_lo  := Max([0, 2*w_xi - n, Ceiling((2*w_xi - t) / 2)]);
    for s in [s_lo..w_xi] do
        log_ps := LogBinom(w_xi, s) + LogBinom(n - w_xi, w_xi - s) - log_denom;
        c_min  := Max(0, 2*w_xi - t - s);
        if c_min eq 0 then
            ps := Exp(log_ps);
            total +:= ps;
            exp_c +:= ps * (RR!s) * p;
        else
            for c in [c_min..s] do
                log_pc := LogBinom(s, c) + (RR!c)*log_p + (RR!(s - c))*log_1mp;
                term   := Exp(log_ps + log_pc);
                total +:= term;
                exp_c +:= term * (RR!c);
            end for;
        end if;
    end for;
    return -Log2(total), exp_c / total;
end function;

/* --- Numerical Solvers (Bisection & Golden-Section Maximization) --- */
function Bisect(f, a, b, iters)
    lo := RR!a; hi := RR!b;
    flo := f(lo);
    for i in [1..iters] do
        mid := (lo + hi) / 2;
        fmid := f(mid);
        if (flo le 0 and fmid le 0) or (flo ge 0 and fmid ge 0) then
            lo := mid; flo := fmid;
        else
            hi := mid;
        end if;
    end for;
    return (lo + hi) / 2;
end function;

function Maximize1D(f, a, b, iters)
    phi := (Sqrt(RR!5) - 1) / 2;
    lo := RR!a; hi := RR!b;
    c := hi - phi * (hi - lo);
    d := lo + phi * (hi - lo);
    fc := f(c); fd := f(d);
    for i in [1..iters] do
        if fc lt fd then
            lo := c; c := d; fc := fd;
            d := lo + phi * (hi - lo);
            fd := f(d);
        else
            hi := d; d := c; fd := fc;
            c := hi - phi * (hi - lo);
            fc := f(c);
        end if;
    end for;
    x_opt := (lo + hi) / 2;
    return x_opt, f(x_opt);
end function;

/* ============================================================================
 * PART 1: Universal Rate Cutoff R_min and Critical Budget Thresholds
 * ============================================================================ */
print "==============================================================================";
print "PART 1: Universal Rate Cutoff R_min (Theorem 1, Part 3)";
print "==============================================================================";
u_star := Bisect(func<u | u^3 + u^2 + u - 1>, RR!0, RR!1, 140);
R_min  := u_star^4;
printf "u* (root of u^3 + u^2 + u = 1) = %o\n", Fmt(u_star, 8);
printf "R_min = (u*)^4                 = %o\n", Fmt(R_min, 8);
printf "Canonical (R=0.75): beta_crit^BM = %o, beta_crit^GS = %o\n\n",
       Fmt(RR!1 / (RR!1 + Sqrt((RR!1 + RR!0.75)/2)), 5),
       Fmt(2 * (RR!1 - (RR!0.75)^(1/4)) / (RR!1 - RR!0.75), 5);

/* ============================================================================
 * PART 2: Optimal Code Rates R*(beta) and Low-Rate Trap (Table 1 & Fig. 1b)
 * ============================================================================ */
print "==============================================================================";
print "PART 2: Optimal Rates R*(beta) in the MDS/RS Limit (Table 1)";
print "==============================================================================";
print "beta   | R_min^GS | R*_BM  | E*_BM  | R*_GS  | E*_GS";
print "-------+----------+--------+--------+--------+-------";
betas_num := [65, 70, 74, 75, 80, 85, 90];
for bn in betas_num do
    beta := RR!bn / 100;
    f_root := func<R | (2*(beta*(1-R)/2) - (beta*(1-R)/2)^2) - (1 - Sqrt(R))>;
    rmin_b := Bisect(f_root, RR!1/100000, RR!999/1000, 120);
    r_bm, e_bm := Maximize1D(func<R | EcollBM(R, beta)>, RR!1/1000, RR!999/1000, 120);
    r_gs, e_gs := Maximize1D(func<R | EcollGS(R, beta)>, rmin_b + RR!1/10000, RR!999/1000, 120);
    printf "%o   | %o   | %o | %o | %o | %o\n",
           Fmt(beta, 2), Fmt(rmin_b, 4), Fmt(r_bm, 4), Fmt(e_bm, 4), Fmt(r_gs, 4), Fmt(e_gs, 4);
end for;

beta_can := RR!23 / 31;
r_gs_can, e_gs_can := Maximize1D(func<R | EcollGS(R, beta_can)>, RR!0.31, RR!0.99, 120);
printf "\nCanonical beta = 23/31 (~0.7419): R*_GS = %o (E*_GS = %o)\n",
       Fmt(r_gs_can, 4), Fmt(e_gs_can, 4);
printf "  At R=0.75: E_GS = %o bits/sym | At R=0.50: E_GS = %o bits/sym\n\n",
       Fmt(EcollGS(RR!0.75, beta_can), 4), Fmt(EcollGS(RR!0.50, beta_can), 4);

/* ============================================================================
 * PART 3: Field-Size Decoupling at Fixed n = 255 (Fig. 2 & Subsec. 4.2)
 * ============================================================================ */
print "==============================================================================";
print "PART 3: Alphabet-Size Decoupling at Fixed n=255, k=192, w_xi=23 (Fig. 2)";
print "==============================================================================";
print "m  | H_pair(BM) | E[C|BM] | H_pair(GS) | E[C|GS] | H_xi  | H_grind(W=70,GS)";
print "---+------------+---------+------------+---------+-------+-----------------";
for m in [8, 9, 10, 12, 14, 16, 20, 24, 32] do
    q := 2^m;
    h_bm, ec_bm := ExactPairCollision(255, q, 23, 31);
    h_gs, ec_gs := ExactPairCollision(255, q, 23, 34);
    h_st := StaticHxi(255, q, 23);
    h_gr := OverThreshGrind(255, q, 70, 34);
    printf "%2o | %10o | %7o | %10o | %7o | %5o | %16o\n",
           m, Fmt(h_bm, 2), Fmt(ec_bm, 2), Fmt(h_gs, 2), Fmt(ec_gs, 2),
           Fmt(h_st, 1), Fmt(h_gr, 1);
end for;
print "";

/* ============================================================================
 * PART 4: Finite-Length Profiles, Asymptotic Accuracy & Saddle Point c* (Table 2)
 * ============================================================================ */
print "==============================================================================";
print "PART 4: Finite-Length Profiles & Low-Rate Trap Verification (Table 2 & Fig. 3)";
print "==============================================================================";
print "Profile              | [n,k,d]          | (w_xi,wa)|  H_xi  | BM(asymp)     | GS(asymp)     | E[C|GS] c*_GS";
print "---------------------+------------------+----------+--------+---------------+---------------+--------------";

profiles := [
    <"Compact Hybrid",     8,  192,  23>,
    <"Tight Hybrid",       8,  192,  25>,
    <"Low-Rate Trap m=8",  8,  128,  47>,
    <"Mid-Tier",           9,  384,  47>,
    <"Mid-Tier Tight",     9,  384,  51>,
    <"Standalone 128b",   10,  768,  95>,
    <"Standalone GS-128", 10,  768, 100>,
    <"Low-Rate Trap m=10",10,  512, 191>,
    <"High-Margin",       11, 1536, 191>,
    <"Split-Limb",        12, 3072, 383>
];

for prof in profiles do
    name := prof[1]; m := prof[2]; k := prof[3]; w_xi := prof[4];
    q := 2^m; n := q - 1; d := n - k + 1;
    aleph := (d - 1) div 2;
    J := Floor(RR!n - Sqrt(RR!(n * (k - 1))));
    w_attr := aleph - w_xi;

    h_st := StaticHxi(n, q, w_xi);
    bm, ec_bm := ExactPairCollision(n, q, w_xi, aleph);
    gs, ec_gs := ExactPairCollision(n, q, w_xi, J);

    tau := RR!w_xi / n;
    th_bm := RR!aleph / n;
    th_gs := RR!J / n;
    a_bm := (RR!n) * Ihyp(2*tau - th_bm, tau);
    a_gs := (RR!n) * Ihyp(2*tau - th_gs, tau);
    c_star_gs := (1 - th_gs) * (2*tau - th_gs)^2 / ((th_gs - tau)^2);

    code_str := Sprintf("[%o,%o,%o]", n, k, d);
    w_str    := Sprintf("(%o,%o)", w_xi, w_attr);
    bm_str   := Fmt(bm, 2) cat "(" cat Fmt(a_bm, 2) cat ")";
    gs_str   := Fmt(gs, 2) cat "(" cat Fmt(a_gs, 2) cat ")";
    printf "%-20o | %-16o | %-8o | %6o | %-13o | %-13o | %5o   %5o\n",
           name, code_str, w_str, Fmt(h_st, 1), bm_str, gs_str, Fmt(ec_gs, 2), Fmt(c_star_gs, 2);
end for;
print "";

/* ============================================================================
 * PART 5: Live Algebraic RS[255, 192, 64]_{256} Masking & BM Decoding Test
 * ============================================================================ */
print "==============================================================================";
printf "PART 5: Live RS[255, 192, 64]_{256} Masking & Blind Decoding (%o frames)\n", N_LIVE;
print "==============================================================================";
F256<theta> := GF(256);
F_nz := [x : x in F256 | x ne 0];
C_RS := ReedSolomonCode(F256, 64);
G_RS := GeneratorMatrix(C_RS);
V255 := AmbientSpace(C_RS);
K192 := VectorSpace(F256, Dimension(C_RS));

function SampleSphereVec(n, w)
    v := Zero(V255);
    idx := [1..n];
    for i in [1..w] do
        j := Random(i, n);
        tmp := idx[i]; idx[i] := idx[j]; idx[j] := tmp;
        v[idx[i]] := F_nz[Random(1, 255)];
    end for;
    return v;
end function;

mask_static := Random(V255); // represents m + (w || r_S)
N_LIVE := 55000;
ok_in_cohort   := 0;
ok_xsession_bm := 0; // true bounded-distance decodings: dec in C_RS and wt(sum - dec) <= 31
ok_xsession_eq := 0; // exact recovery of cw + prev_cw

prev_frame := Zero(V255);
prev_cw    := Zero(V255);
for trial in [1..N_LIVE] do
    msg  := Random(K192);
    cw   := msg * G_RS;
    xi   := SampleSphereVec(255, 23);
    c_tx := cw + mask_static + xi;

    // Authorized receiver with worst-case fuzzy discrepancy d_H = 8
    delta_r := SampleSphereVec(255, 8);
    rx := c_tx + mask_static + delta_r;
    ok, dec := Decode(C_RS, rx);
    if ok and (dec in C_RS) and (Weight(rx - dec) le 31) and (dec eq cw) then
        ok_in_cohort +:= 1;
    end if;

    // Cross-session linear addition attack (c'_1 + c'_2)
    if trial gt 1 then
        sum_frame := c_tx + prev_frame;
        ok_x, dec_x := Decode(C_RS, sum_frame);
        if ok_x and (dec_x in C_RS) and (Weight(sum_frame - dec_x) le 31) then
            ok_xsession_bm +:= 1;
            if dec_x eq (cw + prev_cw) then
                ok_xsession_eq +:= 1;
            end if;
        end if;
    end if;
    prev_frame := c_tx;
    prev_cw    := cw;
end for;
printf "Authorized in-cohort decodings (d_H = 8, w_xi = 23):     %o / %o\n", ok_in_cohort, N_LIVE;
printf "Cross-session valid BM decodings (wt(sum - dec) <= 31):  %o / %o\n", ok_xsession_bm, N_LIVE - 1;
printf "Cross-session exact recovery (dec == cw_1 + cw_2):       %o / %o\n\n", ok_xsession_eq, N_LIVE - 1;

/* ============================================================================
 * PART 6: Empirical Monte Carlo Demonstration of the Low-Rate Trap (GF(256))
 *   RS[255, 192, 64] (R=0.75, w_xi=23, J=34) vs RS[255, 128, 128] (R=0.50, w_xi=47, J=75)
 * ============================================================================ */
print "==============================================================================";
print "PART 6: Empirical Monte Carlo Simulation of the Low-Rate Trap (30,000 trials)";
print "==============================================================================";

procedure RunMonteCarlo(label, n, q, w_xi, aleph, J, N_trials)
    idx1 := [1..n];
    idx2 := [1..n];
    val1 := [0 : i in [1..n]];

    sum_wt  := 0;
    min_wt  := 2*w_xi;
    bm_hits := 0;
    gs_hits := 0;

    for tr in [1..N_trials] do
        // Partial Fisher-Yates for xi_1
        for i in [1..w_xi] do
            j := Random(i, n);
            tmp := idx1[i]; idx1[i] := idx1[j]; idx1[j] := tmp;
            val1[idx1[i]] := Random(1, q - 1);
        end for;

        // Partial Fisher-Yates for xi_2 and overlap/cancellation count
        s_cnt := 0;
        c_cnt := 0;
        for i in [1..w_xi] do
            j := Random(i, n);
            tmp := idx2[i]; idx2[i] := idx2[j]; idx2[j] := tmp;
            pos := idx2[i];
            if val1[pos] ne 0 then
                s_cnt +:= 1;
                if val1[pos] eq Random(1, q - 1) then
                    c_cnt +:= 1;
                end if;
            end if;
        end for;

        // Clear val1 entries cleanly
        for i in [1..w_xi] do
            val1[idx1[i]] := 0;
        end for;

        wt := 2*w_xi - s_cnt - c_cnt;
        sum_wt +:= wt;
        if wt lt min_wt then min_wt := wt; end if;
        if wt le aleph then bm_hits +:= 1; end if;
        if wt le J then gs_hits +:= 1; end if;
    end for;

    emp_mean  := RR!sum_wt / N_trials;
    theo_mean := RR!(2*w_xi) - (RR!(w_xi^2) / n) * (RR!q / (q - 1));
    printf "%-24o | emp E[wt]=%o (theo %o) | min wt=%2o | BM(<=%2o)=%o | GS(<=%2o)=%o\n",
           label, Fmt(emp_mean, 2), Fmt(theo_mean, 2), min_wt, aleph, bm_hits, J, gs_hits;
end procedure;

N_MC := 30000;
RunMonteCarlo("R=0.75 [255,192,64]", 255, 256, 23, 31, 34, N_MC);
RunMonteCarlo("R=0.50 [255,128,128]", 255, 256, 47, 63, 75, N_MC);

printf "\nTotal execution time: %o s\n", Fmt(Cputime(t_start), 3);
