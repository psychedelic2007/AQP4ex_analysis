
/* ---- profile-likelihood CI for the global omega (added; FEL site tests skipped) ---- */
prof.lf = fel.final_partitioned_mg_results[terms.likelihood_function];
prof.omega_id = ((fel.final_partitioned_mg_results[terms.global])["non-synonymous/synonymous rate ratio for *test*"])[terms.id];
fprintf (stdout, "\nPROFILE FEL_MLE ", Eval(prof.omega_id), " ", fel.final_partitioned_mg_results[terms.fit.log_likelihood], "\n");
Optimize (prof.res0, ^prof.lf);
prof.mle  = Eval(prof.omega_id);
prof.Lmax = prof.res0[1][0];
fprintf (stdout, "PROFILE MLE ", prof.mle, " ", prof.Lmax, "\n");
prof.target = prof.Lmax - 1.920729410347062;
Export (prof.saved, ^prof.lf);

function prof.eval (w) {
    ExecuteCommands (prof.saved);
    ExecuteCommands (prof.omega_id + ":=" + w);
    Optimize (prof.r, ^prof.lf);
    fprintf (stdout, "PROFILE POINT ", w, " ", prof.r[1][0], "\n");
    return prof.r[1][0];
}

function prof.bound (inside, step) {
    outside = inside + step;
    while (prof.eval (outside) > prof.target) {
        inside = outside; step = step * 2; outside = inside + step;
        if (outside <= 0) { outside = 1e-6; break; }
    }
    for (it = 0; it < 12; it += 1) {
        mid = (inside + outside) / 2;
        if (prof.eval (mid) > prof.target) { inside = mid; } else { outside = mid; }
    }
    return (inside + outside) / 2;
}

prof.lo = prof.bound (prof.mle, -0.25 * prof.mle);
prof.hi = prof.bound (prof.mle,  0.25 * prof.mle);
fprintf (stdout, "PROFILE CI ", prof.lo, " ", prof.mle, " ", prof.hi, "\n");
return 0;
