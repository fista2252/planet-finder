# Experimental transit diagnostics

Run two optional diagnostics against the original FITS and fixed event windows in a completed workflow report.

```sh
python techniques.py --report runs/example-b/report.json --tic 238187210 --out diagnostics-238187210.json
```

By default, the command reads `.audit-cache/shared/fits/`. Use `--cache PATH` to point at the shared cache for that run. It checks each FITS digest against the report and refuses to replace an existing output. It performs no downloads or period search.

**Exposure-integrated transit shapes.** Fit a box and a trapezoid to the same quality-masked cadences. Both models use an independent linear baseline and a depth coefficient, with center and duration fixed by the report. Ingress and egress each occupy one quarter of the trapezoid duration. Both templates are averaged over the sample exposure. The output gives the two residual sums of squares and depths. A smaller residual is a morphology check. Source blending, grazing geometry and finite samples limit the comparison. This implementation uses its own fixed trapezoid, in contrast to the limb-darkened search described by [Hippke and Heller](https://arxiv.org/abs/1901.02015).

**Control-window noise scaling.** Fit the same templates to up to 32 nonoverlapping windows in the same sector, excluding the neighborhoods of the reported events. Convert the median absolute deviation of their nominal scores to a robust scale. The adjustment only decreases absolute scores and requires at least eight usable controls. The count of comparable control scores helps a reviewer distinguish an unusual event from a noisy sector. Stationarity and independence still need investigation. These scores are descriptive and have no calibrated false-alarm interpretation. See [Carter and Winn](https://arxiv.org/abs/0909.0747) for correlated noise and the limits of simple error models.

Neither diagnostic changes the original report, period acceptance, follow-up queue or transit planner. Missing coverage stays unmeasurable. Validation currently covers known injected shapes, white versus correlated noise, and missing cadence or control data. Before using either diagnostic in automated decisions, reserve independent validation fixtures and compare false alarms, recovery and wrong period counts with the existing baseline. See the [maintenance protocol](maintenance.md).
