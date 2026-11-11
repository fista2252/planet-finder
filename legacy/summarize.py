"""Join pilot outputs to per-target coverage and contamination summaries."""
import csv
import json

from screen import ROOT


def load(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def neighbors(tic, depth):
    path = ROOT / "data/nearby" / f"{tic}.json"
    if not path.exists():
        return None
    rows = json.loads(path.read_text(encoding="utf-8"))
    own = next((x for x in rows if str(x["ID"]) == tic), None)
    if own is None or own["Tmag"] is None:
        return None
    result = []
    for row in rows:
        if str(row["ID"]) == tic or row["Tmag"] is None or row["dstArcSec"] > 42:
            continue
        ratio = 10 ** (-.4*(row["Tmag"]-own["Tmag"]))
        if ratio / (1+ratio) >= depth:
            result.append(row)
    return sorted(result, key=lambda x: x["Tmag"])


def main():
    s, a, v, short = map(load, ("screen.json", "aliases.json", "vetting.json", "short_period.json"))
    with (ROOT/"data/tess_ebs.csv").open(newline="") as f:
        eb = {row["tess_id"] for row in csv.DictReader(f)}
    rows = []
    for measured, bounds, vet, periodic in zip(s,a,v,short,strict=True):
        tic = measured["catalog"]["TICID"]
        assert tic == bounds["tic"] == vet["tic"] == periodic["tic"]
        sector = int(measured["catalog"]["Sector"])
        orig = measured.get("original") or {}
        raw = measured.get("original_raw") or {}
        coverage = [p for p in measured.get("sectors", []) if "n" in p]
        alternate = next((p["peaks"][0] for p in coverage if p["sector"] == sector and p["peaks"]), {})
        candidates = vet["vetting"]
        nearby = neighbors(tic, orig.get("depth", 1))
        short_fit = periodic["short_period"] or {}
        short_base = short_fit.get("training_result") or {}
        short_held = [h["result"] for h in short_fit.get("holdout", []) if h["result"]]
        short_confirm = sum(h["positive_fraction"] >= .75 and h["depth"] >= .35*short_base.get("depth", 1)
                            for h in short_held)
        grid = bounds.get("single_grid", {})
        row = dict(tic=tic, original_sector=sector, nominal_original_score=orig.get("score"),
                   original_btjd=orig.get("tc"), original_depth=orig.get("depth"), original_duration_days=orig.get("duration"),
                   nominal_raw_original_score=raw.get("score"), nominal_alternate_original_score=alternate.get("score"),
                   valid_later_sectors=len({p["sector"] for p in coverage if p["sector"] > sector}),
                   valid_other_sectors=len({p["sector"] for p in coverage if p["sector"] != sector}),
                   product_errors=sum("error" in p for p in measured.get("sectors", [])),
                   matching_later_peaks=bounds.get("compatible_peaks"),
                   period_grid_excluded=grid.get("excluded"), period_grid_tested=grid.get("tested"),
                   pairs_checked=len(candidates), pairs_with_surviving_aliases=sum(p["remaining_aliases"]>0 for p in candidates),
                   best_additional_sectors=max((len({h["sector"] for h in a["additional"]})
                                               for p in candidates for a in p["supported_aliases"]),default=0),
                   strong_short_period_holdout_count=short_confirm, short_period_days=short_fit.get("period"),
                   eb_catalog_match=tic in eb, plausible_neighbor_count=None if nearby is None else len(nearby),
                   brightest_plausible_neighbor=None if not nearby else nearby[0]["ID"],
                   brightest_plausible_neighbor_separation_arcsec=None if not nearby else nearby[0]["dstArcSec"],
                   nearby_eb_ids=";".join(str(n["ID"]) for n in nearby or [] if str(n["ID"]) in eb))
        rows.append(row)
    with (ROOT/"pilot_summary.csv").open("w",encoding="utf-8",newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
