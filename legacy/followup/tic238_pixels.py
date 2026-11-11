"""Compare full-cutout deficit images to target and catalog-neighbor models."""
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def measure(result):
    image=np.array(result["missing"])
    yy,xx=np.indices(image.shape)
    # Fit the sampled deficit over the full 7x7 stamp, with a free flat offset.
    # This is a simple model comparison, without a calibrated TESS PRF.
    guesses=[]
    for star in result["catalog_pixels"]:
        best=None
        for width in np.linspace(.35,1.4,36):
            model=np.exp(-((xx-star["x"])**2+(yy-star["y"])**2)/(2*width**2)).ravel()
            mat=np.column_stack((np.ones(model.size),model))
            fit=np.linalg.lstsq(mat,image.ravel(),rcond=None)[0]
            if fit[1]<0:
                continue
            sse=float(np.square(image.ravel()-mat@fit).sum())
            if best is None or sse<best["sse"]:
                best=dict(tic=star["tic"],sse=sse,width=width,Tmag=star["Tmag"],peak_flux=float(fit[1]))
        if best:
            guesses.append(best)
    return dict(sector=result["sector"],ranked=sorted(guesses,key=lambda r:r["sse"]),
                peak_flux=result["peak_flux"],control_peak_max=max(result["null_max_abs_pixel"]))


if __name__=="__main__":
    rows=json.loads((HERE/"tic238_pixel_inputs.json").read_text(encoding="utf-8"))
    out=[measure(r) for r in rows if r["tic"]=="238187210"]
    (HERE/"tic238_pixels.json").write_text(json.dumps(out,indent=2))
