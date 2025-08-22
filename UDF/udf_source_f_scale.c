/* udf_source_f_scale.c
   UDF to apply a volumetric energy source in solid cells that emulates scaling
   the convective heat transfer by factor f_ref computed externally (MATLAB).
   The UDF reads results/f_ref.txt for:
     f_ref <value>
     h_base <value>    (optional - typical h in W/m2K, if not present a default is used)
     Tfluid <value>    (optional - ambient/fluid reference temperature, default 300K)
   The source returned is S = k * (T_cell - Tfluid), with
     k = (f_ref - 1) * h_base / delta
   where delta ~ characteristic length ~ cbrt(cell_volume). Units: W/m3.
   Author: Daniel Arteaga
   Date: 2025
*/

#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static real f_ref_global = 1.0;
static real h_base_global = 40.0;    /* default nominal h [W/m2K] if not provided */
static real Tfluid_global = 293.0;   /* default fluid reference temperature [K] */
static int params_read = 0;

/* Read parameters from results/f_ref.txt in working directory */
DEFINE_INIT(read_fref_file, domain)
{
    const char *fname = "results/f_ref.txt";
    FILE *fp = fopen(fname, "r");
    if (fp != NULL)
    {
        char key[128];
        double val;
        /* parse simple "key value" pairs */
        while (fscanf(fp, " %127s %lf", key, &val) == 2)
        {
            if (strcmp(key, "f_ref") == 0)
            {
                f_ref_global = (real) val;
            }
            else if (strcmp(key, "h_base") == 0)
            {
                h_base_global = (real) val;
            }
            else if (strcmp(key, "Tfluid") == 0 || strcmp(key, "T_inf") == 0)
            {
                Tfluid_global = (real) val;
            }
        }
        fclose(fp);
        params_read = 1;
        Message("UDF init: read %s -> f_ref=%g, h_base=%g, Tfluid=%g\n",
                fname, (double)f_ref_global, (double)h_base_global, (double)Tfluid_global);
    }
    else
    {
        params_read = 0;
        Message("UDF init: file '%s' not found. Using defaults f_ref=%g, h_base=%g, Tfluid=%g\n",
                fname, (double)f_ref_global, (double)h_base_global, (double)Tfluid_global);
    }
}

/* DEFINE_SOURCE to be hooked to the ENERGY equation in the solid cell zone(s).
   eqn: equation index (ignored, required by macro). */
DEFINE_SOURCE(solid_energy_source, c, t, dS, eqn)
{
    real Vcell = C_VOLUME(c,t);           /* cell volume [m^3] */
    real Tcell = C_T(c,t);                /* cell temperature [K] */
    real delta;                           /* characteristic length [m] */
    real k_surf;                          /* coefficient (W/m2K) * (f_ref -1) */
    real k_vol;                           /* volumetric coefficient (W/m3K) */
    real src;                             /* volumetric source [W/m3] */

    /* characteristic length: cube root of cell volume (approx thickness) */
    if (Vcell > 0.0)
        delta = cbrt(Vcell);
    else
        delta = 1e-3; /* fallback */

    /* effective surface coefficient from f_ref */
    k_surf = (f_ref_global - 1.0) * h_base_global; /* W/m2K */
    /* convert to volumetric coefficient (approx) */
    k_vol = k_surf / delta; /* W/m3K  (approximation) */

    /* source (volumetric): S = k_vol * (Tcell - Tfluid) */
    src = k_vol * (Tcell - Tfluid_global);

    /* linearization term dS/dT for solver (Jacobian) */
    dS[eqn] = k_vol;

    return src;
}