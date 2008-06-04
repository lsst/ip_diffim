// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of Kriging and support functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */

#include <vario.h>  
#include <glvars.h>

using namespace std;

lsst::ip::diffim::Variogram::Variogram()
    :    
    lsst::mwi::data::LsstBase(typeid(this)) {
    
    /* vario.c */
    init_variogram(v);
}

lsst::ip::diffim::Variogram::fillVariogram() {
}

lsst::ip::diffim::Variogram::calcVariogram() {
    /* sem.c */
    assert(v);
    
    d1 = d[v->id1];
    d2 = d[v->id2];
    if ((d1->dummy || d2->dummy)) {
        ErrMsg(ER_READ, "could not read sample variogram");
        v->ev->cloud = 0;
        v->ev->recalc = 0;
        return 0;
    }
    if (d1->sel == NULL)
        select_at(d1, NULL); /* global selection (sel = list) */
    if (d2->sel == NULL)
        select_at(d2, NULL);
    
    if (v->ev->evt == CROSSVARIOGRAM && v->ev->is_asym == -1) {
        /* v's first time */
        if (coordinates_are_equal(d[v->id1], d[v->id2]))
            v->ev->pseudo = 0;
        else
            v->ev->pseudo = 1;
        if (gl_sym_ev == 0)
            v->ev->is_asym = v->ev->pseudo;
        /* pseudo: always, else: only if set */
        else
            v->ev->is_asym = 0;
    }
    if (gl_zero_est == ZERO_DEFAULT) { /* choose a suitable default */
        if (is_covariogram(v))
            v->ev->zero = ZERO_SPECIAL;
        else { /* v is variogram */
            if (v->ev->pseudo)
                v->ev->zero = ZERO_SPECIAL;
            else
                v->ev->zero = ZERO_INCLUDE;
        }
    } else
        v->ev->zero = zero_int2enum(gl_zero_est);
    
    assert(v->ev->zero != ZERO_DEFAULT);
    
    fill_cutoff_width(d1, v);
    
    if (v->ev->map && v->ev->S_grid == NULL)
        return -1;
    
    v->ev->cloud = (v->ev->iwidth <= 0.0);
    if (v->ev->cloud &&
        (d[v->id1]->n_sel >= MAX_NH ||  d[v->id2]->n_sel >= MAX_NH))
        pr_warning("observation numbers in cloud will be wrong");
    set_direction_values(gl_alpha, gl_beta, gl_tol_hor, gl_tol_ver);
    
    v->ev->is_directional = is_directional(v);
    if (v->ev->recalc) {
        switch (v->ev->evt) {
        case SEMIVARIOGRAM:
            semivariogram(d[v->id1], v->ev);
            break;
        case CROSSVARIOGRAM:
            cross_variogram(d[v->id1], d[v->id2], v->ev);
            break;
        case COVARIOGRAM:
            v->ev->is_asym = gl_sym_ev;
            covariogram(d[v->id1], v->ev);
            break;
        case CROSSCOVARIOGRAM:
            cross_covariogram(d[v->id1], d[v->id2], v->ev);
            break;
        case NOTSPECIFIED:
        default:
            assert(0); /* aborts */
            break;
        }
    }
    if (v->ev->map && !v->ev->S_grid)
        ev2map(v);

    return 0;
}


lsst::ip::diffim::Variogram::doVariogram(int nvars, METHOD m) {
    if (nvars == 0)
        return;
    
    /* NOTE; LTI is in utils.h */
    for (int i = 0; i < nvars; i++) {
        for (int j = i; j >= 0; j--) {
            vp = get_vgm(LTI(i,j));
            vp->id1 = j;
            vp->id2 = i;
            if (m == COV)
                vp->ev->evt = (i != j) ? CROSSCOVARIOGRAM : COVARIOGRAM;
            else
                vp->ev->evt = (i != j) ? CROSSVARIOGRAM : SEMIVARIOGRAM;
            if (vp->fname != NULL || o_filename != NULL) {
                calc_variogram(vp, vp->fname ? vp->fname : o_filename);
                if (vp->n_models > 0 && gl_fit) {
                    vp->ev->fit = fit_int2enum(gl_fit);
                    if (fit_variogram(vp))
                        pr_warning("error during variogram fit");
                    else
                        logprint_variogram(vp, 1);
                }
            }
        }
    }
}

lsst::ip::diffim::Variogram::SpatialInterpolator()
    :    
    lsst::mwi::data::LsstBase(typeid(this))
{}



