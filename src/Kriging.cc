// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of Kriging and support functions; currently wrapper around gstat
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */

/*
#include <defs.h>
#include <utils.h>
#include <data.h>
#include <vario.h>  
#include <glvars.h>
*/
#include <lsst/ip/diffim/Kriging.h>

using namespace std;

lsst::ip::diffim::Variogram::Variogram()
    :    
    lsst::daf::data::LsstBase(typeid(this)) {

    /* VARIOGRAM */
    v  = NULL;
    /* DATA */
    *d = NULL;
    d1 = NULL;
    d2 = NULL;
    
    /* vario.c */
    init_variogram(v);
}

void lsst::ip::diffim::Variogram::fillVariogram(std::vector<double> x, 
                                                std::vector<double> y,
                                                std::vector<double> z,
                                                std::vector<double> values,
                                                std::vector<double> variance) {
    /* data.c : read_table */
    static DPOINT current;
    static int sizeof_currentX = 0;
    int colmax = 0;

    *d = init_one_data(*d);
    *d.minX = *d.maxX = 0.0; 
    d->minY = d->maxY = 0.0;
    d->minZ = d->maxZ = 0.0;

    d->colnx = 1;
    d->colny = 2;
    d->colnz = 3;
    d->colnvalue = 4;
    d->colnvariance = 5;


    current.u.stratum = 0;
    if (sizeof_currentX == 0)
        current.X = NULL;
    colmax = MAX(d->colnx, MAX(d->colny, MAX(d->colnz, d->colnvalue)));
    colmax = MAX(colmax, MAX(d->colnvariance, d->colns));
    colmax = MAX(colmax, d->coln_id);

    if (d->colnx > 0)
        d->mode = d->mode | X_BIT_SET;
    if (d->colny > 0)
        d->mode = d->mode | Y_BIT_SET;
    if (d->colnz > 0)
        d->mode = d->mode | Z_BIT_SET;
    if (d->colns > 0) 
        d->mode = d->mode | S_BIT_SET;
    if (d->colnvalue > 0) 
        d->mode = d->mode | V_BIT_SET;
    setup_polynomial_X(d);
    
    /* pretend we read in from ascii table with x y z value */
    d->type = data_types[DATA_ASCII_TABLE];

    /* do this by hand? */
    mk_var_names(d);

    d->n_list = d->n_max = 0; 
    SET_POINT(&current);

    if (d->n_X > sizeof_currentX) { /* and therefore > 0: */
        if (sizeof_currentX == 0)
            current.X = (double *) emalloc(d->n_X * sizeof(double));
        else
            current.X = (double *) erealloc(current.X, d->n_X * sizeof(double));
        sizeof_currentX = d->n_X;
    }

    /* data.c : read_data_line */
    for (int i = 0; i < x.size(); i++) {
        current->x = x[i];
        current->y = y[i];
        current->z = z[i];
        current->attr = values[i];
        current->variance = variance[i];
        
        if (d->n_list == 0) {
            d->minX = d->maxX = current->x; 
            d->minY = d->maxY = current->y; 
            d->minZ = d->maxZ = current->z; 
            d->minvariance = d->maxvariance = current->variance;
            d->minstratum = d->maxstratum = current->u.stratum;
        } else {
            d->minX = MIN(d->minX, current->x); 
            d->maxX = MAX(d->maxX, current->x);
            d->minY = MIN(d->minY, current->y); 
            d->maxY = MAX(d->maxY, current->y); 
            d->minZ = MIN(d->minZ, current->z); 
            d->maxZ = MAX(d->maxZ, current->z);
            d->minvariance = MIN(d->minvariance, current->variance);
            d->maxvariance = MAX(d->maxvariance, current->variance);
            d->minstratum = MIN(d->minstratum, current->u.stratum);
            d->maxstratum = MAX(d->maxstratum, current->u.stratum);
        }

        push_point(d, &current);
    }
    /* data.c : read_gstat_data */
    set_norm_fns(d);
    transform_data(d); 
    average_duplicates(d);
    calc_data_mean_std(d);
    correct_strata(d);

}

void lsst::ip::diffim::Variogram::calcVariogram() {
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


void lsst::ip::diffim::Variogram::doVariogram(int nvars, METHOD m) {
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

/*
lsst::ip::diffim::Variogram::SpatialInterpolator()
    :    
    lsst::mwi::data::LsstBase(typeid(this))
{}
*/
