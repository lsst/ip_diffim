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
#include <select.h>
#include <sem.h>
#include <direct.h>
#include <fit.h>

using namespace std;

lsst::ip::diffim::Variogram::Variogram(int nVars)
    :    
    lsst::daf::data::LsstBase(typeid(this)) {

    /* DATA */
    /* data.c */
    data = (DATA **) emalloc(nVars * sizeof(DATA *));
    for (int i = 0; i < nVars; i++)
        data[i] = init_one_data(NULL);

    /* VARIOGRAM */
    /* vario.c */
    int nVgms = (nVars * (nVars + 1))/2;
    vgm = (VARIOGRAM **) emalloc(nVgms * sizeof(VARIOGRAM *));
    for (int i = 0; i < nVgms; i++)
        vgm[i] = NULL;
}

void lsst::ip::diffim::Variogram::fillData(int plane,
                                           std::vector<double> x, 
                                           std::vector<double> y,
                                           std::vector<double> z,
                                           std::vector<double> values,
                                           std::vector<double> variance) {
    /* data.c : read_table */
    static DPOINT current;
    static int sizeof_currentX = 0;
    int colmax = 0;

    data[plane]->minX = data[plane]->maxX = 0.0; 
    data[plane]->minY = data[plane]->maxY = 0.0;
    data[plane]->minZ = data[plane]->maxZ = 0.0;

    data[plane]->colnx = 1;
    data[plane]->colny = 2;
    data[plane]->colnz = 3;
    data[plane]->colnvalue = 4;
    data[plane]->colnvariance = 5;


    current.u.stratum = 0;
    if (sizeof_currentX == 0)
        current.X = NULL;
    colmax = MAX(data[plane]->colnx, MAX(data[plane]->colny, MAX(data[plane]->colnz, data[plane]->colnvalue)));
    colmax = MAX(colmax, MAX(data[plane]->colnvariance, data[plane]->colns));
    colmax = MAX(colmax, data[plane]->coln_id);

    if (data[plane]->colnx > 0)
        data[plane]->mode = data[plane]->mode | X_BIT_SET;
    if (data[plane]->colny > 0)
        data[plane]->mode = data[plane]->mode | Y_BIT_SET;
    if (data[plane]->colnz > 0)
        data[plane]->mode = data[plane]->mode | Z_BIT_SET;
    if (data[plane]->colns > 0) 
        data[plane]->mode = data[plane]->mode | S_BIT_SET;
    if (data[plane]->colnvalue > 0) 
        data[plane]->mode = data[plane]->mode | V_BIT_SET;
    setup_polynomial_X(data[plane]);
    
    /* pretend we read in from ascii table with x y z value */
    data[plane]->type = data_types[DATA_ASCII_TABLE];

    /* do this by hand? */
    mk_var_names(data[plane]);

    data[plane]->n_list = data[plane]->n_max = 0; 
    SET_POINT(&current);

    if (data[plane]->n_X > sizeof_currentX) { /* and therefore > 0: */
        if (sizeof_currentX == 0)
            current.X = (double *) emalloc(data[plane]->n_X * sizeof(double));
        else
            current.X = (double *) erealloc(current.X, data[plane]->n_X * sizeof(double));
        sizeof_currentX = data[plane]->n_X;
    }

    /* data.c : read_data_line */
    for (int i = 0; i < x.size(); i++) {
        current.x = x[i];
        current.y = y[i];
        current.z = z[i];
        current.attr = values[i];
        current.variance = variance[i];
        
        if (data[plane]->n_list == 0) {
            data[plane]->minX = data[plane]->maxX = current.x; 
            data[plane]->minY = data[plane]->maxY = current.y; 
            data[plane]->minZ = data[plane]->maxZ = current.z; 
            data[plane]->minvariance = data[plane]->maxvariance = current.variance;
            data[plane]->minstratum = data[plane]->maxstratum = current.u.stratum;
        } else {
            data[plane]->minX = MIN(data[plane]->minX, current.x); 
            data[plane]->maxX = MAX(data[plane]->maxX, current.x);
            data[plane]->minY = MIN(data[plane]->minY, current.y); 
            data[plane]->maxY = MAX(data[plane]->maxY, current.y); 
            data[plane]->minZ = MIN(data[plane]->minZ, current.z); 
            data[plane]->maxZ = MAX(data[plane]->maxZ, current.z);
            data[plane]->minvariance = MIN(data[plane]->minvariance, current.variance);
            data[plane]->maxvariance = MAX(data[plane]->maxvariance, current.variance);
            data[plane]->minstratum = MIN(data[plane]->minstratum, current.u.stratum);
            data[plane]->maxstratum = MAX(data[plane]->maxstratum, current.u.stratum);
        }

        push_point(data[plane], &current);
    }
    /* data.c : read_gstat_data */
    set_norm_fns(data[plane]);
    transform_data(data[plane]); 
    average_duplicates(data[plane]);
    calc_data_mean_std(data[plane]);
    correct_strata(data[plane]);

}

void lsst::ip::diffim::Variogram::transform_data(DATA *d) {
}
void lsst::ip::diffim::Variogram::mk_var_names(DATA *d) {
}
int lsst::ip::diffim::Variogram::average_duplicates(DATA *d) {
}
void lsst::ip::diffim::Variogram::calc_data_mean_std(DATA *d) {
}
void lsst::ip::diffim::Variogram::correct_strata(DATA *d) {
}
void lsst::ip::diffim::Variogram::ev2map(VARIOGRAM *v) {
}

//int lsst::ip::diffim::Variogram::calcVariogram() {
//    /* sem.c calc_variogram */
//    assert(v);
//    DATA *d1=NULL, *d2=NULL;
//    
//    d1 = d[v->id1];
//    d2 = d[v->id2];
//
//    if ((d1->dummy || d2->dummy)) {
//        //ErrMsg(ER_READ, "could not read sample variogram");
//        v->ev->cloud = 0;
//        v->ev->recalc = 0;
//        return 0;
//    }
//    if (d1->sel == NULL)
//        select_at(d1, NULL); /* global selection (sel = list) */
//    if (d2->sel == NULL)
//        select_at(d2, NULL);
//    
//    if (v->ev->evt == CROSSVARIOGRAM && v->ev->is_asym == -1) {
//        /* v's first time */
//        if (coordinates_are_equal(d[v->id1], d[v->id2]))
//            v->ev->pseudo = 0;
//        else
//            v->ev->pseudo = 1;
//        if (gl_sym_ev == 0)
//            v->ev->is_asym = v->ev->pseudo;
//        /* pseudo: always, else: only if set */
//        else
//            v->ev->is_asym = 0;
//    }
//    if (gl_zero_est == ZERO_DEFAULT) { /* choose a suitable default */
//        if (is_covariogram(v))
//            v->ev->zero = ZERO_SPECIAL;
//        else { /* v is variogram */
//            if (v->ev->pseudo)
//                v->ev->zero = ZERO_SPECIAL;
//            else
//                v->ev->zero = ZERO_INCLUDE;
//        }
//    } else
//        v->ev->zero = zero_int2enum(gl_zero_est);
//    
//    assert(v->ev->zero != ZERO_DEFAULT);
//    
//    fill_cutoff_width(d1, v);
//    
//    if (v->ev->map && v->ev->S_grid == NULL)
//        return -1;
//    
//    v->ev->cloud = (v->ev->iwidth <= 0.0);
//    if (v->ev->cloud &&
//        (d[v->id1]->n_sel >= MAX_NH ||  d[v->id2]->n_sel >= MAX_NH))
//        //pr_warning("observation numbers in cloud will be wrong");
//        ;
//    set_direction_values(gl_alpha, gl_beta, gl_tol_hor, gl_tol_ver);
//    
//    v->ev->is_directional = is_directional(v);
//    if (v->ev->recalc) {
//        switch (v->ev->evt) {
//        case SEMIVARIOGRAM:
//            semivariogram(d[v->id1], v->ev);
//            break;
//        case CROSSVARIOGRAM:
//            cross_variogram(d[v->id1], d[v->id2], v->ev);
//            break;
//        case COVARIOGRAM:
//            v->ev->is_asym = gl_sym_ev;
//            covariogram(d[v->id1], v->ev);
//            break;
//        case CROSSCOVARIOGRAM:
//            cross_covariogram(d[v->id1], d[v->id2], v->ev);
//            break;
//        case NOTSPECIFIED:
//        default:
//            assert(0); /* aborts */
//            break;
//        }
//    }
//    if (v->ev->map && !v->ev->S_grid)
//        ev2map(v);
//
//    return 0;
//}
//
//
//void lsst::ip::diffim::Variogram::doVariogram(int nvars, METHOD m) {
//    /* gstat.c do_variogram */
//    VARIOGRAM *vp = NULL;
//
//    if (nvars == 0)
//        return;
//    
//    /* NOTE; LTI is in utils.h */
//    for (int i = 0; i < nvars; i++) {
//        for (int j = i; j >= 0; j--) {
//            vp = get_vgm(LTI(i,j));
//            vp->id1 = j;
//            vp->id2 = i;
//            if (m == COV)
//                vp->ev->evt = (i != j) ? CROSSCOVARIOGRAM : COVARIOGRAM;
//            else
//                vp->ev->evt = (i != j) ? CROSSVARIOGRAM : SEMIVARIOGRAM;
//            if (vp->fname != NULL || o_filename != NULL) {
//                calc_variogram(vp, vp->fname ? vp->fname : o_filename);
//                if (vp->n_models > 0 && gl_fit) {
//                    vp->ev->fit = fit_int2enum(gl_fit);
//                    if (fit_variogram(vp))
//                        //pr_warning("error during variogram fit");
//                        ;
//                    else
//                        logprint_variogram(vp, 1);
//                }
//            }
//        }
//    }
//}
//
///*
//lsst::ip::diffim::Variogram::SpatialInterpolator()
//    :    
//    lsst::mwi::data::LsstBase(typeid(this))
//{}
//*/
//
