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

#include <math.h>

#include <select.h>
#include <sem.h>
#include <direct.h>
#include <fit.h>
#include <random.h>
#include <mapio.h>

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
    /* data.c */
    Double_index *values;
    double nscore, q;
    FILE *f = NULL;
    int i, j, tie_length, tie_total = 0;
    DPOINT *p;
    
    if (d->log) {
        for (i = 0; i < d->n_list; i++) {
            p = d->list[i];
            if (p->attr <= 0.0) 
                //ErrMsg(ER_IMPOSVAL, "log of non-positive value");
                ;
            p->attr = log(p->attr);
        }
    } else if (! is_mv_double(&(d->Icutoff))) {
        for (i = 0; i < d->n_list; i++) {
            p = d->list[i];
            if (p->attr <= d->Icutoff)
                p->attr = 1.0;
            else
                p->attr = 0.0;
        }
    } else if (d->nscore_table) {
        if (strlen(d->nscore_table) > 0) {
            f = efopen(d->nscore_table, "w");
            fprintf(f, "normal scores for %s\n", d->variable);
            fprintf(f, "2\nobserved value\nnormal score\n");
        }
        /* ACB */
        values = (Double_index *) emalloc(d->n_list * sizeof(Double_index));
        for (i = 0; i < d->n_list; i++) {
            values[i].d = d->list[i]->attr; /* */
            values[i].index = i; /* unique identifiers */
        }
        qsort((void *) values, (size_t) d->n_list, sizeof(Double_index),
              (int CDECL (*)(const void *, const void *)) double_index_cmp);
        for (i = 0; i < d->n_list; i += tie_length) { /* ignore ties: */
            /* assume they're all ties, with min run length 1: */
            tie_length = 1;
            while (i + tie_length < d->n_list && 
                   values[i].d == values[i+tie_length].d)
                tie_length++;
            q = (i + 0.5 * tie_length) / d->n_list;
            nscore = q_normal(q);
            for (j = 0; j < tie_length; j++) {
                if (f)
                    fprintf(f, "%g %g\n", d->list[values[i+j].index]->attr, 
                            nscore);
                /* transform: */
                d->list[values[i+j].index]->attr = nscore;
            }
            tie_total += (tie_length - 1);
        }
        efree(values);
        if (f)
            efclose(f);
        if (tie_total)
            //pr_warning("tied normal scores are assigned to tied data (%.2g%%)", 
            //100.0 * tie_total / d->n_list);
            ;
    }
}

void lsst::ip::diffim::Variogram::mk_var_names(DATA *d) {
    /* data.c */
    char tmp1[100], tmp2[100]; 
    /* make variable names, if not read from EAS header: */
    if (d->variable == NULL || d->variable[0] == '\0')
        sprintf(tmp1, "col[%d]", d->colnvalue);
    else {
        if (strlen(d->variable) > 50)
            d->variable[50] = '\0';
        sprintf(tmp1, "%s", d->variable);
    }
    if (d->log)
        sprintf(tmp2, "log(%s)", tmp1);
    else if (!is_mv_double(&(d->Icutoff)))
        sprintf(tmp2, "I(%s,%g)", tmp1, d->Icutoff);
    else if (d->Category)
        sprintf(tmp2, "I(%s='%s')", tmp1, d->Category);
    else if (d->nscore_table != NULL)
        sprintf(tmp2, "%s (nscores)", tmp1);
    else /* no transform, copy: */
        sprintf(tmp2, "%s", tmp1);
    d->variable = string_dup(tmp2);
    if (d->x_coord == NULL && d->colnx) {
        sprintf(tmp1, "x_%d", d->colnx);
        d->x_coord = string_dup(tmp1);
    }
    if (d->y_coord == NULL && d->colny) {
        sprintf(tmp1, "y_%d", d->colny);
        d->y_coord = string_dup(tmp1);
    }
    if (d->z_coord == NULL && d->colnz)  {
        sprintf(tmp1, "z_%d", d->colnz);
        d->z_coord = string_dup(tmp1);
    }
    if (d->V_coord == NULL && d->colnvariance)  {
        sprintf(tmp1, "var_%d", d->colnvariance);
        d->V_coord = string_dup(tmp1);
    }
    if (d->s_coord == NULL && d->colns) {
        sprintf(tmp1, "stra_%d", d->colns);
        d->s_coord = string_dup(tmp1);
    }
    if (d->id_name == NULL && d->coln_id) {
        sprintf(tmp1, "ID_%d", d->colns);
        d->id_name = string_dup(tmp1);
    }
}

int lsst::ip::diffim::Variogram::average_duplicates(DATA *d) {
    /* data.c */

    /*
     * average duplicate records (equal coordinates) 
     * and correct number of records;
     * new variance = (sum variances)/(n * n)
     */
    int n_tot = 0, n_dupl, i, j;
    double dzero2;
    
    if (! d->average)
        return 0;
    dzero2 = gl_zero * gl_zero;
    for (i = 0; i < d->n_list; i++) {
        n_dupl = 0;
        j = i + 1;
        while (j < d->n_list) { /* can't use a for() loop here */
            if (d->pp_norm2(d->list[i], d->list[j]) < dzero2) {
                d->list[i]->attr += d->list[j]->attr;
                d->list[i]->variance += d->list[j]->variance;
                pop_point(d, j); /* decrements d->n_list */
                n_dupl++;
            } else 
                j++;
        } /* while j */
        if (n_dupl > 0) {
            n_tot += n_dupl;
            d->list[i]->attr /= (n_dupl + 1.0);
            d->list[i]->variance /= SQR(n_dupl + 1.0);
        }
    }
    return d->n_averaged = n_tot;
}

void lsst::ip::diffim::Variogram::calc_data_mean_std(DATA *d) {
    /* data.c */
    
    /* 
     * Calculates fields mean and std of d with mean and standard dev. (/(n-1))
     */
    int i;
    
    if (d->standard == 2) { /* we did this already.. */
        for (i = 0; i < d->n_list; i++)
            d->list[i]->attr *= d->std;
    }
    
    d->mean = 0.0; 
    d->std = 0.0;
    if (d->n_list <= 0) {
        //pr_warning("calc_data_mean_std: n_list <= 0: %d", d->n_list);
        return;
    }
    for (i = 0; i < d->n_list; i++) 
        d->mean += d->list[i]->attr;
    d->mean /= d->n_list;
    if (d->n_list == 1)
        return;
    for (i = 0; i < d->n_list; i++) 
        d->std += SQR(d->list[i]->attr - d->mean);
    d->std = sqrt((d->std)/(d->n_list - 1));
    
    if (d->standard > 0) {
        for (i = 0; i < d->n_list; i++)
            d->list[i]->attr /= d->std;
        d->standard = 2;
    }
    return;
}

void lsst::ip::diffim::Variogram::correct_strata(DATA *d) {
    /* data.c */

    int i;
    
    if (d->colns == 0)
        return;
    for (i = 0; i < d->n_list; i++)
        d->list[i]->u.stratum -= d->minstratum;
}

void lsst::ip::diffim::Variogram::ev2map(VARIOGRAM *v) {
    /* sem.c */
    
    GRIDMAP *m1 = NULL, *m2 = NULL;
    unsigned int row, col, i;
    SAMPLE_VGM *ev;
    
    if (v->fname == NULL)
        return;
    ev = v->ev;
    m1 = map_dup(v->fname, ev->map);
    if (v->fname2 != NULL)
        m2 = map_dup(v->fname2, ev->map);
    for (row = i = 0; row < m1->rows; row++) {
        for (col = 0; col < m1->cols; col++) {
            if (ev->nh[i] > 0)
                map_put_cell(m1, row, col, ev->gamma[i]);
            if (m2 != NULL)
                map_put_cell(m2, row, col, 1.0 * ev->nh[i]);
            i++;
        }
    }
    m1->write(m1);
    if (m2 != NULL)
        m2->write(m2);
    return;
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
