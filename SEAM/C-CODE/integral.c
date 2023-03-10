
#include<stdio.h>
#include<math.h>

#include "cgammalike.c"

actr_params get_actr_params(swift_parameters * params) {
	actr_params aparams = {1, val(params, F), val(params, G), val(params, ans), val(params, mas), val(params, d), val(params, match_penalty)};
	return aparams;
}

typedef struct {
	double ** hist;
	int ar_size;
	int length;
} W_hist;

typedef struct {
	FILE * fixation_sequence;
	FILE * activation_history;
	FILE * events;
} swift_output_files;

typedef struct swift_run swift_run;

struct swift_run {
	int * prev_n_count;
	int * n_count;
	int * N_count;
	int * states;
	int N;
	int is_fitting;
	double t;
	int gaze_word;
	int prev_gaze_word;
	int fix_saccade_target;
	double gaze_letter;
	swift_corpus * corpus;
	swift_parameters * params;
	RANSEED_TYPE seed;
	int s;
	double * aa;
	double * procrate;
	double dist;
	int saccade_target;
	double * view;
	double * border;
	double * ptar;
	double * W;
	struct {
		W_hist glob, lab, nlab, sacc;
	} w_hist_by_stage;
	int * len;
	int is_mislocated, is_refixation;
	int canc;
	double dt;
	struct {
        int * word_waits_for_retrieval;
        double * word_processing_block_times;
        double * retrieval_share;
        double R_count;
        actr_trial * atrial;
        actr_retrieval_result retrieval_result;
        int current_retrieval_id;
        int current_retrieval_trigger;
        double current_retrieval_started;
        double current_retrieval_ends;
        actr_retrieval_item current_retrieval;
        int last_retrieval_at;
	} actr;
	struct {
		void (*transition_rates_calculated)(struct swift_run *);
		void (*state_selected)(struct swift_run *, int state);
		void (*target_selected)(struct swift_run *, int target);
		void (*saccade_executed)(struct swift_run *);
		void (*counters_propagated)(struct swift_run *);
	} handlers;
	swift_output_files files;
};

double sum(double * summands, int N) {
	int i;
	double sum = 0.0;
	for(i = 1; i <= N; i++) sum += summands[i];
	return sum;
}

int sample_boolean(double threshold, RANSEED_TYPE * seed) {
	return ran1(seed) < threshold;
}

int linear_sample_normalized(double * weights, double sum, int N, RANSEED_TYPE * seed) {
	int i = 1;
	double test = 0.0;
	double cutoff = ran1(seed) * sum;
	while((test += weights[i]) < cutoff && i < N) i++;
	return i;
}

int log_linear_sample_normalized(double * weights, double log_sum, int N, RANSEED_TYPE * seed) {
	int i = 1;
	double test = -INFINITY;
	double cutoff = log(ran1(seed)) + log_sum;
	while((test = logaddexp(test, weights[i])) < cutoff && i < N) i++;
	return i;
}

int linear_sample(double * weights, int N, RANSEED_TYPE * seed) {
	return linear_sample_normalized(weights, sum(weights, N), N, seed);
}

int log_linear_sample(double * weights, int N, RANSEED_TYPE * seed) {
	return log_linear_sample_normalized(weights, logsumexp(weights, N), N, seed);
}



#define log_handler(trial, event, ...) ((trial)->handlers.event == NULL ? NULL : (trial)->handlers.event((trial) __VA_OPT__(,) __VA_ARGS__))

#define log_event(trial, fmt, ...) \
	if(trial->files.events != NULL) { \
		_Pragma("omp critical(event_file)") \
		{ \
			fprintf(trial->files.events, "%d\t%.2lf\t", trial->s, trial->t); \
			fprintf(trial->files.events, fmt __VA_OPT__(,) __VA_ARGS__); \
			fputs("\n", trial->files.events); \
		} \
	}

void log_activation_to_file(swift_run * trial) {
	if(trial->files.activation_history != NULL) {
		#pragma omp critical(activation_history_file)
		{
			fprintf(trial->files.activation_history, "%d\t%.2lf\t", trial->s, trial->t);
			for(int i = 1; i <= trial->N + 4; i++) {
				fprintf(trial->files.activation_history, "%d ", trial->n_count[i]);
			}
			fputs("\n", trial->files.activation_history);
		}
	}
}



typedef struct {
	double fixation_location, fixation_duration, saccade_duration;
} swift_likelihood;

void init_w_hist(W_hist * w_hist) {
	if(w_hist->hist != NULL) {
		w_hist->length = 0;
	} else {
		w_hist->length = 0;
		w_hist->ar_size = 100;
		w_hist->hist = matrix(double, 2, w_hist->ar_size);
	}
}

void clear_w_hist(W_hist w_hist) {
	if(w_hist.hist != NULL) free_matrix(double, w_hist.hist);
	w_hist.ar_size = 0;
	w_hist.length = 0;
}

void append_w_hist(W_hist * w_hist, double a, double b) {
	int found = -1, i;
	for(i = 1; i <= w_hist->length; i++) {
		if(b >= w_hist->hist[2][i] - DBL_EPSILON && b <= w_hist->hist[2][i] + DBL_EPSILON) {
			found = i;
			break;
		}
	}
	if(found != -1) {
		w_hist->hist[1][found] += a;
	} else {
		if(w_hist->length >= w_hist->ar_size) {
			w_hist->hist = resize_matrix(double, w_hist->hist, 2, 2, w_hist->ar_size, w_hist->ar_size*2);
			w_hist->ar_size*=2;
		}
		w_hist->length++;
		w_hist->hist[1][w_hist->length] = a;
		w_hist->hist[2][w_hist->length] = b;
	}
}

void transfer_w_hist(W_hist * from, W_hist * to) {
	int i;
	for(i = 1; i <= from->length; i++) {
		append_w_hist(to, from->hist[1][i], from->hist[2][i]);
	}
}

swift_run * new_swift_trial(swift_model * model, int s, unsigned int seed) {
	int i, j;
	int N = nwords(model->corpus, s);
	swift_run * ret = malloc(sizeof(swift_run));
	ret->prev_n_count = vector(int, N+4);
	ret->n_count = vector(int, N+4);
	ret->aa = vector(double, N);
	ret->procrate = vector(double, N);
	ret->N_count = vector(int, N+4);
	ret->ptar = vector(double, N);
	ret->W = vector(double, N+4);
	ret->s = s;
	ret->handlers.transition_rates_calculated = NULL;
	ret->handlers.state_selected = NULL;
	ret->handlers.target_selected = NULL;
	ret->handlers.saccade_executed = NULL;
	ret->handlers.counters_propagated = NULL;
	ret->w_hist_by_stage.glob.hist = NULL;
	ret->w_hist_by_stage.lab.hist = NULL;
	ret->w_hist_by_stage.nlab.hist = NULL;
	ret->w_hist_by_stage.sacc.hist = NULL;
	ret->files.fixation_sequence = NULL;
	ret->files.activation_history = NULL;
	ret->files.events = NULL;
	ret->is_mislocated = 0;
	ret->is_refixation = 0;
	ret->canc = 0;
    ret->N_count[1] = (int) model->params->cord;
    ret->N_count[2] = (int) model->params->lord;
    ret->N_count[3] = (int) model->params->nord;
    ret->N_count[4] = (int) model->params->xord;
    double maxfreq = 0.0;
    for(i=1;i<=nsentences(model->corpus);i++) {
        for(j=1;j<=nwords(model->corpus, i);j++) {
            if(word_prop(model->corpus, i, j, freq) > maxfreq) {
                maxfreq = word_prop(model->corpus, i, j, freq);
            }
        }
    }
    initSeed(seed ? seed : ranint(&model->seed), &ret->seed);
	ret->states = vector(int, N+4);
	ret->states[1] = 1;
	ret->n_count[1] = (int) (model->params->msac0 * ret->N_count[1]);
	ret->t = 0.0;
	ret->dist = 0.0;
	ret->gaze_word = 1;
	ret->gaze_letter = 1.0 + 0.5 * word_prop(model->corpus, s, 1, nl);
	ret->N = N;
	ret->is_fitting = 0;
	double logf, lex;
	for(i = 1; i <= N; i++) {
		ret->states[4+i] = STATE_LEXICAL;
		logf = log(word_prop(model->corpus, s, i, freq))/log(maxfreq);
		if ( logf<0.0 )  logf = 0.0;
		lex = 1.0 - model->params->beta*logf;   // word frequency effect
		ret->N_count[4+i] = (int) ceil(model->params->aord*lex);  // proportion of maximum activation aord
	}

    ret->view = vector(double, N);
    ret->border = vector(double, N);
    ret->len = vector(int, N);
    double b = 0.0;
    for(i = 1; i <= N; b = ret->border[i], i++) {
    	ret->len[i] = word_prop(model->corpus, s, i, nl);
    	ret->view[i] = b + model->params->iws + 0.5 * ret->len[i];
    	ret->border[i] = b + model->params->iws + ret->len[i];
    }
	ret->corpus = model->corpus;
	ret->params = model->params;


	ret->actr.word_waits_for_retrieval = vector(int, N);
	ret->actr.word_processing_block_times = vector(double, N);
	ret->actr.retrieval_share = vector(double, N);
	ret->actr.R_count = 0;
	ret->actr.atrial = actr_new_trial(1);
	ret->actr.current_retrieval_id = 0;
	ret->actr.current_retrieval_trigger = 0;
	ret->actr.last_retrieval_at = 0;
	ret->s = s;

	actr_add_memory_for(ret->actr.atrial, sentence_prop(ret->corpus, s, actr_template.memory_item_count), sentence_prop(ret->corpus, s, actr_template.memory_template), 0, 0, 0);

	return ret;
}

swift_run * clone_swift_trial(swift_run * src, unsigned int seed) {
	swift_run * ret = malloc(sizeof(swift_run));
	*ret = *src;
	ret->w_hist_by_stage.glob.hist = duplicate_matrix(double, src->w_hist_by_stage.glob.hist, 2, src->w_hist_by_stage.glob.ar_size);
	ret->w_hist_by_stage.nlab.hist = duplicate_matrix(double, src->w_hist_by_stage.nlab.hist, 2, src->w_hist_by_stage.nlab.ar_size);
	ret->w_hist_by_stage.lab.hist = duplicate_matrix(double, src->w_hist_by_stage.lab.hist, 2, src->w_hist_by_stage.lab.ar_size);
	ret->w_hist_by_stage.sacc.hist = duplicate_matrix(double, src->w_hist_by_stage.sacc.hist, 2, src->w_hist_by_stage.sacc.ar_size);
	ret->n_count = duplicate_vector(int, src->n_count, src->N + 4);
	ret->prev_n_count = duplicate_vector(int, src->prev_n_count, src->N + 4);
	ret->N_count = duplicate_vector(int, src->N_count, src->N + 4);
	ret->states = duplicate_vector(int, src->states, src->N + 4);
	ret->W = duplicate_vector(double, src->W, src->N + 4);
	ret->ptar = duplicate_vector(double, src->ptar, src->N);
	ret->aa = duplicate_vector(double, src->aa, src->N);
	ret->view = duplicate_vector(double, src->view, src->N);
	ret->border = duplicate_vector(double, src->border, src->N);
	ret->len = duplicate_vector(int, src->len, src->N);
	ret->procrate = duplicate_vector(double, src->procrate, src->N);
	ret->actr.word_waits_for_retrieval = duplicate_vector(int, src->actr.word_waits_for_retrieval, src->N);
	ret->actr.word_processing_block_times = duplicate_vector(double, src->actr.word_processing_block_times, src->N);
	ret->actr.retrieval_share = duplicate_vector(double, src->actr.retrieval_share, src->N);
	ret->actr.atrial = actr_duplicate_trial(src->actr.atrial);
	if(seed) {
		initSeed(seed, &ret->seed);
	}
	return ret;
}

void free_swift_trial(swift_run * trial) {
	clear_w_hist(trial->w_hist_by_stage.glob);
	clear_w_hist(trial->w_hist_by_stage.nlab);
	clear_w_hist(trial->w_hist_by_stage.lab);
	clear_w_hist(trial->w_hist_by_stage.sacc);
	free_vector(int, trial->prev_n_count);
	free_vector(int, trial->n_count);
	free_vector(int, trial->N_count);
	free_vector(int, trial->states);
	free_vector(double, trial->aa);
	free_vector(double, trial->W);
	free_vector(double, trial->procrate);
	free_vector(double, trial->view);
	free_vector(double, trial->border);
	free_vector(int, trial->len);
	free_vector(double, trial->ptar);
	free_vector(int, trial->actr.word_waits_for_retrieval);
	free_vector(double, trial->actr.word_processing_block_times);
	free_vector(double, trial->actr.retrieval_share);
	actr_free_trial(trial->actr.atrial);
	free(trial);
}

#define sq(A) ((A)*(A))

void processing_rate(swift_run * trial, double * procrate)

{
    double     dynfac, deltar, deltal, prate;
    double     ecc, leftact, bl, br, c0;
    int        j, l;

    
	/* dynamical processing span */ 
    dynfac = 1.0 - trial->aa[trial->gaze_word];
	deltar = trial->params->delta0 + trial->params->delta1*dynfac;
	deltal = trial->params->delta0*trial->params->asym;
	
    /* parameters of the inverted quadratic function */
    c0 = 3.0/2.0/(deltal+deltar);
    bl = c0/sq(deltal);
    br = c0/sq(deltar);

    /* calculating lexical processing rates */
    for ( j=1; j<=trial->N; j++ ) {
        procrate[j] = 0.0;
        for ( l=1; l<=trial->len[j]; l++ )  {
            ecc = trial->view[j] - trial->len[j]/2.0 + 1.0*l - 0.5 - trial->gaze_letter;  
            prate = 0.0;
            if ( -deltal<ecc && ecc<0 )  prate = c0 - bl*sq(ecc);
            if ( 0<=ecc && ecc<deltar )  prate = c0 - br*sq(ecc);
            procrate[j] += prate;
        }
        procrate[j] *= exp(-trial->params->eta*log(trial->len[j]));
        
        /* computing activation to the left */
        for ( l=1, leftact=1.0; l<j; l++ )  leftact *= exp(-trial->params->ppf*trial->aa[l]);
        procrate[j] *= leftact;
        
        /* predictability effect */
        //procrate[j] *= exp(trial->params->theta*sq(word_prop(trial->corpus, trial->s, j, pred)));


    }
    
    /* adding decay if as[.]==2 */
    for ( j=1; j<=trial->N; j++ ) {
        if ( trial->states[4+j] == STATE_POSTLEXICAL || trial->states[4+j] == STATE_RETRIEVAL || trial->states[4+j] == STATE_POSTRETRIEVAL )  {
            procrate[j] *= trial->params->proc;
            if ( procrate[j]<trial->params->decay )  procrate[j] = trial->params->decay;
        }
        if ( trial->states[4] && trial->states[4+j]==STATE_LEXICAL )  procrate[j] = 0.0;
    }
}

#undef sq

double compute_inhibition_rate(swift_run * trial, int word) {
	double ifovea, iparafovea, inhib;
	int i;

	/* compute inhibition */
    ifovea = trial->aa[word];
    for(i=word+1, iparafovea = 0.0; i <= trial->N; i++)  {
        iparafovea += trial->aa[i];
    }
    inhib = trial->params->h*ifovea + trial->params->h1*iparafovea;
    if ( inhib<0.0 )  inhib = 0.0;

	return 1.0/(1.0+inhib);
}


void counter_base_rates(swift_run * trial, double * W) {
	double inhibrate;
	double kapparate;

	/* compute inhibition */
	inhibrate = compute_inhibition_rate(trial, trial->gaze_word);

	kapparate = 1.0/(1.0 + trial->params->kappa0 * exp(-trial->params->kappa1*trial->dist*trial->dist));

	double labrate = 1.0;
	if(trial->is_mislocated) {
		labrate = trial->params->misfac;
	} else if(trial->is_refixation) {
		labrate = trial->params->refix;
	}

	W[1] = trial->N_count[1] / (trial->params->msac*100.0) * inhibrate;                          /* rate of random walk for timer */
	W[2] = trial->N_count[2] / (trial->params->tau_l*100.0) * labrate;                  /* rate of random walk for labile sacprog stage */
	W[3] = trial->N_count[3] / (trial->params->tau_n*100.0) * kapparate;        /* ... nonlabile stage */
	W[4] = trial->N_count[4] / (trial->params->tau_ex*100.0);                 /* ... saccade execution */

}

void transition_rates(swift_run * trial) {
	int i;

	counter_base_rates(trial, trial->W);

	trial->W[2] *= trial->states[2]; 
	trial->W[3] *= trial->states[3]; 
	trial->W[4] *= trial->states[4]; 

	// Update ACT-R retrieval latencies
    if ( trial->actr.atrial->n_items > trial->actr.last_retrieval_at && trial->actr.current_retrieval_id ) {
        actr_retrieval_result * result = vector(actr_retrieval_result, trial->actr.atrial->runs);
        actr_retrieval_result ** all_results = matrix(actr_retrieval_result, trial->actr.atrial->runs, trial->actr.atrial->n_items);
        actr_retrieve(get_actr_params(trial->params), (int) trial->t, trial->actr.current_retrieval, trial->actr.atrial->items, trial->actr.atrial->n_items, trial->actr.atrial->moments, trial->actr.atrial->n_moments, trial->actr.atrial->history, trial->actr.atrial->n_history, trial->actr.atrial->runs, result, all_results, &trial->seed, INFINITY);
        // only accept retrieval result for items we haven't looked at so far in the current retrieval
        for( i = trial->actr.last_retrieval_at + 1; i <= trial->actr.atrial->n_items ; i++ ) {
            if( trial->actr.atrial->items[i].trigger == result[1].memory_trigger ) {
                trial->actr.retrieval_result = result[1];
                //current_retrieval_ends = t + actr_mu1 + retrieval_result.latency;
                //word_processing_block_times[current_retrieval_trigger] = current_retrieval_ends;
                break;
            }
        }
        for( i = 1; i <= trial->actr.atrial->n_items; i++ ) {
            int memory_trigger = all_results[1][i].memory_trigger;
            if(memory_trigger && (trial->states[4+memory_trigger] != STATE_TRIGGERRETRIEVAL) && memory_trigger != trial->actr.current_retrieval_trigger) {
                trial->actr.retrieval_share[memory_trigger] = all_results[1][i].noisy_activation;
                log_event(trial, "ract %d BA=%.2lf MS=%.2lf NA=%.2lf RL=%.2lf", memory_trigger, all_results[1][i].base_activation, all_results[1][i].memory_strength, all_results[1][i].noisy_activation, all_results[1][i].latency);
                trial->states[4+memory_trigger] = STATE_TRIGGERRETRIEVAL;
            }
        }
        //printf_vector(retrieval_share, NW, "%.2lf ");
        free_vector(actr_retrieval_result, result);
        free_matrix(actr_retrieval_result, all_results);
        trial->actr.last_retrieval_at = trial->actr.atrial->n_items;
    }


	processing_rate(trial, trial->procrate);

	for ( i=1; i<= trial->N; i++ )  {
		if(trial->actr.word_processing_block_times[i] > trial->t) trial->W[4+i] = 0.0;
		else if(trial->states[4+i] == STATE_COMPLETE) trial->W[4+i] = 0.0;
		else if(trial->states[4+i] == STATE_RETRIEVAL && trial->n_count[4+i] <= trial->N_count[4+i] * trial->params->mu2) trial->W[4+i] = 0.0;
		else if(trial->states[4+i] == STATE_WAITFORRETRIEVAL) trial->W[4+i] = 0.0;
        else if(trial->states[4+i] == STATE_TRIGGERRETRIEVAL) {
			if(trial->actr.current_retrieval_started + trial->params->mu1 < trial->t) {
				trial->W[4+i] = trial->actr.R_count/(1000.0*trial->params->F*exp(-trial->actr.retrieval_share[i]));
			} else {
				trial->W[4+i] = 0.0;
			}
        }
        else if(trial->states[4+i] == STATE_POSTRETRIEVAL) {
            trial->W[4+i] = trial->procrate[i] * trial->params->alpha;
        }
	    else trial->W[4+i] = trial->procrate[i] * trial->params->alpha;
	}

}

void select_state(swift_run * trial, double * dt, int * state) {
	double Wsum = 0.0;
	int i;
	for(i = 1; i <= 4+trial->N; i++) {
		Wsum += trial->W[i];
	}

	*state = linear_sample_normalized(trial->W, Wsum, trial->N+4, &trial->seed);
	*dt = rexp(Wsum, &trial->seed);
}

void select_target(swift_run * trial) {
	selectar(trial->aa,&trial->states[4],trial->gaze_word,trial->N,trial->params->gamma,&trial->seed,1,trial->params->minact,trial->ptar);

	trial->saccade_target = linear_sample_normalized(trial->ptar, 1.0, trial->N, &trial->seed);
	trial->dist = fabs( trial->view[trial->saccade_target] - trial->gaze_letter );

	log_handler(trial, target_selected, trial->saccade_target);
}

void execute_saccade(swift_run * trial) {

	execsacc(trial->params, &trial->gaze_letter,&trial->gaze_word,&trial->saccade_target,trial->view,trial->border,trial->len,trial->N,&trial->seed,0,0.0,0);

	trial->is_mislocated = trial->saccade_target != trial->gaze_word;
	trial->is_refixation = trial->prev_gaze_word == trial->gaze_word;

	log_handler(trial, saccade_executed);
}

double loglik_spat_cond(swift_run * trial, int word, double x) {
	return log(trial->ptar[word]) + execsacc(trial->params, &trial->gaze_letter, &trial->gaze_word, &word, trial->view, trial->border, trial->len, trial->N, &trial->seed, 1, x, 0);
}

double loglik_spat(swift_run * trial, double x) {

	int i;

	double * p = vector(double, trial->N);

	for(i = 1; i <= trial->N; i++) {
		p[i] = loglik_spat_cond(trial, i, x);
	}

	//printf_vector(p, trial->N, " %.1lf");

	double loglik = logsumexp(p, trial->N);

	if(isnan(loglik)) {
		printf("%d %.1lf %.1lf\n", trial->gaze_word, trial->gaze_letter, x);
		printf_vector(trial->ptar, trial->N, " %.1lf");
		printf_vector(p, trial->N, " %.1lf");
	}

	free_vector(double, p);

	return loglik;
}

double loglik_temp(W_hist * w_hist, double t) {
	if(t <= 0.0) {
		return -INFINITY;
	} else {
		return cgammaengine2(t, w_hist->hist[1], w_hist->hist[2], w_hist->length, 1);
	}
}

int counter_direction(swift_run * trial, int state) {
	if(state > 4 && state <= 4 + trial->N) {
		if(trial->states[state] == STATE_POSTLEXICAL) return -1;
		if(trial->states[state] == STATE_POSTRETRIEVAL) return -1;
		if(trial->states[state] == STATE_RETRIEVAL) return -1;
	}
	return 1;
}

void propagate_counters(swift_run * trial, double dt, int state) {

	copy_vector(int, trial->n_count, trial->prev_n_count, trial->N + 4);

	trial->n_count[state] += counter_direction(trial, state);
	trial->t += dt;

	if(trial->n_count[1] >= trial->N_count[1]) {
		/* global timer ended -> (re-)start labile stage and reset global timer */
		trial->n_count[1] = 0;
		trial->n_count[2] = 0;
		if(trial->states[2]) {
			trial->canc++;
			log_event(trial, "sc");
		} else {
			trial->canc = 0;
			trial->states[2] = 1;
			log_event(trial, "lab");
		}
	}

	if(trial->n_count[2] >= trial->N_count[2] && !trial->states[3] && !trial->states[4]) {
		/* labile stage ended and non-labile stage can start -> select target and start non-labile stage */
		trial->n_count[2] = 0;
		trial->n_count[3] = 0;
		trial->states[3] = 1;
		trial->states[2] = 0;
		select_target(trial);
		log_event(trial, "st %d", trial->saccade_target);
		log_event(trial, "nlab");
	}

	if(trial->n_count[3] >= trial->N_count[3]) {
		/* non-labile stage ended -> execute saccade */
		trial->n_count[3] = 0;
		trial->n_count[4] = 0;
		trial->states[4] = 1;
		trial->states[3] = 0;
		log_event(trial, "sb");
	}

	if(trial->n_count[4] >= trial->N_count[4]) {
		trial->prev_gaze_word = trial->gaze_word;
		execute_saccade(trial);
		trial->n_count[4] = 0;
		trial->states[4] = 0;
		if(!trial->is_fitting && (trial->is_mislocated || trial->is_refixation)) {
			trial->n_count[1] = 0;
			trial->n_count[2] = 0;
			trial->states[2] = 1;
		}
		trial->fix_saccade_target = trial->saccade_target;
		log_event(trial, "se");
	}

	if(state > 4 && state <= 4 + trial->N) {

		int i, j = state - 4;

		if(trial->n_count[state] >= trial->N_count[state] && trial->states[state] == STATE_LEXICAL) {
			trial->n_count[state] = trial->N_count[state];
			log_event(trial, "s %d %d %d", j, trial->states[state], STATE_POSTLEXICAL);
			trial->states[state] = STATE_POSTLEXICAL;

            int already_encoded_in_memory = 0, i;
            for(i=1; i<=trial->actr.atrial->n_items; i++) {
                if(trial->actr.atrial->items[i].trigger == j) {
                    already_encoded_in_memory = 1;
                    break;
                }
            }
            if(!already_encoded_in_memory) {
                actr_add_memory_for(trial->actr.atrial, sentence_prop(trial->corpus, trial->s, actr_template.memory_item_count), sentence_prop(trial->corpus, trial->s, actr_template.memory_template), j, trial->t - trial->dt, 1);
				log_event(trial, "encode %d", j);
            }
            actr_prepare_retrievals_for(trial->actr.atrial, sentence_prop(trial->corpus, trial->s, actr_template.retrieval_item_count), sentence_prop(trial->corpus, trial->s, actr_template.retrieval_template), j, trial->t - trial->dt);
            trial->actr.word_waits_for_retrieval[j] = 0;
            for(i=trial->actr.atrial->n_retrievals_complete+1;i<=trial->actr.atrial->n_retrievals;i++) {
                if(j == trial->actr.atrial->retrievals[i].trigger) trial->actr.word_waits_for_retrieval[j]++;
            }
            if(trial->actr.word_waits_for_retrieval[j]) {
				log_event(trial, "s %d %d %d", j, trial->states[state], STATE_WAITFORRETRIEVAL);
                trial->states[state] = STATE_WAITFORRETRIEVAL;
            }
		}
		
		if(trial->n_count[state] <= 0 && trial->states[state] == STATE_POSTLEXICAL) {
			trial->n_count[state] = 0;
			log_event(trial, "s %d %d %d", j, trial->states[state], STATE_COMPLETE);
			trial->states[state] = STATE_COMPLETE;
		}


        if( trial->states[state] == STATE_POSTRETRIEVAL && trial->n_count[state] <= 0) {
			log_event(trial, "s %d %d %d", j, trial->states[state], STATE_COMPLETE);
            trial->states[state] = STATE_COMPLETE;
            trial->n_count[state] = 0;
        }

        if( trial->states[state] == STATE_TRIGGERRETRIEVAL && trial->n_count[state] >= trial->actr.R_count) {
            for(i=1; i<=trial->N; i++) {
                if(trial->states[i+4] == STATE_RETRIEVAL) {
					log_event(trial, "s %d %d %d", i, trial->states[i+4], STATE_POSTLEXICAL);
                    trial->states[i+4] = STATE_POSTLEXICAL;
                    trial->actr.word_processing_block_times[i] = trial->t + trial->params->mu1;
                } else if(trial->states[i+4] == STATE_TRIGGERRETRIEVAL) {
					log_event(trial, "s %d %d %d", i, trial->states[i+4], STATE_POSTRETRIEVAL);
                    trial->states[i+4] = STATE_POSTRETRIEVAL;
                }
            }
            log_event(trial, "match %d %d", trial->actr.current_retrieval_id, j);
            trial->actr.retrieval_result.memory_trigger = j;
            trial->actr.retrieval_result.latency = trial->t - trial->actr.current_retrieval_started;
            actr_end_retrieval(trial->actr.atrial, trial->actr.retrieval_result);
            trial->actr.atrial->n_retrievals_complete++;
            trial->actr.current_retrieval_id = 0;
            trial->actr.current_retrieval_trigger = 0;
            trial->actr.current_retrieval_ends = trial->t;
            trial->actr.last_retrieval_at = 0;
        }

	}

	int i;
	for(i = 1; i <= trial->N; i++) trial->aa[i] = (double) trial->n_count[4+i] / trial->params->aord;


    if(trial->actr.atrial->n_retrievals_complete < trial->actr.atrial->n_retrievals && !trial->actr.current_retrieval_id) {
        trial->actr.current_retrieval_id = trial->actr.atrial->n_retrievals_complete + 1;
        trial->actr.current_retrieval = trial->actr.atrial->retrievals[trial->actr.current_retrieval_id];
        trial->actr.current_retrieval_trigger = trial->actr.current_retrieval.trigger;
        log_event(trial, "retrieve %d %d", trial->actr.current_retrieval_id, trial->actr.current_retrieval_trigger);
        if(trial->states[4+trial->actr.current_retrieval_trigger] == STATE_WAITFORRETRIEVAL) {
            //trial->actr.word_processing_block_times[trial->actr.current_retrieval_trigger] = INFINITY;
            trial->actr.R_count = fmax(1.0, trial->params->aord * trial->params->rfrac);
			log_event(trial, "s %d %d %d", trial->actr.current_retrieval_trigger, trial->states[trial->actr.current_retrieval_trigger+4], STATE_RETRIEVAL);
            trial->states[4+trial->actr.current_retrieval_trigger] = STATE_RETRIEVAL;
            for(i=1;i<=trial->N;i++) {
                trial->actr.retrieval_share[i] = -INFINITY;
            }
        }
        trial->actr.last_retrieval_at = 0;
        trial->actr.current_retrieval_started = trial->t;
        trial->actr.current_retrieval_ends = INFINITY;
    }
}

int check_all_words_processed(swift_run * trial) {
	int i;
	for(i = 1; i <= trial->N; i++) if(trial->states[4+i] != STATE_COMPLETE && trial->states[4+i] != STATE_POSTRETRIEVAL) return 0;
	return 1;
}

int event_all_words_processed(swift_run * trial, va_list args) {
	return check_all_words_processed(trial);
}

int check_state_changed(swift_run * trial, int state, int stage) {
	return trial->states[state] == stage;
}

int event_state_changed(swift_run * trial, va_list args) {
	int state = va_arg(args, int);
	int stage = va_arg(args, int);
	return check_state_changed(trial, state, stage);
}

int check_time_passed(swift_run * trial, double t) {
	return trial->t >= t;
}

int event_time_passed(swift_run * trial, va_list args) {
	double t = va_arg(args, double);
	return check_time_passed(trial, t);
}

int check_saccade_cancelled(swift_run * trial, int ncanc) {
	return trial->canc >= ncanc;
}

int event_saccade_cancelled(swift_run * trial, va_list args) {
	int ncanc = va_arg(args, int);
	return check_saccade_cancelled(trial, ncanc);
}

int event_saccade_cancelled_or_nonlabile_stage_started(swift_run * trial, va_list args) {
	int ncanc = va_arg(args, int);
	return check_saccade_cancelled(trial, ncanc) || check_state_changed(trial, 3, 1);
}

int event_time_passed_or_saccade_started(swift_run * trial, va_list args) {
	double t = va_arg(args, double);
	return check_time_passed(trial, t) || check_state_changed(trial, 4, 1);
}

void log_transition_to_history(struct swift_run * trial, int state) {
	if(state == 1) append_w_hist(&trial->w_hist_by_stage.glob, 1.0, trial->W[1]);
	else if(state == 2) append_w_hist(&trial->w_hist_by_stage.lab, 1.0, trial->W[2]);
	else if(state == 3) append_w_hist(&trial->w_hist_by_stage.nlab, 1.0, trial->W[3]);
	else if(state == 4) append_w_hist(&trial->w_hist_by_stage.sacc, 1.0, trial->W[4]);
}

void swift_run_until_event(swift_run * trial, int check(swift_run * trial, va_list params), ...) {
	int state;
	va_list args;
	/* update transition rates */
	while(1) {
		transition_rates(trial);
		log_handler(trial, transition_rates_calculated);
		/* pass variable arguments as va_list to conditional function */
		va_start(args, check);
		if(check(trial, args)) break;
		va_end(args);
		/* select next transition and timestep */
		select_state(trial, &trial->dt, &state);
		log_handler(trial, state_selected, state);
		/* execute transition */
		propagate_counters(trial, trial->dt, state);
		log_handler(trial, counters_propagated);
		/* log new state */
		/* update transition rates */
	}
}

void loglik_swift(swift_model * model, swift_dataset * dataset, int * trials, int n_trials, swift_likelihood * likelihood) {

	int N = model->params->nsims;
	int i, j, k, l, n, u, v;

	if(trials == NULL) n = dataset->n;
	else n = n_trials;

	double **** ll = vector(double***, n);
	unsigned int *** seeds = vector(unsigned int**, n);

	for(i = 1; i <= n; i++) {
		swift_trial * sequence = &dataset->trials[trials == NULL ? i : trials[i]];
		int nfix = sequence->nfix;
		int nw = nwords(model->corpus, sequence->sentence);
		ll[i] = (double***) array(double, 3, nfix, 3, N);
		seeds[i] = matrix(unsigned int, nfix, N);
		for(j = 1; j <= N; j++) {
			for(u = 1; u <= nfix; u++) {
				seeds[i][u][j] = (unsigned int) ranint(&model->seed);
			}
		}
	}

	#pragma omp parallel for collapse(1) private(i,j,k,l,u,v) schedule(dynamic)
	for(k = 1; k <= n; k++) {
		swift_run ** ensemble = vector(swift_run*, N);
		swift_trial * sequence = &dataset->trials[trials == NULL ? k : trials[k]];
		double t_fix_started, t_sac_started;
		double p_variation;

		swift_run * trial = new_swift_trial(model, sequence->sentence, seeds[k][1][1]);
		trial->is_fitting = 1;
		trial->ptar[sequence->fixations[1].fw] = 1.0;
		trial->saccade_target = 1;
		trial->fix_saccade_target = 1;
		// set an event handler for when a new state is selected -> save to transition log within trial (trial->w_hist_by_stage)
		trial->handlers.state_selected = log_transition_to_history;
		trial->gaze_letter = sequence->fixations[1].fl;
		if(sequence->fixations[1].fw > 1) {
			trial->gaze_letter += trial->border[sequence->fixations[1].fw-1];
		}
		trial->gaze_word = sequence->fixations[1].fw;


		W_hist * w_hist = vector(W_hist, N);

		int * saccade_targets = vector(int, N);
		double * ensemble_logliks = vector(double, N);
		double * p_spat = vector(double, trial->N);

		for(i = 1; i <= sequence->nfix; i++) {


			t_fix_started = trial->t;

			trial->gaze_letter = sequence->fixations[i].fl;
			if(sequence->fixations[i].fw > 1) {
				trial->gaze_letter += trial->border[sequence->fixations[i].fw-1];
			}
			trial->gaze_word = sequence->fixations[i].fw;

			init_w_hist(&trial->w_hist_by_stage.glob);
			init_w_hist(&trial->w_hist_by_stage.lab);
			init_w_hist(&trial->w_hist_by_stage.nlab);
			init_w_hist(&trial->w_hist_by_stage.sacc);

			double upcoming_location, current_location;

			current_location = sequence->fixations[i].fl;
			if(sequence->fixations[i].fw > 1) current_location += trial->border[sequence->fixations[i].fw-1];
			if(i < sequence->nfix) {
				upcoming_location = sequence->fixations[i+1].fl;
				if(sequence->fixations[i+1].fw > 1) upcoming_location += trial->border[sequence->fixations[i+1].fw-1];
			} else {
				upcoming_location = NAN;
			}




			for(j = 1; j <= N; j++) {

				ensemble[j] = j == N ? trial : clone_swift_trial(trial, seeds[k][i][j]);

				if(i > 1) {
					ensemble[j]->saccade_target = saccade_targets[j];


					double mlp = 1.0 - exp(p_spat[trial->gaze_word] - logsumexp(p_spat, ensemble[j]->N));


					ensemble[j]->is_mislocated = sample_boolean(mlp, &ensemble[j]->seed);

					//printf_vector(p_spat, trial->N, "%.2lf ");

					//printf("%d %.2lf\n", ensemble[j]->is_mislocated, mlp);

					ensemble[j]->is_refixation = sequence->fixations[i-1].fw == ensemble[j]->gaze_word;

					if(ensemble[j]->is_mislocated || ensemble[j]->is_refixation) {
						ensemble[j]->n_count[1] = 0;
						ensemble[j]->n_count[2] = 0;
						ensemble[j]->states[2] = 1;
					}
				}


				w_hist[j].hist = NULL;

				init_w_hist(&w_hist[j]);

				swift_run_until_event(ensemble[j], event_state_changed, 2, 1); // first labile stage started

				do {
					transfer_w_hist(&ensemble[j]->w_hist_by_stage.glob, &w_hist[j]);
					init_w_hist(&ensemble[j]->w_hist_by_stage.glob);
					init_w_hist(&ensemble[j]->w_hist_by_stage.lab);
					swift_run_until_event(ensemble[j], event_saccade_cancelled_or_nonlabile_stage_started, ensemble[j]->canc+1);
				} while(!check_state_changed(ensemble[j], 3, 1));
				transfer_w_hist(&ensemble[j]->w_hist_by_stage.lab, &w_hist[j]);

				if(i < sequence->nfix) {
					ll[k][i][1][j] = loglik_spat(ensemble[j], upcoming_location);
				}

				swift_run_until_event(ensemble[j], event_state_changed, 4, 1); // saccade execution (4) started (1)
				transfer_w_hist(&ensemble[j]->w_hist_by_stage.nlab, &w_hist[j]);

				ll[k][i][2][j] = loglik_temp(&w_hist[j], sequence->fixations[i].tfix);
				//trial->t = t_fix_started + sequence->fixations[i].tfix;
				swift_run_until_event(ensemble[j], event_state_changed, 4, 0); // saccade execution (4) finished (0)
				
				if(i < sequence->nfix) {
					ll[k][i][3][j] = loglik_temp(&ensemble[j]->w_hist_by_stage.sacc, sequence->fixations[i].tsac);
				}

				ensemble_logliks[j] = ll[k][i][1][j] + ll[k][i][2][j] + ll[k][i][3][j];

				clear_w_hist(w_hist[j]);
			}
			

			int best_trial = 1;

			for(j = 2; j <= N; j++) {
				if(ensemble_logliks[j] > ensemble_logliks[best_trial]) {
					best_trial = j;
				}
			}

			for(j = 1; j <= N; j++) {
				saccade_targets[j] = ensemble[j]->saccade_target;
				if(j == best_trial) {
					trial = ensemble[j];
					trial->gaze_letter = current_location;
					trial->gaze_word = trial->prev_gaze_word;
					for(u = 1; u <= trial->N; u++) {
						p_spat[u] = loglik_spat_cond(trial, u, upcoming_location);
					}
				} else {
					free_swift_trial(ensemble[j]);
				}
			}


			trial->t = t_fix_started + sequence->fixations[i].tfix + sequence->fixations[i].tsac;
		}
		free_vector(double, p_spat);
		free_vector(double, ensemble_logliks);
		free_vector(int, saccade_targets);
		free_vector(W_hist, w_hist);
		free_swift_trial(trial);
		free_vector(swift_run*, ensemble);
	}

	likelihood->fixation_location = 0.0;
	likelihood->fixation_duration = 0.0;
	likelihood->saccade_duration = 0.0;

	#pragma omp parallel for private(k,i) schedule(guided)
	for(k = 1; k <= n; k++) {
		swift_trial * sequence = &dataset->trials[trials == NULL ? k : trials[k]];
		for(i = 1; i <= sequence->nfix; i++) {
			#pragma omp atomic update
			likelihood->fixation_location += logsumexp(ll[k][i][1], N) - log(N);
			#pragma omp atomic update
			likelihood->fixation_duration += logsumexp(ll[k][i][2], N) - log(N);
			#pragma omp atomic update
			likelihood->saccade_duration += logsumexp(ll[k][i][3], N) - log(N);
		}
		free_array(double, ll[k], 3);
		free_matrix(unsigned int, seeds[k]);
	}

	free_vector(double***, ll);
	free_vector(unsigned int**, seeds);
	
}


void generate_swift_single(swift_model * model, int s, swift_trial * sequence, unsigned int seed, swift_output_files * files) {
	swift_run * trial = new_swift_trial(model, s, seed);

	trial->handlers.counters_propagated = log_activation_to_file;

	if(files != NULL) trial->files = *files;

	int Nfix = 0;
	int ar_size = 20;
	double t_fix_started;

	sequence->fixations = vector(swift_fixation, ar_size);
	sequence->sentence = s;

	sequence->fixations[1].fw = 1;
	sequence->fixations[1].fl = trial->view[1];
	trial->saccade_target = 1;
	trial->fix_saccade_target = 1;

	trial->ptar[trial->saccade_target] = 1.0;


	while(1) {
		t_fix_started = trial->t;
		Nfix++;
		if(Nfix >= ar_size) {
			sequence->fixations = resize_vector(swift_fixation, sequence->fixations, ar_size, ar_size*2);
			ar_size*=2;
		}
		swift_run_until_event(trial, event_state_changed, 4, 1); /* saccade started */
		sequence->fixations[Nfix].tfix = (int) (trial->t - t_fix_started);
		t_fix_started = trial->t;
		if(check_all_words_processed(trial) && Nfix >= 2) break;
		swift_run_until_event(trial, event_state_changed, 4, 0); /* saccade ended */
		sequence->fixations[Nfix+1].fw = trial->gaze_word;
		if(trial->gaze_word > 1) sequence->fixations[Nfix+1].fl = trial->gaze_letter - trial->border[trial->gaze_word-1];
		else sequence->fixations[Nfix+1].fl = trial->gaze_letter;
		sequence->fixations[Nfix].tsac = (int) (trial->t - t_fix_started);
	}


	sequence->nfix = Nfix;
	if(Nfix != ar_size) sequence->fixations = resize_vector(swift_fixation, sequence->fixations, ar_size, Nfix);

	free_swift_trial(trial);

}

void generate_swift_all(swift_model * model, swift_dataset * data, swift_output_files * files) {
	int i, j;
	data->n = nsentences(model->corpus) * model->params->runs;
	unsigned int * seeds = vector(unsigned int, data->n);
	for(i = 1; i <= data->n; i++) {
		seeds[i] = (unsigned int) ranint(&model->seed);
	}
	data->trials = vector(swift_trial, data->n);
	#pragma omp parallel for collapse(2) private(i,j)
	for(i = 1; i <= nsentences(model->corpus); i++) {
		for(j = 1; j <= model->params->runs; j++) {
			generate_swift_single(model, i, &data->trials[(i-1)*model->params->runs+j], seeds[(i-1)*model->params->runs+j], files);
		}
	}
	free_vector(unsigned int, seeds);
}

void generate_swift_some(swift_model * model, int * items, int n_items, swift_dataset * data, swift_output_files * files) {
	if(items == NULL) {
		generate_swift_all(model, data, files);
		return;
	}

	int i, j;
	data->n = n_items * model->params->runs;
	data->trials = vector(swift_trial, data->n);
	unsigned int * seeds = vector(unsigned int, data->n);
	for(i = 1; i <= data->n; i++) {
		seeds[i] = (unsigned int) ranint(&model->seed);
	}
	#pragma omp parallel for collapse(2) private(i,j)
	for(i = 1; i <= n_items; i++) {
		for(j = 1; j <= model->params->runs; j++) {
			generate_swift_single(model, items[i], &data->trials[(i-1)*model->params->runs+j], seeds[(i-1)*model->params->runs+j], files);
		}
	}
	free_vector(unsigned int, seeds);
}
