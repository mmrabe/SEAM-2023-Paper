
#include<stdio.h>
#include<math.h>

typedef struct {
	int * n_count;
	int * N_count;
	int * states;
	int N;
	double t;
	int gaze_word;
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
	int * len;
	int is_mislocated, is_refixation;
	double loglik_temp, loglik_spat;
} swift_run;

typedef struct {
	double spatial, temporal;
} swift_likelihood;

swift_run * new_swift_trial(swift_model * model, int s) {
	int i, j;
	int N = nwords(model->corpus, s);
	swift_run * ret = malloc(sizeof(swift_run));
	ret->n_count = vector(int, N+4);
	ret->aa = vector(double, N);
	ret->procrate = vector(double, N);
	ret->N_count = vector(int, N+4);
	ret->ptar = vector(double, N);
	ret->W = vector(double, N+4);
	ret->is_mislocated = 0;
	ret->is_refixation = 0;
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
	ret->states = vector(int, N+4);
	ret->states[1] = 1;
	ret->n_count[1] = (int) (model->params->msac0 * ret->N_count[1]);
	ret->t = 0.0;
	ret->dist = 0.0;
	ret->gaze_word = 1;
	ret->gaze_letter = 1.0 + 0.5 * word_prop(model->corpus, s, 1, nl);
	ret->N = N;
	initSeed(ranlong(&model->seed), &ret->seed);
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
    	ret->view[i] = b + 1.0 + 0.5 * ret->len[i];
    	ret->border[i] = b + 1.0 + ret->len[i];
    }
	ret->corpus = model->corpus;
	ret->params = model->params;
	return ret;
}

swift_run * clone_swift_trial(swift_run * src) {
	swift_run * ret = malloc(sizeof(swift_run));
	ret->n_count = duplicate_vector(int, src->n_count, src->N + 4);
	ret->N_count = duplicate_vector(int, src->N_count, src->N + 4);
	ret->states = duplicate_vector(int, src->states, src->N + 4);
	ret->W = duplicate_vector(double, src->W, src->N + 4);
	ret->N = src->N;
	ret->t = src->t;
	ret->is_mislocated = src->is_mislocated;
	ret->is_refixation = src->is_refixation;
	ret->loglik_temp = src->loglik_temp;
	ret->loglik_spat = src->loglik_spat;
	ret->ptar = duplicate_vector(double, src->ptar, src->N);
	ret->aa = duplicate_vector(double, src->aa, src->N);
	ret->view = duplicate_vector(double, src->view, src->N);
	ret->border = duplicate_vector(double, src->border, src->N);
	ret->len = duplicate_vector(int, src->len, src->N);
	ret->aa = duplicate_vector(double, src->aa, src->N);
	ret->procrate = duplicate_vector(double, src->procrate, src->N);
	ret->gaze_word = src->gaze_word;
	ret->gaze_letter = src->gaze_letter;
	ret->corpus = src->corpus;
	ret->params = src->params;
	ret->s = src->s;
	ret->saccade_target = src->saccade_target;
	initSeed(ranlong(&src->seed), &ret->seed);
	return ret;
}

void free_swift_trial(swift_run * trial) {
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
        if ( trial->states[4+j] == STATE_POSTLEXICAL )  {
            procrate[j] *= trial->params->proc;
            if ( procrate[j]<trial->params->decay )  procrate[j] = trial->params->decay;
        }
        if ( trial->states[4] && trial->states[4+j]==STATE_LEXICAL )  procrate[j] = 0.0;
    }
}

#undef sq

void transition_rates(swift_run * trial) {
	double ifovea, iparafovea, inhib, inhibrate, kapparate;
	int i;


	/* compute inhibition */
    ifovea = trial->aa[trial->gaze_word];
    for(i=trial->gaze_word+1, iparafovea = 0.0; i <= trial->N; i++)  {
        iparafovea += trial->aa[i];
    }
    inhib = trial->params->h*ifovea + trial->params->h1*iparafovea;
    if ( inhib<0.0 )  inhib = 0.0;
	inhibrate = 1.0/(1.0+inhib);

	kapparate = 1.0/(1.0 + trial->params->kappa0 * exp(-trial->params->kappa1*trial->dist*trial->dist));

	trial->W[1] = trial->N_count[1] / (trial->params->msac*100.0) * inhibrate;                          /* rate of random walk for timer */
	trial->W[2] = trial->N_count[2] / (trial->params->tau_l*100.0) * trial->states[2];                  /* rate of random walk for labile sacprog stage */
	if(trial->is_mislocated) trial->W[2] *= trial->params->misfac;
	else if(trial->is_refixation) trial->W[2] *= trial->params->refix;
	trial->W[3] = trial->N_count[3] / (trial->params->tau_n*100.0) * trial->states[3]*kapparate;        /* ... nonlabile stage */
	trial->W[4] = trial->N_count[4] / (trial->params->tau_ex*100.0) * trial->states[4];                 /* ... saccade execution */

	processing_rate(trial, trial->procrate);

	for ( i=1; i<= trial->N; i++ )  {
		if(trial->states[4+i] == STATE_COMPLETE) trial->W[4+i] = 0.0;
	    else trial->W[4+i] = trial->procrate[i] * trial->params->alpha;
	}

	//printf("W: ");printf_vector(trial->W, trial->N+4, " %.3lf");

}

void select_state(swift_run * trial, double * dt, int * state) {
	double Wsum = 0.0;
	int i;
	for(i = 1; i <= 4+trial->N; i++) {
		Wsum += trial->W[i];
	}
	double test = ran1(&trial->seed) * Wsum;
	for(i = 1; i <= 4+trial->N; i++) {
		if(test < trial->W[i]) break;
		else test -= trial->W[i];
	}
	*dt = rexp(Wsum, &trial->seed);
	*state = i;
}

void select_target(swift_run * trial) {
	selectar(trial->aa,&trial->states[4],trial->gaze_word,trial->N,trial->params->gamma,&trial->seed,1,trial->params->minact,trial->ptar);
	int i;
	double test = ran1(&trial->seed);
	for(i = 1; i <= trial->N; i++) {
		if(test < trial->ptar[i]) break;
		else test -= trial->ptar[i];
	}
	trial->saccade_target = i;
	trial->dist = fabs( trial->view[i] - trial->gaze_letter );
}

void execute_saccade(swift_run * trial) {

	int previous_word = trial->gaze_word;

	execsacc(trial->params, &trial->gaze_letter,&trial->gaze_word,&trial->saccade_target,trial->view,trial->border,trial->len,trial->N,&trial->seed,0,0.0,0);

	trial->is_mislocated = trial->saccade_target != trial->gaze_word;
	trial->is_refixation = previous_word == trial->gaze_word;

}

double loglik_spat(swift_run * trial, double x) {

	int i;

	double * p = vector(double, trial->N);

	for(i = 1; i <= trial->N; i++) {
		double ll = execsacc(trial->params, &trial->gaze_letter, &trial->gaze_word, &i, trial->view, trial->border, trial->len, trial->N, &trial->seed, 1, x, 0);
		p[i] = trial->ptar[i] + ll;
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

double loglik_temp(swift_run * trial, int state, double t) {
	double rv;
	if(t <= 0.0) {
		rv = log(DBL_EPSILON) + t;
	} else {
		rv = gammaloglike(t, trial->N_count[state] - trial->n_count[state], trial->W[state]) + log1p(-DBL_EPSILON);
	}
	return rv;
}

int counter_direction(swift_run * trial, int state) {
	if(state > 4 && state <= 4 + trial->N && trial->states[state] == STATE_POSTLEXICAL) {
		return -1;
	}
	return 1;
}

void propagate_counters(swift_run * trial, double dt, int state) {

	trial->n_count[state] += counter_direction(trial, state);
	trial->t += dt;

	if(trial->n_count[1] >= trial->N_count[1]) {
		/* global timer ended -> (re-)start labile stage and reset global timer */
		trial->n_count[1] = 0;
		trial->n_count[2] = 0;
		trial->states[2] = 1;
	}

	if(trial->n_count[2] >= trial->N_count[2] && !trial->states[3] && !trial->states[4]) {
		/* labile stage ended and non-labile stage can start -> select target and start non-labile stage */
		trial->n_count[2] = 0;
		trial->n_count[3] = 0;
		trial->states[3] = 1;
		trial->states[2] = 0;
		select_target(trial);
	}

	if(trial->n_count[3] >= trial->N_count[3]) {
		/* non-labile stage ended -> execute saccade */
		trial->n_count[3] = 0;
		trial->n_count[4] = 0;
		trial->states[4] = 1;
		trial->states[3] = 0;
	}

	if(trial->n_count[4] >= trial->N_count[4]) {
		execute_saccade(trial);
		trial->n_count[4] = 0;
		trial->states[4] = 0;
		if(trial->is_mislocated || trial->is_refixation) {
			trial->n_count[1] = 0;
			trial->n_count[2] = 0;
			trial->states[2] = 1;
		}
	}

	if(state > 4 && state <= 4 + trial->N) {

		if(trial->n_count[state] >= trial->N_count[state] && trial->states[state] == STATE_LEXICAL) {
			trial->n_count[state] = trial->N_count[state];
			trial->states[state] = STATE_POSTLEXICAL;
		}
		
		if(trial->n_count[state] <= 0 && trial->states[state] == STATE_POSTLEXICAL) {
			trial->n_count[state] = 0;
			trial->states[state] = STATE_COMPLETE;
		}

	}

	int i;
	for(i = 1; i <= trial->N; i++) trial->aa[i] = (double) trial->n_count[4+i] / trial->params->aord;
}

int check_all_words_processed(swift_run * trial) {
	int i;
	for(i = 1; i <= trial->N; i++) if(trial->states[4+i] != STATE_COMPLETE) return 0;
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

void swift_run_until_event(swift_run * trial, int check(swift_run * trial, va_list params), void log_trial(swift_run * trial, double dt, int state), ...) {
	double dt;
	int state;
	va_list args;
	/* update transition rates */
	while(1) {
		transition_rates(trial);
		/* pass variable arguments as va_list to conditional function */
		va_start(args, log_trial);
		if(check(trial, args)) break;
		va_end(args);
		/* select next transition and timestep */
		select_state(trial, &dt, &state);
		/* execute transition */
		propagate_counters(trial, dt, state);
		/* if logging handler was given, log current state */
		if(log_trial != NULL) {
			printf("LOG! %lx\n", log_trial);
			log_trial(trial, dt, state);
		}
		/* update transition rates */
	}
}

/*
void loglik_swift(swift_model * model, swift_dataset * dataset, int * trials, int n_trials, swift_likelihood * likelihood) {

	int N = model->params->nsims;
	int i, j, k, n;

	if(trials == NULL) n = dataset->n;
	else n = n_trials;

	double *** ll = (double***) array(double, 3, n, 2, N);

	likelihood->temporal = 0.0;
	likelihood->spatial = 0.0;

	#pragma omp parallel for collapse(2) private(i,j,k) schedule(dynamic)
	for(k = 1; k <= n; k++) {
		for(j = 1; j <= N; j++) {
			swift_trial * sequence;
			if(trials == NULL) sequence = &dataset->trials[k];
			else sequence = &dataset->trials[trials[k]];
			double t_fix_started;
			double t_nl_started;
			swift_run * trial = new_swift_trial(model, sequence->sentence);
			for(i = 1; i <= sequence->nfix; i++) {
				t_fix_started = trial->t;
				trial->gaze_letter = sequence->fixations[i].fl;
				trial->gaze_word = sequence->fixations[i].fw;
				swift_run_until_event(trial, event_state_changed, NULL, 3, 1); // non-labile saccade program (3) started (1)
				t_nl_started = trial->t;
				ll[k][1][j] += loglik_temp(trial, 3, sequence->fixations[i].tfix - (t_nl_started - t_fix_started));
				swift_run_until_event(trial, event_state_changed, NULL, 4, 1); // saccade execution (4) started (1)
				trial->t = t_fix_started + sequence->fixations[i].tfix;
				if(i < sequence->nfix) {
					double x = sequence->fixations[i+1].fl;
					if(sequence->fixations[i+1].fw > 1) x += trial->border[sequence->fixations[i+1].fw-1];
					ll[k][2][j] += loglik_spat(trial, x);
					#pragma omp critical(iop)
					if(isnan(ll[k][2][j]) || isinf(ll[k][2][j])) {
						write_trial(stdout, *sequence);
						printf("S%d R%d F%d/%d: %lf\n", k, j, i+1, sequence->nfix, loglik_spat(trial, x));
						printf_vector(trial->states, trial->N+4, "% 4d");
						printf_vector(trial->n_count, trial->N+4, "% 4d");
						printf_vector(trial->N_count, trial->N+4, "% 4d");
						printf_vector(trial->ptar, trial->N, " %.3lf");
						exit(1);
					}
				}
				swift_run_until_event(trial, event_state_changed, NULL, 4, 0); // saccade execution (4) finished (0)
				trial->t = t_fix_started + sequence->fixations[i].tfix + sequence->fixations[i].tsac;
			}
			free_swift_trial(trial);
		}
	}

	#pragma omp parallel for private(k)
	for(k = 1; k <= n; k++) {
		//printf("Loglik(%d/%d): %lf, %lf\n", i, Nfix, a, b);
		#pragma omp atomic update
		likelihood->temporal += logsumexp(ll[k][1], N) - log(N);
		#pragma omp atomic update
		likelihood->spatial += logsumexp(ll[k][2], N) - log(N);
	}

	free_array(double, ll, 3);
	
}
*/

void loglik_swift(swift_model * model, swift_dataset * dataset, int * trials, int n_trials, swift_likelihood * likelihood) {

	int N = model->params->nsims;
	int i, j, k, n;

	if(trials == NULL) n = dataset->n;
	else n = n_trials;

	double **** ll = vector(double***, n);

	for(k = 1; k <= n; k++) {
		ll[k] = (double***) array(double, 3, dataset->trials[k].nfix, 2, N);
	}

	likelihood->temporal = 0.0;
	likelihood->spatial = 0.0;

	#pragma omp parallel for collapse(2) private(i,j,k) schedule(dynamic)
	for(k = 1; k <= n; k++) {
		for(j = 1; j <= N; j++) {
			swift_trial * sequence;
			if(trials == NULL) sequence = &dataset->trials[k];
			else sequence = &dataset->trials[trials[k]];
			double t_fix_started;
			double t_nl_started;
			swift_run * trial = new_swift_trial(model, sequence->sentence);
			for(i = 1; i <= sequence->nfix; i++) {
				t_fix_started = trial->t;
				trial->gaze_letter = sequence->fixations[i].fl;
				trial->gaze_word = sequence->fixations[i].fw;
				swift_run_until_event(trial, event_state_changed, NULL, 3, 1); // non-labile saccade program (3) started (1)
				t_nl_started = trial->t;
				printf("Start of fix %.1lf, nl %.1lf, data %.1lf, loglik of %.1lf\n", t_fix_started, t_nl_started, sequence->fixations[i].tfix, sequence->fixations[i].tfix - (t_nl_started - t_fix_started));
				ll[k][i][1][j] = loglik_temp(trial, 3, sequence->fixations[i].tfix - (t_nl_started - t_fix_started));
				swift_run_until_event(trial, event_state_changed, NULL, 4, 1); // saccade execution (4) started (1)
				trial->t = t_fix_started + sequence->fixations[i].tfix;
				if(i < sequence->nfix) {
					double x = sequence->fixations[i+1].fl;
					if(sequence->fixations[i+1].fw > 1) x += trial->border[sequence->fixations[i+1].fw-1];
					ll[k][i+1][2][j] += loglik_spat(trial, x);
					#pragma omp critical(iop)
					if(isnan(ll[k][i+1][2][j]) || isinf(ll[k][i+1][2][j])) {
						write_trial(stdout, *sequence);
						printf("S%d R%d F%d/%d: %lf\n", k, j, i+1, sequence->nfix, loglik_spat(trial, x));
						printf_vector(trial->states, trial->N+4, "% 4d");
						printf_vector(trial->n_count, trial->N+4, "% 4d");
						printf_vector(trial->N_count, trial->N+4, "% 4d");
						printf_vector(trial->ptar, trial->N, " %.3lf");
						exit(1);
					}
				}
				swift_run_until_event(trial, event_state_changed, NULL, 4, 0); // saccade execution (4) finished (0)
				trial->t = t_fix_started + sequence->fixations[i].tfix + sequence->fixations[i].tsac;
			}
			free_swift_trial(trial);
		}
	}

	#pragma omp parallel for private(k,i) schedule(guided)
	for(k = 1; k <= n; k++) {
		swift_trial * sequence;
		if(trials == NULL) sequence = &dataset->trials[k];
		else sequence = &dataset->trials[trials[k]];
		for(i = 1; i <= sequence->nfix; i++) {
			//printf("Loglik(%d/%d): %lf, %lf\n", i, Nfix, a, b);
			#pragma omp atomic update
			likelihood->temporal += logsumexp(ll[k][i][1], N) - log(N);
			if(i > 1) {
				#pragma omp atomic update
				likelihood->spatial += logsumexp(ll[k][i][2], N) - log(N);
			}
		}
		free_array(double, ll[k], 3);
	}

	free_vector(double***, ll);
	
}


void generate_swift_single(swift_model * model, int s, swift_trial * sequence) {
	swift_run * trial = new_swift_trial(model, s);

	int Nfix = 0;
	int ar_size = 20;
	double t_fix_started;

	sequence->fixations = vector(swift_fixation, ar_size);
	sequence->sentence = s;

	sequence->fixations[1].fw = 1;
	sequence->fixations[1].fl = trial->view[1];


	while(1) {
		t_fix_started = trial->t;
		Nfix++;
		if(Nfix >= ar_size) {
			sequence->fixations = resize_vector(swift_fixation, sequence->fixations, ar_size, ar_size*2);
			ar_size*=2;
		}
		swift_run_until_event(trial, event_state_changed, NULL, 4, 1); /* saccade started */
		sequence->fixations[Nfix].tfix = (int) (trial->t - t_fix_started);
		t_fix_started = trial->t;
		if(check_all_words_processed(trial) && Nfix >= 2) break;
		swift_run_until_event(trial, event_state_changed, NULL, 4, 0); /* saccade ended */
		sequence->fixations[Nfix+1].fw = trial->gaze_word;
		if(trial->gaze_word > 1) sequence->fixations[Nfix+1].fl = trial->gaze_letter - trial->border[trial->gaze_word-1];
		else sequence->fixations[Nfix+1].fl = trial->gaze_letter;
		sequence->fixations[Nfix].tsac = (int) (trial->t - t_fix_started);
	}


	sequence->nfix = Nfix;
	if(Nfix != ar_size) sequence->fixations = resize_vector(swift_fixation, sequence->fixations, ar_size, Nfix);

	free_swift_trial(trial);

}

void generate_swift_all(swift_model * model, swift_dataset * data) {
	int i, j;
	data->n = nsentences(model->corpus) * model->params->runs;
	data->trials = vector(swift_trial, data->n);
	#pragma omp parallel for collapse(2) private(i,j)
	for(i = 1; i <= nsentences(model->corpus); i++) {
		for(j = 1; j <= model->params->runs; j++) {
			generate_swift_single(model, i, &data->trials[(i-1)*model->params->runs+j]);
		}
	}
}

void generate_swift_some(swift_model * model, int * items, int n_items, swift_dataset * data) {
	int i, j;
	data->n = n_items * model->params->runs;
	data->trials = vector(swift_trial, data->n);
	#pragma omp parallel for collapse(2) private(i,j)
	for(i = 1; i <= n_items; i++) {
		for(j = 1; j <= model->params->runs; j++) {
			generate_swift_single(model, items[i], &data->trials[(i-1)*model->params->runs+j]);
		}
	}
}

void test_swift(swift_model * model, int s) {

	swift_dataset data;

	data.name = NULL;
	generate_swift_all(model, &data);


	write_dataset(stdout, data);

/*

	double ll_temp, ll_spat;

	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 3.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 2.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	
	model->params->msac = 1.5;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 3.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 4.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 5.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 6.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);

	model->params->msac = 7.0;
	loglik_swift(model, &data, &ll_temp, &ll_spat);
	printf("Loglik (msac=%.1lf): %lf, %lf\n", model->params->msac, ll_temp, ll_spat);
	*/

	clear_dataset(data);



}

