
#ifndef __ACTR_LOGIC

#define __ACTR_LOGIC

#include "array.c"
#include "logsumexp.c"

double rlogistic(double mu, double s, RANSEED_TYPE * seed) {
	double p = ran1(seed);
	return mu + s * (log(p) - log1p(-p));
}


// ACT-R stuff starts here


#define actr_get_value(type, obj) *((type*) obj.value)
#define actr_get_label(obj) ((char*) ((size_t) &((actr_feature_type_enum_extra*) actr_feature_types[obj.type].extra)->values + (size_t) ((size_t*) &((actr_feature_type_enum_extra*) actr_feature_types[obj.type].extra)->values)[actr_get_value(int, obj)-1]))
#define actr_get_type_name(obj) (actr_feature_types[obj.type].name)
#define actr_get_type_class(obj) (actr_feature_types[obj.type].class)
#define actr_is_nil(obj) (actr_feature_types[obj.type].legal_nil && obj.is_nil)
#define actr_is_null(obj) (actr_feature_types[obj.type].legal_null && obj.is_null)

typedef enum {
	ACTR_FEATURE_INT, ACTR_FEATURE_ENUM, ACTR_FEATURE_STRING, ACTR_FEATURE_BOOL
} actr_feature_class;

typedef struct {
	actr_feature_class class;
	char name[50];
	int legal_null;
	int legal_nil;
	void * extra;
} actr_feature_type;

typedef struct {
	int nvalues;
	char values[];
} actr_feature_type_enum_extra;

actr_feature_type * actr_feature_types = NULL;
int actr_n_feature_types = 0;

typedef struct {
	int type;
	int is_nil;
	void * value;
} actr_feature;

typedef struct {
	int type;
	int is_null;
	void * value;
} actr_cue;

typedef struct {
	char name[128];
	int creation;
	int trigger;
	int nfeatures;
	actr_feature * features;
} actr_memory_item;

typedef struct {
	int moment;
	int trigger;
	int ncues;
	actr_cue * cues;
} actr_retrieval_item;

int actr_find_feature_type(char * name) {
	int i;
	for(i = 1; i <= actr_n_feature_types; i++) {
		if(!strcmp(actr_feature_types[i].name, name)) {
			return i;
		}
	}
	return -1;
}

void actr_free_feature_types() {
	int i, j;
	for(i = 1; i <= actr_n_feature_types; i++) {
		if(actr_feature_types[i].class == ACTR_FEATURE_ENUM) {
			actr_feature_type_enum_extra * extra = (actr_feature_type_enum_extra *) &actr_feature_types[i].extra;
			// for(j = 1; j <= extra->nvalues; j++) {
			// 	free(extra->values[j]);
			// }
			free_vector(char**, extra->values);
			free(extra);
		} else if(actr_feature_types[i].class == ACTR_FEATURE_STRING) {
			free(actr_feature_types[i].extra);
		}
	}
}

void actr_clear_feature(actr_feature x) {
	if(x.value != NULL) free(x.value);
}

void actr_clear_cue(actr_cue x) {
	if(x.value != NULL) free(x.value);
}

void actr_free_feature(actr_feature* x) {
	actr_clear_feature(*x);
	free(x);
}

void actr_free_cue(actr_cue* x) {
	actr_clear_cue(*x);
	free(x);
}

void actr_cleanup() {
	if(actr_n_feature_types > 0) {
		actr_n_feature_types = 0;
		free_vector(actr_feature_type, actr_feature_types);
	}
}

int actr_register_feature_type(actr_feature_type feature_type) {
	if(actr_find_feature_type(feature_type.name) != -1) {
		error(1, "Feature with name “%s” already registered!", feature_type.name);
		return 0;
	}
	#pragma omp critical(actr_feature_types)
	{
		if(actr_feature_types==NULL) {
			actr_feature_types = vector(actr_feature_type, 1);
		} else {
			actr_feature_types = resize_array(actr_feature_type, actr_feature_types, 1, actr_n_feature_types, actr_n_feature_types+1);
		}
		actr_feature_types[++actr_n_feature_types] = feature_type;
	}
	return actr_n_feature_types;
}

int actr_register_bool_feature_type(char * name, int legal_null, int legal_nil) {
	actr_feature_type feature_type;
	feature_type.class = ACTR_FEATURE_BOOL;
	strcpy(feature_type.name, name);
	feature_type.legal_null = legal_null;
	feature_type.legal_nil = legal_nil;
	feature_type.extra = NULL;
	return actr_register_feature_type(feature_type);
}

int actr_register_int_feature_type(char * name, int legal_null, int legal_nil) {
	actr_feature_type feature_type;
	feature_type.class = ACTR_FEATURE_INT;
	strcpy(feature_type.name, name);
	feature_type.legal_null = legal_null;
	feature_type.legal_nil = legal_nil;
	feature_type.extra = NULL;
	return actr_register_feature_type(feature_type);
}

int actr_register_enum_feature_type(char * name, int legal_null, int legal_nil, int nvalues, char ** values) {
	actr_feature_type feature_type;
	feature_type.class = ACTR_FEATURE_ENUM;
	strcpy(feature_type.name, name);
	feature_type.legal_null = legal_null;
	feature_type.legal_nil = legal_nil;
	int total_length = 0;
	int i, o;
	for(i=1;i<=nvalues;i++) {
		total_length += strlen(values[i])+1;
	}
	int id = actr_register_feature_type(feature_type);
	actr_feature_type_enum_extra * extra = malloc(sizeof(actr_feature_type_enum_extra) + sizeof(size_t) * nvalues + total_length);
	size_t * extra_values = (size_t*) &extra->values;
	extra->nvalues = nvalues;
	for(i=0,o=0;i<nvalues;o+=strlen(values[i+1])+1,i++) {
		extra_values[i] = sizeof(char*) * nvalues + o;
		char * addr = (char*) ((size_t) extra_values + (size_t) extra_values[i]);
		strcpy(addr, values[i+1]);
		// printf("%d/%d %lx %lx %s %lx %s\n", i+1	, nvalues, (size_t) extra_values[i], (size_t) addr, addr, (size_t) values[i], values[i]);
	}
	actr_feature_types[id].extra = (void*) extra;
	// printf("Extra addr: %lx\n", extra);
	return id;
}

int actr_register_string_feature_type(char * name, int legal_null, int legal_nil) {
	actr_feature_type feature_type;
	feature_type.class = ACTR_FEATURE_STRING;
	strcpy(feature_type.name, name);
	feature_type.legal_null = legal_null;
	feature_type.legal_nil = legal_nil;
	feature_type.extra = NULL;
	return actr_register_feature_type(feature_type);
}

int actr_compare(actr_feature feature, actr_cue cue) {
	int class = actr_feature_types[feature.type].class;
	if(feature.type != cue.type) {
		warn("Cannot compare feature and cue of different types (%s != %s)!", actr_get_type_name(feature), actr_get_type_name(cue));
		return 0;
	} else if(feature.is_nil) {
		return -100;
	} else if(cue.is_null) {
		return 100;
	} else if(class == ACTR_FEATURE_INT || class == ACTR_FEATURE_ENUM || class == ACTR_FEATURE_BOOL) {
		return actr_get_value(int, feature) - actr_get_value(int, cue);
	} else if(class == ACTR_FEATURE_STRING) {
		return strcmp(actr_get_value(char*, feature), actr_get_value(char*, cue));
	} else {
		stop(1, "Cannot compare unknown types");
	}
}

int actr_get_feature_type_id_for_name(char * name) {
	int i;
	for(i=0;i<actr_n_feature_types;i++) {
		if(!strcmp(actr_feature_types[i].name, name)) {
			return i;
		}
	}
	return -1;
}

int actr_find_feature(int feature_type, actr_memory_item memory_item) {
	int i;
	for(i=1;i<=memory_item.nfeatures;i++) {
		if(memory_item.features[i].type == feature_type) {
			return i;
		}
	}
	return -1;
}

int actr_find_cue(int feature_type, actr_retrieval_item retrieval_item) {
	int i;
	for(i=1;i<=retrieval_item.ncues;i++) {
		if(retrieval_item.cues[i].type == feature_type) {
			return i;
		}
	}
	return -1;
}

int actr_find_and_compare(int feature_type, actr_memory_item memory_item, actr_retrieval_item retrieval_item) {
	int a = actr_find_feature(feature_type, memory_item);
	int b = actr_find_cue(feature_type, retrieval_item);
	if(a == -1) {
		//warn("Memory item doesn’t have feature %s!", actr_feature_types[feature_type].name);
		return -1;
	}

	if(b == -1) {
		//warn("Retrieval item doesn’t have cue %s!", actr_feature_types[feature_type].name);
		return 1;
	}
	return actr_compare(memory_item.features[a], retrieval_item.cues[b]);
}

void actr_parse_value(int feature_type_id, char * value, int * is_nil, void ** out_value) {
	actr_feature_type * ftype = &actr_feature_types[feature_type_id];
	if((ftype)->legal_nil && !strcmp("nil", value)) {
		*is_nil = 1;
		*out_value = NULL;
	} else if((ftype)->legal_null && !strcmp("NULL", value)) {
		*is_nil = 1;
		*out_value = NULL;
	} else if((ftype)->class == ACTR_FEATURE_INT) {
		*is_nil = 0;
		*out_value = malloc(sizeof(int));
		if(sscanf(value, "%d", (int*) *out_value) != 1) {
			warn("“%s” is not a valid integer value for %s!", value, actr_feature_types[feature_type_id].name);
			*((int*) *out_value) = 0;
		}
	} else if((ftype)->class == ACTR_FEATURE_BOOL) {
		*is_nil = 0;
		*out_value = malloc(sizeof(int));
		if(!strcmp(value, "true") || !strcmp(value, "TRUE") || !strcmp(value, "yes") || !strcmp(value, "1")) {
			*((int*) *out_value) = 1;
		} else if(!strcmp(value, "false") || !strcmp(value, "FALSE") || !strcmp(value, "no") || !strcmp(value, "0")) {
			*((int*) *out_value) = 0;
		} else {
			warn("“%s” is not a valid boolean value for %s! Must be true, TRUE, yes, 1, false, FALSE, no, or 0.", value, actr_feature_types[feature_type_id].name);
			*((int*) *out_value) = 0;
		}
	} else if((ftype)->class == ACTR_FEATURE_ENUM) {
		*is_nil = 0;
		*out_value = malloc(sizeof(int));
		int i, j = 0;
		actr_feature_type_enum_extra * extra = (actr_feature_type_enum_extra*) (ftype)->extra;
		size_t * extra_values = (size_t *) &extra->values;
		// printf("Extra addr %s: %lx\n", name, extra);
		for(i = 0; i < extra->nvalues; i++) {
			char * addr = (char*) ((size_t) extra_values + (size_t) extra_values[i]);
			// printf("%d/%d %lx %lx\n", i, extra->nvalues, (size_t) extra_values[i], (size_t) addr);
			//printf("%d/%d %lx %s\n", i, extra->nvalues, (size_t) addr, addr);
			if(!strcmp(value, addr)) {
				j = i+1;
				break;
			}
		}
		if(!j) {
			warn("“%s” is not a valid enum value for %s!", value, actr_feature_types[feature_type_id].name);
			*((int*) *out_value) = 0;
		} else {
			*((int*) *out_value) = j;
		}
	} else if((ftype)->class == ACTR_FEATURE_STRING) {
		*is_nil = 0;
		*out_value = malloc(sizeof(char*) + strlen(value) + 1);
		((size_t*) *out_value)[0] = (size_t) &((size_t*) *out_value)[1];
		strcpy((char*) &((size_t*) *out_value)[1], value);
	} else {
		stop(1, "Don’t know how to parse feature “%s”!", (ftype)->name);
	}
}

void actr_parse_type_and_value(char * name, char * value, int * type, int * is_nil, void ** out_value) {
	int feature_type_id = actr_find_feature_type(name);
	if(feature_type_id == -1) {
		stop(1, "Feature “%s” has not yet been registered!", name);
	}
	*type = feature_type_id;
	actr_parse_value(feature_type_id, value, is_nil, out_value);
}

void actr_parse_cue_type_and_value(char * name, char * value, actr_cue * out) {
	actr_parse_type_and_value(name, value, &out->type, &out->is_null, &out->value);
}

void actr_parse_cue_value(char * value, actr_cue * out) {
	actr_parse_value(out->type, value, &out->is_null, &out->value);
}

void actr_parse_feature_type_and_value(char * name, char * value, actr_feature * out) {
	actr_parse_type_and_value(name, value, &out->type, &out->is_nil, &out->value);
}

void actr_parse_feature_value(char * value, actr_feature * out) {
	actr_parse_value(out->type, value, &out->is_nil, &out->value);
}


void actr_clear_memory_item(actr_memory_item x) {
	int i;
	for(i=1;i<=actr_n_feature_types;i++) {
		actr_clear_feature(x.features[i]);
	}
	free_vector(actr_feature, x.features);
}

void actr_free_memory_item(actr_memory_item * x) {
	actr_clear_memory_item(*x);
	free(x);
}

void actr_clear_retrieval_item(actr_retrieval_item x) {
	int i;
	for(i=1;i<=actr_n_feature_types;i++) {
		actr_clear_cue(x.cues[i]);
	}
	free_vector(actr_cue, x.cues);
}

void actr_free_retrieval_item(actr_retrieval_item * x) {
	actr_clear_retrieval_item(*x);
	free(x);
}

int actr_read_memory_items(FILE * f, actr_memory_item ** read_items, int * items) {
    long fpos = ftell(f);
    char buf[256], buf2[256];
    *items = 0;
    int i, j, is_space = 0, prev_is_space = 0, c;
    // iterate through file byte by byte to count words in first line
    while(!feof(f) && (c=fgetc(f)) != '\n') {
    	is_space = (c == ' ' || c == '\t');
    	if(prev_is_space && !is_space) {
    		(*items)++;
    	}
    	prev_is_space = is_space;
    }
    // go back to start of first line
    fseek(f, fpos, SEEK_SET);
    *read_items = vector(actr_memory_item, *items);
    for(i=1;i<=*items;i++) {
    	(*read_items)[i].features = vector(actr_feature, actr_n_feature_types);
    	(*read_items)[i].nfeatures = 0;
    }
    j = 0;
    while(!feof(f)) {
    	if(fscanf(f, "%s", buf) != 1) {
    		break;
    	}
    	// buf now holds the feature name
    	for(i=1;i<=*items;i++) {
	    	if(!strcmp("name", buf)) {
	    		fscanf(f, "%s", (*read_items)[i].name);
	    	} else if(!strcmp("created", buf)) {
	    		fscanf(f, "%d", &((*read_items)[i].creation));
	    	} else {
	    		int feature_type = actr_find_feature_type(buf);
	    		if(feature_type != -1) {
	    			j=++(*read_items)[i].nfeatures;
		    		(*read_items)[i].features[j].type = feature_type;
			    	if(fscanf(f, "%s", buf2) == 0) {
			    		stop(1, "Unexpected IO error when reading entry %d for feature %s!", i, buf);
			    		break;
			    	}
			    	// buf2 now holds the feature value
			    	actr_parse_feature_value(buf2, &((*read_items)[i].features[j]));
	    		} else {
	    			stop(1, "There is no registered feature type “%s”!", buf);
	    		}
	    	}
    	}
    }
    return 0;
}

int actr_read_retrieval_items(FILE * f, actr_retrieval_item ** read_items, int * items) {
    long fpos = ftell(f);
    char buf[256], buf2[256];
    *items = 0;
    int i, j, is_space = 0, prev_is_space = 0, c;
    // iterate through file byte by byte to count words in first line
    while(!feof(f) && (c=fgetc(f)) != '\n') {
    	is_space = (c == ' ' || c == '\t');
    	if(prev_is_space && !is_space) {
    		(*items)++;
    	}
    	prev_is_space = is_space;
    }
    // go back to start of first line
    fseek(f, fpos, SEEK_SET);
    *read_items = vector(actr_retrieval_item, *items);
    for(i=1;i<=*items;i++) {
    	(*read_items)[i].cues = vector(actr_cue, actr_n_feature_types);
    	(*read_items)[i].ncues = 0;
    }
    while(!feof(f)) {
    	if(fscanf(f, "%s", buf) != 1) {
    		break;
    	}
    	// buf now holds the feature name
    	for(i=1;i<=*items;i++) {
	    	if(!strcmp("moment", buf)) {
	    		fscanf(f, "%d", &(*read_items)[i].moment);
	    	} else {
	    		j=++(*read_items)[i].ncues;
	    		int feature_type = actr_find_feature_type(buf);
	    		(*read_items)[i].cues[j].type = feature_type;
		    	if(fscanf(f, "%s", buf2) == 0) {
		    		stop(1, "Unexpected IO error when reading entry %d for feature %s!", i, buf);
		    		break;
		    	}
		    	// buf2 now holds the feature value
		    	actr_parse_cue_value(buf2, &(*read_items)[i].cues[j]);
	    	}
    	}
    }
    return 0;
}

void actr_read_feature_types(FILE * f) {
	char buf1[256], buf2[256], buf3[256], buf4[256];
	int legal_nil, legal_null;
	while(!feof(f)) {
		if(fscanf(f, "%s", buf1) != 1) {
			break;
		}
		if(fscanf(f, "%s%s%s", buf2, buf3, buf4) != 3) {
			stop(1, "Feature type declaration for `%s` must have format: type name nil:[yes|no] null:[yes|no] [#values value1 value2 ...]", buf1);
		}
		if(!strcmp(buf3, "nil:yes")) {
			legal_nil = 1;
		} else if(!strcmp(buf3, "nil:no")) {
			legal_nil = 0;
		} else {
			stop(1, "Illegal value for yes/no nil: %s", buf3);
		}
		if(!strcmp(buf4, "null:yes")) {
			legal_null = 1;
		} else if(!strcmp(buf4, "null:no")) {
			legal_null = 0;
		} else {
			stop(1, "Illegal value for yes/no null: %s", buf4);
		}
		if(!strcmp(buf1, "int")) {
			actr_register_int_feature_type(buf2, legal_null, legal_nil);
		} else if(!strcmp(buf1, "bool")) {
			actr_register_bool_feature_type(buf2, legal_null, legal_nil);
		} else if(!strcmp(buf1, "string")) {
			actr_register_string_feature_type(buf2, legal_null, legal_nil);
		} else if(!strcmp(buf1, "enum")) {
			int nvalues, i;
			if(fscanf(f, "%d", &nvalues) != 1) {
				stop(1, "Illegal integer value for number of enum values!");
			}
			const size_t valbufsize = 256;
			char valbuf[nvalues*valbufsize];
			char ** ivalues = vector(char*, nvalues);
			for(i=1;i<=nvalues;i++) {
				ivalues[i] = &valbuf[(i-1)*valbufsize];
				fscanf(f, "%s", ivalues[i]);
			}
			actr_register_enum_feature_type(buf2, legal_null, legal_nil, nvalues, ivalues);
			free_vector(char*, ivalues);
		} else  {
			stop(1, "Illegal feature type class “%s” (must be int, enum or string).", buf1);
		}
	}
}

typedef struct {
	int runs;
	double F; // latency factor
	double G; // total source activation
	double ans; // activation noise parameter for logistics distribution
	double mas; // fan parameter
	double d; // base level decay parameter
	double match_penalty; // match penalty
} actr_params;

void actr_compute_base_level_activations(actr_params params, double moment, int * moments, int n_moments, int n_trials, int n_items, int ** history, double ** base_levels) {
	int i, j, k;

	// time since last retrieval for each retrieval (converted from milliseconds) [actr.r, l. 29]
	double * tj = vector(double, n_moments);
	int * is_past = vector(int, n_moments);
	int num_past_retrievals = 0;
	for(i=1;i<=n_moments;i++) {
		tj[i] = (moment - moments[i]) / 1000.0;
		is_past[i] = tj[i] > 0.0;
		num_past_retrievals++;
	}

	// the decay function tj^-d [act.r, l. 37]
	double ** decay = matrix(double, n_trials, n_moments);
	for(i=1;i<=n_trials;i++) {
		for(j=1;j<=n_moments;j++) {
			if(is_past[j]) {
				decay[i][j] = pow(tj[j], -params.d);
			} else {
				decay[i][j] = 0.0;
			}
		}
	}
	double ** activations = matrix(double, n_trials, n_moments);

	// compute base level activation for each item (for each trial) [actr.r, l. 40]
	for(i=1;i<=n_items;i++) {
		// activations = (past == i) * decay
		for(j=1;j<=n_trials;j++) {
			for(k=1;k<=n_moments;k++) {
				activations[j][k] = decay[j][k] * is_past[k] * (history[j][k] == i);
			}
		}

		for(j=1;j<=n_trials;j++) {
			double b = 0.0;
			for(k=1;k<=n_moments;k++) {
				b += activations[j][k];
			}
			b = log(b);
			if(isinf(b)) b = 0.0;
			base_levels[i][j] = b;
		}
	}

	free_matrix(double, activations);
	free_matrix(double, decay);
	free_vector(int, is_past);
	free_vector(double, tj);
}

void actr_compute_overall_base_level_activation(actr_params params, double moment, int * moments, int n_moments, double * base_level) {
	int i, j, k;

	*base_level = 0.0;

	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		if(moments[j] < moment) {
			*base_level += exp(-(moment - moments[j])/1000.0 * moment_decay);
		}
	}

}


void actr_compute_overall_base_level_activation_by_word(actr_params params, double moment, int * moments, int n_moments, int ** history, int n_words, actr_memory_item * items, double * base_level) {
	int i, j, k;

	for(i=0;i<=n_words;i++) {
		base_level[i] = 0.0;
	}

	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		int belonging_memory_item = history[1][j];
		if(moments[j] < moment && belonging_memory_item) {
			base_level[items[belonging_memory_item].trigger] += exp(-(moment - moments[j])/1000.0 * moment_decay);
		}
	}

}

void actr_compute_overall_base_level_activation_by_memory_item(actr_params params, double moment, int * moments, int n_moments, int ** history, int n_items, actr_memory_item * items, double * base_level) {
	int i, j, k;

	for(i=1;i<=n_items;i++) {
		base_level[i] = 0.0;
	}

	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		int belonging_memory_item = history[1][j];
		if(moments[j] < moment && belonging_memory_item) {
			base_level[belonging_memory_item] += exp(-(moment - moments[j])/1000.0 * moment_decay);
		}
	}

}



void actr_compute_avg_overall_base_level_activation(actr_params params, double moment1, double moment2, int * moments, int n_moments, double * base_level) {
	int i, j, k;

	if(fabs(moment2 - moment1) < 0.01) {
		actr_compute_overall_base_level_activation(params, moment1, moments, n_moments, base_level);
		return;
	}

	*base_level = 0.0;

	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		if(moments[j] < moment2) {
			*base_level += exp(-(moment2 - moments[j])/1000.0 * moment_decay) / moment_decay * 1000.0;
			if(moments[j] <= moment1) {
				*base_level -= exp(-(moment1 - moments[j])/1000.0 * moment_decay) / moment_decay * 1000.0;
			} else {
				*base_level -= 1.0 / moment_decay;
			}
		}
	}

	*base_level /= (moment1 - moment2);

}


void actr_compute_avg_overall_base_level_activation_by_word(actr_params params, double moment1, double moment2, int * moments, int n_moments, int ** history, int n_words, actr_memory_item * items, double * base_level) {
	int i, j, k;


	if(fabs(moment2 - moment1) < 0.01) {
		actr_compute_overall_base_level_activation_by_word(params, moment1, moments, n_moments, history, n_words, items, base_level);
		return;
	}

	for(i=0;i<=n_words;i++) {
		base_level[i] = 0.0;
	}


	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		int belonging_memory_item = history[1][j];
		//printf("%d ", belonging_memory_item);
		if(moments[j] < moment2 && belonging_memory_item) {
			double act = exp(-(moment2 - moments[j]) / 1000.0 * moment_decay) / moment_decay * 1000.0;
			if(moments[j] <= moment1) {
				act -= exp(-(moment1 - moments[j]) / 1000.0 * moment_decay) / moment_decay * 1000.0;
			} else {
				act -= 1.0 / moment_decay;
			}
			base_level[items[belonging_memory_item].trigger] += act / (moment1 - moment2);
		}
	}
	//printf("\n");

}


void actr_compute_avg_overall_base_level_activation_by_memory_item(actr_params params, double moment1, double moment2, int * moments, int n_moments, int ** history, int n_items, actr_memory_item * items, double * base_level) {
	int i, j, k;


	if(fabs(moment2 - moment1) < 0.01) {
		actr_compute_overall_base_level_activation_by_memory_item(params, moment1, moments, n_moments, history, n_items, items, base_level);
		return;
	}

	for(i=1;i<=n_items;i++) {
		base_level[i] = 0.0;
	}


	// the decay function tj^-d [act.r, l. 37]
	for(j=1;j<=n_moments;j++) {
		double moment_decay = params.d;
		int belonging_memory_item = history[1][j];
		//printf("%d ", belonging_memory_item);
		if(moments[j] < moment2 && belonging_memory_item) {
			double act = exp(-(moment2 - moments[j]) / 1000.0 * moment_decay) / moment_decay * 1000.0;
			if(moments[j] <= moment1) {
				act -= exp(-(moment1 - moments[j]) / 1000.0 * moment_decay) / moment_decay * 1000.0;
			} else {
				act -= 1.0 / moment_decay;
			}
			base_level[belonging_memory_item] += act / (moment1 - moment2);
		}
	}
	//printf("\n");

}

typedef struct {
	int item;
	int retrieval_trigger;
	int memory_trigger;
	int match;
	double latency;
	double activation;
	double base_activation, memory_strength, penalty;
	double noisy_activation;
} actr_retrieval_result;



void actr_retrieve(actr_params params, int moment, actr_retrieval_item cues, actr_memory_item * items, int n_items, int * moments, int n_moments, int ** history, int n_history, int n_trials, actr_retrieval_result * result, actr_retrieval_result ** all_results, RANSEED_TYPE * seed, double timeout) {
	int i, j, k;

	double ** base_levels = matrix(double, n_items, n_trials);
	// [actr.r, l. 60]
	int * is_retrieval_cue = vector(int, actr_n_feature_types), num_cues = 0;
	for(i=1;i<=actr_n_feature_types;i++) {
		j = actr_find_cue(i, cues);
		if(j == -1 || cues.cues[j].is_null) {
			is_retrieval_cue[i] = 0;
		} else {
			is_retrieval_cue[i] = 1;
			num_cues++;
		}
	}

	// compute base level activations [actr.r, l. 62]

	actr_compute_base_level_activations(params, moment, moments, n_moments, n_trials, n_items, history, base_levels);

	//printf_matrix(base_levels, n_items, n_trials, "%8.3lf");

	int cat_id = actr_get_feature_type_id_for_name("cat");
	if(cat_id == -1) {
		stop(1, "cat feature not defined!");
	}

	// compute the match between items and retrieval cues (a boolean matrix) [actr.r, l. 65]

	int num_features = actr_n_feature_types;
	int ** match = matrix(int, n_items, num_features);
	int ** exists = matrix(int, n_items, num_features);
	int ** matches_category = matrix(int, n_items, n_trials);
	for(i=1;i<=n_items;i++) {
		for(j=1;j<=num_features;j++) {
			match[i][j] = actr_find_and_compare(j, items[i], cues) == 0;
			// checks which items exist at this moment [actr.r, l. 73]
			exists[i][j] = 1; // items[i].creation <= moment;
		}
		// checks which items match category cue
		for(j=1;j<=n_trials;j++) {
			matches_category[i][j] = match[i][cat_id];
		}
	}

	// compute fan for each feature: number of existing items matching each feature [actr.r, l. 81]
	double * strength = vector(double, num_features);
	for(j=1;j<=num_features;j++) {
		double fan = 0.0;
		for(i=1;i<=n_items;i++) {
			fan += match[i][j] && exists[i][j];
		}
		strength[j] = params.mas - log(fan);
		if(0 && isinf(strength[j])) {
			printf("Fan %d: %lf %lf\n", j, fan, log(fan));
			for(i=1;i<=n_items;i++) {
				printf("%d %d %d>%d, ", match[i][j], exists[i][j], items[i].creation, moment);
			}
			printf("\n");
		}
	}

	// printf_vector(strength, num_features, "%10.2lf");
	// printf_vector(is_retrieval_cue, num_features, "%d ");

	// compute source activation available for each cue (source spread over
    // cues) and multiply by the fan (S * W in act-r equation). [actr.r, l. 86]

	double * cue_weights = vector(double, num_features), sum_cue_weights = 0.0;
	for(i=1;i<=num_features;i++) {
		cue_weights[i] = 1.0;
		sum_cue_weights += cue_weights[i];
	}


	double ** sw = matrix(double, n_items, num_features);
	double W;
	for(i=1;i<=n_items;i++) {
		for(j=1;j<=num_features;j++) {
			W = params.G * cue_weights[j] / sum_cue_weights;
			sw[i][j] = strength[j] * W;
		}
	}

	// printf_matrix(sw, n_items, num_features, "%10.2lf");

    // compute extra activation for each item; sum must ignore NA's because
    // those correspond to features that are not retrieval cues. [actr.r, l. 100]
	// compute mismatch penalty [actr.r, l. 104]
	double * extra = vector(double, n_items);
	double * penalty = vector(double, n_items);
	for(i=1;i<=n_items;i++) {
		extra[i] = 0.0;
		penalty[i] = 0.0;
		for(j=1;j<=num_features;j++) {
			if(is_retrieval_cue[j] && match[i][j]) {
				extra[i] += sw[i][j];
			}
			if(!match[i][j] && is_retrieval_cue[j]) {
				penalty[i] += params.match_penalty;
			}
		}
	}



	//printf_vector(penalty, n_items, "%10.3lf");

	// compute activation boost/penalty [actr.r, l. 116]
	double ** activation = matrix(double, n_items, n_trials);
	double ** noisy_activation = matrix(double, n_items, n_trials);
	for(i=1;i<=n_items;i++) {
		double activation_adjustment = extra[i] + penalty[i];
		for(j=1;j<=n_trials;j++) {
			activation[i][j] = activation_adjustment + base_levels[i][j];
			noisy_activation[i][j] = rlogistic(activation[i][j], params.ans, seed);
		}
	}

	//printf_matrix(noisy_activation, n_items, n_trials, "%10.3lf");

	// make items that don't exist yet, or that don't match category cues,  have activation of -999 [actr.r, l. 132]
	double ** retrieval_latency = matrix(double, n_items, n_trials);

	for(i=1;i<=n_items;i++) {
		for(j=1;j<=n_trials;j++) {
			retrieval_latency[i][j] = params.F * exp(-noisy_activation[i][j]) * 1000.0;
		}
	}

	for(j=1;j<=n_trials;j++) {
		int winner = 0;
		for(i=1;i<=n_items;i++) {
			int exists = items[i].creation <= moment;
			if(exists && matches_category[i][j] && retrieval_latency[i][j] <= timeout && (winner == 0 || retrieval_latency[i][j] < retrieval_latency[winner][j])) {
				winner = i;
			}
			if(all_results) {
				all_results[j][i].retrieval_trigger = cues.trigger;
				all_results[j][i].item = i;
				all_results[j][i].memory_trigger = items[i].trigger;
				all_results[j][i].match = matches_category[i][j];
				all_results[j][i].latency = retrieval_latency[i][j];
				all_results[j][i].activation = activation[i][j];
				all_results[j][i].base_activation = base_levels[i][j];
				all_results[j][i].memory_strength = extra[i];
				if(isinf(extra[i])) {
					printf("Item %d extra: %lf\n", i, extra[i]);
				}
				all_results[j][i].penalty = penalty[i];
				all_results[j][i].noisy_activation = noisy_activation[i][j];
			}
		}

		result[j].retrieval_trigger = cues.trigger;

		if(winner) {
			result[j].item = winner;
			result[j].memory_trigger = items[winner].trigger;
			result[j].match = 1;
			result[j].latency = retrieval_latency[winner][j];
			result[j].activation = activation[winner][j];
			result[j].base_activation = base_levels[winner][j];
			result[j].memory_strength = extra[winner];
			result[j].penalty = penalty[winner];
			result[j].noisy_activation = noisy_activation[winner][j];
		} else {
			result[j].item = 0;
			result[j].memory_trigger = 0;
			result[j].match = 0;
			result[j].latency = timeout;
			result[j].activation = NAN;
			result[j].base_activation = NAN;
			result[j].memory_strength = NAN;
			result[j].penalty = NAN;
			result[j].noisy_activation = NAN;
		}
	}

    // free memory
    free_matrix(double, retrieval_latency);
	free_vector(int, is_retrieval_cue);
	free_matrix(double, noisy_activation);
	free_matrix(double, activation);
	free_vector(double, penalty);
	free_vector(double, extra);
	free_matrix(double, sw);
	free_vector(double, cue_weights);
	free_vector(double, strength);
	free_matrix(int, matches_category);
	free_matrix(int, exists);
	free_matrix(int, match);


	free_matrix(double, base_levels);
}



#define actr_fprintf_val(feature) \
    if(actr_get_type_class(feature) == ACTR_FEATURE_INT) { \
		fprintf(f, "%d", actr_get_value(int, feature)); \
	} else if(actr_get_type_class(feature) == ACTR_FEATURE_BOOL) { \
		fprintf(f, actr_get_value(int, feature) ? "yes" : "no"); \
	} else if(actr_get_type_class(feature) == ACTR_FEATURE_ENUM) { \
		fprintf(f, "%s(%d)", actr_get_label(feature), actr_get_value(int, feature)); \
	} else if(actr_get_type_class(feature) == ACTR_FEATURE_STRING) { \
		fprintf(f, "“%s”", actr_get_value(char*, feature)); \
	} else { \
		fprintf(f, "?"); \
	}

void actr_fprintf_feature_value(FILE * f, actr_feature feature) {
	if(actr_is_nil(feature)) {
		fprintf(f, "nil");
	} else {
		actr_fprintf_val(feature);
	}
}

void actr_fprintf_cue_value(FILE * f, actr_cue cue) {
	if(actr_is_null(cue)) {
		fprintf(f, "NULL");
	} else {
		actr_fprintf_val(cue);
	}
}

#undef actr_fprintf_val

void actr_fprintf_features(FILE * f, actr_memory_item item) {
	int j, k;
	for(j=1;j<=actr_n_feature_types;j++) {
		if(j>1) fprintf(f, ", ");
		fprintf(f, "%s=", actr_feature_types[j].name);
		k = actr_find_feature(j, item);
		if(k == -1) {
			fprintf(f, "NA");
		} else {
			actr_fprintf_feature_value(f, item.features[k]);
		}
	}
}

void actr_fprintf_memory_item(FILE * f, actr_memory_item item) {
	fprintf(f, "%s: created=%d, ", item.name, item.creation);
	actr_fprintf_features(f, item);
}

void actr_fprintf_cues(FILE * f, actr_retrieval_item item) {
	int j, k;
	for(j=1;j<=actr_n_feature_types;j++) {
		if(j>1) fprintf(f, ", ");
		fprintf(f, "%s=", actr_feature_types[j].name);
		k = actr_find_cue(j, item);
		if(k == -1) {
			fprintf(f, "NA");
		} else {
			actr_fprintf_cue_value(f, item.cues[k]);
		}
	}
}

void actr_fprintf_retrieval_item(FILE * f, actr_retrieval_item item) {
	fprintf(f, "moment=%d, ", item.moment);
	actr_fprintf_cues(f, item);
}

typedef struct {
	int runs;
	actr_memory_item * items;
	int n_items;
	actr_retrieval_item * retrievals;
	actr_retrieval_result ** retrieval_results;
	int n_retrievals;
	int n_retrievals_complete;
	int n_retrieval_results;
	int ** history;
	int n_history;
	int * moments;
	int n_moments;
} actr_trial;

#ifdef CHECK_MEMORY
#define actr_new_trial(runs) track_array_memory(actr_new_trial_(runs), __FILE__, __LINE__)
#define actr_free_trial(trial) { actr_free_trial_(trial); untrack_array_memory(trial, __FILE__, __LINE__); }
#else
#define actr_new_trial(runs) actr_new_trial_(runs)
#define actr_free_trial(trial) actr_free_trial_(trial)
#endif

actr_trial * actr_new_trial_(int runs) {
	actr_trial * trial = malloc(sizeof(actr_trial));
	trial->runs = runs;
	trial->n_items = 0;
	trial->items = NULL;
	trial->n_retrievals = 0;
	trial->n_retrievals_complete = 0;
	trial->retrievals = NULL;
	trial->n_retrieval_results = 0;
	trial->retrieval_results = NULL;
	trial->n_history = 0;
	trial->history = NULL;
	trial->n_moments = 0;
	trial->moments = NULL;
	return trial;
}

void actr_free_trial_(actr_trial * trial) {
	if(trial->items != NULL) {
		free_vector(actr_memory_item, trial->items);
	}
	if(trial->retrievals != NULL) {
		free_vector(actr_retrieval_item, trial->retrievals);
	}
	if(trial->retrieval_results != NULL) {
		free_matrix(actr_retrieval_result, trial->retrieval_results);
	}
	if(trial->history != NULL) {
		free_matrix(int, trial->history);
	}
	if(trial->moments != NULL) {
		free_vector(int, trial->moments);
	}
	free(trial);
}

actr_trial * actr_duplicate_trial(actr_trial * src) {
	actr_trial * ret = actr_new_trial(src->runs);
	*ret = *src;
	ret->items = duplicate_vector(actr_memory_item, src->items, ret->n_items);
	ret->retrievals = duplicate_vector(actr_retrieval_item, src->retrievals, ret->n_retrievals);
	ret->retrieval_results = duplicate_matrix(actr_retrieval_result, src->retrieval_results, ret->runs, ret->n_retrieval_results);
	ret->history = duplicate_matrix(int, src->history, ret->runs, ret->n_history);
	ret->moments = duplicate_vector(int, src->moments, ret->n_moments);
	return ret;
}

void actr_add_moment(actr_trial * trial, int moment) {
	if(trial->moments == NULL) {
		trial->moments = vector(int, 1);
	} else {
		trial->moments = resize_vector(int, trial->moments, trial->n_moments, trial->n_moments+1);
	}
	trial->moments[++trial->n_moments] = moment;
}

void actr_add_history_all(actr_trial * trial, int item) {
	int i;
	if(trial->history == NULL) {
		trial->history = matrix(int, trial->runs, 1);
	} else {
		trial->history = resize_matrix(int, trial->history, trial->runs, trial->runs, trial->n_history, trial->n_history+1);
	}
	trial->n_history++;
	for(i=1;i<=trial->runs;i++) {
		trial->history[i][trial->n_history] = item;
	}
}

void actr_add_history_vec(actr_trial * trial, int * item) {
	int i;
	if(trial->history == NULL) {
		trial->history = matrix(int, trial->runs, 1);
	} else {
		trial->history = resize_matrix(int, trial->history, trial->runs, trial->runs, trial->n_history, trial->n_history+1);
	}
	trial->n_history++;
	for(i=1;i<=trial->runs;i++) {
		trial->history[i][trial->n_history] = item[i];
	}
}

void actr_add_memory_item(actr_trial * trial, actr_memory_item item) {
	// printf("Adding ");
	// actr_fprintf_memory_item(stdout, item);
	// printf("...\n");
	if(trial->items == NULL) {
		trial->items = vector(actr_memory_item, 1);
	} else {
		trial->items = resize_vector(actr_memory_item, trial->items, trial->n_items, trial->n_items+1);
	}
	actr_add_moment(trial, item.creation);
	trial->n_items++;
	actr_add_history_all(trial, trial->n_items);
	trial->items[trial->n_items] = item;
}

void actr_begin_retrieval(actr_trial * trial, actr_retrieval_item item) {
	// printf("Start retrieval for ");
	// actr_fprintf_retrieval_item(stdout, item);
	// printf("...\n");
	if(trial->retrievals == NULL) {
		trial->retrievals = vector(actr_retrieval_item, 1);
	} else {
		trial->retrievals = resize_vector(actr_retrieval_item, trial->retrievals, trial->n_retrievals, trial->n_retrievals+1);
	}
	actr_add_moment(trial, item.moment);
	trial->retrievals[++trial->n_retrievals] = item;
}

void actr_end_retrieval(actr_trial * trial, actr_retrieval_result result) {
	int i;
	if(trial->retrieval_results == NULL) {
		trial->retrieval_results = matrix(actr_retrieval_result, trial->runs, 1);
	} else {
		trial->retrieval_results = resize_matrix(actr_retrieval_result, trial->retrieval_results, trial->runs, trial->runs, trial->n_retrieval_results, trial->n_retrieval_results+1);
	}
	trial->n_retrieval_results++;
	for(i=1; i<=trial->runs; i++) {
		trial->retrieval_results[i][trial->n_retrieval_results] = result;
	}
	actr_add_history_all(trial, result.item);
}

void actr_end_retrieval_vec(actr_trial * trial, actr_retrieval_result * result) {
	int i;
	int * winners = vector(int, trial->runs);
	if(trial->retrieval_results == NULL) {
		trial->retrieval_results = matrix(actr_retrieval_result, trial->runs, 1);
	} else {
		trial->retrieval_results = resize_matrix(actr_retrieval_result, trial->retrieval_results, trial->runs, trial->runs, trial->n_retrieval_results, trial->n_retrieval_results+1);
	}
	trial->n_retrieval_results++;
	for(i=1; i<=trial->runs; i++) {
		winners[i] = result[i].item;
		trial->retrieval_results[i][trial->n_retrieval_results] = result[i];
	}
	actr_add_history_vec(trial, winners);
	free_vector(int, winners);
}

void actr_add_memory_for(actr_trial * trial, int n_items, actr_memory_item * items, int k, int moment, int add_fail_dummy) {
	int i;
	int item_found = 0;
	actr_memory_item item, dummy;
	dummy.features = NULL;
	dummy.nfeatures = 0;
	for(i=1;i<=n_items;i++) {
		if(items[i].creation == k) {
			item = items[i];
			item.trigger = k;
			item.creation = moment;
			actr_add_memory_item(trial, item);
			item_found = 1;
		}
	}
	if(!item_found && add_fail_dummy) {
		dummy.trigger = k;
		dummy.creation = moment;
		actr_add_memory_item(trial, dummy);
	}
}

void actr_prepare_retrievals_for(actr_trial * trial, int n_retrievals, actr_retrieval_item * items, int k, int moment) {
	int i;
	actr_retrieval_item item;
	for(i=1;i<=n_retrievals;i++) {
		if(items[i].moment == k) {
			item = items[i];
			item.trigger = k;
			item.moment = moment;
			actr_begin_retrieval(trial, item);
		}
	}
}

void actr_perform_retrievals(actr_params params, RANSEED_TYPE * seed, actr_trial * trial, actr_retrieval_result *** output, actr_retrieval_result **** all_output, int * retrieved_items, double timeout) {
	int i, j, n, m, k;
	actr_retrieval_result * result = vector(actr_retrieval_result, params.runs);
	actr_retrieval_result ** all_results = matrix(actr_retrieval_result, params.runs, trial->n_items);
	n = trial->n_retrievals - trial->n_retrievals_complete;
	if(retrieved_items != NULL) {
		*retrieved_items = n;
	}
	if(output != NULL) {
		*output = matrix(actr_retrieval_result, params.runs, n);
	}
	if(all_output != NULL) {
		*all_output = (actr_retrieval_result***) array(actr_retrieval_result, params.runs, n, trial->n_items);
	}
	for(i=trial->n_retrievals_complete+1,m=0;i<=trial->n_retrievals;i++) {
		m++;
		actr_retrieve(params, trial->retrievals[i].moment, trial->retrievals[i], trial->items, trial->n_items, trial->moments, trial->n_moments, trial->history, trial->n_history, trial->runs, result, all_results, seed, timeout);
		actr_end_retrieval_vec(trial, result);
		// printf("Retrieved ");
		// actr_fprintf_retrieval_item(stdout, item);
		// printf(" in");
		if(output != NULL || all_output != NULL) {
			for(j=1;j<=params.runs;j++) {
				if(output != NULL) (*output)[j][m] = result[j];
				if(all_output != NULL) for(k=1;k<=trial->n_items;k++) {
					(*all_output)[j][m][k] = all_results[j][k];
				}
			}
		}
		// printf(".\n");
	}
	trial->n_retrievals_complete = trial->n_retrievals;
	free_vector(actr_retrieval_result, result);
	free_matrix(actr_retrieval_result, all_results);
}

void actr_serial_retrievals_for(actr_params params, RANSEED_TYPE * seed, actr_trial * trial, int n_retrievals, actr_retrieval_item * items, int tk, int moment, actr_retrieval_result *** output, actr_retrieval_result **** all_output, int * retrieved_items, double timeout) {
	if(params.runs != 1) stop(1, "Serial retrievals only work with actr runs = 1!");
	int i, j, n, m, k;
	actr_retrieval_item item;
	actr_retrieval_result * result = vector(actr_retrieval_result, params.runs);
	actr_retrieval_result ** all_results = matrix(actr_retrieval_result, params.runs, trial->n_items);
	for(i=1,n=0;i<=n_retrievals;i++) {
		if(items[i].moment == tk) {
			n++;
		}
	}
	if(retrieved_items != NULL) {
		*retrieved_items = n;
	}
	if(output != NULL) {
		*output = matrix(actr_retrieval_result, params.runs, n);
	}
	if(all_output != NULL) {
		*all_output = (actr_retrieval_result***) array(actr_retrieval_result, 3, params.runs, n, trial->n_items);
	}
	double * t = vector(double, params.runs);
	for(i=1;i<=params.runs;i++) {
		t[i] = moment;
	}
	for(i=1,m=0;i<=n_retrievals;i++) {
		if(items[i].moment == tk) {
			m++;
			item = items[i];
			// ONLY WORKS WITH RUNS=1 AT THE MOMENT BECAUSE actr_begin_retrieval is not vectorized!
			item.trigger = tk;
			item.moment = (int) t[1];
			actr_begin_retrieval(trial, item);
			actr_retrieve(params, item.moment, item, trial->items, trial->n_items, trial->moments, trial->n_moments, trial->history, trial->n_history, trial->runs, result, all_results, seed, timeout);
			actr_end_retrieval_vec(trial, result);
			if(output != NULL || all_output != NULL) {
				for(j=1;j<=params.runs;j++) {
					if(output != NULL) (*output)[j][m] = result[j];
					if(all_output != NULL) for(k=1;k<=trial->n_items;k++) {
						(*all_output)[j][m][k] = all_results[j][k];
					}
				}
			}
		}
	}
	free_vector(actr_retrieval_result, result);
	free_matrix(actr_retrieval_result, all_results);
	free_vector(double, t);
}

#endif
