/**
 * @file sjqbn_util.h
 *
**/

#ifndef SJQBN_UTIL_H
#define SJQBN_UTIL_H

#define SJQBN_DATASET_MAX 10

typedef struct sjqbn_properties_t sjqbn_properties_t;

/** The SJQBN a dataset's working structure. */
typedef struct sjqbn_dataset_t {
	/** tracking netcdf id **/
        int ncid;

	/** Number of x(lon) points */
	int nx;
	/** Number of y(lat) points */
	int ny;
	/** Number of z(dep) points */
	int nz;

	/** list of longitudes **/
	float *longitudes;
	/** list of latitudes **/
	float *latitudes;
	/** list of depths **/
	float *depths;

	int vp_varid;
	int vs_varid;
	int rho_varid;

	int elems;
	float *vp_buffer;
	float *vs_buffer;
	float *rho_buffer; 

/* flag to show if data i read in memory */
        int in_memory;

} sjqbn_dataset_t;

typedef struct sjqbn_pt_info_t {
        float lon;
        float lat;
        float dep;
        int lon_idx;
        int lat_idx;
        int dep_idx;
        float lon_percent;
        float lat_percent;
        float dep_percent;
} sjqbn_pt_info_t;


/* utilitie functions */
sjqbn_dataset_t *make_a_sjqbn_dataset(char *datadir, char *datafile, int tooBig);
int free_sjqbn_dataset(sjqbn_dataset_t *data);

int get_one_property(sjqbn_dataset_t *dataset, sjqbn_pt_info_t *pt, sjqbn_properties_t *data);
void get_interp_property(sjqbn_dataset_t *dataset, sjqbn_pt_info_t *pt, sjqbn_properties_t *data);

#endif

