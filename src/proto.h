/*! \file proto.h
 *  \brief this file contains all function prototypes of the code
 */

#ifndef ALLVARS_H
#include "allvars.h"
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

void   allocate_commbuffers(void);
void   allocate_memory(void);
void   begrun(void);
int    blockpresent(enum iofields blocknr);
void   catch_abort(int sig);
void   catch_fatal(int sig);
void   check_omega(void);
void   close_outputfiles(void);
int    compare_key(const void *a, const void *b);
void   compute_accelerations(int mode);
void   compute_global_quantities_of_system(void);
void   compute_potential(void);
int    dens_compare_key(const void *a, const void *b);

void density(void);
void density_evaluate(int, int);

void   distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);
double dmax(double, double);
double dmin(double, double);
void   do_box_wrapping(void);

void   domain_Decomposition(void); 
int    domain_compare_key(const void *a, const void *b);
int    domain_compare_key(const void *a, const void *b);
int    domain_compare_toplist(const void *a, const void *b);
void   domain_countToGo(void);
void   domain_decompose(void);
void   domain_determineTopTree(void);
void   domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count);
void   domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void   domain_findExtent(void);
int    domain_findSplit(int cpustart, int ncpu, int first, int last);
void   domain_shiftSplit(void);
void   domain_sumCost(void);
void   domain_topsplit(int node, peanokey startkey);
void   domain_topsplit_local(int node, peanokey startkey);

void   dump_particles(void);
void   empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);
void   endrun(int);

void   ewald_corr(double dx, double dy, double dz, double *fper);
void   ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void   ewald_init(void);
double ewald_pot_corr(double dx, double dy, double dz);
double ewald_psi(double x[3]);

void   fill_Tab_IO_Labels(void);
void   fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
void   find_dt_displacement_constraint(double hfac);
int    find_files(char *fname);

void   force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount, int *nextfree);
void   force_exchange_pseudodata(void);
void   force_flag_localnodes(void);
void   force_insert_pseudo_particles(void);
void   force_setupnonrecursive(int no);
void   force_treeallocate(int maxnodes, int maxpart); 
int    force_treebuild(int npart);
int    force_treebuild_single(int npart);
int    force_treeevaluate(int target, int mode, double *ewaldcountsum);
int    force_treeevaluate_direct(int target, int mode);
int    force_treeevaluate_ewald_correction(int target, int mode, double pos_x, double pos_y, double pos_z, double aold);
void   force_treeevaluate_potential(int target, int type);
void   force_treeevaluate_potential_shortrange(int target, int mode);
int    force_treeevaluate_shortrange(int target, int mode);
void   force_treefree(void);
void   force_treeupdate_pseudos(void);
void   force_update_hmax(void);
void   force_update_len(void);
void   force_update_node(int no, int flag);
void   force_update_node_hmax_local(void);
void   force_update_node_hmax_toptree(void);
void   force_update_node_len_local(void);
void   force_update_node_len_toptree(void);
void   force_update_node_recursive(int no, int sib, int father);
void   force_update_pseudoparticles(void);
void   force_update_size_of_parent_node(int no);

void   free_memory(void);

int    get_bytes_per_blockelement(enum iofields blocknr);
void   get_dataset_name(enum iofields blocknr, char *buf);
int    get_datatype_in_block(enum iofields blocknr);
int    get_particles_in_block(enum iofields blocknr, int *typelist);
double get_random_number(int id);
int    get_timestep(int p, double *a, int flag);
int    get_values_per_blockelement(enum iofields blocknr);

int    imax(int, int);
int    imin(int, int);

void   init(void);
void   init_peano_map(void);

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

int    ngb_clear_buf(FLOAT searchcenter[3], FLOAT hguess, int numngb);
void   ngb_treeallocate(int npart);
void   ngb_treebuild(void);
int    ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode);
int    ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hguess, int *startnode);
void   ngb_treefree(void);
void   ngb_treesearch(int);
void   ngb_treesearch_pairs(int);
void   ngb_update_nodes(void);

void   open_outputfiles(void);

peanokey peano_hilbert_key(int x, int y, int z, int bits);
void   peano_hilbert_order(void);

double pow(double, double);  /* on some old DEC Alphas, the correct prototype for pow() is missing, even when math.h is included */

void   read_file(char *fname, int readTask, int lastTask);
void   read_header_attributes_in_hdf5(char *fname);
void   read_ic(char *fname);
int    read_outputlist(char *fname);
void   read_parameter_file(char *fname);
void   readjust_timebase(double TimeMax_old, double TimeMax_new);

void   reorder_gas(void);
void   reorder_particles(void);
void   restart(int mod);

double second(void);

void   set_random_numbers(void);
void   set_softenings(void);
void   set_units(void);

void   setup_smoothinglengths(void);
void   statistics(void);
void   terminate_processes(void);
double timediff(double t0, double t1);

/* error_check.c */
void error_check_header(void);

/* interpolate.c */
float interpolate(GRID *** grid, float x, float y, float z, int i, int j, int k);
float integrate_cell(int i, int j, int k, int mode);

//array2minc.c
void array2mnc(float ***array, /*char* format_filename,*/ int *size, char *name,
          float min, float max, int is_signed);

/*
  WARNING : Some of these functions are defined in both feedback_w_kernal.c
  and feedback.c. FEEDBACK_W_SPH_KERNAL determines which file controls the
 function's behaviour.
*/
void feedback(void);
void get_number_of_gas_nearby(void);
long get_number_gas_part_nearby(FLOAT radius, int target, int mode);
int ngb_treefind_gas_in_box(FLOAT searchcenter[3], FLOAT radius,int *startnode);

int gaspartsearch_compare_key(const void *a, const void *b);
int yielddata_compare_key(const void *a, const void *b);

void update_star_hsml(void);
void get_feedback_normalization(void);
void distribute_metals_and_energy_w_kernal(void);
long distribute_feedback_to_nearby_gas_w_kernal(FLOAT radius, int target,
int mode);
FLOAT get_feedback_normalization_constant(FLOAT h, double * gasdensity, 
double * gaspressure, int target, int mode);
int feedbacknorm_compare_key(const void *a, const void *b);

// simplify.c
void simplify_particle_data(void);

// write.c
void write_data(char *);
