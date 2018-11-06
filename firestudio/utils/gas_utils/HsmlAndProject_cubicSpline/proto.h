
//#include <malloc/malloc.h>

void determine_hsml(void);
void make_map(void);
//int findHsmlAndProject(int argc, void *argv[]);
void peano_hilbert_order(void);

float ngb_treefind(float xyz[3], int desngb, float hguess);

size_t tree_treeallocate(int maxnodes, int maxpart);
void tree_treefree(void);

void endrun(int ierr);
void free_particle_data(void);


int get_next_file(int begsnap, int begfilenr, int endsnap, int endfilenr, 
		  int *snapshot, int *filenr);

int compare_grp_particles(const void *a, const void *b);
int compare_dens(const void *a, const void *b);
int compare_grp_particles(const void *a, const void *b);
int compare_energy(const void *a, const void *b);

void unbind(int lev, int head, int len);
double tree_treeevaluate_potential(int target);

void save_hsml(void);

void check(int i, double h);

void process_group(int gr);

int load_hash_table(void);
void load_group_catalogue(void);
void mark_hash_cells(void);
void load_particle_data(void);
void find_group_indices(void);
int id_sort_compare_id(const void *a, const void *b);

void determine_densities(void);


int tree_treebuild(void);

void tree_update_node_recursive(int no, int sib, int father);


void *mymalloc(size_t n);
void myfree(void *p);

void set_units(void);




int ngb_compare_key(const void *a, const void *b);
float ngb_treefind(float xyz[3], int desngb, float hguess);
int ngb_treefind_variable(float searchcenter[3], float hguess);


void read_parameter_file(char *);
