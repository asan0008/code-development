// This code computes AS forces which are called from MMTK_brownian.c (updated on 10/08/2017 by Aritra)
// Neighborlist for Mixed Flow is updated on 01/11/12 by A. Jain
#include "functionDef.h"
#include "sharedmacros.h"
#include <sys/time.h>
#include "stdio.h"
#include <math.h>

#define Nbdef 100000

#define vector_substract(v1, v2, v3) {                                  \
    int i_;                                                             \
    for (i_=0;i_<3;i_++) {                                              \
      (v1)[i_] = (v2)[i_] - (v3)[i_];                                   \
      }                                                                 \
}
/* 
void rotateCounterClockwiseXY(double *x1, double *y1, double *angle)
{
  double xrot1;
  double yrot1;
  xrot1 = *x1 * cos(*angle) - *y1 * sin(*angle);
  yrot1 = *x1 * sin(*angle) + *y1 * cos(*angle);

  *x1 = xrot1;
  *y1 = yrot1;
}

double  rotateClockwiseXY(double *x, double *y, double *angle)
{
  double xrot;
  double yrot;
  xrot = *x*cos(*angle) + *y*sin(*angle);
  yrot = -*x*sin(*angle) + *y*cos(*angle);

  *x = xrot;
  *y = yrot;
}
*/
////////////////////////////////////////////////////////////////////////////
// Subboxes Beginning
////////////////////////////////////////////////////////////////////////////

// Free subbox objects.
// A Common interface to take care of freeing all the mallocs spaces when a new object is created.
void
subbox_deallocAS(SubBoxListObjectAS *self)
{
  free(self->box_numberAS);
  free(self->boxesAS);
  free(self);
}

/*
 * Create a subbox object.
 *
 * Subboxes store the atoms in the simulation domain in sub_boxes.
 * They are used to compute the nearest neighbour list for short range
 * interactions.  The description of variables is given in the header
 * file where this is defined.
 */
SubBoxListObjectAS *
subbox_newAS(void)
{
  SubBoxListObjectAS *self;
  self = (SubBoxListObjectAS  *) malloc(sizeof(SubBoxListObjectAS));
  //freed in subbox_dealloc() from the function that calls subbox_init()

  self->universe_spec = NULL;
  self->box_numberAS = NULL;
  self->box_atomsAS = NULL;
  self->boxesAS = NULL;
  self->nboxesAS = 0;
  return self;
}

/*
 * Initialise the subboxes in a Universe.
 *
 * This Creates a pointer to SBLO (nblist), and initialises it.  If it
 * already exists, then it is de allocated and fresh nblist is created
 * with the specifications given in other arguments.  The code is
 * adapted from MMTK.nblist_update().  We use the similar notation as in
 * that code for convenience.
 *
 * The return value is also the same as the first input argument. This
 * is so that it can be assigned in the calling function (for which
 * the compiler will otherwise warn "may be used uninitialised").  It
 * is required as an input argument in cases where a current object
 * needs a reset.
 */
SubBoxListObjectAS *
subbox_initAS(SubBoxListObjectAS *nblistAS,
    int natoms, double *radii, int periodicity, int Nbpc, double cutoffAS1, double cutoffAS2, int subbox_per_rcAS,
    PyUniverseSpecObject *universe_spec)
{
  if (nblistAS != NULL) {
    subbox_deallocAS(nblistAS);
  }
  nblistAS = subbox_newAS();
  if (nblistAS == NULL) return NULL;
  nblistAS->cutoffAS1 = cutoffAS1;
  nblistAS->cutoffAS2 = cutoffAS2;
  nblistAS->universe_spec = universe_spec;

  double *geometry_data = universe_spec->geometry_data;

  vector3 box1AS, box2AS;
  double box_sizeAS[3];
  int i, ix, iy, iz, minx, miny, minz, maxx, maxy, maxz;

  double bead_rad = radii[0]; //*Non-dimensional bead radius*  This line needs to be modified for a polydisperse system. Choose bead_rad to be the the maximum of all the radii for polydisperse system.
  //nblistEV->box_countEV[0] = nblistEV->box_countEV[1] = nblistEV->box_countEV[2] = 0;
  if (nblistAS->box_numberAS == NULL) {
    nblistAS->box_numberAS = (int *)malloc(2*natoms*sizeof(int));
    //freed in subbox_dealloc()
    if (nblistAS->box_numberAS == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    /* Both box_number and box_atoms require natoms int each.  A
       total of two*natoms is allocated above to box_number.
       box_number uses the first natoms, and the second natoms is
       used by box_atoms.  Here we assign the begining of box_atoms
       pointer to a shifted reference from box_number pointer.
       */
    nblistAS->box_atomsAS = nblistAS->box_numberAS + natoms;
  }
  double L10x = geometry_data[0]; double L20x = geometry_data[1]; double L30x = geometry_data[2];
  double L10y = geometry_data[3]; double L20y = geometry_data[4]; double L30y = geometry_data[5];
  double L10z = geometry_data[6]; double L20z = geometry_data[7]; double L30z = geometry_data[8];
  double Len1, Len2, Len3;
  // Magnitude of the three sides
  Len1 = sqrt((L10x*L10x)+(L10y*L10y)+(L10z*L10z));
  Len2 = sqrt((L20x*L20x)+(L20y*L20y)+(L20z*L20z));
  Len3 = sqrt((L30x*L30x)+(L30y*L30y)+(L30z*L30z));

  // Calculating the number of subboxes
  nblistAS->box_countAS[0] = ceil((subbox_per_rcAS*(Len1)/nblistAS->cutoffAS2) +1);
  nblistAS->box_countAS[1] = ceil((subbox_per_rcAS*(Len2)/nblistAS->cutoffAS2) +1);
  nblistAS->box_countAS[2] = ceil((subbox_per_rcAS*(Len3)/nblistAS->cutoffAS2) +1);

  // Calculating the sizes of subboxes
  box_sizeAS[0] = (Len1)/nblistAS->box_countAS[0];
  box_sizeAS[1] = (Len2)/nblistAS->box_countAS[1];
  box_sizeAS[2] = (Len3)/nblistAS->box_countAS[2];
  
 /* 
  print_location();
  debug_vector("box sizeEV",box_sizeEV);
  printf("box countEV = %lu, %lu, %lu\n",
      nblistEV->box_countEV[0],
      nblistEV->box_countEV[1],
      nblistEV->box_countEV[2]
      );
*/
  // Making sure that the size of subbox must not be less than 2a
  while (box_sizeAS[0] <= (2.0*bead_rad)) {
        nblistAS->box_countAS[0] = nblistAS->box_countAS[0] - 1;
        box_sizeAS[0] = (Len1)/nblistAS->box_countAS[0];
 }

  while (box_sizeAS[1] <= (2.0*bead_rad)) {
        nblistAS->box_countAS[1] = nblistAS->box_countAS[1] - 1;
        box_sizeAS[1] = (Len2)/nblistAS->box_countAS[1];
 }

  while (box_sizeAS[2] <= (2.0*bead_rad)) {
        nblistAS->box_countAS[2] = nblistAS->box_countAS[2] - 1;
        box_sizeAS[2] = (Len3)/nblistAS->box_countAS[2];
 }
// Making sure that the number of subboxes are odd numbers (this condition cannot be violated as it will results in conflicts while calculating intercations)
if (nblistAS->box_countAS[0]%2 == 0) {
        nblistAS->box_countAS[0] = nblistAS->box_countAS[0] + 1;
        box_sizeAS[0] = (Len1)/nblistAS->box_countAS[0];
}

if (nblistAS->box_countAS[1]%2 == 0) {
        nblistAS->box_countAS[1] = nblistAS->box_countAS[1] + 1;
        box_sizeAS[1] = (Len2)/nblistAS->box_countAS[1];
}

if (nblistAS->box_countAS[2]%2 == 0) {
        nblistAS->box_countAS[2] = nblistAS->box_countAS[2] + 1;
        box_sizeAS[2] = (Len3)/nblistAS->box_countAS[2];
}
  // Total number of boxes
  nblistAS->nboxesAS =
    nblistAS->box_countAS[0]*nblistAS->box_countAS[1]*nblistAS->box_countAS[2];
/*
  print_location();
  debug_vector("box sizeEV",box_sizeEV);
  printf("box countEV = %lu, %lu, %lu\n",
      nblistEV->box_countEV[0],
      nblistEV->box_countEV[1],
      nblistEV->box_countEV[2]
      );
*/
  nblistAS->boxesAS = malloc((nblistAS->nboxesAS)*sizeof(nbbox));
  //freed in subbox_dealloc()
  if (nblistAS->boxesAS == NULL) {
    PyErr_NoMemory();
    return 0;
  }

//printf("bc[1] = %d \t bs[1] = %f \n", nblistAS->box_countAS[1], box_sizeAS[1]);

  int abs_minx, abs_miny, abs_minz;
  // Rotating the box
  double theta_t = atan(L10y/L10x);
  if (theta_t < 0.0) {theta_t = (M_PI) + theta_t;}
  rotateClockwiseXY(&L10x, &L10y, &theta_t);
  rotateClockwiseXY(&L20x, &L20y, &theta_t);
  double lx, ly, delx, dely_vert, delx_subbox;
  lx = L10x;
  ly = L20y;
  delx = L20x;
  dely_vert = ly/nblistAS->box_countAS[1];
  delx_subbox = delx/nblistAS->box_countAS[1];

  /* Each box stores its own coordinates in ix,iy,iz. Boxes are
     numbered linearly starting in the x-direction, followed by y,
     z.
     */
  i = 0;
  for (iz = 0; iz < nblistAS->box_countAS[2]; iz++)
    for (iy = 0; iy < nblistAS->box_countAS[1]; iy++)
      for (ix = 0; ix < nblistAS->box_countAS[0]; ix++) {
        nblistAS->boxesAS[i].ix = ix;
        nblistAS->boxesAS[i].iy = iy;
        nblistAS->boxesAS[i].iz = iz;
        i++;
      }

  /* neighbors, is a linear list of relative integer coordinates of
   * a neighbouring box from any box, that lie within the cutoff
   * distance. The first "neighbors" is the box itself with relative
   * coordinate 0,0,0. A box is not considered a neighbour only if
   * the minimum distance separating them is greater then the cutoff.
   */
  nblistAS->neighborsAS[0][0] = 0;
  nblistAS->neighborsAS[0][1] = 0;
  nblistAS->neighborsAS[0][2] = 0;
  i = 1;
//printf("nblistEV->cutoffEV = %f \t dely_vert = %f \t fabs(delx_subbox) = %f \t box_sizeEV[0] = %f \n", nblistEV->cutoffEV, dely_vert, fabs(delx_subbox), box_sizeEV[0]);
  //Determination of the limits of the neighbor-list
  minx = -((nblistAS->cutoffAS2+((1.0 + (nblistAS->cutoffAS2/dely_vert))*fabs(delx_subbox)))/box_sizeAS[0]);
  miny = -(nblistAS->cutoffAS2/dely_vert) - 1;
  minz = -(nblistAS->cutoffAS2/box_sizeAS[2]) - 1;

  maxx = ((nblistAS->cutoffAS2+((1.0 + (nblistAS->cutoffAS2/dely_vert))*fabs(delx_subbox)))/box_sizeAS[0]);
  maxy = (nblistAS->cutoffAS2/dely_vert) + 1;
  maxz = (nblistAS->cutoffAS2/box_sizeAS[2]) + 1;
  // Rotating back...
  rotateCounterClockwiseXY(&L10x, &L10y, &theta_t);
  rotateCounterClockwiseXY(&L20x, &L20y, &theta_t);
//printf("minx = %d \t maxx = %d \t miny = %d \t maxy = %d \t minz = %d \t maxz = %d \n", minx, maxx, miny, maxy, minz, maxz);
  // Making sure that the limits determined for the neighbor-list are within rc
  if ((maxx*box_sizeAS[0]) < nblistAS->cutoffAS2)
  {maxx++;}
  if ((maxy*box_sizeAS[1]) < nblistAS->cutoffAS2)
  {maxy++;}

  if ((maxz*box_sizeAS[2]) < nblistAS->cutoffAS2)
  {maxz++;}

  if ((minx*box_sizeAS[0]) > -nblistAS->cutoffAS2)
  {minx--;}

  if ((miny*box_sizeAS[1]) > -nblistAS->cutoffAS2)
  {miny--;}

  if ((minz*box_sizeAS[2]) > -nblistAS->cutoffAS2)
  {minz--;}

//printf("minx = %d \t maxx = %d \t miny = %d \t maxy = %d \t minz = %d \t maxz = %d \n", minx, maxx, miny, maxy, minz, maxz);

  // Making sure that the limits are equally distributed at below and at above 0,0,0
  if ( (minx+maxx) != 0 ){
        abs_minx = abs(minx);
        if (abs_minx > maxx) {maxx = abs_minx;}
        else { minx = -maxx;}
}
  if ( (miny+maxy) != 0 ){
        abs_miny = abs(miny);
        if (abs_miny > maxy) {maxy = abs_miny;}
        else { miny = -maxy;}
}
  if ( (minz+maxz) != 0 ){
        abs_minz = abs(minz);
        if (abs_minz > maxz) {maxz = abs_minz;}
        else { minz = -maxz;}
}

//printf("minx = %d \t maxx = %d \t miny = %d \t maxy = %d \t minz = %d \t maxz = %d \n", minx, maxx, miny, maxy, minz, maxz);

// When the fineness parameter is 0, we should have the minimum image convention
if (subbox_per_rcAS == 0){
        minx = miny = minz = -1;
        maxx = maxy = maxz = 1;
        }

  //printf("minx = %d \t maxx = %d \t miny = %d \t maxy = %d \t minz = %d \t maxz = %d \n", minx, maxx, miny, maxy, minz, maxz);
 // Generating the neighbor-list 
 for (ix = minx; ix <= maxx; ix++)  {
    for (iy = miny; iy <= maxy; iy++)  {
      for (iz = minz; iz <= maxz; iz++)  {
        if (!(ix == 0 && iy == 0 && iz == 0)) {
          //i = ((ix-1)*box_count*box_count) + ((iy-1)*box_count) + iz;
          //double dx = ((abs(ix)-1.)*(box_sizeEV[0])) - (f_rc*sign*(delx_subbox)*(abs(iy)));
         // double dy = (abs(iy)-1.)*dely_vert;
          //double dx = (abs(ix)-1.)*box_sizeEV[0];
          //double dy = (abs(iy)-1.)*box_sizeEV[1];
          //double dz = (abs(iz)-1.)*box_sizeEV[2];
          //if (dx < 0.) dx = 0.;
          //if (dy < 0.) dy = 0.;
          //if (dz < 0.) dz = 0.;
          //if (dx*dx+dy*dy+dz*dz <= sqr(nblistEV->cutoffEV)) { // This "if" condition will further refine the neighbor-list but for the time being it is postponed.The formulaes written above for dx, dy and dz need to corrected.
            nblistAS->neighborsAS[i][0] = ix;
            nblistAS->neighborsAS[i][1] = iy;
            nblistAS->neighborsAS[i][2] = iz;
            i++;
        // }
        }
      }
    }
  }

  nblistAS->nneighborsAS = i; // i is incremented by one till end of loop
  //debug
  //printf("nboxesEV = %lu, nneighboursEV = %d\n",nblistEV->nboxesEV, nblistEV->nneighborsEV);
  return nblistAS;
}

/* Update the subbox numbers for the atom positions.
 *
 * Each atom is assigned box numbers, and the boxes contain a list of
 * atoms. This function must be called everytime there is a change in
 * the position of any atom.  This function is adapted from
 * nblist_update of MMTK, with the initialisation separated from the
 * finding of the box numbers.
 */
int
subbox_updateAS(SubBoxListObjectAS *nblistAS, int natoms, vector3 *x )
{
  vector3 box1AS, box2AS;
  double box_sizeAS[3];
  int *p;
  int i, n, mu;
  // atom_subset variables, not used here (adapted from MMTK)
  int n_sub;
  long subset[1];
  double *geometry_data = nblistAS->universe_spec->geometry_data;


  if (nblistAS == NULL) {
    //PyErr_SetString(PyExc_TypeError, "subbox is NULL");
    Py_FatalError("subboxAS is NULL");
    return -1;
  }

double L1x = geometry_data[0]; double L2x = geometry_data[1]; double L3x = geometry_data[2];
double L1y = geometry_data[3]; double L2y = geometry_data[4]; double L3y = geometry_data[5];
double L1z = geometry_data[6]; double L2z = geometry_data[7]; double L3z = geometry_data[8];

double Len1, Len2, Len3;
Len1 = sqrt((L1x*L1x)+(L1y*L1y)+(L1z*L1z));
Len2 = sqrt((L2x*L2x)+(L2y*L2y)+(L2z*L2z));
Len3 = sqrt((L3x*L3x)+(L3y*L3y)+(L3z*L3z));
//printf("Before rotation: L1x = %f \t L1y = %f \t L2x = %f \t L2y = %f \n", L1x, L1y, L2x, L2y);
// Rotating the box (see http://www.sciencedirect.com/science/article/pii/S0010465598001787)
// This rotatioal part is irrelevant for PSF
double theta_t = atan(L1y/L1x);
if (theta_t < 0.0) {theta_t = (M_PI) + theta_t;}
rotateClockwiseXY(&L1x, &L1y, &theta_t);
rotateClockwiseXY(&L2x, &L2y, &theta_t);
//printf("angle = %f \n", theta_t);
//printf("After rotation: L1x = %f \t L1y = %f \t L2x = %f \t L2y = %f \n", L1x, L1y, L2x, L2y);
double lx, ly, delx, dely_vert, delx_subbox;
lx = L1x;
ly = L2y;
delx = L2x;
dely_vert = ly/nblistAS->box_countAS[1];
delx_subbox = delx/nblistAS->box_countAS[1];

// Setting the limits of the box
box1AS[2] = -Len3/2.0;
box2AS[2] = Len3/2.0;
box_sizeAS[0] = (Len1)/nblistAS->box_countAS[0];
box_sizeAS[1] = (Len2)/nblistAS->box_countAS[1];
box_sizeAS[2] = (Len3)/nblistAS->box_countAS[2];
//printf("bc0 = %d \t bc1 = %d \t bc2 = %d \n", nblistEV->box_countEV[0], nblistEV->box_countEV[1], nblistEV->box_countEV[2]);
  /* Initialise/Reset ... */
  for (i = 0; i < nblistAS->nboxesAS; i++) {
    nblistAS->boxesAS[i].n = 0; // # atoms in a box
    nblistAS->boxesAS[i].i = 0; // atom counter
  }

  n_sub = 0; // atom_subset not implemented
  subset[0] = 0;

  n = (n_sub == 0) ? natoms : n_sub;
  //print_location();
  for (mu = 0; mu < n; mu++) {
    int box = 0;
    int ai, nb;
    int ix3, iy3, iz3;
    double xai0;
    ai = mu;
    rotateClockwiseXY(&x[ai][0], &x[ai][1], &theta_t);
    xai0 = x[ai][0] - ((x[ai][1]/dely_vert)*delx_subbox);
    box = (int)((xai0+(lx/2.0))/box_sizeAS[0]);
    if (box == nblistAS->box_countAS[0]) box--;
    nb = (int)((x[ai][1]+(ly/2.0))/dely_vert);
    if (nb == nblistAS->box_countAS[1]) nb--;
    box += nblistAS->box_countAS[0]*nb;
    nb = (int)((x[ai][2]-box1AS[2])/box_sizeAS[2]);
    if (nb == nblistAS->box_countAS[2]) nb--;
    box += nblistAS->box_countAS[0]*nblistAS->box_countAS[1]*nb;
    nblistAS->box_numberAS[ai] = box; // box number for each atom
    nblistAS->boxesAS[box].n++;
    rotateCounterClockwiseXY(&x[ai][0], &x[ai][1], &theta_t);
  }
// Rotating back...
rotateCounterClockwiseXY(&L1x, &L1y, &theta_t);
rotateCounterClockwiseXY(&L2x, &L2y, &theta_t);
  /* To understand this, see the header file where the the variables
     are declared */
  p = nblistAS->box_atomsAS; // p points to box_atoms[0]
  for (i = 0; i < nblistAS->nboxesAS; i++) {
    nblistAS->boxesAS[i].atoms = p;
    nblistAS->boxesAS[i].i = 0;
    p += nblistAS->boxesAS[i].n; // the next p is after n atoms from the last
    /* This means that for boxes with atoms, boxes[i].atom[0] to
       boxes[i].atom[n-1], will "belong" to i-th box. A box with no
       atoms is also assigned the same pointer location for .atoms as
       end of atoms of the previous box. But this does not matter,
       because when the atoms in a box are found out (for caclculating
       interactions), only those boxes with n>0 are checked.  While assigning
       atoms to this box, a counter i (which is a member of nbbox) is used
       */
  }
  for (mu = 0; mu < n; mu++) {
    int ai = (n_sub == 0) ? mu : subset[mu];
    /* find the box for each atom, shifted from the origin by box_number*/
    nbbox *boxAS = nblistAS->boxesAS + nblistAS->box_numberAS[ai]; // allocate a pointer to the (boxesAS) + box number of the atom 
    /* add atom to the box*/
    boxAS->atoms[boxAS->i++] = ai;  // increment atom counter inside the box
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////
// subboxes end
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
#define compute_AS_force1() {                                             \
  double R8, R14, lj2_2, lj2_6, C;                                        \
  vector_copy(rij1, rij);                                                 \
  lj2_2 = Lja2*Lja2;                                                      \
  lj2_6 = lj2_2*lj2_2*lj2_2;                                              \
  R8 = R2*R2*R2*R2;                                                       \
  R14 = R8*R2*R2*R2;                                                      \
  C = 48.0*((sqr(lj2_6)/R14) - (0.5*lj2_6/R8))/Lja1;                      \
  vector_scale(rij, C);                                                   \
  vector_copy(Fij, rij);                                                  \
  vector_add(FAS[mu], rij, -1.0);                                         \
  vector_add(FAS[nu], rij, 1.0);                                          \
  ASTensor[0][0] = ASTensor[0][0] + (rij1[0]*Fij[0]);                     \
  ASTensor[1][1] = ASTensor[1][1] + (rij1[1]*Fij[1]);                     \
  ASTensor[0][1] = ASTensor[0][1] + (rij1[0]*Fij[1]);                     \
} 

#define compute_AS_force2() {                                             \
  double Rmin8, Rmin14, lj2_2, lj2_6, C;                                  \
  vector_copy(rij1, rij);                                                 \
  lj2_2 = Lja2*Lja2;                                                      \
  lj2_6 = lj2_2*lj2_2*lj2_2;                                              \
  Rmin8 = Rmin2*Rmin2*Rmin2*Rmin2;                                        \
  Rmin14 = Rmin8*Rmin2*Rmin2*Rmin2;                                       \
  C = 48.0*((sqr(lj2_6)/Rmin14) - (0.5*lj2_6/Rmin8))/Lja1;                \
  vector_scale(rij, C);                                                   \
  vector_copy(Fij, rij);                                                  \
  vector_add(FAS[mu], rij, -1.0);                                         \
  vector_add(FAS[nu], rij, 1.0);                                          \
  ASTensor[0][0] = ASTensor[0][0] + (rij1[0]*Fij[0]);                     \
  ASTensor[1][1] = ASTensor[1][1] + (rij1[1]*Fij[1]);                     \
  ASTensor[0][1] = ASTensor[0][1] + (rij1[0]*Fij[1]);                     \
}


// Attractive component of the LJ interaction for association
#define compute_AS_force3() {                                             \
  double C, lj2_2;                                                        \
  vector_copy(rij1, rij);                                                 \
  lj2_2 = Lja2*Lja2; 						          \
  C = (phi/lj2_2)*alphaAS*sin(alphaAS*(R2/lj2_2)+betaAS);                 \
  vector_scale(rij, C);                                                   \
  vector_copy(Fij, rij);                                                  \
  vector_add(FAS[mu], rij, -1.0);                                         \
  vector_add(FAS[nu], rij, 1.0);                                          \
  ASTensor[0][0] = ASTensor[0][0] + (rij1[0]*Fij[0]);                     \
  ASTensor[1][1] = ASTensor[1][1] + (rij1[1]*Fij[1]);                     \
  ASTensor[0][1] = ASTensor[0][1] + (rij1[0]*Fij[1]);                     \
}

// Attractive component of the LJ interaction for association to define the quality of solvent
#define compute_AS_force4() {                                             \
  double C, lj2_2, phi_solvent;                                           \
  vector_copy(rij1, rij);                                                 \
  lj2_2 = Lja2*Lja2;                                                      \
  phi_solvent = 0.45;                                                     \
  C = (phi_solvent/lj2_2)*alphaAS*sin(alphaAS*(R2/lj2_2)+betaAS);                 \
  vector_scale(rij, C);                                                   \
  vector_copy(Fij, rij);                                                  \
  vector_add(FAS[mu], rij, -1.0);                                         \
  vector_add(FAS[nu], rij, 1.0);                                          \
  ASTensor[0][0] = ASTensor[0][0] + (rij1[0]*Fij[0]);                     \
  ASTensor[1][1] = ASTensor[1][1] + (rij1[1]*Fij[1]);                     \
  ASTensor[0][1] = ASTensor[0][1] + (rij1[0]*Fij[1]);                     \
}

ASConstsObject * 
as_force(SubBoxListObjectAS *nblistAS,
    int Nb,             // Number of beads
    double Lja1,         // LJ paramter   (kBT/epsilon)
    double Lja2,         // LJ paramter   (sigma/lk) 
    int AStype,         // 1 for LJ, 0 for NO ASSOCIATION
    int periodicity,     // position of the associating beads on the chain
    int Nbpc,           // number of beads per chain
    double phi,         // LJ well depth
    double alphaAS,         // LJ parameter, Soddermann, Duenweg, 2001
    double betaAS,    // LJ parameter, Soddermann, Duenweg, 2001
    double cutoffAS1,  // cutoff for repulsive force
    double cutoffAS2,  // cutoff for attractive force
    vector3 *x, vector3 *FAS, vector3 *ASTensor,
    PyUniverseSpecObject *universe_spec)
{
  int mu,nu;
  int Nc_mu, Nc_nu;
  int ibox1,ibox2;
  double *geometry_data = universe_spec->geometry_data;
  double Rmin, Rmin2; 
  double R2;
  vector3 rij, Fij, rij1;
  ASTensor[0][0] = 0.0; ASTensor[0][1] = 0.0; ASTensor[1][1] = 0.0; // Tau[0] = txx, Tau[1] = tyy, Tau[2] = txy
  int first_atom, last_atom;
  first_atom = 0;
  last_atom = Nb;
  Rmin = 0.7*Lja2;  // Min cap for the repulsive part of the interaction
  Rmin2 = Rmin*Rmin;   
 // printf("Rmin = %lf \t Rmin2 = %lf\t cutoffAS1 = %lf \n", Rmin, sqr(Rmin), cutoffAS1);
  double xlimit, ylimit, zlimit;
  int bx, by, bz;
  first_atom = 0;
  last_atom = Nb;
  int pair[Nb]; // defining the functionality of the stickers
  int ix, iy, iz;
  int ix2, iy2, iz2;
  double arg1x1, arg1x2, arg1y1, arg1y2, arg1z1, arg1z2;
  int arg2x1, arg2x2, arg2y1, arg2y2, arg2z1, arg2z2;
  double boxvecx, boxvecy, boxvecz;
  double box_sizeAS[3];
  double L1x = geometry_data[0]; double L2x = geometry_data[1]; double L3x = geometry_data[2];
  double L1y = geometry_data[3]; double L2y = geometry_data[4]; double L3y = geometry_data[5];
  double L1z = geometry_data[6]; double L2z = geometry_data[7]; double L3z = geometry_data[8];
  double Len1, Len2, Len3;
  Len1 = sqrt((L1x*L1x)+(L1y*L1y)+(L1z*L1z));
  Len2 = sqrt((L2x*L2x)+(L2y*L2y)+(L2z*L2z));
  Len3 = sqrt((L3x*L3x)+(L3y*L3y)+(L3z*L3z));

  box_sizeAS[0] = (Len1)/nblistAS->box_countAS[0];
  box_sizeAS[1] = (Len2)/nblistAS->box_countAS[1];
  box_sizeAS[2] = (Len3)/nblistAS->box_countAS[2];

  for (mu = first_atom; mu < last_atom; mu++){
      pair[mu] = 0;  // 0 indicates no pairing        
  }
  // Rotating the box
  double theta_t = atan(L1y/L1x);
  //printf("theta1 = %lf \n", 180*theta_t/M_PI);
  if (theta_t < 0.0) {theta_t = (M_PI) + theta_t;}
  //printf("theta2 = %lf \n", 180*theta_t/M_PI);
  rotateClockwiseXY(&L1x, &L1y, &theta_t);
  rotateClockwiseXY(&L2x, &L2y, &theta_t);
  double lx, ly, delx, dely_vert, delx_subbox;
  lx = L1x;
  ly = L2y;
  delx = L2x;
  dely_vert = ly/nblistAS->box_countAS[1];
  delx_subbox = delx/nblistAS->box_countAS[1];

  // Initialization of AS Force
  for (mu = first_atom; mu < last_atom; mu++) {
    FAS[mu][0] = 0.0;
    FAS[mu][1] = 0.0;
    FAS[mu][2] = 0.0;
  }

// Rotating all the beads
for (mu = first_atom; mu < last_atom; mu++) {
        rotateClockwiseXY(&x[mu][0], &x[mu][1], &theta_t);
}
clock_t total1;
clock_t total2;
double t_elapsed;

//printf("Nneighbors = %d \n", nblistEV->nneighborsEV);
  // we start at the second atom as we only sum over mu > nu  /********/ mu starts from "first_atom" instead of "first_atom + 1" as nu>mu  
  for (mu = first_atom; mu < last_atom; mu++) {
    ibox1 = nblistAS->box_numberAS[mu];
    nbbox *box1 = &nblistAS->boxesAS[ibox1];
    bx = box1->ix - ((nblistAS->box_countAS[0] - 1)/2);
    by = box1->iy - ((nblistAS->box_countAS[1] - 1)/2);
    bz = box1->iz - ((nblistAS->box_countAS[2] - 1)/2);
    int ineighbor;
    for (ineighbor = 0; ineighbor < nblistAS->nneighborsAS; ineighbor++) {
      //box1->ix, box1->iy and box1->iz are according to 0:box_count numbering system (not minx:maxx)
      //bx, by and bz are according to minx:maxx system
      ix = nblistAS->neighborsAS[ineighbor][0] + bx;
      iy = nblistAS->neighborsAS[ineighbor][1] + by;
      iz = nblistAS->neighborsAS[ineighbor][2] + bz;
      // ix1, iy1, iz1 are original locations according to (0,0) being the center of the box
      int ix1 = ix;
      int iy1 = iy;
      int iz1 = iz;
      nbbox *box2;
      int j;
      xlimit = (ix*box_sizeAS[0]);
      ylimit = (iy*box_sizeAS[1]);
      zlimit = (iz*box_sizeAS[2]);
      //ix, iy and iz are now the locations in the box (folded into the box) according to (0,0) being the center of the box
      // Bringing subboxes back into the original box if they are outside the original box
if (universe_spec->is_periodic) {
arg1x1 = (-(Len1-box_sizeAS[0])/2.0);   // subtracting the size of the box to measure the distance by centering a box at the origin
arg2x1 = (-(nblistAS->box_countAS[0]-1)/2);
if ((xlimit < arg1x1) && (ix!=arg2x1)) {ix = ix + nblistAS->box_countAS[0];}

arg1y1 = (-(Len2-box_sizeAS[1])/2.0);
arg2y1 = (-(nblistAS->box_countAS[1]-1)/2);
if ((ylimit < arg1y1) && (iy!=arg2y1)) {iy = iy + nblistAS->box_countAS[1];}

arg1z1 = (-(Len3-box_sizeAS[2])/2.0);
arg2z1 = (-(nblistAS->box_countAS[2]-1)/2);
if ((zlimit < arg1z1) && (iz!=arg2z1)) {iz = iz + nblistAS->box_countAS[2];}

arg1x2 = ((Len1-box_sizeAS[0])/2.0);
arg2x2 = ((nblistAS->box_countAS[0]-1)/2);
if ((xlimit > arg1x2) && (ix!=arg2x2)) {ix = ix - nblistAS->box_countAS[0];}

arg1y2 = ((Len2-box_sizeAS[1])/2.0);
arg2y2 = ((nblistAS->box_countAS[1]-1)/2);
if ((ylimit > arg1y2) && (iy!=arg2y2)) {iy = iy - nblistAS->box_countAS[1];}

arg1z2 = ((Len3-box_sizeAS[2])/2.0);
arg2z2 = ((nblistAS->box_countAS[2]-1)/2);
if ((zlimit > arg1z2) && (iz!=arg2z2)) {iz = iz - nblistAS->box_countAS[2];}
      }

      //ix2, iy2, iz2 are locations accorrding to (0,0) being at corner and not at the center of the box
      ix2 = ix + ((nblistAS->box_countAS[0]-1)/2);
      iy2 = iy + ((nblistAS->box_countAS[1]-1)/2);
      iz2 = iz + ((nblistAS->box_countAS[2]-1)/2);

      ibox2 = ix2 + nblistAS->box_countAS[0]*(iy2 + nblistAS->box_countAS[1]*iz2);
      box2 = &nblistAS->boxesAS[ibox2];
      //printf("box2->n = %d \n", box2->n);
     /* double cutoffsq;
      if (EVtype == 3) // for WCA: hard coded cutoff
        cutoffsq = sqr(lj2) * pow(2., 1./3.);*/
      
      for (j = 0; j < box2->n; j++) { // only for non-empty boxes
        nu = box2->atoms[j];
        if (nu>mu) {
          boxvecx = ((ix1 - ix)*box_sizeAS[0]) + (delx_subbox*(iy1 - iy));
          boxvecy = (iy1 - iy)*dely_vert;
          boxvecz = (iz1 - iz)*box_sizeAS[2];
          vector3 boxvec = {boxvecx, boxvecy, boxvecz};
          vector_substract(rij, x[nu], x[mu]);
          vector_add(rij, boxvec, 1.0);  // this is for minimum image convention
          R2 = vector_length_sq(rij);
          rotateCounterClockwiseXY(&rij[0], &rij[1], &theta_t);
          Nc_mu = (int)(mu/Nbpc);  // chain number for bead no. mu
          Nc_nu = (int)(nu/Nbpc);  // chain number for bead no. nu
          //printf("mu = %d\t Nc_mu = %d\n", mu, Nc_mu);
          //printf("nu = %d\t Nc_nu = %d\t R2 = %lf\n", nu, Nc_nu, R2);

          if (AStype == 1) {  // Lennard-Jones capped at Rmin
           if ((((mu-4)+Nc_mu) % periodicity) == 0 && (((nu-4)+Nc_nu) % periodicity) == 0){
             if (R2 <= sqr(nblistAS->cutoffAS2)) {
                 if (R2 < sqr(Rmin)){ 
                     //printf ("R2 = %lf\n", R2);
                     compute_AS_force2();
                  } 
                 else {
                     if (R2 >= sqr(Rmin) && R2 <= sqr(nblistAS->cutoffAS1)){ 
                            compute_AS_force1();
                         }
                    /*   
                    else if (R2 > sqr(nblistAS->cutoffAS1) && ((mu - Nc_mu) % periodicity) == 0 && ((nu - Nc_nu) % periodicity) == 0){
                           //printf("Association!!\n");
                           //printf("mu = %d\t Nc_mu = %d\n", mu, Nc_mu);
                           //printf("nu = %d\t Nc_nu = %d\n", nu, Nc_nu);
                           compute_AS_force3(); 
                     }
                   */
                   /*Multi-sticker system*/
                     else if (R2 > sqr(nblistAS->cutoffAS1)){ 
                         //if ((((mu-4)+Nc_mu) % periodicity) == 0 && (((nu-4)+Nc_nu) % periodicity) == 0)
                           /*printf("Association!!\n");
                           printf("mu = %d\t Nc_mu = %d\n", mu, Nc_mu);
                           printf("nu = %d\t Nc_nu = %d\n", nu, Nc_nu);
                           printf("pair[%d] = %d\n", mu, pair[mu]);
                           printf("pair[%d] = %d\n", nu, pair[nu]);*/
                             if (pair[mu] < 1 && pair[nu] < 1){
                                 compute_AS_force3(); 
                                 pair[mu] = pair[mu] + 1;
                                 pair[nu] = pair[nu] + 1;
                                 /*printf("pair[%d] = %d\n", mu, pair[mu]);
                                 printf("pair[%d] = %d\n", nu, pair[nu]);*/
                             }
                        /*
                        else {
                            compute_AS_force4(); 
                         }
                        */      
                         }
                      }   
                  } // R2 <= sqr(cutoffAS2)
               }
             } // AStype == 1
          } // if(nu>=mu)
       } // nu loop closes here
    } // ineighbor loop closes here
  } // mu loop closes here
// Rotating back...
for (mu = first_atom; mu < last_atom; mu++) {
        rotateCounterClockwiseXY(&x[mu][0], &x[mu][1], &theta_t);
}
//printf("txyEV = %f \n", Tau[0][2]);
  return 0;
}
