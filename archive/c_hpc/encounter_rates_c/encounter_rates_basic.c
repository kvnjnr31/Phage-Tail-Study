#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

const int P = 10000000;//1000000; //100000000;
const int B = 1;
const int nsteps = 400;//500; //300000;
double PIref = 3.141592653589793238;

void boundarys(double x_p_new, double y_p_new, double z_p_new, double x_p_old, double y_p_old, double z_p_old,
               double x_left, double x_right, double y_left, double y_right, double z_upper, double z_lower, double R,
               double* x_new1, double* y_new1, double* z_new1, double* x_old1, double* y_old1, double* z_old1);

int main(int argc, char *argv[]){

//Parameters----------------------------
double x_old, x_new, y_old, y_new, z_old, z_new, wt, wx, wy, wz, k1, k2, k3, Ux, Uy, Uz, U1, D0, nu, tmin, Tr, theta, alpha, tau, vx, vy, vz, phage;
double bac, dist, omega, Nenc, Nact, v_comp_x, v_comp_y, v_comp_z, x_comp_x, x_comp_y, x_comp_z, r, dot, theta_x, theta_z, x1, y1, z1, x2, y2, z2;
double data_count, veff_x, veff_y, veff_z, A1, B1, C1, tc, rbound, indx;
double v=30, Rb=1, Rp=1, dt=0.01, Time=0, LX = 100, M = 1;
double x_left=-LX/2, x_right=LX/2, y_left=-LX/2, y_right=LX/2, z_lower=-LX/2, z_upper=LX/2, V = LX*LX*LX, Rcrit=Rb+Rp;

int i, p, b, n, m, k, Mf, proc, tag=314, master=0;                                                                                                                      //No. of bacteria and phage
double MatP[P][8], MatB[B][10], my_kernel[nsteps], kernel_old[nsteps], kernel_new[nsteps], en_indx[P];    //array of particle positions
MPI_Status status;

// Initialize the MPI environment
MPI_Init(&argc, &argv);

// Get the number of processes
int comm_sz;
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

// Get the rank of the process
int my_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

// Seed random numbers to rank
srand(time(NULL) + my_rank);

//Arguments
omega = strtol(argv[2], NULL, 10);   //Number of Omega's
D0 = 2.54;                           //Diffusion Constant
nu = 0.82;                            //Diffusion Exponent
tmin = 0.001;                        //power-law gamma parameter

Mf = M/comm_sz;                      //each worker does a chunk of iterations

printf("omega is: %f, D0 is: %f,  nu is : %f \n", omega, D0, nu);

//open save locations
FILE *pha_x = fopen("pha_x.dat", "wb");
FILE *pha_y = fopen("pha_y.dat", "wb");
FILE *pha_z = fopen("pha_z.dat", "wb");

FILE *bac_x = fopen("bac_x.dat", "wb");
FILE *bac_y = fopen("bac_y.dat", "wb");
FILE *bac_z = fopen("bac_z.dat", "wb");

FILE *accum = fopen("accum.dat", "wb");

//***************************************************************************
//                            Main Loop
//***************************************************************************

for(m=0;m<Mf;m++){ //iterations of experiment

//initialize old kernel
for(n=0;n<nsteps;n++){
    kernel_old[n] = 0;
}
//initialize phage activation index.
for(p=0;p<P;p++){
    en_indx[p]=1;
}

Time=0;

//***************************************************************************
//                            Initialize System
//***************************************************************************
//initialize bacteria
for(b=0;b<B;b++){

//random locations
double x_old=(x_right-x_left-2*Rb)*((double)rand()/(RAND_MAX)) + x_left + Rb;
double y_old=(y_right-y_left-2*Rb)*((double)rand()/(RAND_MAX)) + y_left + Rb;
double z_old=(z_upper-z_lower-2*Rb)*((double)rand()/(RAND_MAX)) + z_lower + Rb;

//store the bacteria tracks
MatB[b][1] = 0*x_old; MatB[b][2] = 0*y_old; MatB[b][3] = 0*z_old;  //pos-new
MatB[b][4] = 0*x_old; MatB[b][5] = 0*y_old; MatB[b][6] = 0*z_old;  //pos-old

//store the bacteria velocities
MatB[b][7] = 0; MatB[b][8] = 0; MatB[b][9] = 0;  //vel-new

//transition timer
MatB[b][10]=0;
}

//initialize phage
for(p=0;p<P;p++){

r = 0;

while (r<(Rb+Rp)){
//random locations
double x_old=(x_right-x_left-2*Rp)*((double)rand()/(RAND_MAX)) + x_left + Rp;
double y_old=(y_right-y_left-2*Rp)*((double)rand()/(RAND_MAX)) + y_left + Rp;
double z_old=(z_upper-z_lower-2*Rp)*((double)rand()/(RAND_MAX)) + z_lower + Rp;

//store the phage tracks
MatP[p][1] = x_old; MatP[p][2] = y_old; MatP[p][3] = z_old;  //pos-new
MatP[p][4] = x_old; MatP[p][5] = y_old; MatP[p][6] = z_old;  //pos-old

r = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));
}

//transition timer
MatP[p][7]=0;

//contact timer
MatP[p][8]=0;

}

//***************************************************************************
//                     Simulate system over nsteps
//***************************************************************************
data_count = 1;                      //number of data points
n = 0;
indx=0;

while((n<nsteps)&&(data_count<10000000)){

//Phage Dynamics
for(p=0;p<P;p++){


    //check for overlaps
    r = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));  //radius asociated with each point.

    //Compute hard core interactions
    /*
    if (r <(Rb+Rp)){

    //Compute effective phage velocities
    veff_x = (MatP[p][1] - MatP[p][4])/dt;
    veff_y = (MatP[p][2] - MatP[p][5])/dt;
    veff_z = (MatP[p][3] - MatP[p][6])/dt;

    A1 = pow((MatB[0][7] - veff_x),2) + pow((MatB[0][8] - veff_y),2) + pow((MatB[0][9] - veff_z),2);
    C1 = (pow((MatB[0][4] - MatP[p][4]),2) + pow((MatB[0][5] - MatP[p][5]),2) + pow((MatB[0][6] - MatP[p][6]),2)-(Rb+Rp)*(Rb+Rp))/A1;
    B1 = (2*(MatB[0][7] - veff_x)*(MatB[0][4] - MatP[p][4]) + 2*(MatB[0][8] - veff_y)*(MatB[0][5] - MatP[p][5]) + 2*(MatB[0][9] - veff_z)*(MatB[0][6] - MatP[p][6]))/A1;

    tc = (-B1-sqrt(B1*B1-4*C1))/2; //time of collision

    //update collision points
    MatP[p][4] = MatP[p][4] + veff_x*tc;
    MatP[p][5] = MatP[p][5] + veff_y*tc;
    MatP[p][6] = MatP[p][6] + veff_z*tc;

    //update new points via impusle
    MatP[p][1] = MatP[p][4] + (dt-tc)*2*(MatB[0][7] - veff_x)*((MatB[0][4]+MatB[b][7]*tc)-MatP[p][4])/r;
	
    MatP[p][2] = MatP[p][5] + (dt-tc)*2*(MatB[0][8] - veff_y)*((MatB[0][5]+MatB[b][8]*tc)-MatP[p][5])/r;
    MatP[p][3] = MatP[p][6] + (dt-tc)*2*(MatB[0][9] - veff_z)*((MatB[0][6]+MatB[b][9]*tc)-MatP[p][6])/r;

    r = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));

    }//end potential

/*
r = sqrt(pow((MatP[p][1]-MatB[0][1]),2) + pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));  //radius associated with each point.

//--------------------Location density of phage-- projections relative to bacteria--------------------------

if((r<5)&&(n>0)){

theta_z = -atan(MatB[0][9]/MatB[0][8]);

//rotation in x-coord
x1 = (MatP[p][1]-MatB[0][1]);
y1 = (MatP[p][2]-MatB[0][2])*cos(theta_z) - (MatP[p][3]-MatB[0][3])*sin(theta_z);
z1 = (MatP[p][2]-MatB[0][2])*sin(theta_z) + (MatP[p][3]-MatB[0][3])*cos(theta_z);

theta_x = acos(MatB[0][7]/v);  //angle of v w.r.t x-axis.

//rotation in z-coord
x2 = x1*cos(theta_x)-y1*sin(theta_x);
y2 = x1*sin(theta_x)+y1*cos(theta_x);
z2 = z1;

//save snapshots
fprintf(pha_x, "%d %f \n", n, x2);
fprintf(pha_y, "%d %f \n", n, y2);
fprintf(pha_z, "%d %f \n", n, z2);

//printf("The snapshot is %d data is : %f \n", n, data_count);

data_count = data_count + 1;  //index saved data
}
*/
//----------------------------collect "time in contact" of surface particles------------------------------
/*
rbound = sqrt(pow(MatB[0][1],2) + pow(MatB[0][2],2) + pow(MatB[0][3],2)); //distance of the bacteria from origin.

if ((r<=1.2)&&(rbound < LX/2-1.2)){
MatP[p][8] = MatP[p][8] + dt;  //time in contact
}

if ((r>1.2)&&(rbound < LX/2-1.2)){

if (MatP[p][8]>0){
indx = indx + 1;
fprintf(accum, "%f %f\n", indx, MatP[p][8]);
}

MatP[p][8] = 0;

}
*/
    wx = 0;
    wy = 0;
    wz = 0;

    k1 = 1;
    k2 = 1;
    k3 = 1;

    if (Time >= MatP[p][7]){

    //generate random noise
        Ux = ((double)rand()/(RAND_MAX))-0.5;
        wx = sqrt(24*D0*dt)*Ux;
        Uy = ((double)rand()/(RAND_MAX))-0.5;
        wy = sqrt(24*D0*dt)*Uy;
        Uz = ((double)rand()/(RAND_MAX))-0.5;
        wz = sqrt(24*D0*dt)*Uz;

        k1 = 1;
        k2 = 1;
        k3 = 1;

    //generate waiting time
    U1 = (double)rand()/(RAND_MAX);
    tau = tmin*pow(1/U1,(1/nu));
    MatP[p][7] = Time + tau; //store waiting time

    }//end update

    //compute velocity field!!
    r = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));

    dot = MatB[0][7]*(MatP[p][1]-MatB[0][1]) + MatB[0][8]*(MatP[p][2]-MatB[0][2]) +MatB[0][9]*(MatP[p][3]-MatB[0][3]);  //dot product at each pt.(scalar).

    v_comp_x = (3*Rb/(4*r))*(1+(Rb*Rb)/(3*r*r))*MatB[0][7];
    v_comp_y = (3*Rb/(4*r))*(1+(Rb*Rb)/(3*r*r))*MatB[0][8];
    v_comp_z = (3*Rb/(4*r))*(1+(Rb*Rb)/(3*r*r))*MatB[0][9];

    x_comp_x = (3*Rb/(4*r))*(1/(r*r) - (Rb*Rb)/(pow(r,4)))*(MatP[p][1]-MatB[0][1])*dot;
    x_comp_y = (3*Rb/(4*r))*(1/(r*r) - (Rb*Rb)/(pow(r,4)))*(MatP[p][2]-MatB[0][2])*dot;
    x_comp_z = (3*Rb/(4*r))*(1/(r*r) - (Rb*Rb)/(pow(r,4)))*(MatP[p][3]-MatB[0][3])*dot;

    //printf("Greetings from processor %d my vel is : %f \n", my_rank, v_comp_x);

if ((kernel_new[n-1]>100)&&(en_indx[p]==0)){
//save snapshots
fprintf(pha_x, "%d %f \n", n, MatP[p][1]);
fprintf(pha_y, "%d %f \n", n, MatP[p][2]);
fprintf(pha_z, "%d %f \n", n, MatP[p][3]);
}
    //Langevin x-dir
        x_old = MatP[p][1];
        x_new = x_old + en_indx[p]*(v_comp_x+x_comp_x)*dt + en_indx[p]*wx;

    //Langevin y-dir
        y_old = MatP[p][2];
        y_new = y_old + en_indx[p]*(v_comp_y+x_comp_y)*dt + en_indx[p]*wy;

    //Langevin z-dir
        z_old = MatP[p][3];
        z_new = z_old + en_indx[p]*(v_comp_z+x_comp_z)*dt + en_indx[p]*wz;


     //apply B.C's
    boundarys(x_new, y_new, z_new, x_old, y_old, z_old, x_left, x_right, y_left, y_right, z_upper, z_lower, Rp, &x_new, &y_new, &z_new, &x_old, &y_old, &z_old);

    //update the phage tracks
    MatP[p][1] = x_new; MatP[p][2] = y_new; MatP[p][3] = z_new;  //pos-new
    MatP[p][4] = x_old; MatP[p][5] = y_old; MatP[p][6] = z_old;  //pos-old

}//end phage

//printf("The data count is %f the time is : %f \n", data_count, Time);

if (kernel_new[n-1]>100){
// open save locations
fprintf(bac_x, "%d %f\n", n, MatB[0][1]);
fprintf(bac_y, "%d %f\n", n, MatB[0][2]);
fprintf(bac_z, "%d %f\n", n, MatB[0][3]);
}

//Bacteria Dynamics
for(b=0;b<B;b++){

    if (Time>=MatB[b][10]){

    //generate random run angles
        U1 = (double)rand()/(RAND_MAX);
        theta = 2*PIref*U1;
        U1 = (double)rand()/(RAND_MAX);
        alpha = 2*PIref*U1;

    //generate run time
        U1 = (double)rand()/(RAND_MAX);
        tau = -log(U1)/omega;
        MatB[b][10] = Time + tau;

    //compute velocity's
        vx=v*cos(alpha)*cos(theta);
        vy=v*cos(alpha)*sin(theta);
        vz=v*sin(alpha);
    //store the bacteria velocities
        MatB[b][7] = vx; MatB[b][8] = vy; MatB[b][9] = vz;  //vel-new

    }//end update


//Langevin x-dir
    x_old = MatB[b][1];
    x_new = MatB[b][7]*dt + x_old;

//Langevin y-dir
    y_old = MatB[b][2];
    y_new = MatB[b][8]*dt + y_old;

//Langevin z-dir
    z_old = MatB[b][3];
    z_new = MatB[b][9]*dt + z_old;


//apply B.C's
    boundarys(x_new, y_new, z_new, x_old, y_old, z_old, x_left, x_right, y_left, y_right, z_upper, z_lower, Rb, &x_new, &y_new, &z_new, &x_old, &y_old, &z_old);

//update the bacteria tracks
MatB[b][1] = x_new; MatB[b][2] = y_new; MatB[b][3] = z_new;  //pos-new
MatB[b][4] = x_old; MatB[b][5] = y_old; MatB[b][6] = z_old;  //pos-old

}//end bacteria

//*********************************Encounter Function***********************************

Nenc = 0; //set number of encounters to zero.
Nact = 0; //set number of active particles to zero.

//count number of active particles.
for(p=0;p<P;p++){
    Nact = Nact + en_indx[p];
}

//loop through phage
for(p=0;p<P;p++){

    dist = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));

    if(dist<=Rcrit){

        MatP[p][8] = MatP[p][8] + dt;  //time in contact

        //printf("timer is: %f, \n", MatP[p][8]);

        //capture time
        //U1 = (double)rand()/(RAND_MAX);
        //tau = -log(U1)/(0.01);


        indx = indx + 1;  //index accumulated phage


        //if (0<MatP[p][8]){
            Nenc = Nenc + en_indx[p];   //conditional encounter index
            en_indx[p] = 0;   // deactivate
        //}

    }

   if((fabs(MatP[p][1]-MatB[0][1])>=x_right)||(fabs(MatP[p][2]-MatB[0][2])>=y_right)||(fabs(MatP[p][3]-MatB[0][3])>=z_upper)){
        en_indx[p] = 1;   //activate
        MatP[p][8] = 0;  //time in contact
    }
} //loop on phage

kernel_new[n] = (Nenc/(M*Nact))*(V/dt);  //update kernel
kernel_old[n] = kernel_new[n];

fprintf(accum, "%f %f\n", n, indx);

//**************************End Encounter Function*************************//

Time = Time + dt; //update time
n = n + 1; //update steps
data_count = 0; //reset data count!!

}// end time

//printf("The snapshot is %d data count is : %f \n", n, data_count);
}//end M

//close save locations
fclose(pha_x);
fclose(pha_y);
fclose(pha_z);

fclose(bac_x);
fclose(bac_y);
fclose(bac_z);

fclose(accum);

if(my_rank==0){ // I am the master process and will gather results from the workers//

double kernel[nsteps];

for(n=0;n<nsteps;n++){

    kernel[n] = kernel_new[n];  //receive kernel_new from workers and update the final kernel.
    kernel_old[n] = 0;
    }//time

// open save locations
FILE *en_mat = fopen("en_mat.dat", "wb");

    for (proc=1;proc<comm_sz;proc++){

    //receive worker encounters
    MPI_Recv(&kernel_new,1,MPI_DOUBLE,proc,tag,MPI_COMM_WORLD,&status);

    for(n=0;n<nsteps;n++){

    kernel[n] = kernel_new[n] + kernel_old[n];  //receive kernel_new from workers and update the final kernel.
    kernel_old[n] = kernel[n];  //update old kernel.

    }//time

    printf("Master has received encounters from : %d\n", proc);
    }//end receive

    // Write to file
    for(n=0;n<nsteps;n++){
    fprintf(en_mat, "%d %f\n", n, kernel[n]);
    }

//close save locations
fclose(en_mat);
}
else{ //workers send encounters to the master

double tot_en=0;
    MPI_Send(&kernel_new,1,MPI_DOUBLE,master,tag,MPI_COMM_WORLD);

        for(n=0;n<nsteps;n++){
        tot_en = tot_en + kernel_new[n];
        }
    printf("Greetings from processor %d my encounters is : %f \n", my_rank, tot_en);

}

MPI_Finalize(); //let MPI finish up
return 0;
} // end main

void boundarys(double x_p_new, double y_p_new, double z_p_new, double x_p_old, double y_p_old, double z_p_old, double x_left, double x_right, double y_left, double y_right, double z_upper, double z_lower, double R, double* x_new_p, double* y_new_p, double* z_new_p, double* x_old_p, double* y_old_p, double* z_old_p){
// Apply random boundary conditions to bacteria
// R = radius of particle

    if(abs(x_p_new)>=x_right){
    x_p_new = x_p_new - (2*x_right)*x_p_new/fabs(x_p_new);
    x_p_old = x_p_old - (2*x_right)*x_p_old/fabs(x_p_old);
    }
    else if(abs(y_p_new)>=y_right){
    y_p_new = y_p_new - (2*y_right)*y_p_new/fabs(y_p_new);
    y_p_old = y_p_old - (2*y_right)*y_p_old/fabs(y_p_old);
    }
    else if(abs(z_p_new)>=z_upper){
    z_p_new = z_p_new - (2*z_upper)*z_p_new/fabs(z_p_new);
    z_p_old = z_p_old - (2*z_upper)*z_p_old/fabs(z_p_old);
    }
    *x_new_p = x_p_new;
    *y_new_p = y_p_new;
    *z_new_p = z_p_new;

    *x_old_p = x_p_old;
    *y_old_p = y_p_old;
    *z_old_p = z_p_old;

};// boundary

