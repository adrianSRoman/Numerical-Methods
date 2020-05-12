#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>
#include<sys/param.h>

#define BUF1 32
#define RK_COUNT 2

int main()
{
   /* Declare variables */
   int     n_dt, n;
   double  fn, t0, tn, t_plus, y0, yn, ym, ym2, ym3, kappa, eta, c_rk, loc_err, ac_err;
   double  t_start, t_stop, dt, dt_tar;
   double  k1_m, k2_m, k1_m2, k2_m2, k1_m3, k2_m3; // Runge Kutta k-values
   /* This is for opening a file -- leave it as is unless you know what you are doing */
   char    *rundat = "Dat-adapt";
   char    datfile[BUF1];
   FILE    *f_dat;

   /* Declaration of right-hand side of differential equation */
   extern double func();

   /* Open the data file Dat */
   sprintf(datfile, "%s", rundat);
   if ((f_dat = fopen(datfile, "w")) == NULL) {
      fprintf(stdout, "Can't open Dat file\n");
      fflush(stdout);
      exit(1);
   }

   /* Set system parameters */
   eta = 1.1;
   kappa = 1.0;
   double t_0 = 6.58211952e-12;
   double mic_amp = 10e-6;
   t_start = 0.0;                /* initial time */
   t_stop  = 100.0;//100.0;              /* final time */
   n_dt    = 1000;//100000;             /* number of time steps */
   dt = (t_stop-t_start)/n_dt;   /* time step */


   ac_err = (1/t_stop)*10e-2; //         /* acceptable local error */

   /* Initial condition */
   t0 = t_start; /* initial time */
   y0 = 0.4;     /* initial value */

   /* now dump the initial value and the function to the Dat file */
   fprintf(f_dat, "%8d %.15e %.15e\n", 0, t0, y0); fflush(f_dat);

   c_rk = 0.5;
   /* Prepare for time-step loop */
   yn = y0;
   fn = func(t0, y0, kappa, eta);
   tn = t0; /* initial time */
   while (tn <= t_stop)
   {  /* enter the discrete-time loop of the RK2 map */

      // 1. Calculate ym(n+1) = y(n) + ... 𝚫t_a
      k1_m = func(tn, yn, kappa, eta);
      k2_m = func(tn + c_rk*dt, yn + c_rk*dt*k1_m, kappa, eta);
      ym = yn + (dt/(2.0*c_rk))*((2.0*c_rk - 1.0)*k1_m + k2_m); /* RK2 step */

      // 2. Calculate ym2(n+2) = y(n + 1) + ... 2(ξ𝚫t_a)
      t_plus = tn + 2*dt; // 
      k1_m2 = func(t_plus, ym, kappa, eta);
      k2_m2 = func(t_plus + c_rk*dt, ym + c_rk*dt*k1_m2, kappa, eta);
      ym2 = ym + (dt/(2.0*c_rk))*((2.0*c_rk - 1.0)*k1_m2 + k2_m2);

      // 3. Calculate ym(n+2) = y(n + 2) + ... (ξ2𝚫t_a)
      k1_m3 = func(t_plus, yn, kappa, eta);
      k2_m3 = func(t_plus + c_rk*2*dt, yn + c_rk*2*dt*k1_m3, kappa, eta);
      ym3 = yn + (2*dt/(2.0*c_rk))*((2.0*c_rk - 1.0)*k1_m3 + k2_m3); /* RK2 step doubled */

      // 4. Calculate estimated local error
      loc_err = fabs((ym3 - ym2)/(pow(2.0, 3.0) - 2.0));
      // fprintf(stdout, "Local Error: %.15e\n", loc_err);

      // 5. Calculate 𝚫t - target
      dt_tar = pow(dt, 3/2.0) * pow(((ac_err)/(loc_err)), 0.5); /* find target 𝚫t */
      // fprintf(stdout, "actual dt: %.15e, target: %.15e\n", dt, dt_tar);

      if (dt < dt_tar) /* Accept ym2 */
      {
         tn = tn + 2*dt;
         yn = ym2 - 2*loc_err; /* prepare for next step */
      }/* else: reject ym2 */
      if (dt > 0.1)
      {
         dt = 0.1;
      }

      dt = 0.9*dt_tar; /* adjust time step */
      fprintf(f_dat, "%8d %.15e %.15e\n", n+1, tn+dt, yn); fflush(f_dat);
   }
   fclose(f_dat);
   /* Print final values to screen */
   fprintf(stdout, "%8d %.15e %.15e\n", n, tn+dt, ym);
   fflush(stdout);
}

/* the right-hand side function f */
double func(t,y,c, eta)
double t, y, c, eta;
{
   double f;

   f = -c*sin(y)+eta;

   return f;
}


