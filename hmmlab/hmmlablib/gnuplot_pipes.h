/*
    This file is part of HMMLab.

    HMMLab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HMMLab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>



#define GP_MAX_TMP_FILES        64
#define GP_TMP_NAME_SIZE        512

#define GP_CMD_SIZE     16384

typedef struct _GNUPLOT_CTRL_ {
    /* command file handling */
    FILE*     gnucmd ;

    /* Plotting options */
    int       nplots ;      /* Number of active plots at the moment */
    char      pstyle[32] ;  /* Current plotting style */

    /* temporary files opened */
    char      to_delete[GP_MAX_TMP_FILES][GP_TMP_NAME_SIZE] ;
    int       ntmp ;


} gnuplot_ctrl ;


#ifndef _ECLIPSE_TYPES_H_
typedef struct _DPOINT_ {
    double  x ;
    double  y ;
} dpoint ;
#endif

/*---------------------------------------------------------------------------
   Function :   check_X_display()
   In       :   void
   Out      :   int 0 if display Ok, 1 if not
   Job      :   checks out the DISPLAY environment variable to see if it
                is set or not.
   Notice   :
 ---------------------------------------------------------------------------*/

int check_X_display(void) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_init()
 * In       :   void
 * Out      :   gnuplot control handle: pointer to a gnuplot_ctrl struct
 * Job      :   fork the process to launch a gnuplot session, open up data
 *              and command files for input.
 * Notice   :
 *--------------------------------------------------------------------------*/

gnuplot_ctrl* gnuplot_init(void) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_close()
 * In       :   gnuplot control handle
 * Out      :   void
 * Job      :   close a gnuplot session previously opened by gnuplot_init()
 * Notice   :   kills the child PID.
 *--------------------------------------------------------------------------*/

void gnuplot_close(gnuplot_ctrl* handle) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_cmd()
 * In       :   gnuplot control handle, character string
 * Out      :   void
 * Job      :   sends a character string to gnuplot, as if typed by user in
 *              an interactive session.
 * Notice   :   no control is done on the ouput of the command.
 *--------------------------------------------------------------------------*/

void gnuplot_cmd(gnuplot_ctrl* handle, char* cmd, ...) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_setstyle()
 * In       :   plotting style: a character string taken from the
 *              following set:
 *              "lines"
 *              "points"
 *              "linespoints"
 *              "impulses"
 *              "dots"
 *              "steps"
 *              "errorbars"
 *              "boxes"
 *              "boxeserrorbars"
 * Out      :   void
 * Job      :   set a plotting style for a given gnuplot session
 * Notice   :
 *--------------------------------------------------------------------------*/

void gnuplot_setstyle(gnuplot_ctrl* h, char* plot_style) ;



/*----------------------------------------------------------------------------
 * Function :   gnuplot_set_xlabel()
 * In       :   handle to gnuplot_ctrl, x label
 * Out      :   void
 * Job      :   sets the x label to be printed out with this handle
 * Notice   :
 *--------------------------------------------------------------------------*/
void gnuplot_set_xlabel(gnuplot_ctrl* h, char* label) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_set_ylabel()
 * In       :   handle to gnuplot_ctrl, y label
 * Out      :   void
 * Job      :   sets the y label to be printed out with this handle
 * Notice   :
 *--------------------------------------------------------------------------*/
void gnuplot_set_ylabel(gnuplot_ctrl* h, char* label) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_resetplot()
 * In       :   handle to gnuplot_ctrl
 * Out      :   void
 * Job      :   reset a plot (next plot with erase previous ones)
 * Notice   :
 *--------------------------------------------------------------------------*/

void gnuplot_resetplot(gnuplot_ctrl* h) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_plot1d_var1()
 * In       :   handle to gnuplot_ctrl, list of doubles, number of
 *              doubles in this list, title for the output graph.
 * Out      :   void
 * Job      :   plots out a 2d graph from a list of doubles. The
 *              x-coordinate is the double index, the y-coordinate is
 *              the double itself.
 * Notice   :
 * Example  :

    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot1d_var1(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;

 *--------------------------------------------------------------------------*/

void gnuplot_plot1d_var1(
    gnuplot_ctrl*       handle,
    double*             d,
    int                 n_point,
    char*               title
) ;

/*----------------------------------------------------------------------------
 * Function :   gnuplot_plot1d_var2()
 * In       :   handle to gnuplot_ctrl, list of dpoints, number of
 *              points in this list, title for the output graph.
 * Out      :   void
 * Job      :   plots out a 2d graph from a list of dpoints.
 * Notice   :
 * Example  :

    gnuplot_ctrl    *h ;
    dpoint          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i].x = (double)(i)/10.0 ;
        d[i].y = d[i].x * d[i].x ;
    }
    gnuplot_plot1d_var2(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;

 *--------------------------------------------------------------------------*/


void gnuplot_plot1d_var2(
    gnuplot_ctrl*       handle,
    dpoint*             d,
    int                 n_points,
    char*               title
) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_plot_slope()
 * In       :   handle to gnuplot_ctrl, slope, y-intercept, graph title.
 * Out      :   void
 * Job      :   Plot out a slope on a gnuplot session
 * Notice   :   The provided slope has an equation of the form:
 *              y = a * x + b
 *              a and b are doubles.
 * Example  :

    gnuplot_ctrl    *   h ;
    double              a, b ;

    h = gnuplot_init() ;
    gnuplot_plot_slope(h, 1.0, 0.0, "unity slope") ;
    sleep(2) ;
    gnuplot_close(h) ;

 *--------------------------------------------------------------------------*/


void gnuplot_plot_slope(
    gnuplot_ctrl*       handle,
    double              a,
    double              b,
    char*               title
) ;


/*----------------------------------------------------------------------------
 * Function :   gnuplot_plot_equation()
 * In       :   handle to gnuplot_ctrl, character string containing the
 *              equation to plot, graph title.
 * Out      :   void
 * Job      :   Plots out a curve of given equation. The general form of
 *              the equation is y = f(x)
 * Notice   :
 * Example  :

        gnuplot_ctrl    *h ;
        char            eq[80] ;

        h = gnuplot_init() ;
        strcpy(eq, "sin(x) * cos(2*x)") ;
        gnuplot_plot_equation(h, eq, "sine wave") ;
        gnuplot_close(h) ;

 *--------------------------------------------------------------------------*/

void gnuplot_plot_equation(
    gnuplot_ctrl*       h,
    char*               equation,
    char*               title
) ;

#endif
/*----------------------------- end of file --------------------------------*/
