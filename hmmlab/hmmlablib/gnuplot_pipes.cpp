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

#include "gnuplot_pipes.h"

#define GP_TITLE_SIZE   80
#define GP_EQ_SIZE      512
#define GP_LINE_SIZE    80

int check_X_display(void)
{
    char*       display ;

    display = getenv("DISPLAY") ;
    if(display == NULL) {
        fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
        return 1 ;
    } else {
        return 0 ;
    }
}

gnuplot_ctrl* gnuplot_init(void)
{
    gnuplot_ctrl*   handle ;

    if(check_X_display()) {
        return NULL ;
    }

    handle = (gnuplot_ctrl*) malloc(sizeof(gnuplot_ctrl)) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;

    handle->gnucmd = popen("gnuplot", "w") ;
    if(handle->gnucmd == NULL) {
        fprintf(stderr, "cannot find gnuplot program: check $path\n") ;
        free(handle) ;
        return NULL ;
    }
    return handle;
}

void gnuplot_close(gnuplot_ctrl* handle)
{
    int     i ;
    if(check_X_display()) {
        return ;
    }
    if(handle->ntmp) {
        for(i = 0 ; i < handle->ntmp ; i++) {
            remove(handle->to_delete[i]) ;
        }
    }
    if(pclose(handle->gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    free(handle) ;
    return ;
}

void gnuplot_cmd(
    gnuplot_ctrl*   handle,
    char*           cmd,
    ...)
{
    assert(strlen(cmd) < GP_LINE_SIZE);
    va_list ap ;
    char    local_cmd[GP_CMD_SIZE];

    va_start(ap, cmd);
    vsprintf(local_cmd, cmd, ap);
    va_end(ap);

    strcat(local_cmd, "\n");

    fputs(local_cmd, handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}

void gnuplot_setstyle(gnuplot_ctrl* h, char* plot_style)
{
    if(strcmp(plot_style, "lines") &&
       strcmp(plot_style, "points") &&
       strcmp(plot_style, "linespoints") &&
       strcmp(plot_style, "impulses") &&
       strcmp(plot_style, "dots") &&
       strcmp(plot_style, "steps") &&
       strcmp(plot_style, "errorbars") &&
       strcmp(plot_style, "boxes") &&
       strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        (void)strcpy(h->pstyle, "points") ;
    } else {
        (void)strcpy(h->pstyle, plot_style) ;
    }
    return ;
}

void gnuplot_set_xlabel(gnuplot_ctrl* h, char* label)
{
    char    cmd[GP_CMD_SIZE] ;

    (void)sprintf(cmd, "set xlabel \"%s\"", label) ;
    gnuplot_cmd(h, cmd) ;
    return ;
}

void gnuplot_set_ylabel(gnuplot_ctrl* h, char* label)
{
    char    cmd[GP_CMD_SIZE] ;

    (void)sprintf(cmd, "set ylabel \"%s\"", label) ;
    gnuplot_cmd(h, cmd) ;
    return ;
}

void gnuplot_resetplot(gnuplot_ctrl* h)
{
    int     i ;
    if(h->ntmp) {
        for(i = 0 ; i < h->ntmp ; i++) {
            remove(h->to_delete[i]) ;
        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}

void gnuplot_plot1d_var1(
    gnuplot_ctrl*       handle,
    double*             d,
    int                 n_point,
    char*               title
)
{
    int         i ;
    FILE*       tmp ;
    char*       name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_LINE_SIZE] ;

    /* can we open one more temporary file? */
    if(handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }

    /* Open temporary file for output   */
    if((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }
    if((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Store file name in array for future deletion */
    (void)strcpy(handle->to_delete[handle->ntmp], name) ;
    handle->ntmp ++ ;

    /* Write data to this file  */
    for(i = 0 ; i < n_point ; i++) {
        (void)fprintf(tmp, "%f\n", d[i]) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;

    /* Command to be sent to gnuplot    */
    if(handle->nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }

    if(title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, handle->pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, handle->pstyle) ;
    }

    /* send command to gnuplot  */
    gnuplot_cmd(handle, line) ;
    handle->nplots++ ;
    return ;
}

void gnuplot_plot1d_var2(
    gnuplot_ctrl*       handle,
    dpoint*             d,
    int                 n_points,
    char*               title
)
{
    int         i ;
    FILE*       tmp ;
    char*       name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_LINE_SIZE] ;

    /* can we open one more temporary file? */
    if(handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }

    /* Open temporary file for output   */
    if((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }
    if((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Store file name in array for future deletion */
    (void)strcpy(handle->to_delete[handle->ntmp], name) ;
    handle->ntmp ++ ;

    /* Write data to this file  */
    for(i = 0 ; i < n_points ; i++) {
        (void)fprintf(tmp, "%f %f\n", d[i].x, d[i].y) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;

    /* Command to be sent to gnuplot    */
    if(handle->nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }

    if(title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, handle->pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, handle->pstyle) ;
    }

    /* send command to gnuplot  */
    gnuplot_cmd(handle, line) ;
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_slope(
    gnuplot_ctrl*       handle,
    double              a,
    double              b,
    char*               title
)
{
    char    stitle[GP_TITLE_SIZE] ;
    char    cmd[GP_CMD_SIZE] ;

    if(title == NULL) {
        (void)strcpy(stitle, "no title") ;
    } else {
        (void)strcpy(stitle, title) ;
    }

    if(handle->nplots > 0) {
        (void)sprintf(cmd, "replot %f * x + %f title \"%s\" with %s",
                      a, b, title, handle->pstyle) ;
    } else {
        (void)sprintf(cmd, "plot %f * x + %f title \"%s\" with %s",
                      a, b, title, handle->pstyle) ;
    }
    gnuplot_cmd(handle, cmd) ;
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_equation(
    gnuplot_ctrl*       h,
    char*               equation,
    char*               title
)
{
    char    cmd[GP_CMD_SIZE];
    char    plot_str[GP_EQ_SIZE] ;
    char    title_str[GP_TITLE_SIZE] ;

    if(title == NULL) {
        (void)strcpy(title_str, "no title") ;
    } else {
        (void)strcpy(title_str, title) ;
    }
    if(h->nplots > 0) {
        (void)strcpy(plot_str, "replot") ;
    } else {
        (void)strcpy(plot_str, "plot") ;
    }

    (void)sprintf(cmd, "%s %s title \"%s\" with %s",
                  plot_str, equation, title_str, h->pstyle) ;
    gnuplot_cmd(h, cmd) ;
    h->nplots++ ;
    return ;
}
