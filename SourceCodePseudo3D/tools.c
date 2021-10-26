#include <stdio.h>
#include <math.h>

double ShiftXY(double x1, double x2) {

    if (x2 > x1) {
        return x2 - 1000.0;
    }
    else if (x2 < x1) {
        return x2 + 1000.0;
    }
    return 0;
}

int trueorfalse(double x1, double y1, double x2, double y2, double a1, double a2, double d, double as, double al) {

    double abs_angle;
    double abs_dista;
    abs_angle = fabs(a1 - a2);
    if (abs_angle < as || abs_angle > al) {
        abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);
        if (abs_dista < d) {
            return 1;
        }
        else if (abs_dista > d*9){
            if (fabs(x1-x2)>900){
                x2 = ShiftXY(x1, x2);
                }
            if (fabs(y1-y2)>900){
                y2 = ShiftXY(y1, y2);
                }

            abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);
            
            if (abs_dista < d) {
                return 1;
            }
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
    return 0;
}


int trueorfalsetrans(double x1, double y1, double x2, double y2, double a1, double a2, double d, double as, double al) {

    double abs_angle0;
    double abs_angle360;
    double abs_dista;
    const double PI = acos(-1);

    abs_angle0 = fabs(a1 - (a2+PI/2));
    abs_angle360 = fabs(fabs(a1 - (a2 + PI / 2))-2*PI);

    if (abs_angle0 < as || abs_angle360 < as) {
        abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);
        if (abs_dista < d) {
            return 1;
        }
        else if (abs_dista > d * 9) {
            if (fabs(x1 - x2) > 900) {
                x2 = ShiftXY(x1, x2);
            }
            if (fabs(y1 - y2) > 900) {
                y2 = ShiftXY(y1, y2);
            }
            abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);
            if (abs_dista < d) {
                return 1;
            }
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
    return 0;
}

int collision(double x1, double y1, double x2, double y2, double d) {

    double abs_dista;
 
    abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);
    if (abs_dista < d) {
        return 1;
    }
    else if (abs_dista > d * 9) {
        if (fabs(x1 - x2) > 900) {
            x2 = ShiftXY(x1, x2);
        }
        if (fabs(y1 - y2) > 900) {
            y2 = ShiftXY(y1, y2);
        }

        abs_dista = pow(x1 - x2, 2) + pow(y1 - y2, 2);

        if (abs_dista < d) {
            return 1;
        }
    }
    else {
        return 0;
    }
    return 0;
}



int main() {
    /*int b; */
    /* x1, y1, x2, y2, a1, a2, d, as, al */
    /*b = trueorfalse(2.5, 500, 999.6, 500, 1, 1, 9, 0.5, 5.0);*/
    /*printf("%i", b);*/
    return 0;
}