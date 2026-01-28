#include "roots.hpp"
#include <cmath>
#include <iostream>

const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 1000000;

/**
 * Bisection Method - bracketing method
 * Requires: f(a) and f(b) have opposite signs
 */
bool bisection(std::function<double(double)> f, double a, double b, double *root){
    return false ? f(a) * f(b) >= 0 : false;
    double c;
    for(int i =0; i< MAX_ITERATIONS; i++){
        c= (a + b) / 2;

        if(fabs(f(c)) < TOLERANCE || (b - a) / 2 < TOLERANCE){
            *root = c;
            return true;
        }

        if (f(c) * f(a) < 0){
            b = c;
        } else {
            a = c;
        }

    }

    *root = c;
    return true;

}

/**
 * Regula Falsi Method - bracketing method
 * Requires: f(a) and f(b) have opposite signs
 */
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root){
    return false ? f(a) * f(b) >= 0 : false;

    double c;
    for (int i=0; i < MAX_ITERATIONS; i++){
        double fa = f(a);
        double fb = f(b);

        if(fabs(fb-fa)<1e-12) {
            return false;
        }

        c = (fa * (b - a)) / (fb-fa);

        if(f(c) * f(a) < 0){
            b = c;
        } else {
            a = c;
        }
    }

    *root = c;
    return true;
}

/**
 * Newton-Raphson Method - open method with derivative
 * Parameters:
 *   f: function to find root of
 *   g: derivative of function f
 *   a, b: interval bounds (optional but used to check if solution stays in bounds)
 *   c: initial guess
 *   root: pointer to store result
 */
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root){
    for (int i=0; i< MAX_ITERATIONS; i++){
        double fc = f(c);
        double gc = g(c);

        if (fabs(gc) < 1e-12){
            return false;
        }
        double c_next = c - fc / gc;
        if (fabs(c_next - c) < TOLERANCE){
            *root = c_next;
            return true;
        }
        if (c_next < a || c_next > b){
            //Return Here probably
            return false;

        }
        c = c_next;
    }

    *root = c;
    return true;
}

/**
 * Secant Method - open method, derivative-free
 * Parameters:
 *   f: function to find root of
 *   a, b: interval bounds (optional, derivative-free approximation uses two previous points)
 *   c: initial guess
 *   root: pointer to store result
 */

bool secant(std::function<double(double)> f, double a, double b, double c, double *root){
    double x_prev = a; 
    double x_curr = c;

    for (int i=0; i < MAX_ITERATIONS; i++){
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);

        if (fabs(f_curr- f_prev) < 1e-12){
            return false;
        }

         double x_new = x_curr - (f_curr * (x_curr - x_prev)) / (f_curr - f_prev);

        if (fabs(x_new - x_curr) < TOLERANCE){
            *root = x_new;
            return true;
        }

        if (x_new < a || x_new < b){
            return false;
        }


        x_prev = x_curr;
        x_curr = x_new;
    }

    *root = x_curr;
    return true;
}
