#include <cmath>
#include <iostream>
#include "roots.hpp"

const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 1000000;

/**
 * Bisection Method - bracketing method
 * Requires: f(a) and f(b) have opposite signs
 */
bool bisection(std::function<double(double)> f, double a, double b, double *root)
{

    if(!root) return false;

    double fa = f(a);
    double fb = f(b);

    if (fa == 0.0) { *root = a; return true;}

    if (fb == 0.0) { *root = b; return true;}

    if (fa * fb > 0) return false;


    double c;
    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        c = 0.5 * (a + b);
        double fc = f(c);

        if (std::fabs(fc) < TOLERANCE || std::fabs(b - a) < TOLERANCE)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    *root = 0.5 * (a + b);
    return true;
}

/**
 * Regula Falsi Method - bracketing method
 * Requires: f(a) and f(b) have opposite signs
 */
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root)
{

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) return false; // No sign change, no guarantee of root
    
    double c = a;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        if (fabs(fb - fa) < 1e-12) return false;

        c = a - (fa * (b - a)) / (fb - fa);

        if (c <= a || c >= b) {
            c = 0.5 * (a + b);  // Fallback to bisection
        }

        double fc = f(c);
        
        if (fabs(fc) < TOLERANCE || fabs(b - a) < TOLERANCE) {
            *root = c;
            return true;
        }

        if (fc * fa < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return false;
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
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root)
{
    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        double fc = f(c);
        double gc = g(c);

        if (fabs(gc) < 1e-12) return false;
        
        double c_next = c - fc / gc;

        if (c_next < a || c_next > b) return false;

        if (fabs(c_next - c) < TOLERANCE)
        {
            *root = c_next;
            return true;
        }

        c = c_next;
    }

    return false;
}

/**
 * Secant Method - open method, derivative-free
 * Parameters:
 *   f: function to find root of
 *   a, b: interval bounds (optional, derivative-free approximation uses two previous points)
 *   c: initial guess
 *   root: pointer to store result
 */

bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{
    double x_prev = a;
    double x_curr = c;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);

        if (std::fabs(f_curr - f_prev) < 1e-12) return false;
        

        double x_new = x_curr - (f_curr * (x_curr - x_prev)) / (f_curr - f_prev);

        if (x_new < a || x_new > b) {
            x_new = 0.5 * (a + b);  // Fallback to bisection
        }
        
        double f_new = f(x_new);

        if (std::fabs(x_new - x_curr) < TOLERANCE || std::fabs(f_new) < TOLERANCE)
        {
            *root = x_new;
            return true;
        }


        x_prev = x_curr;
        x_curr = x_new;
    }
    *root = x_curr;
    return false;
}

// ctest --test-dir build -C Debug --output-on-failure