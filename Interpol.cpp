#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>

using namespace std;

class Cinterp {

public:
    Cinterp(const string& filename, int pointsNumber);
    ~Cinterp();
    void interpolateAndPlot();

private:
    vector<double> x_values, y_values;  // Vectors to store the x and y data points
    int dim, pointsNumber;     // Number of data points and number of interpolation points
    double xmin, xmax, ymin, ymax;  // Bounds of the data

    void readData(const string& filename);  // Reads data from the input file
    void generateLinearInterpolationData();  // Generates data for linear interpolation
    void generateCubicSplineData();  // Generates data for cubic spline interpolation
    void generateLagrangeInterpolationData();  // Generates data for Lagrange interpolation
    void generateGnuplotScript();  // Generates a Gnuplot script to plot the data
};

// Constructor initializes the class and reads the data from the given file
Cinterp::Cinterp(const string& filename, int pointsNumber) : pointsNumber(pointsNumber) {
    readData(filename);
    xmin = x_values.front();
    xmax = x_values.back();
}
// Destructor
Cinterp::~Cinterp(){

}
// Function to read data from a file and store in vectors x_values and y_values
void Cinterp::readData(const string& filename) {
    ifstream fichier(filename);
    string ligne;

    if (!fichier.is_open()) {  // Error handling if the file cannot be opened
        cerr << "Can't open file" << endl;
        exit(EXIT_FAILURE);
    }

    fichier >> dim;  // Read the number of data points
    x_values.resize(dim);
    y_values.resize(dim);

    ymin = numeric_limits<double>::max();  // Initialize ymin and ymax
    ymax = numeric_limits<double>::lowest();
    getline(fichier, ligne);  // Ignore the first line after the number of points

    // Loop to read each (x, y) data point from the file
    for (int i = 0; i < dim; ++i) {
        getline(fichier, ligne);
        stringstream ss(ligne);
        char virg;
        ss >> x_values[i] >> virg >> y_values[i];
        if (ymin > y_values[i]) ymin = y_values[i];  // Update ymin
        if (ymax < y_values[i]) ymax = y_values[i];  // Update ymax
    }
    fichier.close();
}

// Function to generate linear interpolation data and save to a file
void Cinterp::generateLinearInterpolationData() {
    ofstream linearData("linear_data.txt");

    // Interpolate at evenly spaced points between xmin and xmax
    for (int i = 0; i < pointsNumber; ++i) {
        double x = xmin + i * (xmax - xmin) / (pointsNumber - 1);
        int j = 0;
        // Find the interval where the x value falls
        while (j < dim - 1 && x_values[j + 1] < x) ++j;
        double t = (x - x_values[j]) / (x_values[j + 1] - x_values[j]);
        // Linear interpolation formula
        double y = (1 - t) * y_values[j] + t * y_values[j + 1];
        linearData << x << " " << y << endl;  // Write the interpolated data to a file
    }
    linearData.close();
}

// Function to generate cubic spline interpolation data and save to a file
void Cinterp::generateCubicSplineData() {
    vector<double> a(dim), b(dim), c(dim), d(dim);
    vector<double> h(dim - 1), alpha(dim - 1), l(dim), mu(dim), z(dim);

    // Calculate coefficients for the cubic spline
    for (int i = 0; i < dim - 1; ++i) {
        h[i] = x_values[i + 1] - x_values[i];
        if (i == 0) {
            alpha[i] = (3.0 / h[i]) * (y_values[i + 1] - y_values[i]);
        } else {
            alpha[i] = (3.0 / h[i]) * (y_values[i + 1] - y_values[i]) - (3.0 / h[i - 1]) * (y_values[i] - y_values[i - 1]);
        }
    }

    // Forward pass to solve the tridiagonal system
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;
    l[dim - 1] = 1;
    z[dim - 1] = 0;
    c[dim - 1] = 0;
    for (int i = 1; i < dim - 1; ++i) {
        l[i] = 2 * (x_values[i + 1] - x_values[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    // Backward pass to solve for c, b, d coefficients
    l[dim - 1] = 1;
    c[dim - 1] = 0;
    for (int j = dim - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y_values[j + 1] - y_values[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    // Interpolate at evenly spaced points using the cubic spline formula
    ofstream splineData("spline_data.txt");
    for (int i = 0; i < pointsNumber; ++i) {
        double x = xmin + i * (xmax - xmin) / (pointsNumber - 1);
        int j = 0;
        while (j < dim - 1 && x_values[j + 1] < x) ++j;
        double dx = x - x_values[j];
        double y = y_values[j] + b[j] * dx + c[j] * dx * dx + d[j] * dx * dx * dx;
        splineData << x << " " << y << endl;
    }
    splineData.close();
}

// Function to generate Lagrange interpolation data and save to a file
void Cinterp::generateLagrangeInterpolationData() {
    ofstream lagrangeData("lagrange_data.txt");

    // Lambda function to calculate the Lagrange interpolation polynomial
    auto lagrange = [this](double x) {
        double result = 0;
        for (int i = 0; i < dim; ++i) {
            double term = y_values[i];
            for (int j = 0; j < dim; ++j) {
                if (i != j) {
                    term *= (x - x_values[j]) / (x_values[i] - x_values[j]);
                }
            }
            result += term;
        }
        return result;
    };

    // Interpolate at evenly spaced points using the Lagrange formula
    for (int i = 0; i < pointsNumber; ++i) {
        double x = xmin + i * (xmax - xmin) / (pointsNumber - 1);
        double y = lagrange(x);
        lagrangeData << x << " " << y << endl;
    }
    lagrangeData.close();
}

// Function to generate a Gnuplot script to plot the data and interpolations
void Cinterp::generateGnuplotScript() {
    ofstream gnuplotScript("plot.gp");
    gnuplotScript << "set terminal png size 800,600\n";
    gnuplotScript << "set output 'interpolation.png'\n";
    gnuplotScript << "set title 'Interpolation'\n";
    gnuplotScript << "set xlabel 'X'\n";
    gnuplotScript << "set ylabel 'Y'\n";
    gnuplotScript << "set grid\n";
    // Plot original data, linear interpolation, cubic spline, and Lagrange interpolation
    gnuplotScript << "plot 'data.txt' using 1:2 with points pt 7 title 'Data', "
                  << "'linear_data.txt' with lines title 'Linear Interpolation', "
                  << "'spline_data.txt' with lines title 'Cubic Spline', "
                  << "'lagrange_data.txt' with lines title 'Lagrange Interpolation', "
                  << "'data.txt' using 1:2 smooth csplines title 'Gnuplot Cubic Spline'\n";
    gnuplotScript.close();
}

// Function to perform the interpolation and plot the results using Gnuplot
void Cinterp::interpolateAndPlot() {
    generateLinearInterpolationData();  // Generate linear interpolation data
    generateCubicSplineData();  // Generate cubic spline data
    generateLagrangeInterpolationData();  // Generate Lagrange interpolation data
    generateGnuplotScript();  // Generate the Gnuplot script
    system("gnuplot plot.gp");  // Execute the Gnuplot script to generate the plot
}

// Main function
int main() {
    
        Cinterp interp("data.txt", 1000);  // Initialize with input data file and number of interpolation points
        interp.interpolateAndPlot();

    return 0;
}